#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <dlfcn.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define FREE_ON_SUCCESS 1
#define KEEP_ON_SUCCESS 2
#define MATCH 0
#define NO_DIFFERENCE 0
// SEXP stands for S-Expression, and it can be any R data object (or function)
// in this program, we'll only use data from R
// and SEXP is so weird to read, so...
typedef SEXP Rdata;

typedef int(*jacp)(double, const double [], double *, void *);
typedef enum {DIAG,MATVEC} tf_t;

typedef struct {
  tf_t type;
  int length_b;
  double *A;
  double *b;
  gsl_vector *y; /* result vector */
} affine_tf;
  
typedef struct {
  char *name;
  int nt;
  double *time;
  affine_tf *state;
  affine_tf *par;
} event_t;

Rdata from_list(Rdata List, const char *name){
  assert(isVector(List));
  int i;
  int N=length(List);
  Rdata names = GET_NAMES(List);
  Rdata E;
  for (i=0;i<length(names);i++){
    if (strcmp(CHAR(STRING_ELT(names,i)),name) == MATCH){
      E = VECTOR_ELT(List,i);
    }
  }
  return E;  
}

/* creates an affine transformation struct from R objects, Rdata
   objects need to be kept alive as A and b arfe taken from R via
   pointers. Some memory is allocated to store an intermediate result. */
affine_tf* /* an affine transformation (linear with offset) map: x -> A*x+b */
affine_transformation(
 Rdata A,/*a series of matrices, possibly just a set of diagonals*/
 Rdata b)/*a series of offsets*/
{
  assert(IS_NUMERIC(A));
  assert(IS_NUMERIC(b));
  Rdata dA=GET_DIM(A);
  Rdata db=GET_DIM(b);
  int n=length(dA);
  int l=length(db);
  int *dim=INTEGER(dA);
  assert(dim[0] == INTEGER(db)[0]);
  affine_tf *L=malloc(sizeof(affine_tf));
  if (n==3 && dim[0]==dim[1]){
    L->type=MATVEC;
    assert(l==n-1);
  } else if (n==2) {
    L->type=DIAG;
    assert(l==n);
  }
  L->A=REAL(A);
  L->b=REAL(b);
  L->length_b=dim[0];
  L->y=gsl_vector_alloc(dim[0]);
  return L;
}

/* This function applies a affine transformation to the vector z. We
   assume that A and b and z in the input are all of sufficient size:
   n×n for matrices and n for vectors. A can be n-sized as well, if A
   is a diagonal matrix, then we only store the diagonal. n is stored
   in the transformation structure as length_b.  */
int /* the returned status of the gsl operations */
apply_tf(affine_tf *L, /* a transformation struct: A and b are cast to gsl_vectors here*/
	 double *z)/* an array of size n, it is updated using L */{
  assert(L);
  assert(z);
  int n=L->length_b;
  int status=GSL_SUCCESS;
  gsl_vector_view x=gsl_vector_view_array(z,n);
  gsl_vector_view b=gsl_vector_view_array(L->b,n);
  if (L->type == MATVEC){
    /* y <- A*x + b */
    gsl_matrix_view A=gsl_matrix_view_array(L->A,n,n);
    status|=gsl_blas_dgemv(CblasNoTrans, 1.0, &(A.matrix), &(x.vector), 0.0, L->y);
  } else if (L->type == DIAG){
    /* y <- diag(A)*x + b */
    gsl_vector_view a=gsl_vector_view_array(L->A,n);
    status|=gsl_vector_memcpy(L->y,&(x.vector));
    status|=gsl_vector_mul(L->y,&(a.vector));
  } else {
    fprintf(stderr,"[%s] unknown transformation type %i.\n",__func__,L->type);
    abort();
  }
  assert(status==GSL_SUCCESS);
  gsl_vector_add(L->y,&(b.vector));
  gsl_vector_memcpy(&(x.vector),L->y);
  return status;
}

event_t* event_from_R(Rdata E){
  int j;
  event_t *event=malloc(sizeof(event_t));
  
  event->state=affine_transformation(from_list(E,"A"),from_list(E,"a"));
  event->par=affine_transformation(from_list(E,"B"),from_list(E,"b"));
  
  Rdata time = from_list(E,"time");
  int lt=length(time);
  event->time = REAL(time);
  event->nt=lt;
  printf("[%s] event times:",__func__);
  for (j=0;j<lt;j++) printf("\t%g",event->time[j]);
  putchar('\n');
  return event;
}

/*This function allocates memory and concatenates two strings in that
  memory. It is used to make function names (this is for loading the model
  functions by name from a shared library). The model function names
  have this pattern: `MODEL_vf`, `MODEL_jac`, `MODEL_jacp`.*/
char* /* string with model_name and suffix (free after loading the function) */
model_function(const char *model_name, /* the base name of the model */
 const char *suffix) /* suffix, usually `"_vf"` or `"_jac"` */
{
  assert(model_name);
  size_t size=(strlen(model_name)+strlen(suffix)+1);
  assert(size);
  char *f=malloc(sizeof(char)*size);
  assert(f);
  strcat(strcpy(f,model_name),suffix);
  fprintf(stderr,"[%s] «%s»\n",__func__,f); fflush(stderr);
  return f;
}


/* Loads a function from an `.so` file, using `dlsym()`. If dlsym
   fails to find the file `abort()` is called. Optionally, this
   function frees the storage assosiated with the name of the
   function.*/
void *load_or_exit(void *lib, /* file pointer, previously opened via `dlopen()` */
 char *name, /* function to be loaded from file */
 int opt) /* whether to call `free()` on `name` (either: `KEEP_ON_SUCCESS` or `FREE_ON_SUCCESS`). */
{
  assert(lib && name);
  void *symbol=dlsym(lib,name);
  if (symbol) {
    fprintf(stderr,"[%s] «%s» loaded successfully.\n",__func__,name);
    if (opt==FREE_ON_SUCCESS){
      free(name);
    }
  }else{
    fprintf(stderr,"[%s] %s\n",__func__,dlerror());
    abort();
  } 
  return symbol;
}


/* Loads the ODE system from an `.so` file, the file is given by name,
   the returned structure is intended for the `gsl_odeiv2` library of
   solvers. The jacobian dfdx is loaded alongside the right hand side;
   `gsl_odeiv2_system` has no slot for the parameter derivative `dfdp`,
   this is returned as a function pointer instead. */
gsl_odeiv2_system /* the system structure, see gsl documentation. */
load_system(const char *model_name, /* the file-name will be constructed from this name, possibly from @link first_so@ */
 size_t n, /* number of state variables */
 double *p, /* default parameter vector */
 jacp *dfdp) /* [output] additional return value: a pointer to the parameter derivative (matrix) function. */
{
  char *so=model_function(model_name,".so");
  char *local_so = model_function("./",so);
  void *lib=dlopen(local_so,RTLD_LAZY);
  free(local_so);
  void *f,*dfdy;
  //dfdp=malloc(sizeof(jacp*));
  char *symbol_name; // symbol name in .so
  if (lib){
    symbol_name=model_function(model_name,"_vf");
    f=load_or_exit(lib,symbol_name,FREE_ON_SUCCESS);
    symbol_name=model_function(model_name,"_jac");
    dfdy=load_or_exit(lib,symbol_name,FREE_ON_SUCCESS);
    symbol_name=model_function(model_name,"_jacp");
    *dfdp=load_or_exit(lib,symbol_name,FREE_ON_SUCCESS);
  } else {
    fprintf(stderr,"[%s] library «%s» could not be loaded: %s\n",__func__,so,dlerror());
    abort();
  }
  gsl_odeiv2_system sys={f,dfdy,n,p};
  free(so);
  fprintf(stderr,"[%s] ode system created.\n",__func__); fflush(stderr);
  return sys;
}



/* Intergrates the system `sys` using the specified `driver` and
   simulation instructions `sim` (an array of structs, one element per
   simulation). The results are saved to an hdf5 file and also printed
   to standard output. */
int /* error code if any */
simulate_timeseries(const gsl_odeiv2_system sys, /* the system to integrate */
 gsl_odeiv2_driver* driver, /* the driver that is used to integrate `sys` */
 const gsl_vector *y0, /* initial value */
 const gsl_vector *tspan, /* a vector of time-points */
 const event_t *event, /*a struct array with scheduled events */
 gsl_matrix *Yout) /* time span vector: initial time, increment, final time */
{
  assert(driver);
  assert(tspan);
  int nt=tspan->size;
  int ny=(int) sys.dimension;
  gsl_vector *y=gsl_vector_alloc(ny);
  gsl_vector_memcpy(y,y0);
  
  gsl_vector_view Yout_row;
  int i=0;
  int j;
  double t,tf,te;
  double dfdt[ny];
  int status;
  /* initialize time-point 0 values */
  t=gsl_vector_get(tspan,0);
  Yout_row = gsl_matrix_row(Yout,0);
  gsl_vector_memcpy(&(Yout_row.vector),y0);
  
  for (j=1, i=0; j<nt; j++){
    tf=gsl_vector_get(tspan,j);
    if (i<event->nt && event->time[i] < tf) {
      te=event->time[i];
      printf("[%s] a scheduled event occurs at t=%g.\n",__func__,te);
      status=gsl_odeiv2_driver_apply(driver, &t, tf, y->data);
      assert(status==GSL_SUCCESS);
      apply_tf(event->state,y->data);
      apply_tf(event->par,(double*) sys.params);
      status=gsl_odeiv2_driver_reset(driver);
      assert(status==GSL_SUCCESS);
      i++;
    }
    //printf("[%s] t=%f -> %f\n",__func__,t,tf);
    status=gsl_odeiv2_driver_apply(driver, &t, tf, y->data);
    //report any error codes to the R user
    switch (status){
    case GSL_EMAXITER:
      error("[%s] time_point %li: maximum number of steps reached.\n\t\tfinal time: %.10g (short of %.10g)",__func__,j,t,tf);
      break;
    case GSL_ENOPROG:
      error("[%s] time_point %li: step size dropeed below set minimum.\n\t\tfinal time: %.10g (short of %.10g)",__func__,j,t,tf);
      break;
    case GSL_EBADFUNC:
      error("[%s] time_point %li: bad function.\n\t\tfinal time: %.10g (short of %.10g)",__func__,j,t,tf);
      break;
    case GSL_SUCCESS:
      Yout_row = gsl_matrix_row(Yout,j);
      gsl_vector_memcpy(&(Yout_row.vector),y);
      break;
    default:
      error("[%s] unhandled error code: %#x\n",__func__,status);
      abort();
    }
    if (status!=GSL_SUCCESS) break;
  }
  gsl_odeiv2_driver_reset(driver);
  return status;
}





/* This prgram loads an ODE model, specified for the `gsl_odeiv2`
   library. It simulates the model for each column of initial
   conditions y0 and parameters p.
```
   y'=f(y,t;p)
   y0=y(t[0])
```
   if y0 and p are matrices with m columns (both), then the model is
   simulated m times, in parallel if possible.
*/
Rdata /* the trajectories as a three dimensional array */
r_gsl_odeiv2(
 Rdata model_name, /* a string */
 Rdata tspan, /* a vector of output times with tspan[0] crresponding to y0 */
 Rdata y0, /* initial conditions at tspan[0], can be a matrix*/
 Rdata p, /* parameter matrix */
 Rdata event){ 
  int i,j;
  double abs_tol=1e-6,rel_tol=1e-5,h=1e-3;
  assert(IS_CHARACTER(model_name));
  assert(IS_NUMERIC(tspan));
  assert(IS_NUMERIC(y0));
  assert(IS_NUMERIC(p));

  /* We will solve the ODE several times: p may be a matrix, with each
     column providing a different parameterisation of the initial
     value problem. We solve the ode one for each parameter vector
     (a column).

     But, R stores matrices column-wise while GSL stores them row-wise.
     To compensate, we switch the roles here.
  */
  size_t nt=length(tspan);
  size_t ny;
  size_t np,N;
  Rdata experiment_names;
  int ExperimentsAreNamed = 0;
  if (isMatrix(p)){
    np=nrows(p);
    N=ncols(p);
    assert(isMatrix(y0));
    ny=nrows(y0);
    assert(ncols(y0)==N);
    experiment_names = GET_COLNAMES(p);
    if (experiment_names != R_NilValue) ExperimentsAreNamed = 1;
  } else if (isVector(p)) {
    np=length(p);
    N=1;
    assert(isVector(y0));
    ny=length(y0);
  }
  printf("[%s] ny=%li, np=%li, N=%li.\n",__func__,ny,np,N);
  gsl_vector_view t=gsl_vector_view_array(REAL(tspan),nt);
  gsl_matrix_view initial_value = gsl_matrix_view_array(REAL(y0),N,ny);
  gsl_matrix_view ode_parameter = gsl_matrix_view_array(REAL(p),N,np);
  
  // load system from file
  jacp dfdp;
  gsl_odeiv2_system sys = load_system(CHAR(STRING_ELT(model_name,0)), ny, REAL(p), &dfdp);
  printf("[%s] system dimension: %li\n",__func__,sys.dimension);

  // check whether events are happening during integration:
  int l,lt;
  double *event_time;
  Rdata E;
  Rdata temp;
  event_t **ev=calloc(N,sizeof(event_t*));
  if (event != R_NilValue && event){
    l=length(event);
    assert(l==N || ExperimentsAreNamed);
    for (i=0;i<l;i++){
      ev[i] = event_from_R(VECTOR_ELT(event, i));
    }
  }
    
  /* Rdata dM=GET_DIM(M); */
  /* assert(IS_NUMERIC(dM) && IS_INTEGER(dM)); */
  /* int M_nd = length(dM); */
  /* printf("[%s] M has %i indices:\n",__func__,M_nd); */
  /* for (i=0;i<M_nd;i++) printf("\t%i",INTEGER(dM)[i]); */
  /* putchar('\n'); */

  /* Rdata dK=GET_DIM(K); */
  /* assert(IS_NUMERIC(dK) && IS_INTEGER(dK)); */
  /* int K_nd = length(dK); */
  /* printf("[%s] K has %i indices:\n",__func__,K_nd); */
  /* for (i=0;i<K_nd;i++) printf("\t%i",INTEGER(dK)[i]); */
  /* putchar('\n'); */

  const gsl_odeiv2_step_type * T=gsl_odeiv2_step_msbdf;
  gsl_odeiv2_driver *driver;

  Rdata Y = PROTECT(alloc3DArray(REALSXP,ny,nt,N));
  double *ydata;
  gsl_vector_view iv_row;
  size_t nyt=nt*ny;
  gsl_matrix_view y;
  int status;
  ydata = REAL(Y);

#pragma omp parallel for private(driver,y,iv_row) firstprivate(sys,ydata,initial_value,t)
  for (i=0;i<N;i++){
    driver=gsl_odeiv2_driver_alloc_y_new(&sys,T,h,abs_tol,rel_tol);
    printf("[%s] solving %i of %li.\n",__func__,i,N);
    y=gsl_matrix_view_array(&(ydata[i*nyt]),nt,ny);
    sys.params = gsl_matrix_ptr(&(ode_parameter.matrix),i,0);
    iv_row=gsl_matrix_row(&(initial_value.matrix),i);
    status=simulate_timeseries
      (sys,
       driver,
       &(iv_row.vector),
       &(t.vector),
       ev[i],
       &(y.matrix));
    assert(status==GSL_SUCCESS);
    gsl_odeiv2_driver_free(driver);
  }
  UNPROTECT(1);
  return Y;
}
