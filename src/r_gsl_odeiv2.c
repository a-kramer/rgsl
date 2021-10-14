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

#define FREE_ON_SUCCESS 1
#define KEEP_ON_SUCCESS 2

// SEXP stands for S-Expression, and it can be any R data object (or function)
// in this program, we'll only use data from R
// and SEXP is so weird to read, so...
typedef SEXP Rdata;

typedef int(*jacp)(double, const double [], double *, void *);


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
 gsl_matrix *Yout) /* time span vector: initial time, increment, final time */
{
  assert(driver);
  assert(tspan);
  int nt=tspan->size;
  int ny=(int) sys.dimension;
  gsl_vector *y=gsl_vector_alloc(ny);
  gsl_vector_memcpy(y,y0);
  
  gsl_vector_view Yout_row;
  int j;
  double t,tf,dt;
  double dfdt[ny];
  int status;
  /* initialize time-point 0 values */
  t=gsl_vector_get(tspan,0);
  Yout_row = gsl_matrix_row(Yout,0);
  gsl_vector_memcpy(&(Yout_row.vector),y0);
  
  for (j=1; j<nt; j++){
    tf=gsl_vector_get(tspan,j);
    printf("[%s] t=%f -> %f\n",__func__,t,tf);
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
  return status;
}

/* This prgram loads an ODE model, specified for the `gsl_odeiv2`
   library. The initial value problems are specified in an hdf5 file,
   intended for use in systems biology applications. So, some of the
   terminology in the expected hdf5 _data_ file is vaguely related to
   biological systems. All command line arguments are optional and
   names of files are guessed, based on the contents of the current
   working directoy. The hdf5 fiel is expected to have a group called
   "data", this group shall contain hdf5 DATASETS with ATTRIBUTES that
   describe initial value problems suitable to replicate these
   datasets: InitialValue, time, and index. Currenty, the model is
   parameterized using the value of the DATASET *mu*, in the GROUP
   called *prior*. 

   The possible command line options are documented in the [README.md](../README.md) .
*/
Rdata /* `EXIT_SUCESS` if all files are found and integration succeeds, default `abort()` signal otherwise.*/
r_gsl_odeiv2(Rdata model_name, Rdata tspan, Rdata y0, Rdata p){ 
  int i=0;
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
  size_t ny=length(y0);
  size_t np,N;
  if (isMatrix(p)){
    np=nrows(p);
    N=ncols(p);
  } else if (isVector(p)) {
    np=length(p);
    N=1;
  }
  gsl_vector_view t=gsl_vector_view_array(REAL(tspan),nt);
  gsl_vector_view initial_value = gsl_vector_view_array(REAL(y0),ny);
  gsl_matrix_view ode_parameter = gsl_matrix_view_array(REAL(p),N,np);
  
  // load system from file and test it
  jacp dfdp;
  gsl_odeiv2_system sys = load_system(CHAR(STRING_ELT(model_name,0)), ny, REAL(p), &dfdp);
  
  //test_evaluation(sys,dfdp,&(initial_value.vector),&(ode_parameter.vector));

  const gsl_odeiv2_step_type * T=gsl_odeiv2_step_msbdf;
  gsl_odeiv2_driver *driver=gsl_odeiv2_driver_alloc_y_new(&sys,T,h,abs_tol,rel_tol);
  // most CPU work happens here:
  Rdata Y = PROTECT(alloc3DArray(REALSXP,ny,nt,N));
  double *ydata;
  size_t nyt=nt*ny;
  gsl_matrix_view y;
  int status;
  for (i=0;i<N;i++){
    printf("[%s] solving %i of %li.\n",__func__,i,N);
    ydata = REAL(Y);
    y=gsl_matrix_view_array(&(ydata[i*nyt]),nt,ny);
    sys.params = gsl_matrix_ptr(&(ode_parameter.matrix),i,0);
    status=simulate_timeseries
      (sys,
       driver,
       &(initial_value.vector),
       &(t.vector),
       &(y.matrix));
    assert(status==GSL_SUCCESS);
  }
  gsl_odeiv2_driver_free(driver);
  UNPROTECT(1);
  return Y;
}
