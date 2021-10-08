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

typedef SEXP Rdata;

typedef struct {
  gsl_vector *u;
  gsl_vector *y0;
  gsl_vector *t;
  char *name;
} simulation_t;

typedef int(*jacp)(double, const double [], double *, void *);


/*This function allocates memory and concatenates two strings in that
  memory. It is used to make function names (this is for loading the model
  functions by name from a shared library). The model function names
  have this pattern: `MODEL_vf`, `MODEL_jac`, `MODEL_jacp`.*/
char* /* string with model_name and suffix (free after loading the function) */
model_function(char *model_name, /* the base name of the model */
 char *suffix) /* suffix, usually `"_vf"` or `"_jac"` */
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
load_system(char *model_name, /* the file-name will be constructed from this name, possibly from @link first_so@ */
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
solution_t** /* a structure that contains the trajectories y(t) as well as jacobians Jy, Jp and the y-sensitivity (dy/dp).*/
simulate_timeseries(gsl_odeiv2_system sys, /* the system to integrate */
 gsl_odeiv2_driver* driver, /* the driver that is used to integrate `sys` */
 simulation_t *sim, /*`N` simulation instructions (e.g.: initial value y0)*/
 hsize_t N, /* number of simulations to perform, length of `sim` */
 double *tspan, /* time span vector: initial time, increment, final time */
 gsl_vector *u, /* input vector (a pointer that can be used to change `sys`)*/ 
 gsl_vector *par, /* parameters of the model, a pointer that can be used to change `sys`.*/
 jacp dfdp, /* a function that returns the parameter derivative of the model's right hand side function.*/
 hid_t h5f) /* an hdf5 file opened for writing.*/
{
  assert(driver && sim && N>0);
  assert(tspan);
  assert(sim);
  size_t d=sim[0].y0->size;
  size_t p=par->size;
  gsl_vector *y=gsl_vector_alloc(d);
  double *y_ptr;
  size_t i,j,k;
  solution_t **solution=malloc(sizeof(solution_t*)*N);
  int index[3]={0,0,0};
  double t,tf,dt;
  double *Jy,*Jp,*f;
  double dfdt[d];
  size_t n;
  hid_t g_id;
  int status;
  for (i=0;i<N;i++){
    assert(sim[i].name && strlen(sim[i].name));
    if (h5f) g_id = H5Gcreate2(h5f, sim[i].name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    gsl_odeiv2_driver_reset(driver);
    // fix possible issues with simulation time specification:
    fix_tspan_if_necessary(sim[i].t,tspan);    
    gsl_vector_memcpy(y,sim[i].y0);
    //gsl_vector_add_constant(y,1e-5);
    gsl_vector_memcpy(u,sim[i].u);
     t=tspan[0];
    dt=tspan[1];
    // print table header:
    printf("# Simulation %li: \n",i);
    printf("#t\t");
    for (k=0;k<y->size;k++)
      printf("y%li\t",k);
    printf("\n");
    n=(size_t) ((tspan[2]-tspan[0])/tspan[1]);
    solution[i]=solution_alloc(d,p,n);
    for (j=0;j<n;j++){
      tf=t+dt;
      status=gsl_odeiv2_driver_apply(driver, &t, tf, y->data);
      //report any error codes to the user
      switch (status){
      case GSL_EMAXITER:
	fprintf(stderr,"[%s] simulation %li, time_point %li: maximum number of steps reached.\n",__func__,i,j);
	fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
	break;
      case GSL_ENOPROG:
	fprintf(stderr,"[%s] simulation %li, time_point %li: step size dropeed below set minimum.\n",__func__,i,j);
	fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
	break;
      case GSL_EBADFUNC:
	fprintf(stderr,"[%s] simulation %li, time_point %li: bad function.\n",__func__,i,j);
	fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
	break;
      case GSL_SUCCESS:
	printf("%g\t",t);
	for (k=0;k<y->size;k++)
	  printf("%g\t",gsl_vector_get(y,k));
	printf("\n");
	index[2]=j;
	y_ptr = ndarray_ptr(solution[i]->y,&(index[1]));
	f = ndarray_ptr(solution[i]->f,&(index[1]));
	Jy = ndarray_ptr(solution[i]->Jy,index);
	Jp = ndarray_ptr(solution[i]->Jp,index);
	sys.function(t, y->data, f, par->data);
	sys.jacobian(t, y->data, Jy, dfdt, par->data);
	dfdp(t, y->data, Jp, par->data);
	memcpy(y_ptr,y->data,sizeof(double)*(y->size));
	*ndarray_ptr(solution[i]->t,&(index[2]))=t;
	break;
      default:
	fprintf(stderr,"[%s] unhandled error code: %#x\n",__func__,status);
	abort();
      }
      if (status!=GSL_SUCCESS) break;
    }
    printf("\n\n");
  }
  return solution;
}


/* Intergrates the system `sys` using the specified `driver`, but
   using the solvers steps within the boundaries of the time span. The
   simulation instructions are stored in `sim` (an array of structs,
   one element per simulation). */
solution_t** /* a structure that contains the trajectories y(t).*/
simulate_evolve(gsl_odeiv2_system sys, /* the system to integrate */
 gsl_odeiv2_driver* driver, /* the driver that is used to integrate `sys` */
 simulation_t *sim, /*`N` simulation instructions (e.g.: initial value y0)*/
 hsize_t N, /* number of simulations to perform, length of `sim` */
 double *tspan, /* time span vector: initial time, increment, final time */
 gsl_vector *u, /* input vector (a pointer that can be used to change `sys`)*/ 
 gsl_vector *par, /* parameters of the model, a pointer that can be used to change `sys`.*/
 jacp dfdp, /* a function that returns the parameter derivative of the model's right hand side function.*/
 hid_t h5f) /* an hdf5 file opened for writing.*/
{
  assert(driver && sim && N>0);
  assert(tspan);
  assert(sim);
  size_t d=sim[0].y0->size;
  size_t p=par->size;
  gsl_vector *y=gsl_vector_alloc(d);
  double h=1e-2;
  size_t i,k;
  solution_t **solution=malloc(sizeof(solution_t*)*N);
  int index[3]={0,0,0};
  double t,tf;
  double *Jy,*Jp,*y_ptr,*f;
  double dfdt[d];
  size_t n,n_max=100;
  hid_t g_id;
  int status;
  for (i=0;i<N;i++){// simulations
    assert(sim[i].name && strlen(sim[i].name));
    if (h5f) g_id = H5Gcreate2(h5f, sim[i].name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    gsl_odeiv2_driver_reset(driver);
    // fix possible issues with simulation time specification:
    fix_tspan_if_necessary(sim[i].t,tspan);    
    gsl_vector_memcpy(y,sim[i].y0);
    //gsl_vector_add_constant(y,1e-5);
    gsl_vector_memcpy(u,sim[i].u);
    t=tspan[0];
    tf=tspan[2];
    // print table header:
    printf("# Simulation %li: \n",i);
    printf("#t\t");
    for (k=0;k<y->size;k++)
      printf("y%li\t",k);
    printf("\n");
    // estimate number of points:
    n=0;
    solution[i]=solution_alloc(d,p,n_max);
    while (t<tf){
      if (n==n_max){
	n_max=n_max+100;
	solution_resize(solution[i],n_max);
      }
      
      printf("%g\t",t);
      for (k=0;k<y->size;k++)
	printf("%g\t",gsl_vector_get(y,k));
      printf("\n");
      index[2]=n;
      y_ptr = ndarray_ptr(solution[i]->y,&(index[1]));
      f = ndarray_ptr(solution[i]->f,&(index[1]));
      Jy = ndarray_ptr(solution[i]->Jy,index);
      Jp = ndarray_ptr(solution[i]->Jp,index);
      sys.function(t, y->data, f, par->data);
      sys.jacobian(t, y->data, Jy, dfdt, par->data);
      dfdp(t, y->data, Jp, par->data);
      
      memcpy(y_ptr,y->data,sizeof(double)*(y->size));
      *ndarray_ptr(solution[i]->t,&(index[2]))=t;

      status = gsl_odeiv2_evolve_apply(driver->e, driver->c, driver->s, &sys, &t, tf, &h, y->data);
      //report any error codes to the user
      switch (status){
      case GSL_EMAXITER:
	fprintf(stderr,"[%s] simulation %li, time_point %li: maximum number of steps reached.\n",__func__,i,n);
	fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
	break;
      case GSL_ENOPROG:
	fprintf(stderr,"[%s] simulation %li, time_point %li: step size dropeed below set minimum.\n",__func__,i,n);
	fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
	break;
      case GSL_EBADFUNC:
	fprintf(stderr,"[%s] simulation %li, time_point %li: bad function.\n",__func__,i,n);
	fprintf(stderr,"\t\tfinal time: %.10g (short of %.10g)",t,tf);
	break;
      case GSL_SUCCESS:
	n++;
	break;
      default:
	fprintf(stderr,"[%s] unhandled error code: %#x\n",__func__,status);
	abort();
      }
      if (status!=GSL_SUCCESS) break;
    }
    solution_resize(solution[i],n);
  }
  return solution;
}


/* Interprets a string as a range specification from three values:
   "initial increment final", the values can be separated by spaces or
   colons `:`. */
double* /* a three element array of `double`s */
read_tspan(char *val)/*a string of the form "a:b:c" or "a b c". */{
  int j;
  char *s,*p;
  double *t=calloc(3,sizeof(double));
  if (val){
    s=val;
    p=val;
    for (j=0;j<3;j++){
      if (p[0]==':') p++;
      s=p;
      t[j]=strtod(s,&p);
      if (s==p) break;
    }
  }
  return t;
}

void test_evaluation(gsl_odeiv2_system sys, jacp dfdp, gsl_vector *y0, gsl_vector *par){
  //test evaluation:
  size_t d=sys.dimension;
  int jacp_size=d*par->size;
  double *Jy=malloc(sizeof(double)*(d*d));
  double *Jp=malloc(sizeof(double)*jacp_size);
  double *f=malloc(sizeof(double)*d);
  int i;
  // f
  sys.function(0,y0->data,f,par->data);
  fprintf(stderr,"[%s] test evaluation of right hand side function (f):\n",__func__);
  for (i=0;i<d;i++) fprintf(stderr,"%g ",f[i]);
  fprintf(stderr,"\n\n");
  fflush(stderr);

  //jacp
  dfdp(0,y0->data,Jp,par->data);
  fprintf(stderr,"[%s] test evaluation of jacp (df/dp):\n",__func__);
  for (i=0;i<jacp_size;i++) fprintf(stderr,"%g ",Jp[i]);
  fprintf(stderr,"\n (that was a flat %li × %li matrix)\n",d,par->size);
  fflush(stderr);
  // jac

  sys.jacobian(0,y0->data,Jy,f,par->data);
  fprintf(stderr,"[%s] test evaluation of jacobian (df/dy):\n",__func__);
  for (i=0;i<d*d;i++) fprintf(stderr,"%g ",Jy[i]);
  fprintf(stderr,"\n (that was a flat %li × %li matrix)\n",d,d);
  fflush(stderr);
  free(Jy);
  free(Jp);
  free(f);


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
r_gsl_odeiv(Rdata model_so, Rdata tspan, Rdata y0, Rdata p){
  
  int i=0;
  double abs_tol=1e-6,rel_tol=1e-5,h=1e-3;
  double *t=NULL;
  assert(model_name);
  
  assert(h5file);
  hid_t h5f_id=H5Fopen(h5file,H5F_ACC_RDONLY,H5P_DEFAULT);
  assert(h5f_id);
  hid_t prior=H5Gopen2(h5f_id,"prior",H5P_DEFAULT);
  gsl_vector *mu=h5_to_gsl(prior,"mu",NULL);
 
  fprintf(stderr,"[%s] prior: ",__func__);
  for (i=0;i<mu->size;i++) fprintf(stderr,"%g ",gsl_vector_get(mu,i));
  fprintf(stderr,"\n");
  
  hsize_t N;
  hid_t data=H5Gopen2(h5f_id,"data",H5P_DEFAULT);
  fprintf(stderr,"[%s] data id: %ld\n",__func__,data);
  simulation_t *sim = sim_from_h5(data,&N);
  size_t d=sim[0].y0->size; // get from initial conditions

  size_t nu=sim[0].u->size;
  gsl_vector *par=gsl_vector_alloc(mu->size + nu);
  gsl_vector_view p=gsl_vector_subvector(par,0,mu->size);
  gsl_vector_view u=gsl_vector_subvector(par,mu->size,nu);
  // initialize:
  gsl_vector_memcpy(&(p.vector),mu);
  gsl_vector_memcpy(&(u.vector),sim[0].u);
  
  for (i=0;i<mu->size;i++) par->data[i]=exp(par->data[i]);

  // load system from file and test it
  jacp dfdp;
  gsl_odeiv2_system sys = load_system(model_name, d, par->data, &dfdp);
  fprintf(stderr,
	  "[%s] ode system dim: %li and %li parameters (%li inputs) with parameters:\n",
	  __func__,
	  sys.dimension,
	  par->size,
	  nu);
  for (i=0;i<mu->size;i++) fprintf(stderr,"%g ",((double*) sys.params)[i]);
  fprintf(stderr,"(that is exp(mu))\n");
  
  test_evaluation(sys,dfdp,sim[0].y0,par);

  const gsl_odeiv2_step_type * T=gsl_odeiv2_step_msbdf;
  gsl_odeiv2_driver *driver=gsl_odeiv2_driver_alloc_y_new(&sys,T,h,abs_tol,rel_tol);
  fprintf(stderr,"[%s] driver allocated with tolearances: abs_tol=%g and rel_tol=%g\n",
	 __func__,abs_tol,rel_tol);
  char *h5out_name = model_function(model_name,"_out.h5");
  fprintf(stderr,"[%s] output file: %s\n",__func__,h5out_name); 
  hid_t h5f = H5Fcreate(h5out_name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  assert(h5f);

  // the most CPU work happens here:
  solution_t **solution=simulate_evolve(sys,driver,sim,N,t,&(u.vector),par,dfdp,h5f);
  gsl_odeiv2_driver_free(driver);
  H5Fclose(h5f);
  return EXIT_SUCCESS;
}
