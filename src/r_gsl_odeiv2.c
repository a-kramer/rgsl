#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <dlfcn.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#define FREE_ON_SUCCESS 1
#define KEEP_ON_SUCCESS 2
#define MATCH 0
#define NO_DIFFERENCE 0
/* SEXP stands for S-Expression, and it can be any R data object (or
 * function) in this program, we'll only use data from R. However SEXP is
 * a bad type name; I am compelled to pronounce it in my head ...
 */
typedef SEXP Rdata;

/* these two function types could be different, but currently aren't */
typedef int(*jacp)(double, const double [], double *, void *);
typedef int(*func)(double, const double [], double *, void *);
typedef int(*def_par)(double, void *);
typedef int (*vf) (double t, const double y[], double dydt[], void * params);
typedef enum {SCALE,DIAG,MATVEC} tf_t;

typedef struct {
	tf_t type;
	int length_b;
	int l;
	double *A;
	double *b;
	gsl_vector *y; /* result vector */
} affine_tf;

typedef struct {
	int nt;
	double *time;
	affine_tf *state;
	affine_tf *par;
} event_t;

/* finds named item in List, `name` can be a space separated list of possible names */
int in_list(Rdata List, const char *name){
	if (!isVector(List)) return -1;
	int i;
	int N=length(List);
	int l=strlen(name);
	char *str=malloc(l+1);
	char *context=NULL;
	*(((char*) memcpy(str,name,l))+l)='\0';
	char *t=strtok_r(str," ",&context);
	while (t){
		for (i=0;i<N;i++){
			if (strcmp(CHAR(STRING_ELT(List,i)),t) == MATCH){
				free(str);
				return i;
			}
		}
		t=strtok_r(NULL," ",&context);
	}
	free(str);
	return -1;
}

Rdata from_list(Rdata List, const char *name){
	if (!isVector(List)) return R_NilValue;
	Rdata names = GET_NAMES(List);
	int i=in_list(names,name);
	Rdata E=R_NilValue;
	if (i>=0){
		E=VECTOR_ELT(List,i);
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
	if (A == R_NilValue || b == R_NilValue) return NULL;

	Rdata dA=GET_DIM(A);
	Rdata db=GET_DIM(b);
	assert(length(dA)==length(db));
	int n=length(dA);
	int *dim=INTEGER(dA);
	int *dim_b=INTEGER(db);
	int j;
	if (dim[0] != dim_b[0]) return NULL;
	affine_tf *L=malloc(sizeof(affine_tf));
	if (n==3) {
		L->l=dim[2];
	} else {
		L->l=1;
	}

	if (dim[0]==dim[1]){
		L->type=MATVEC;
	} else if (dim[1] == 1) {
		L->type=DIAG;
	} else {
		fprintf(stderr,"[%s] A and b have weird dimensionality: A is ",__func__);
		for (j=0;j<n;j++) printf("%i%s",dim[j],j==n-1?"×":" ");
		fprintf(stderr," and b is ");
		for (j=0;j<n;j++) printf("%i%s",dim_b[j],j==n-1?"×":" ");
		fprintf(stderr,"\nThey should be both three dimensional (length(dim(A))==3).\n");
		return NULL;
	}
	L->A=REAL(A);
	L->b=REAL(b);
	L->length_b=dim[0];
	L->y=gsl_vector_alloc(dim[0]);
	return L;
}

void free_tf(affine_tf *L){
	if (L){
		gsl_vector_free(L->y);
		free(L);
	}
}

/* This function applies a affine transformation to the vector z. We
	 assume that A and b and z in the input are all of sufficient size:
	 n×n for matrices and n for vectors. A can be n-sized as well, if A
	 is a diagonal matrix, then we only store the diagonal. n is stored
	 in the transformation structure as length_b. */
int /* the returned status of the gsl operations */
apply_tf(affine_tf *L, /* a transformation struct: A and b are cast to gsl_vectors here */
	 double *z,/* an array of size n, it is updated using L */
	 int t_index)/* if A and b are each a series of matrices, pick the one with this offset */
{
	if (!L) return GSL_SUCCESS; /* nothing to be done */
	if (!z || !(t_index >=0 && L->l > 0)) {
		fprintf(stderr,"[%s] (%p) z must be a pointer to a double array. t_index must be a valid index (0<=%i<%i), L must be applicable to z\n",__func__,z,t_index,L->l);
		return GSL_EINVAL;
	}
	int n=L->length_b;
	int i=t_index % (L->l);
	int status=GSL_SUCCESS;
	gsl_vector_view x=gsl_vector_view_array(z,n);
	gsl_vector_view b=gsl_vector_view_array((L->b)+i*n,n);
	gsl_vector_view a;
	gsl_matrix_view A;
	switch (L->type){
	case SCALE:
		status|=gsl_vector_memcpy(L->y,&(x.vector));
		status|=gsl_vector_scale(L->y,L->A[i]);
		break;
	case DIAG:
		/* y <- diag(A)*x + b */
		a=gsl_vector_view_array((L->A)+i*n,n);
		status|=gsl_vector_memcpy(L->y,&(x.vector));
		status|=gsl_vector_mul(L->y,&(a.vector));
		break;
	case MATVEC:
		/* y <- A*x + b */
		A=gsl_matrix_view_array((L->A)+i*n*n,n,n);
		status|=gsl_blas_dgemv(CblasNoTrans, 1.0, &(A.matrix), &(x.vector), 0.0, L->y);
		break;
	default:
		fprintf(stderr,"[%s] unknown transformation type %i.\n",__func__,L->type);
		return GSL_EINVAL;
	}
	if (status==GSL_SUCCESS){
		gsl_vector_add(L->y,&(b.vector));
		gsl_vector_memcpy(&(x.vector),L->y);
	}
	return status;
}

event_t* event_from_R(Rdata E){
	int j;
	Rdata tf=from_list(E,"tf");
	Rdata time = from_list(E,"time");
	if (tf == R_NilValue || time == R_NilValue) return NULL;
	Rdata state_tf=from_list(tf,"state");
	Rdata param_tf=from_list(tf,"param");
	event_t *event=malloc(sizeof(event_t));
	event->state=affine_transformation(from_list(state_tf,"A"),from_list(state_tf,"b"));
	event->par=affine_transformation(from_list(param_tf,"A"),from_list(param_tf,"b"));

	int lt=length(time);
	event->time = REAL(AS_NUMERIC(time));
	event->nt=lt;
#ifdef DEBUG_PRINT
	printf("[%s] event times:",__func__);
	for (j=0;j<lt;j++) printf("\t%g",event->time[j]);
	putchar('\n');
#endif
	return event;
}

/* this function takes the address of an event structure pointer, clears the
	 memory and changes the pointer to NULL, so that the event cannot be
	 accessed after being freed (except through a different pointer). */
void event_free(event_t **ev){
	if (ev && *ev){
		free_tf((*ev)->state);
		free_tf((*ev)->par);
		free(*ev);
		*ev=NULL;
	}
}

/*This function allocates memory and concatenates two strings in that
	memory. It is used to make function names (this is for loading the model
	functions by name from a shared library). The model function names
	have this pattern: `MODEL_vf`, `MODEL_jac`, `MODEL_jacp`.*/
char* /* string with model_name and suffix (free after loading the function) */
model_function(const char *model_name, /* the base name of the model */
 const char *suffix) /* suffix, usually `"_vf"` or `"_jac"` */
{
	if (!model_name) return NULL;
	size_t size=(strlen(model_name)+strlen(suffix)+1);
	char *f=malloc(sizeof(char)*size);
	strcat(strcpy(f,model_name),suffix);
#ifdef DEBUG_PRINT
	fprintf(stdout,"[%s] «%s»\n",__func__,f);
#endif
	return f;
}


/* Loads a function from an `.so` file, using `dlsym()`.
	 Optionally, this
	 function frees the storage assosiated with the name of the
	 function.*/
void *load_or_warn(void *lib, /* file pointer, previously opened via `dlopen()` */
 char *name, /* function to be loaded from file */
 int opt) /* whether to call `free()` on `name` (either: `KEEP_ON_SUCCESS` or `FREE_ON_SUCCESS`). */
{
  if (!lib || !name) return NULL;
	void *symbol=dlsym(lib,name);
	if (symbol) {
#ifdef DEBUG_PRINT
		fprintf(stderr,"[%s] «%s» loaded successfully.\n",__func__,name);
#endif
		if (opt==FREE_ON_SUCCESS){
			free(name);
		}
	}else{
		fprintf(stderr,"[%s] loading of «%s» failed: slerror was «%s»\n",__func__,name,dlerror());
		return NULL;
	}
	return symbol;
}


/* Loads the ODE system from an `.so` file, the file is given by name,
	 the returned structure is intended for the `gsl_odeiv2` library of
	 solvers. The jacobian dfdx is loaded alongside the right hand side;
	 `gsl_odeiv2_system` has no slot for the parameter derivative `dfdp`,
	 this is returned as a function pointer instead. */
gsl_odeiv2_system /* the system structure, see gsl documentation. */
load_system(
 const char *model_name, /* the name of the model, function names will be inferred from that: model_name_vf, model_name_jac, etc. */
 const char *model_so, /* the path to the shared library that contains the model. */
 size_t n, /* number of state variables, 0 for auto-detection */
 double *p, /* default parameter vector */
 jacp *dfdp, /*[out] additional return value: a pointer to the parameter derivative (matrix) function.*/
 func *F, /* [out] observables for this model (output functions) */
 def_par *DefParFunc) /* [out] function that returns default parameters */
{
	void *lib=dlopen(model_so,RTLD_LAZY);
	void *dfdy;
	gsl_odeiv2_system sys={NULL,NULL,0,NULL};
	vf f;
	//dfdp=malloc(sizeof(jacp*));
	char *symbol_name; // symbol name in .so
	if (lib){
		symbol_name=model_function(model_name,"_vf");
		f=(vf) load_or_warn(lib,symbol_name,FREE_ON_SUCCESS);
		symbol_name=model_function(model_name,"_jac");
		dfdy=load_or_warn(lib,symbol_name,FREE_ON_SUCCESS);
		if (dfdp){
			symbol_name=model_function(model_name,"_jacp");
			*dfdp = (jacp) load_or_warn(lib,symbol_name,FREE_ON_SUCCESS);
		}
		symbol_name=model_function(model_name,"_func");
		if (symbol_name && F){
			*F = (func) load_or_warn(lib,symbol_name,FREE_ON_SUCCESS);
		}
		symbol_name=model_function(model_name,"_default");
		if (symbol_name && DefParFunc){
			*DefParFunc = load_or_warn(lib,symbol_name,FREE_ON_SUCCESS);
		}
	} else {
		fprintf(stderr,"[%s] library «%s» could not be loaded: %s\n",__func__,model_so,dlerror());
		return sys;
	}
	if (!n) n=f(0,NULL,NULL,NULL);
#ifdef DEBUG_PRINT
	fprintf(stderr,"[%s] vf function returns %li components.\n",__func__,n);
#endif
	sys.function=f;
	sys.jacobian=dfdy;
	sys.dimension=n;
	sys.params=p;
#ifdef DEBUG_PRINT
	fprintf(stderr,"[%s] ode system created.\n",__func__); fflush(stderr);
#endif
	return sys;
}

/*This function responds to the status returned by the gsl solvers.*/
void check_status(
		int status, /*the returned value from gsl_odeiv2_driver_apply and similar functions*/
		double current_t, /* the time at which integration stopped*/
		double target_t, /* the time we tried to reach*/
		int iteration)/* the iteration at which the error happened */
{
	int j=iteration;
	double t=current_t;
	double tf=target_t;
	switch (status){
	case GSL_EMAXITER:
		error("[%s] time_point %li: maximum number of steps reached.\n\t\tfinal time: %.10g (short of %.10g)",__func__,j,t,tf);
		break;
	case GSL_ENOPROG:
		error("[%s] time_point %li: step size dropped below set minimum.\n\t\tfinal time: %.10g (short of %.10g)",__func__,j,t,tf);
		break;
	case GSL_EBADFUNC:
		error("[%s] time_point %li: bad function.\n\t\tfinal time: %.10g (short of %.10g)",__func__,j,t,tf);
		break;
	}
}

/* appends the given state to the allocated result buffers */
int append_state(){

}

/* Intergrates the system `sys` using the specified `driver` and
   simulation instructions `sim` (an array of structs, one element per
   simulation). The solver picks the output times */
int /* error code if any, otherwise GSL_SUCCESS */
simulate_timeseries_fine_alloc(const gsl_odeiv2_system sys, /* the system to integrate */
	gsl_odeiv2_driver* driver, /* the driver that is used to integrate `sys` */
	double t0, /* the initial time: y(t0) = y0 */
	const gsl_vector *y0, /* initial value */
	double tf, /* final time of entire simulation */
	const event_t *event, /*a struct array with scheduled events */
	double *Yout, /* return value (a matrix: ny×M, dense in ny), allocated here */
	double *tout) /* return time vector (starts with t0) with length M, allocated here */
{
	gsl_set_error_handler_off();
	int ny=(int) sys.dimension;
	gsl_vector *y=gsl_vector_alloc(ny);
	int i=0,j;
	double t=t0;
	double te;
	int status=GSL_SUCCESS;
	size_t MaxColY=500;
	size_t m=0;
	tout=malloc((sizeof *tout)*MaxColY);
	Yout=malloc((sizeof *Yout)*ny*MaxColY);
	/* initialize t0 values */
	gsl_vector_memcpy(y,y0);

	for (j = 0, i = 0; j < nt; j++){
		// reallocate when we run out of columns to put results in
		if (m==MaxColY) {
			MaxColY += 100;
			tout = realloc(tout,(sizeof *tout) * MaxColY);
			Yout = realloc(Yout,(sizeof *Yout) * ny * MaxColY);
		}
		tout[m]=t;
		
		if (event && i<event->nt && event->time[i] < tf) {
			te=event->time[i];
			while (t<te){
				
				status=gsl_odeiv2_evolve_apply(driver->e, driver->c, driver->s, &sys, &t, te, &(driver->h), y->data);
				if (status!=GSL_SUCCESS){
					// output is initial conditions, with error code from solver
					tout = realloc(tout,(sizeof *tout));
					Yout = realloc(Yout,(sizeof *Yout)*ny);
					return(status);
				} else {
					memcpy(Yout[m*(y->size)],y->data,(sizeof *Yout)*(y->size));
				}
			}
			apply_tf(event->state,y->data,i);
			apply_tf(event->par,(double*) sys.params,i);
			status=gsl_odeiv2_driver_reset(driver);
			if (status!=GSL_SUCCESS){
				return(status);
			}
			i++;
		}
		if (tf>t) status=gsl_odeiv2_evolve_apply(driver->e, driver->c, driver->step, &sys, &t, te, &h, y->data);
		//report any error codes to the R user
		check_status(status,t,tf,j);
		if(status==GSL_SUCCESS){
			Yout_row = gsl_matrix_row(Yout,j);
			gsl_vector_memcpy(&(Yout_row.vector),y);
		} else {
			return(status);
		}
	}
	gsl_odeiv2_driver_reset(driver);
	return status;
}

/* Intergrates the system `sys` using the specified `driver` and
   simulation instructions `sim` (an array of structs, one element per
   simulation). The results are saved to an hdf5 file and also printed
   to standard output. */
int /* error code if any, otherwise GSL_SUCCESS */
simulate_timeseries(const gsl_odeiv2_system sys, /* the system to integrate */
	gsl_odeiv2_driver* driver, /* the driver that is used to integrate `sys` */
	double t0, /* the initial time: y(t0) = y0 */
	const gsl_vector *y0, /* initial value */
	const gsl_vector *time, /* a vector of time-points */
	const event_t *event, /*a struct array with scheduled events */
	gsl_matrix *Yout) /* (OUT) return vaule, pre-allocated */
{
	gsl_set_error_handler_off();
	int nt=time->size;
	int ny=(int) sys.dimension;
	gsl_vector *y=gsl_vector_alloc(ny);
	gsl_vector_memcpy(y,y0);

	gsl_vector_view Yout_row;
	int i=0,j;
	double t=t0;
	double tf,te;
	int status=GSL_SUCCESS;

	/* initialize t0 values */
	gsl_vector_memcpy(y,y0);

	for (j=0, i=0; j<nt; j++){
		tf=gsl_vector_get(time,j);
		if (event && i<event->nt && event->time[i] < tf) {
			te=event->time[i];
			status=gsl_odeiv2_driver_apply(driver, &t, te, y->data);
			if (status!=GSL_SUCCESS){
#ifdef DEBUG_PRINT
				fprintf(stderr,"[%s] before event %i gsl_odeiv2_driver_apply produced an error: %s.\n",__func__,i,gsl_strerror(status));
#endif
				return(status);
			}
			apply_tf(event->state,y->data,i);
			apply_tf(event->par,(double*) sys.params,i);
			status=gsl_odeiv2_driver_reset(driver);
			if (status!=GSL_SUCCESS){
#ifdef DEBUG_PRINT
				fprintf(stderr,"[%s] resetting the system after event %i produced an error: %s.\n",__func__,i,gsl_strerror(status));
#endif
				return(status);
			}
			i++;
		}
		if (tf>t) status=gsl_odeiv2_driver_apply(driver, &t, tf, y->data);
		//report any error codes to the R user
		check_status(status,t,tf,j);
		if(status==GSL_SUCCESS){
			Yout_row = gsl_matrix_row(Yout,j);
			gsl_vector_memcpy(&(Yout_row.vector),y);
		} else {
			return(status);
		}
	}
	gsl_odeiv2_driver_reset(driver);
	return status;
}

/* This prgram loads an ODE model, as needed for `gsl_odeiv2`.
	 It simulates the model for each column of initial
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
 Rdata modelName, /* a string */
 Rdata tspan, /* a vector of output times with tspan[0] crresponding to y0 */
 Rdata t0, /* initial time y(t0) = y0 */
 Rdata y0, /* initial conditions at tspan[0], can be a matrix*/
 Rdata p, /* parameter matrix */
 Rdata event, /* list of events */
 Rdata absolute_tolerance, /* absolute tolerance for GSL's solver */
 Rdata relative_tolerance, /* relative tolerance for GSL's solver */
 Rdata initial_step_size) /* initial guess for the step size */
{
	gsl_set_error_handler_off();
	int i,j;
	double abs_tol=asReal(absolute_tolerance);
	double rel_tol=asReal(relative_tolerance);
	double h=asReal(initial_step_size);
	size_t nt=length(tspan);
	size_t ny,np,N;
	Rdata experiment_names;
	int ExperimentsAreNamed = 0;
	const char* model_so=CHAR(asChar(getAttrib(modelName,install("comment"))));
	const char* model_name=CHAR(STRING_ELT(modelName,0));
	if (isMatrix(p)){
		np=nrows(p);
		N=ncols(p);
		ny=nrows(y0);
		experiment_names = GET_COLNAMES(p);
		if (experiment_names != R_NilValue) ExperimentsAreNamed = 1;
	} else if (isVector(p)) {
		np=length(p);
		N=1;
		ny=length(y0);
	}

	gsl_vector_view t=gsl_vector_view_array(REAL(AS_NUMERIC(tspan)),nt);
	gsl_matrix_view initial_value = gsl_matrix_view_array(REAL(AS_NUMERIC(y0)),N,ny);
	gsl_matrix_view ode_parameter = gsl_matrix_view_array(REAL(AS_NUMERIC(p)),N,np);

	jacp dfdp;
	gsl_odeiv2_system sys = load_system(model_name, model_so, ny, REAL(AS_NUMERIC(p)), &dfdp, NULL,NULL);
	if (sys.dimension == 0) {
		fprintf(stderr,"system dimension is «%li», nothing left to do.\n",sys.dimension);
		return R_NilValue;
	}

	int l;
	Rdata event_names;
	event_t **ev=calloc(N,sizeof(event_t*));
	if (event && event != R_NilValue){
		l=length(event);
		if (ExperimentsAreNamed){
			event_names=GET_NAMES(event);
			for (i=0;i<l;i++){
				j=in_list(experiment_names,CHAR(VECTOR_ELT(event_names,i)));
				if (j>=0) ev[j] = event_from_R(VECTOR_ELT(event, i));
			}
		} else {
			for (i=0;i<l;i++){
				ev[i] = event_from_R(VECTOR_ELT(event, i));
			}
		}
	}
	const gsl_odeiv2_step_type * T=gsl_odeiv2_step_msbdf;
	gsl_odeiv2_driver *driver;
	Rdata Y = PROTECT(alloc3DArray(REALSXP,ny,nt,N));
	for (i=0;i<ny*nt*N;i++) REAL(Y)[i]=NA_REAL;

	double *ydata;
	gsl_vector_view iv_row;
	size_t nyt=nt*ny;
	gsl_matrix_view y;
	int status;
	ydata = REAL(AS_NUMERIC(Y));
	for (i=0;i<N;i++){
		driver=gsl_odeiv2_driver_alloc_y_new(&sys,T,h,abs_tol,rel_tol);
#ifdef DEBUG_PRINT
		printf("[%s] solving %i of %li.\n",__func__,i,N);
#endif
		y=gsl_matrix_view_array(&(ydata[i*nyt]),nt,ny);
		sys.params = gsl_matrix_ptr(&(ode_parameter.matrix),i,0);
		iv_row=gsl_matrix_row(&(initial_value.matrix),i);
		status=simulate_timeseries
			(sys,
			 driver,
			 REAL(AS_NUMERIC(t0))[0],
			 &(iv_row.vector),
			 &(t.vector),
			 ev[i],
			 &(y.matrix));
		if (status!=GSL_SUCCESS){
			fprintf(stderr,"[%s] simulation error: %s\n",__func__,gsl_strerror(status));
			break;
		}
		gsl_odeiv2_driver_free(driver);
		event_free(&(ev[i]));
	}
	free(ev);
	UNPROTECT(1);
	return Y;
}

void set_names(Rdata list, const char *names[], size_t n)
{
	int i;
	Rdata rnames=PROTECT(allocVector(STRSXP, n));
	Rdata elt;
	for (i=0;i<n;i++){
		elt=PROTECT(mkChar(names[i]));
		SET_STRING_ELT(rnames,i,elt);
	}
	SET_NAMES(list,rnames);
	UNPROTECT(1+n);
}

/* This prgram loads an ODE model, specified for `gsl_odeiv2`. It
	 simulates the model for each entry in a list of named items, each
	 describing a single initial value problem (y0, t, parameters p, events).
	 ```
	 y'=f(y,t;p) y0=y(t[0])
	 ```
*/
Rdata /* the trajectories as a list (same size as experiments) */
r_gsl_odeiv2_simulate(
 Rdata modelName, /* a string */
 Rdata experiments, /* a list of simulation experiments */
 Rdata absolute_tolerance, /* absolute tolerance for GSL's solver */
 Rdata relative_tolerance, /* relative tolerance for GSL's solver */
 Rdata initial_step_size) /* initial guess for the step size */
{
	gsl_set_error_handler_off();
	const char* model_so=CHAR(asChar(getAttrib(modelName,install("comment"))));
	const char* model_name=CHAR(STRING_ELT(modelName,0));
	int i,j,status;
	double *f;
	double abs_tol=asReal(absolute_tolerance);
	double rel_tol=asReal(relative_tolerance);
	double h=asReal(initial_step_size);
	int N=GET_LENGTH(experiments);
	const gsl_odeiv2_step_type * T=gsl_odeiv2_step_msbdf;
	gsl_odeiv2_driver *driver;
	Rdata res_list = PROTECT(NEW_LIST(N)); /* use VECTOR_ELT and SET_VECTOR_ELT */
	SET_NAMES(res_list,GET_NAMES(experiments));
	const char *yf_names[2]={"state","func"};
	Rdata iv, t, t0, F, Y, yf_list;
	gsl_vector_view initial_value, time;
	gsl_matrix_view y;
	size_t ny, nt, nf;
	event_t *ev=NULL;
	jacp dfdp=NULL;
	func observable=NULL;
	gsl_odeiv2_system sys = load_system(model_name, model_so, 0, NULL, &dfdp, &observable, NULL);
	if (sys.dimension == 0) {
		UNPROTECT(1); /* res_list */
		fprintf(stderr,"[%s] system dimension is «%li».\n",__func__,sys.dimension);
		return R_NilValue;
	}
	for (i=0;i<N;i++){
		driver=gsl_odeiv2_driver_alloc_y_new(&sys,T,h,abs_tol,rel_tol);
		iv = from_list(VECTOR_ELT(experiments,i),"initial_value initialState");
		t = from_list(VECTOR_ELT(experiments,i),"time outputTimes");
		ev = event_from_R(from_list(VECTOR_ELT(experiments,i),"events scheduledEvents"));
		t0 = from_list(VECTOR_ELT(experiments,i),"initialTime t0 T0");
		ny=length(iv);
		nt=length(t);
		if (ny>0 && nt>0 && ny==sys.dimension){
			initial_value=gsl_vector_view_array(REAL(AS_NUMERIC(iv)),ny);
			time=gsl_vector_view_array(REAL(AS_NUMERIC(t)),nt);
			sys.params = REAL(from_list(VECTOR_ELT(experiments,i),"parameters param par"));
			assert(sys.params);

			Y=PROTECT(allocMatrix(REALSXP,ny,nt));
			for (j=0;j<ny*nt;j++) REAL(Y)[j]=NA_REAL;
			y=gsl_matrix_view_array(REAL(AS_NUMERIC(Y)),nt,ny);

			if (observable){
				nf=observable(0,NULL,NULL,NULL);
				F=PROTECT(allocMatrix(REALSXP,nf,nt));
			}
			status=simulate_timeseries(sys,
				driver,
				REAL(AS_NUMERIC(t0))[0],
				&(initial_value.vector),
				&(time.vector),ev,&(y.matrix)
			);
			if (observable && status==GSL_SUCCESS) {
				for (j=0;j<nt;j++){
					f=&(REAL(AS_NUMERIC(F))[j*nf]);
					observable(gsl_vector_get(&(time.vector),j),gsl_matrix_ptr(&(y.matrix),j,0),f,sys.params);
				}
			}
			yf_list=PROTECT(NEW_LIST(2));
			SET_VECTOR_ELT(yf_list,0,Y);
			SET_VECTOR_ELT(yf_list,1,F);
			set_names(yf_list,yf_names,2);
			SET_VECTOR_ELT(res_list,i,yf_list);
			gsl_odeiv2_driver_free(driver);
			event_free(&ev);
			UNPROTECT(1); /* yf_list */
			if (observable) UNPROTECT(1); /* F */
			UNPROTECT(1); /* Y */
			if (status!=GSL_SUCCESS) break;
		}
	}
	UNPROTECT(1); /* res_list */
	return res_list;
}


/* This prgram loads an ODE model, with functions compatible with `gsl_odeiv2`
	(see odeiv2 documentation on the GSL webpage). It
	 simulates the model for each entry in a list of named items, each
	 describing a single initial value problem (`y0`, `t`, parameters `p`, events).
	 ```
	 y'=f(y,t;p) y0=y(t[0])
	 ```
	 This function is parallel in the list of supplied experiments.
*/
Rdata /* the trajectories as a list (same size as experiments) */
r_gsl_odeiv2_outer(
 Rdata modelName, /* a string */
 Rdata experiments, /* a list of simulation experiments */
 Rdata parameters, /* a matrix of parameterization columns*/
 Rdata absolute_tolerance, /* absolute tolerance for GSL's solver */
 Rdata relative_tolerance, /* relative tolerance for GSL's solver */
 Rdata initial_step_size) /* initial guess for the step size */
{
	gsl_set_error_handler_off();
	const char* model_so=CHAR(asChar(getAttrib(modelName,install("comment"))));
	const char* model_name=CHAR(STRING_ELT(modelName,0));
	int i,j,k,l,status=GSL_SUCCESS;
	double abs_tol=asReal(absolute_tolerance);
	double rel_tol=asReal(relative_tolerance);
	double h=asReal(initial_step_size);
	int N=GET_LENGTH(experiments);
	size_t np=nrows(parameters);
	size_t M=ncols(parameters);
	const gsl_odeiv2_step_type * T=gsl_odeiv2_step_msbdf;
	gsl_odeiv2_driver *driver;
	Rdata res_list = PROTECT(NEW_LIST(N)); /* use VECTOR_ELT and SET_VECTOR_ELT */
	SET_NAMES(res_list,GET_NAMES(experiments));
	Rdata yf_list, input, Y, F, iv, t, t0;
	const char *yf_names[2]={"state","func"};
	gsl_vector_view initial_value, time;
	gsl_matrix_view y;
	size_t ny, nt, nf, nu;
	event_t *ev=NULL;
	double *p, *f;
	jacp dfdp=NULL;
	func observable=NULL;
	def_par default_parameters=NULL;
	gsl_odeiv2_system sys = load_system(model_name, model_so, 0, NULL, &dfdp, &observable, &default_parameters);
	int np_model = default_parameters(0,NULL);
	if (sys.dimension == 0) {
		UNPROTECT(1); /* res_list */
		fprintf(stderr,"[%s] system dimension: «%li».\n",__func__,sys.dimension);
		return R_NilValue;
	}
	for (i=0; i<N; i++){
		driver=gsl_odeiv2_driver_alloc_y_new(&sys,T,h,abs_tol,rel_tol);
		iv = from_list(VECTOR_ELT(experiments,i),"initial_value initialState");
		t = from_list(VECTOR_ELT(experiments,i),"time outputTimes");
		t0 = from_list(VECTOR_ELT(experiments,i),"initialTime t0 T0");
		ev = event_from_R(from_list(VECTOR_ELT(experiments,i),"events scheduledEvents"));
		input = from_list(VECTOR_ELT(experiments,i),"input");
		nu=(input && input!=R_NilValue)?length(input):0;
		p=malloc(sizeof(double)*np_model);
		if (np_model!= np+nu) {
			fprintf(stderr,"[%s] number of parameters given (np=%li + nu=%li) does not agree with the supplied model (%i).\n",__func__,np,nu,np_model);
		}
		if (p && nu) memcpy(p+np,REAL(AS_NUMERIC(input)),nu*sizeof(double));
		ny=length(iv);
		nt=length(t);
		if (ny>0 && nt>0 && ny==sys.dimension){
			initial_value=gsl_vector_view_array(REAL(AS_NUMERIC(iv)),ny);
			time=gsl_vector_view_array(REAL(AS_NUMERIC(t)),nt);
			Y=PROTECT(alloc3DArray(REALSXP,ny,nt,M));
			for (j=0; j<ny*nt*M; j++) REAL(Y)[j]=NA_REAL; /* initialize to NA */
			if (observable){
				nf=observable(0,NULL,NULL,NULL);
				F=PROTECT(alloc3DArray(REALSXP,nf,nt,M));
				for (j=0;j<nf*nt*M;j++) REAL(F)[j]=NA_REAL; /* initialize to NA */
			}
			for (k=0;k<M;k++){
				y=gsl_matrix_view_array(REAL(AS_NUMERIC(Y))+(nt*ny*k),nt,ny);
				memcpy(p,REAL(AS_NUMERIC(parameters))+np*k,np*sizeof(double));
				sys.params = p;
				status=simulate_timeseries(
					sys,
					driver,
					REAL(AS_NUMERIC(t0))[0],
					&(initial_value.vector),
					&(time.vector),
					ev,
					&(y.matrix)
				);
				if (observable && status==GSL_SUCCESS) {
					for (j=0;j<nt;j++){
						f=REAL(F)+(0+j*nf+k*nf*nt);
						observable(gsl_vector_get(&(time.vector),j),gsl_matrix_ptr(&(y.matrix),j,0),f,sys.params);
					}
				}
			}
		}
		yf_list=PROTECT(NEW_LIST(2));
		SET_VECTOR_ELT(yf_list,0,Y);
		SET_VECTOR_ELT(yf_list,1,F);
		set_names(yf_list,yf_names,2);
		SET_VECTOR_ELT(res_list,i,yf_list);
		event_free(&ev);
		gsl_odeiv2_driver_free(driver);
		UNPROTECT(1); /* yf_list */
		if (observable) UNPROTECT(1); /* F */
		UNPROTECT(1); /* Y */
#ifdef DEBUG_PRINT
		if (status!=GSL_SUCCESS){
			fprintf(stderr,"[%s] parameter set lead to solver errors (%s) in experiment %i/%i, values:\n",__func__,gsl_strerror(status),i,N);
			for (j=0; j<np_model; j++) fprintf(stderr,"%g%s",p[j],(j==np_model-1?"\n":", "));
		}
#endif
		free(p);
	}
	UNPROTECT(1); /* res_list */
	return res_list;
}
