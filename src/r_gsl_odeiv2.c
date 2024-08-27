#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#include <math.h>
#include <dlfcn.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#define MATCH 0
#define NO_DIFFERENCE 0
#define RCOND_LIMIT 1e-10
/* SEXP stands for S-Expression, and it can be any R data object (or
 * function) in this program, we'll only use data from R. However SEXP is
 * a bad type name; I am compelled to pronounce it in my head ...
 */
typedef SEXP Rdata;

/* an ODE model y'=f(t,y,p) will have a name, and several functions, prefixed with
 * that name, e.g.: ${name}_vf().
 * The possible suffixes are:
 * _vf          the ODE right-hand-side function f (vector-field)
 * _jac         the y-jacobian df/dy (with respect to y)
 * _jacp        the p-jacobian df/dp (with respect to p)
 * _init        sets the initial values of the state variables y
 * _default     sets the default parameter values p
 * _func        output functions F (observable values)
 * _funcJac     output function y-jacobian dF/dy
 * _funcJacp    output function p-jacobian dF/dp
 * all of these functions return integer error codes and GSL_SUCCESS on success (0).
 *
 * We load these functions by name from a shared library (.so) and
 * give them new, generic names, where we replace the model's name with "ODE".
 *
 */
int (*ODE_vf)(double t, const double y_[], double f_[], void *par);
int (*ODE_jac)(double t, const double y_[], double *jac_, double *dfdt_, void *par);
int (*ODE_jacp)(double t, const double y_[], double *jacp_, double *dfdt_, void *par);
int (*ODE_func)(double t, const double y_[], double *func_, void *par);
int (*ODE_funcJac)(double t, const double y_[], double *funcJac_, void *par);
int (*ODE_funcJacp)(double t, const double y_[], double *funcJacp_, void *par);
int (*ODE_default)(double t, void *par);
int (*ODE_init)(double t, double *y_, void *par);
int (*ODE_event)(double t, double *y_, void *par, int EventLabel, double dose);

void transition_matrix(gsl_matrix *Ji, gsl_matrix *Jf, double ti, double tf, gsl_matrix *phi);

typedef enum {SCALE,DIAG,MATVEC} tf_t;

typedef struct {
	tf_t type;
	int length_b;
	int l;
	double *A;
	double *b;
	gsl_vector *y; /* result vector */
} affine_tf;

enum event_type {affine_tf_event, model_func_event};
struct event {
	enum event_type type;
	int *label;
	int nt;
	int nL;
	int nDose;
	double *time;
	double *dose;
	affine_tf *state;
	affine_tf *par;
};

double sec(clock_t c){
	double t=(double) c;
	return t/CLOCKS_PER_SEC;
}

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
		fprintf(stderr,"[%s] (%p) z must be a pointer to a double array. t_index must be a valid index (0<=%i<%i), L must be applicable to z\n",__func__,(void*) z,t_index,L->l);
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

struct event* event_from_R(Rdata E){
	int j;
	if (E == R_NilValue) return NULL;
	if (!IS_VECTOR(E)) return NULL;
	Rdata tf=from_list(E,"tf");
	Rdata time = from_list(E,"time");
	Rdata label = from_list(E,"label");
	Rdata dose = from_list(E,"dose");
	if (time == R_NilValue) return NULL;
	struct event *event=malloc(sizeof(struct event));
	if (label != R_NilValue){
		event->nL = length(label);
		event->label = INTEGER(AS_INTEGER(label));
		event->dose = REAL(AS_NUMERIC(dose));
		event->nDose = length(dose);
		event->type = model_func_event;
		//fprintf(stderr,"[%s] event of type «model_func_event» %i loaded.\n",__func__,event->type);
	} else if (tf != R_NilValue){
		Rdata state_tf=from_list(tf,"state");
		Rdata param_tf=from_list(tf,"param");
		event->state=affine_transformation(from_list(state_tf,"A"),from_list(state_tf,"b"));
		event->par=affine_transformation(from_list(param_tf,"A"),from_list(param_tf,"b"));
		event->type=affine_tf_event;
		//fprintf(stderr,"[%s] event of type «affine_tf_event» %i loaded.\n",__func__,event->type);
	} else {
		fprintf(stderr,"[%s] unknown type of event.\n",__func__);
	}
	int lt=length(time);
	if (event->type == model_func_event && lt != event->nL){
		fprintf(stderr,"[%s] event time vector and event label vector must be the same length for this type of event. (%i != %i)\n",__func__,lt,event->nL);
	}
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
void event_free(struct event **ev){
	if (ev && *ev){
		switch ((*ev)->type){
			case affine_tf_event:
			free_tf((*ev)->state);
			free_tf((*ev)->par);
			break;
			case model_func_event:
			/* for now, no action is needed */
			break;
		}
		free(*ev);
		*ev=NULL;
	}
}

/* Loads a function from an `.so` file, using `dlsym()`.
	 Optionally, this
	 function frees the storage assosiated with the name of the
	 function.*/
void *load_or_warn(void *lib, /* file pointer, previously opened via `dlopen()` */
 char *name) /* function to be loaded from file */
{
  if (!lib || !name) return NULL;
	void *symbol=dlsym(lib,name);
	if (symbol) {
#ifdef DEBUG_PRINT
		fprintf(stderr,"[%s] «%s» loaded successfully.\n",__func__,name);
#endif
	} else {
#ifdef DEBUG_PRINT
		fprintf(stderr,"[%s] loading of «%s» failed: slerror was «%s»\n",__func__,name,dlerror());
#endif
		return NULL;
	}
	return symbol;
}


/* Loads the ODE system from an `.so` file, the file is given by name,
	 the returned structure is intended for the `gsl_odeiv2` library of
	 solvers. The jacobian dfdx is loaded alongside the right hand side;
	 `gsl_odeiv2_system`. Other model functions (_jacp, _func, _init,
	 _default, funcJac, funcJacp) are assigned to global variables. */
gsl_odeiv2_system /* the system structure, see gsl documentation. */
load_system(
 const char *model_name, /* the name of the model, function names will be inferred from that: model_name_vf, model_name_jac, etc. */
 const char *model_so) /* the path to the shared library that contains the model. */
{
	void *lib=dlopen(model_so,RTLD_LAZY);
	gsl_odeiv2_system sys={NULL,NULL,0,NULL};
	size_t n,l;
	size_t m=strlen(model_name);
	char *symbol_name=malloc(m+32); // symbol name in .so
	char *suffix=mempcpy(symbol_name,model_name,m);
	*suffix='\0';

	if (lib){
		*((char*) mempcpy(suffix,"_vf",3))='\0';
		//printf("[%s] loading «%s» from «%s».\n",__func__,symbol_name,model_so);
		if ((ODE_vf=load_or_warn(lib,symbol_name))==NULL){
			fprintf(stderr,"[%s] loading «%s» is required.\n",__func__,symbol_name);
			free(symbol_name);
			return sys;
		}
		*((char*) mempcpy(suffix,"_jac",4))='\0';
		//printf("[%s] loading «%s» from «%s».\n",__func__,symbol_name,model_so);
		if ((ODE_jac=load_or_warn(lib,symbol_name))==NULL){
			fprintf(stderr,"[%s] loading «%s» is required.\n",__func__,symbol_name);
			free(symbol_name);
			return sys;
		}
		*((char*) mempcpy(suffix,"_jacp",5))='\0';
		//printf("[%s] loading «%s» from «%s».\n",__func__,symbol_name,model_so);
		if ((ODE_jacp = load_or_warn(lib,symbol_name))==NULL){
#ifdef DEBUG_PRINT
			fprintf(stderr,"[%s] not having «%s» is OK if no sensitivities are needed.\n",__func__,symbol_name);
#endif
		}
		*((char*) mempcpy(suffix,"_func",5))='\0';
		//printf("[%s] loading «%s» from «%s».\n",__func__,symbol_name,model_so);
		if ((ODE_func = load_or_warn(lib,symbol_name))==NULL){
#ifdef DEBUG_PRINT
			fprintf(stderr,"[%s] not having «%s» is OK in some cases, output values will not be calculated.\n",__func__,symbol_name);
#endif
		}

		*((char*) mempcpy(suffix,"_funcJac",8))='\0';
		//printf("[%s] loading «%s» from «%s».\n",__func__,symbol_name,model_so);
		if ((ODE_funcJac = load_or_warn(lib,symbol_name))==NULL){
#ifdef DEBUG_PRINT
			fprintf(stderr,"[%s] not having «%s» is OK, but output sensitivities will not be calculated.\n",__func__,symbol_name);
#endif
		}

		*((char*) mempcpy(suffix,"_funcJacp",9))='\0';
		//printf("[%s] loading «%s» from «%s».\n",__func__,symbol_name,model_so);
		if ((ODE_funcJacp = load_or_warn(lib,symbol_name))==NULL){
#ifdef DEBUG_PRINT
			fprintf(stderr,"[%s] not having «%s» is OK, but output sensitivities will not be calculated.\n",__func__,symbol_name);
#endif
		}
		*((char*) mempcpy(suffix,"_default",8))='\0';
		//printf("[%s] loading «%s» from «%s».\n",__func__,symbol_name,model_so);
		if ((ODE_default = load_or_warn(lib,symbol_name))==NULL) {
			fprintf(stderr,"[%s] loading «%s» is required.\n",__func__,symbol_name);
			free(symbol_name);
			return sys;
		}

		*((char*) mempcpy(suffix,"_init",5))='\0';
		//printf("[%s] loading «%s» from «%s».\n",__func__,symbol_name,model_so);
		if ((ODE_init = load_or_warn(lib,symbol_name))==NULL){
#ifdef DEBUG_PRINT
			fprintf(stderr,"[%s] «%s» is optional.\n",__func__,symbol_name);
#endif
		}
		*((char*) mempcpy(suffix,"_event",6))='\0';
		//printf("[%s] loading «%s» from «%s».\n",__func__,symbol_name,model_so);
		if ((ODE_event = load_or_warn(lib,symbol_name))==NULL){
#ifdef DEBUG_PRINT
			fprintf(stderr,"[%s] «%s» is optional.\n",__func__,symbol_name);
#endif
		}
	} else {
		fprintf(stderr,"[%s] library «%s» could not be loaded: %s\n",__func__,model_so,dlerror());
		return sys;
	}
	n=ODE_vf(0,NULL,NULL,NULL);
	l=ODE_default(0,NULL);
#ifdef DEBUG_PRINT
	fprintf(stderr,"[%s] vf function returns %li components.\n",__func__,n);
#endif
	sys.function=ODE_vf;
	sys.jacobian=ODE_jac;
	sys.dimension=n;
	double *p=malloc(sizeof(double)*l);
	ODE_default(0,p);
	sys.params=p;
#ifdef DEBUG_PRINT
	fprintf(stderr,"[%s] ode system created.\n",__func__); fflush(stderr);
#endif
	free(symbol_name);
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

/* Intergrates the system `sys` using the specified `driver` and
   simulation instructions. */
int /* error code if any, otherwise GSL_SUCCESS */
simulate_timeseries(const gsl_odeiv2_system sys, /* the system to integrate */
	gsl_odeiv2_driver* driver, /* the driver that is used to integrate `sys` */
	double t0, /* the initial time: y(t0) = y0 */
	const gsl_vector *y0, /* initial value */
	const gsl_vector *time, /* a vector of time-points */
	const struct event *event, /*a struct array with scheduled events */
	gsl_matrix *Yout) /* (OUT) return vaule, pre-allocated */
{
	gsl_set_error_handler_off();
	int nt=time->size;
	int ny=(int) sys.dimension;
	gsl_vector *y=gsl_vector_alloc(ny);
	gsl_vector_memcpy(y,y0);

	gsl_vector_view Yout_row;
	int i=0,j,nL,nDose;
	double t=t0;
	double tf,te;
	int status=GSL_SUCCESS;

	/* initialize t0 values */
	gsl_vector_memcpy(y,y0);

	for (j=0, i=0; j<nt; j++){
		tf=gsl_vector_get(time,j);
		while (event && i<event->nt && event->time[i] < tf) {
			te=event->time[i];
			status=gsl_odeiv2_driver_apply(driver, &t, te, y->data);
			if (status!=GSL_SUCCESS){
#ifdef DEBUG_PRINT
				fprintf(stderr,"[%s] before event %i gsl_odeiv2_driver_apply produced an error: %s.\n",__func__,i,gsl_strerror(status));
#endif
				return(status);
			}
			//fprintf(stderr,"[%s] event type: %i.\n",__func__,event->type);
			switch (event->type){
			case affine_tf_event:
				apply_tf(event->state,y->data,i);
				apply_tf(event->par,(double*) sys.params,i);
				break;
			case model_func_event:
				nL = event->nL;
				nDose = event->nDose;
				ODE_event(te,y->data,sys.params,event->label[i % nL],event->dose[i % nDose]);
				break;
			}
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
	gsl_vector_free(y);
	gsl_odeiv2_driver_reset(driver);
	return status;
}

/*  QR witch column pivoting is abbreviated as QRPT */
struct pbsrMem {
	double t;
	gsl_matrix *A;     // jacobian: df/dy
	gsl_matrix *B;     // parameter-jacobian: df/dp
	gsl_matrix *Sy;    // sensitivity of y: dy/dp
	gsl_matrix *Sf;    // sensitivity of output functions: dF/dp
	gsl_matrix *phi;   // transfer matrix
	gsl_matrix *temp;  // temporary storage sized like B and Sy
	gsl_vector *dfdt;
	gsl_matrix *FA;
	gsl_matrix *FB;
};

struct pbsrMem pbsr_mem_alloc(int ny, int np, int nf){
	struct pbsrMem m;
	m.temp = gsl_matrix_alloc(ny,np);
	m.phi = gsl_matrix_alloc(ny,ny);
	m.A = gsl_matrix_alloc(ny,ny);
	m.B = gsl_matrix_alloc(ny,np);
	m.Sy = gsl_matrix_alloc(ny,np);
	m.Sf = gsl_matrix_alloc(nf,np);
	m.dfdt = gsl_vector_alloc(ny);
	m.FA = gsl_matrix_alloc(nf,ny);
	m.FB = gsl_matrix_alloc(nf,np);
	return m;
}

void pbsr_mem_free(struct pbsrMem m){
	gsl_matrix_free(m.temp);
	gsl_matrix_free(m.phi);
	gsl_matrix_free(m.A);
	gsl_matrix_free(m.B);
	gsl_matrix_free(m.Sy);
	gsl_matrix_free(m.Sf);
	gsl_vector_free(m.dfdt);
	gsl_matrix_free(m.FA);
	gsl_matrix_free(m.FB);
}

struct QRPT {
	gsl_matrix *QR;    // stores the QR decomposed factors for QR = AP
	gsl_vector *tau;   // QRPT tau vector
	gsl_vector *work;  // QRPT work vector for rcond
	gsl_vector *norm;  // workspace
	gsl_permutation *p;// permutation for WRPT
	int signum;        // QRPT sign
};

struct QRPT qrtp_alloc(int ny){
	struct QRPT m;
	m.QR = gsl_matrix_alloc(ny,ny);
	m.tau = gsl_vector_alloc(ny);
	m.work = gsl_vector_alloc(3*ny);
	m.p = gsl_permutation_alloc(ny);
	m.norm = gsl_vector_alloc(ny);
	return m;
}

void qrtp_free(struct QRPT m){
	gsl_matrix_free(m.QR);
	gsl_vector_free(m.tau);
	gsl_vector_free(m.work);
	gsl_permutation_free(m.p);
	gsl_vector_free(m.norm);
}

int qrtp(struct QRPT qr, const gsl_matrix *A){
	int status;
	gsl_matrix_memcpy(qr.QR,A);
	status = gsl_linalg_QRPT_decomp(qr.QR, qr.tau, qr.p, &(qr.signum), qr.norm);
	return status;
}

double rcond(struct QRPT qr){
	double r;
	int status;
	status = gsl_linalg_QRPT_rcond(qr.QR, &r, qr.work);
	if (status != GSL_SUCCESS) {
		fprintf(stderr,"[%s] rcond (%g) failed: %s\n",__func__,r,gsl_strerror(status));
	}
	return r;
}

struct eigen {
	gsl_eigen_nonsymm_workspace *w;
	gsl_vector_complex *ev;
	gsl_matrix *T;
};

struct eigen eigen_alloc(const int ny){
	struct eigen E;
	E.w = gsl_eigen_nonsymm_alloc(ny);
	E.T = gsl_matrix_alloc(ny,ny);
	E.ev= gsl_vector_complex_alloc(ny);
	gsl_eigen_nonsymm_params(0, 1, E.w);
	return E;
}

void eigen_free(struct eigen E){
	gsl_eigen_nonsymm_free(E.w);
	gsl_matrix_free(E.T);
	gsl_vector_complex_free(E.ev);
}

int stable(const gsl_matrix *A, struct eigen E){
	gsl_matrix_memcpy(E.T,A);
	gsl_eigen_nonsymm(E.T, E.ev, E.w);
	gsl_vector_view re_ev = gsl_vector_complex_real(E.ev);
	return gsl_vector_isneg(&re_ev.vector);
}

/* advances the state y from `left` state to `right` state and also
	 advances the sensitivity approximation with it. the structs letf
	 and right have pre-allocated space for jacobians and other matrices
	 (including the snesitivity matrix). The structs pbsrMem and QRPT
	 are structs of pointers to allocated memory.*/
int advance(const gsl_odeiv2_system sys, gsl_odeiv2_driver *driver, gsl_vector *y, struct pbsrMem left, struct pbsrMem right, struct QRPT qr, struct eigen ev){
	int status=GSL_SUCCESS;
	const double ti=left.t;  // initial time point
	const double tf=right.t; // final time point
	double t=left.t;         // t is advanced from `ti` to `tf`
	double delta_t;
	double *p=sys.params;
	double h=1e-2;
	int k,l=left.B->size2;
	gsl_vector_view col;
	while (left.t < right.t){
		ODE_jac(t,y->data,(left.A)->data,(left.dfdt)->data,p);
		ODE_jacp(t,y->data,(left.B)->data,(left.dfdt)->data,p);

		left.t = t;
		status = gsl_odeiv2_evolve_apply(driver->e, driver->c, driver->s, &sys, &t, tf, &h, y->data); // t is moved forward unless it reached tf
		right.t = t;
		delta_t = right.t - left.t;

		ODE_jac(t,y->data,(right.A)->data,(right.dfdt)->data,p);
		ODE_jacp(t,y->data,(right.B)->data,(right.dfdt)->data,p);

		if (status!=GSL_SUCCESS) {
			check_status(status,t,right.t,-1);
			return(status);
		}
		/* here we branch out based on the invertibility of the jacobian
		   and proximity to a stable steady state */
		if (qrtp(qr,left.A)==GSL_SUCCESS && rcond(qr) > RCOND_LIMIT && stable(left.A,ev)){
			for (k=0;k<l;k++){
				col=gsl_matrix_column(left.B,k);
				gsl_linalg_QRPT_svx(qr.QR, qr.tau, qr.p, &(col.vector));  /* B <- A\B*/
			}
			gsl_matrix_scale(left.A, delta_t);                          /* A <- (df/dy)*(t-t0)*/
			gsl_linalg_exponential_ss(left.A, ev.T, GSL_PREC_SINGLE);   /* T <- exp(A*(t-t0))*/

			gsl_matrix_memcpy(left.temp,left.Sy);
			gsl_matrix_add(left.temp,left.B);

			gsl_matrix_memcpy(right.Sy,left.B);
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ev.T, left.temp, -1.0, right.Sy);
			/* M.Sy is now the piece-wise-constant approximation of the sensitivity */
		} else {
			transition_matrix(left.A,right.A,left.t,right.t,left.phi);
			transition_matrix(right.A,right.A,right.t,left.t,right.phi);
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, right.phi, right.B, 1.0, left.A);
			// A = jac + PHI_bwd * jacp
			gsl_matrix_scale(left.A,0.5*delta_t);
			gsl_matrix_add(left.A,left.Sy); // A = Sy + (jac + PHI_bwd * jacp)*dt/2
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, left.phi, left.A, 0.0, right.Sy);
			// Sy = PHI_fwd (Sy + (jac + PHI_bwd * jacp)*dt/2)
		}
	}
	return status;
}

/* Intergrates the system `sys` using the specified `driver` and
   simulation instructions. But here, we use the solver's step
   function directly, so that we can also step forward the sensitivity
   approximation.*/
int /* error code (GSL_SUCCESS) */
simulate_timeseries_and_sensitivity(
	const gsl_odeiv2_system sys, /* the system to integrate */
	gsl_odeiv2_driver* driver, /* the driver that is used to integrate `sys` */
	const double t0, /* the initial time: y(t0) = y0 */
	const gsl_vector *y0, /* initial value */
	const gsl_vector *time, /* a vector of time-points */
	struct pbsrMem left,
	struct pbsrMem right,
	struct QRPT qr,
	struct eigen ev,
	const struct event *event, /*a struct array with scheduled events */
	double *sens,     /* (OUT) return value, sensitivity */
	double *func_sens,/* (OUT) return value, function sensitivity */
	gsl_matrix *Yout) /* (OUT) return value, state */
{
	gsl_set_error_handler_off();
	int nt=time->size;
	struct pbsrMem temp;
	int ny=(int) sys.dimension;
	gsl_vector *y=gsl_vector_alloc(ny);
	gsl_vector *z=gsl_vector_alloc(ny);
	gsl_vector_memcpy(y,y0);
	double *p=sys.params;
	double h=1e-2;
	gsl_vector_view Yout_row;
	gsl_matrix_view Sy;
	gsl_matrix_view Sf;
	int i=0,j;
	double t=t0;
	double ti,tf,te,td;
	int status=GSL_SUCCESS;
	// gsl_vector *dfdt=gsl_vector_alloc(ny);
	int np = ODE_default(0.0,NULL);
	int nf = ODE_func(0.0,NULL,NULL,NULL);
	/* initialize t0 values */
	gsl_vector_memcpy(y,y0);
	gsl_matrix_set_zero(left.Sy);
	for (j=0, i=0; j<nt; j++){
		td=gsl_vector_get(time,j); // time of data observation
		left.t = t;
		if (event && i<event->nt && event->time[i] < td) {
			right.t = event->time[i]; // time of next event
			status = advance(sys, driver, y, left, right, qr, ev);
			t = right.t;
			// swap left and right (it's a bunch of pointers)
			temp = left; left = right; right = temp;
			apply_tf(event->state,y->data,i);
			apply_tf(event->par,p,i);
			status=gsl_odeiv2_driver_reset(driver);
			if (status!=GSL_SUCCESS){
#ifdef DEBUG_PRINT
				fprintf(stderr,"[%s] resetting the system after event %i produced an error: %s.\n",__func__,i,gsl_strerror(status));
#endif
				return(status);
			}
			i++;
		}
		right.t = td;
		if (td>t) status=advance(sys, driver, y, left, right, qr, ev);
		t = right.t;
		if(status==GSL_SUCCESS){
			ODE_funcJac(t,y->data,(right.FA)->data,p);
			ODE_funcJacp(t,y->data,(right.Sf)->data,p);
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, right.FA, right.Sy, 1.0, right.Sf);

			Sy = gsl_matrix_view_array(sens+j*ny*np,np,ny);
			Sf = gsl_matrix_view_array(func_sens+j*nf*np,np,nf);
			gsl_matrix_transpose_memcpy(&Sy.matrix, right.Sy);
			//report any error codes to the R user
			Yout_row = gsl_matrix_row(Yout,j);
			gsl_vector_memcpy(&(Yout_row.vector),y);
		} else {
			check_status(status,t,td,j);
			return(status);
		}
	}
	gsl_odeiv2_driver_reset(driver);
	return status;
}

// allocates a matrix of parameters, onw row per experiment;
gsl_matrix* params_alloc(Rdata experiments){
	int i,N;
	gsl_matrix *P=NULL;
	int np,nu;
	gsl_vector *p=NULL;
	Rdata input;
	double* u;
	gsl_vector_view row;
	if (ODE_default) {
		np=ODE_default(0.0,NULL);
		p=gsl_vector_alloc(np);
		ODE_default(0.0,p->data);
	}
	if (experiments != R_NilValue && IS_VECTOR(experiments)){
		N=length(experiments);
		P=gsl_matrix_alloc(N,p->size);
		for (i=0;i<N;i++){
			input = from_list(VECTOR_ELT(experiments,i),"input Input u inputs Inputs");
			nu = (input != R_NilValue) ? length(input) : 0;
			row = gsl_matrix_row(P,i);
			gsl_vector_memcpy(&row.vector,p);
			if (nu>0 && np>=nu){
				u = gsl_matrix_ptr(P,i,np-nu);
				memcpy(u, REAL(input), nu*sizeof(double));
			}
		}
	} else {
		P = gsl_matrix_alloc(1,p->size);
		row = gsl_matrix_row(P,0);
		gsl_vector_memcpy(&row.vector,p);
	}
	if (p) gsl_vector_free(p);
	return P;
}

// allocates a matrix of initial values, one row per experiment
gsl_matrix* init_alloc(const gsl_odeiv2_system sys, gsl_matrix *p){
	int ny = sys.dimension;
	int i,N = p->size1;
	gsl_matrix *Y0 = gsl_matrix_alloc(N,ny);
	double t = 0.0;
	gsl_matrix_set_zero(Y0);
	if (ODE_init){
		for (i=0;i<N;i++){
			ODE_init(t,gsl_matrix_ptr(Y0,i,0),gsl_matrix_ptr(p,i,0));
		}
	}
	return Y0;
}

void update_initial_values(gsl_vector *y0, gsl_odeiv2_system sys, Rdata iv){
	gsl_vector_view r_init;
	if (!y0) return;
	if (!sys.params) return;
	if (iv != R_NilValue && IS_VECTOR(iv) && length(iv) == y0->size){
		r_init = gsl_vector_view_array(REAL(AS_NUMERIC(iv)),length(iv));
		gsl_vector_memcpy(y0,&(r_init.vector));
	} else if (ODE_init!=NULL && y0->size == sys.dimension) {
		ODE_init(0.0,y0->data,(double*) sys.params);
	}
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

struct sensApproxMem {
	gsl_matrix *A;
	gsl_matrix *LU;
	gsl_matrix *B;
	gsl_matrix *E;
	gsl_matrix *S_AB;
	gsl_permutation *P;
	gsl_matrix *FA;
	gsl_matrix *Sy;
	gsl_matrix *Sf;
};

/* n: number of state variables (ODE_vf);
 * l: number of parameters (ODE_default);
 * f: number of output values (ODE_func);
 */
struct sensApproxMem sensApproxMemAlloc(size_t n, size_t l, size_t f){
	struct sensApproxMem M; /* contains a bunch of pointers */
	M.A=gsl_matrix_alloc(n,n);
	M.LU=gsl_matrix_alloc(n,n);
	M.B=gsl_matrix_alloc(n,l);
	M.E=gsl_matrix_alloc(n,n);
	M.S_AB=gsl_matrix_alloc(n,l);
	M.P=gsl_permutation_alloc(n);
	M.FA=gsl_matrix_alloc(f,n);
	M.Sy=gsl_matrix_alloc(n,l);
	M.Sf=gsl_matrix_alloc(f,l);
	return M; /* bunch of pointers, by value */
}

void sensApproxMemFree(struct sensApproxMem M){
	if (M.A) gsl_matrix_free(M.A);
	if (M.LU) gsl_matrix_free(M.LU);
	if (M.B) gsl_matrix_free(M.B);
	if (M.E) gsl_matrix_free(M.E);
	if (M.S_AB) gsl_matrix_free(M.S_AB);
	if (M.P) gsl_permutation_free(M.P);
	if (M.FA) gsl_matrix_free(M.FA);
	if (M.Sy) gsl_matrix_free(M.Sy);
	if (M.Sf) gsl_matrix_free(M.Sf);
}

void transition_matrix(gsl_matrix *Ji, /* the jacobian at t=ti */
 gsl_matrix *Jf, /* the jacobian at t=tf */
 double ti, /* initial time of the interval (left bondary) */
 double tf,/* final time of the interval (right boundary) */
 gsl_matrix *phi) /* return buffer */
{
  double s=0.5*(tf-ti);
  size_t ny=Ji->size1;
  int i,n=3;
  gsl_matrix *V=gsl_matrix_alloc(ny,ny);
  gsl_matrix *W=gsl_matrix_alloc(ny,ny);
  gsl_matrix *I1=gsl_matrix_alloc(ny,ny);  // I1(tf;ti) = 0.5*(tf-ti)*(Jf+Ji)
  gsl_matrix_set_identity(phi);
  // I1
  gsl_matrix_memcpy(I1,Jf);
  gsl_matrix_add(I1,Ji);
  gsl_matrix_scale(I1,s);
  gsl_matrix_set_identity(W);
  for (i=0;i<n;i++){
    gsl_matrix_set_identity(V);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, s, W, Jf, 1.0, V);
    gsl_matrix_memcpy(W,V);
  }
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, W, I1, 1.0, phi);

  gsl_matrix_free(I1);
  gsl_matrix_free(V);
  gsl_matrix_free(W);
}

int sensitivityApproximation(double t0, gsl_vector *t, gsl_vector *p, gsl_matrix *Y, double *dYdp, double *dFdp, struct sensApproxMem M)
{
	int j,k;
	int m=t->size;
	int n=Y->size2;
	int l=p->size;
	int f=ODE_func(0,NULL,NULL,NULL);
	double tj,delta_t;
	gsl_vector_view col, diag;
	gsl_matrix *A=M.A;
	gsl_matrix *LU=M.LU;
	gsl_matrix *B=M.B;
	gsl_matrix *E=M.E;
	gsl_matrix *S_AB=M.S_AB;
	gsl_permutation *P=M.P;
	gsl_matrix *FA=M.FA;
	gsl_matrix_view SY0, SY, SF; /* sensitivity matrices (array-views) */
	int sign;
	clock_t ct1,ct2,ct3,ct4;
	clock_t ct_expm, ct_svx, ct_LU, ct_dgemm[2];
	clock_t total_ct_expm=0, total_ct_svx=0, total_ct_LU=0, total_ct_dgemm[2]={0,0};
	const double *y;
	if (m != Y->size1) fprintf(stderr,"[%s] t has length %i, but Y has %li rows.\n",__func__,m,Y->size1);

	for (j=0;j<m;j++){
		ct1=clock();
		tj=gsl_vector_get(t,j);
		delta_t=tj-t0;
		t0=tj;
		y=gsl_matrix_ptr(Y,j,0);
		ODE_jac(tj,y,A->data,NULL,p->data);                                           /* A <- df/dy */
		ODE_jacp(tj,y,B->data,NULL,p->data);                                          /* B <- df/dp */
		gsl_matrix_memcpy(LU,A);                                                      /* LU <- A */
		// dirty hack .. maybe this will work better?
		diag = gsl_matrix_diagonal(A);
		gsl_vector_add_constant(&(diag.vector),-1e-9); // all eigenvalues should be negative, if some are close to zero, we force them
		// .. dirty hack
		ct_LU = clock();
		gsl_linalg_LU_decomp(LU, P, &sign);                                           /* make P*A = L*U */
		ct_LU = clock() - ct_LU;
		ct_svx=clock();
		for (k=0;k<l;k++){
			col=gsl_matrix_column(B,k);
			gsl_linalg_LU_svx(LU, P, &(col.vector));                                    /* B <- A\B*/
		}
		ct_svx=clock()-ct_svx;
		gsl_matrix_scale(A,delta_t);                                                  /* A <- (df/dy)*(t-t0)*/
		ct2=clock();
		ct_expm=clock();
		gsl_linalg_exponential_ss(A,E,GSL_PREC_SINGLE);                               /* E <- exp(A*(t-t0))*/
		ct_expm=clock()-ct_expm;
		ct3=clock();

		if (j==0){
			gsl_matrix_memcpy(S_AB,B);
		} else {
			gsl_matrix_memcpy(S_AB,M.Sy);
			gsl_matrix_add(S_AB,B);
		}
		gsl_matrix_memcpy(M.Sy,B);
		ct_dgemm[0]=clock();
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, E, S_AB, -1.0, M.Sy);         /* S is now the piece-wise-constant approximation of the sensitivity */
		ct_dgemm[0]=clock()-ct_dgemm[0];
		/* write result to output buffer: */
		SY=gsl_matrix_view_array(dYdp+n*l*j,l,n);
		gsl_matrix_transpose_memcpy(&(SY.matrix),M.Sy);

		ODE_funcJac(tj,y,FA->data,p->data);
		ODE_funcJacp(tj,y,M.Sf->data,p->data);

		ct_dgemm[1]=clock();
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, FA, M.Sy, 1.0, M.Sf);
		SF=gsl_matrix_view_array(dFdp+f*l*j,l,f);
		gsl_matrix_transpose_memcpy(&(SF.matrix),M.Sf);
		ct_dgemm[1]=clock()-ct_dgemm[1];
		ct4=clock();
#ifdef DEBUG_PRINT
		fprintf(stderr,"[%s] timePoint %i of %i: %g s (%g s for matrix exponential).\n",__func__,j,m,sec(ct4-ct1),sec(ct3-ct2));
#endif
		total_ct_expm+=ct_expm;
		total_ct_svx+=ct_svx;
		total_ct_LU+=ct_LU;
		total_ct_dgemm[0]+=ct_dgemm[0];
		total_ct_dgemm[1]+=ct_dgemm[1];
	}
	/* stats: */
	//fprintf(stderr,"%li\t%li\t%li\t%li\t",total_ct_expm, total_ct_svx, total_ct_LU, total_ct_dgemm[0]+total_ct_dgemm[1]);
	//fprintf(stderr,"%g\t%g\t%g\t%g\n",sec(total_ct_expm), sec(total_ct_svx), sec(total_ct_LU), sec(total_ct_dgemm[0])+sec(total_ct_dgemm[1]));
	return GSL_SUCCESS;
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
r_gsl_odeiv2_outer_state_only(
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
	int i,j,k;
	double abs_tol=asReal(absolute_tolerance);
	double rel_tol=asReal(relative_tolerance);
	double h=asReal(initial_step_size);
	int N=GET_LENGTH(experiments);
	size_t M=ncols(parameters);
	const gsl_odeiv2_step_type * T=gsl_odeiv2_step_msbdf;

	Rdata res_list = PROTECT(NEW_LIST(N)); /* use VECTOR_ELT and SET_VECTOR_ELT */
	SET_NAMES(res_list,GET_NAMES(experiments));
	Rdata yf_list, Y, iv, t;
	double t0;
	const char *yf_names[1]={"state"};
	gsl_vector_view initial_value, time;
	gsl_matrix_view y;
	size_t ny, nt;
	struct event *ev=NULL;
	gsl_odeiv2_system sys = load_system(model_name, model_so); /* also sets ODE_*() functions */
	if (ODE_default==NULL || ODE_init==NULL){
		fprintf(stderr,"[%s] loading model has failed.\n",__func__);
		UNPROTECT(1); /* res_list */
		return R_NilValue;
	}
	gsl_matrix *P = params_alloc(experiments);
	gsl_matrix *Y0 = init_alloc(sys,P);
	gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys,T,h,abs_tol,rel_tol);
	for (i=0; i<N; i++){
		//fprintf(stderr,"[%s] experiment %i of %i.\n",__func__,i,N); fflush(stderr);
		iv = from_list(VECTOR_ELT(experiments,i),"initial_value initialState");
		t = from_list(VECTOR_ELT(experiments,i),"time outputTimes");
		t0 = REAL(AS_NUMERIC(from_list(VECTOR_ELT(experiments,i),"initialTime t0 T0")))[0];
		ev = event_from_R(from_list(VECTOR_ELT(experiments,i),"event events scheduledEvents scheduled_events"));
		sys.params = gsl_matrix_ptr(P,i,0);
		initial_value = gsl_matrix_row(Y0,i);
		nt = length(t);
		ny = Y0->size2;
		time=gsl_vector_view_array(REAL(AS_NUMERIC(t)),nt);
		Y=PROTECT(alloc3DArray(REALSXP,ny,nt,M));
		for (j=0; j<ny*nt*M; j++) REAL(Y)[j]=NA_REAL; /* initialize to NA */
		for (k=0;k<M;k++){
			y=gsl_matrix_view_array(REAL(AS_NUMERIC(Y))+(nt*ny*k),nt,ny);
			memcpy((double*) sys.params, REAL(AS_NUMERIC(parameters))+nrows(parameters)*k, nrows(parameters)*sizeof(double));
			update_initial_values(&initial_value.vector, sys, iv); // in case system parameters have an effect on them
			simulate_timeseries(
				sys,
				driver,
				t0,
				&(initial_value.vector),
				&(time.vector),
				ev,
				&(y.matrix)
			);
		}
		yf_list=PROTECT(NEW_LIST(1));
		SET_VECTOR_ELT(yf_list,0,Y);
		set_names(yf_list,yf_names,1);
		SET_VECTOR_ELT(res_list,i,yf_list);
		event_free(&ev);

		UNPROTECT(1); /* yf_list */
		UNPROTECT(1); /* Y */
	} // experiments: 0 to N-1
	gsl_matrix_free(P);
	gsl_matrix_free(Y0);
	gsl_odeiv2_driver_free(driver);
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
r_gsl_odeiv2_outer_func(
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
	int i,j,k,status=GSL_SUCCESS;
	double abs_tol=asReal(absolute_tolerance);
	double rel_tol=asReal(relative_tolerance);
	double h=asReal(initial_step_size);
	int N=GET_LENGTH(experiments);
	size_t M=ncols(parameters);
	const gsl_odeiv2_step_type * T=gsl_odeiv2_step_msbdf;

	Rdata res_list = PROTECT(NEW_LIST(N)); /* use VECTOR_ELT and SET_VECTOR_ELT */
	SET_NAMES(res_list,GET_NAMES(experiments));
	Rdata yf_list, Y, F, iv, t;
	double t0;
	const char *yf_names[2]={"state","func"};
	gsl_vector_view initial_value, time;
	gsl_matrix_view y;
	size_t nt;
	struct event *ev=NULL;
	double *f;
	gsl_odeiv2_system sys = load_system(model_name, model_so); /* also sets ODE_*() functions */
	if (ODE_default==NULL || ODE_init==NULL || ODE_func==NULL){
		fprintf(stderr,"[%s] loading model has failed.\n",__func__);
		UNPROTECT(1); /* res_list */
		return R_NilValue;
	}
	size_t nf = ODE_func(0,NULL,NULL,NULL);
	size_t ny = sys.dimension;
	gsl_matrix *P = params_alloc(experiments);
	gsl_matrix *Y0 = init_alloc(sys,P);
	gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys,T,h,abs_tol,rel_tol);

	for (i=0; i<N; i++){
		//fprintf(stderr,"[%s] experiment %i of %i.\n",__func__,i,N); fflush(stderr);
		iv = from_list(VECTOR_ELT(experiments,i),"initial_value initialState initialValue intialValues");
		t = from_list(VECTOR_ELT(experiments,i),"time outputTimes t");
		t0 = REAL(AS_NUMERIC(from_list(VECTOR_ELT(experiments,i),"initial_time initialTime t0 T0")))[0];
		ev = event_from_R(from_list(VECTOR_ELT(experiments,i),"event events scheduledEvents scheduledEvent scheduled_events"));
		sys.params = gsl_matrix_ptr(P,i,0);
		initial_value = gsl_matrix_row(Y0,i);
		nt = length(t);
		time=gsl_vector_view_array(REAL(AS_NUMERIC(t)),nt);
		Y=PROTECT(alloc3DArray(REALSXP,ny,nt,M));
		for (j=0; j<ny*nt*M; j++) REAL(Y)[j]=NA_REAL; /* initialize to NA */

		F=PROTECT(alloc3DArray(REALSXP,nf,nt,M));
		for (j=0;j<nf*nt*M;j++) REAL(F)[j]=NA_REAL;   /* initialize to NA */
		for (k=0;k<M;k++){
			y=gsl_matrix_view_array(REAL(AS_NUMERIC(Y))+(nt*ny*k),nt,ny);
			memcpy((double*) sys.params, REAL(AS_NUMERIC(parameters))+nrows(parameters)*k, nrows(parameters)*sizeof(double));
			update_initial_values(&initial_value.vector,sys,iv);
			status=simulate_timeseries(
				sys,
				driver,
				t0,
				&(initial_value.vector),
				&(time.vector),
				ev,
				&(y.matrix)
			);
			if (status==GSL_SUCCESS) {
				for (j=0;j<nt;j++){
					f=REAL(F)+(0+j*nf+k*nf*nt);
					ODE_func(gsl_vector_get(&(time.vector),j),gsl_matrix_ptr(&(y.matrix),j,0),f,sys.params);
				}
			}
		}
		yf_list=PROTECT(NEW_LIST(2));
		SET_VECTOR_ELT(yf_list,0,Y);
		SET_VECTOR_ELT(yf_list,1,F);
		set_names(yf_list,yf_names,2);
		SET_VECTOR_ELT(res_list,i,yf_list);
		event_free(&ev);
		UNPROTECT(1); /* yf_list */
		UNPROTECT(1); /* F */
		UNPROTECT(1); /* Y */
	} // experiments: 0 to N-1
	gsl_matrix_free(P);
	fprintf(stderr,"[%s] %s final free() statement.\n",__func__,model_name); fflush(stderr);
	//gsl_matrix_free(Y0);
	gsl_odeiv2_driver_free(driver);
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
r_gsl_odeiv2_outer_sens(
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
	int i,j,k,status=GSL_SUCCESS;
	double abs_tol=asReal(absolute_tolerance);
	double rel_tol=asReal(relative_tolerance);
	double h=asReal(initial_step_size);
	int N=GET_LENGTH(experiments);
	size_t M=ncols(parameters);
	const gsl_odeiv2_step_type * T=gsl_odeiv2_step_msbdf;

	Rdata res_list = PROTECT(NEW_LIST(N)); /* use VECTOR_ELT and SET_VECTOR_ELT */
	SET_NAMES(res_list,GET_NAMES(experiments));
	Rdata yf_list, Y, F, iv, t;
	double t0;
	const char *yf_names[4]={"state","func","stateSensitivity","funcSensitivity"};
	gsl_vector_view initial_value, time;
	gsl_matrix_view y;
	size_t nt;
	struct event *ev=NULL;
	double *f;
	gsl_vector_view p;
	Rdata SY, SF;
	Rdata sy_k, sf_k;

	gsl_odeiv2_system sys = load_system(model_name, model_so); /* also sets ODE_*() functions */
	if (sys.dimension == 0 || ODE_default==NULL || ODE_init==NULL || ODE_func==NULL || ODE_funcJac==NULL || ODE_funcJacp==NULL){
		fprintf(stderr,"[%s] loading model has failed (system dimension: «%li»).\n",__func__,sys.dimension);
		UNPROTECT(1); /* res_list */
		return R_NilValue;
	}
	gsl_matrix *P = params_alloc(experiments);
	gsl_matrix *Y0 = init_alloc(sys,P);
	int nf = ODE_func(0,NULL,NULL,NULL);
	int ny = sys.dimension;
	int np = P->size2;
	struct sensApproxMem saMem = sensApproxMemAlloc(ny,np,nf);
	//fprintf(stderr,"expm\tsvx\tLU\tdgemm\texpm/s\tsvx/s\tLU/s\tdgemm/s\n"); // sensitivityApproximation
	//fprintf(stderr,"----\t---\t--\t-----\t------\t-----\t----\t-------\n"); // will print times
	gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys,T,h,abs_tol,rel_tol);
	for (i=0; i<N; i++){
		iv = from_list(VECTOR_ELT(experiments,i),"initial_value initialState initialValue initialValues");
		t = from_list(VECTOR_ELT(experiments,i),"time outputTimes t");
		t0 = REAL(AS_NUMERIC(from_list(VECTOR_ELT(experiments,i),"intial_time initialTime t0 T0")))[0];
		ev = event_from_R(from_list(VECTOR_ELT(experiments,i),"event events scheduledEvents scheduledEvent scheduled_event"));
		initial_value = gsl_matrix_row(Y0,0);
		sys.params = gsl_matrix_ptr(P,i,0);
		nt = length(t);
		time=gsl_vector_view_array(REAL(AS_NUMERIC(t)),nt);
		ny = Y0->size2;
		Y=PROTECT(alloc3DArray(REALSXP,ny,nt,M));
		for (j=0; j<ny*nt*M; j++) REAL(Y)[j]=NA_REAL; /* initialize to NA */
		F=PROTECT(alloc3DArray(REALSXP,nf,nt,M));
		SY=PROTECT(NEW_LIST(M));
		SF=PROTECT(NEW_LIST(M));
		for (j=0;j<nf*nt*M;j++) REAL(F)[j]=NA_REAL;   /* initialize to NA */
		for (k=0;k<M;k++){
			y=gsl_matrix_view_array(REAL(AS_NUMERIC(Y))+(nt*ny*k),nt,ny);
			memcpy((double*) sys.params, REAL(AS_NUMERIC(parameters))+nrows(parameters)*k, nrows(parameters)*sizeof(double));
			update_initial_values(&initial_value.vector,sys,iv);
			status=simulate_timeseries(
				sys,
				driver,
				t0,
				&(initial_value.vector),
				&(time.vector),
				ev,
				&(y.matrix)
			);
			sy_k=PROTECT(alloc3DArray(REALSXP,ny,np,nt));
			sf_k=PROTECT(alloc3DArray(REALSXP,nf,np,nt));
			if (status==GSL_SUCCESS) {
				p = gsl_matrix_row(P,i);
				for (j=0;j<nt;j++){
					f=REAL(F)+(0+j*nf+k*nf*nt);
					ODE_func(gsl_vector_get(&(time.vector),j),gsl_matrix_ptr(&(y.matrix),j,0),f,sys.params);
				}
				sensitivityApproximation(t0,&(time.vector),&(p.vector),&(y.matrix),REAL(sy_k),REAL(sf_k),saMem);
			}
			SET_VECTOR_ELT(SY,k,sy_k);
			SET_VECTOR_ELT(SF,k,sf_k);
			UNPROTECT(2); // sy_k, sf_k;
		}
		yf_list=PROTECT(NEW_LIST(4));
		SET_VECTOR_ELT(yf_list,0,Y);
		SET_VECTOR_ELT(yf_list,1,F);
		SET_VECTOR_ELT(yf_list,2,SY);
		SET_VECTOR_ELT(yf_list,3,SF);
		set_names(yf_list,yf_names,4);
		SET_VECTOR_ELT(res_list,i,yf_list);
		event_free(&ev);

		UNPROTECT(1); /* yf_list */
		UNPROTECT(2); /* SY, SF */
		UNPROTECT(1); /* F */
		UNPROTECT(1); /* Y */
#ifdef DEBUG_PRINT
		if (status!=GSL_SUCCESS){
			fprintf(stderr,"[%s] parameter set lead to solver errors (%s) in experiment %i/%i, values:\n",__func__,gsl_strerror(status),i,N);
			for (j=0; j<np_model; j++) fprintf(stderr,"%g%s",p[j],(j==np_model-1?"\n":", "));
		}
#endif
	} // experiments: 0 to N-1
  sensApproxMemFree(saMem); /* frees the memory of temporary matrices of sensitivity approximation */
	UNPROTECT(1); /* res_list */
	gsl_odeiv2_driver_free(driver);
	gsl_matrix_free(P);
	gsl_matrix_free(Y0);
	return res_list;
}
