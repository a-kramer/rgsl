#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <dlfcn.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#define BATCH_SIZE 100
#define SEED 12653416L

#define max(a,b) (((a)>(b))?(a):(b))
/* SEXP stands for S-Expression, and it can be any R data object (or
 * function) in this program, we'll only use data from R. However SEXP is
 * a bad type name; I am compelled to pronounce it in my head ...
 */
typedef SEXP Rdata;

double frac(int a, int b){
	double y;
	if (b!=0) y=((double) a)/((double) b);
	else y=NAN;
	return y;
}

Rdata estimate_sigma(Rdata r_sample){
	int i,j,k,l;
	int n=ncols(r_sample);
	int m=nrows(r_sample);
	double *sample=REAL(r_sample);
	Rdata r_sigma=PROTECT(NEW_NUMERIC(1));
	int N=max(10,n/BATCH_SIZE); // batches
	gsl_vector_view sample_column,xcol;
	gsl_vector *diff=gsl_vector_alloc(m);
	gsl_vector *g=gsl_vector_alloc(m);
	double sigma=0.1;
	double d;
	gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(r,SEED);
	int a;
	for (j=0;j<N;j++){
		i=gsl_rng_uniform(r)*n;
		xcol=gsl_vector_view_array(sample+i*m,m);
		a=0;
		d=0;
		for (k=0;k<BATCH_SIZE;k++){
			i=gsl_rng_uniform(r)*n;
			sample_column=gsl_vector_view_array(sample+i*m,m);
			for (l=0;l<m;l++) gsl_vector_set(g,l,gsl_ran_gaussian(r,sigma));
			gsl_vector_add(g,&(sample_column.vector));
			gsl_vector_memcpy(diff,&(sample_column.vector));
			gsl_vector_sub(diff,&(xcol.vector)); // diff = x-sample[,j]
			gsl_blas_ddot(diff,diff,&d);  // d = (x-sample[,j])^T * (x-sample[,j])
			d/=gsl_pow_2(sigma);
			a+=(sqrt(d)<6);
		}
		sigma*=2*(1.0-0.9*sqrt(sqrt(frac(a,BATCH_SIZE))));
	}
	gsl_rng_free(r);
	gsl_vector_free(diff);
	gsl_vector_free(g);
	REAL(r_sigma)[0]=sigma;
	UNPROTECT(1);
	return r_sigma;
}

/*! kde calculates density estimates

		The sample is a matrix of continuous observations. The evaluation
		will happen at every point (vector) in the matrix x.

		@param sample a matrix with sampled points from the target
		distribution, as columns.

		@param x a point at which the kde is evaluated. If this is a matrix
		then the return value will be a vector of densities for each volumn.

		@return probability density estimates
 */
Rdata dkde(Rdata sample, Rdata x, Rdata r_sigma){
	int i,j,k,l;
	int ridx;
	int n=ncols(sample);
	int m=nrows(sample);
	Rdata D=PROTECT(NEW_NUMERIC(n));
	double d,E;
	double sigma = REAL(r_sigma)[0];//estimate_sigma(REAL(sample),n,m);
	double C=gsl_pow_int(sqrt(2*M_PI),m) * gsl_pow_int(sigma,m);
	gsl_vector_view sample_column,xcol;
	gsl_vector *diff=gsl_vector_alloc(m);
	gsl_vector *g=gsl_vector_alloc(m);
	gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(r,SEED);
	int N=max(10,n/BATCH_SIZE); // batches
	double f,h;
	double abs_tol=1e-14;
	double rel_tol=1e-14;

	for (i=0;i<ncols(x);i++){
		xcol=gsl_vector_view_array(REAL(x)+i*m,m);
		for (j=0;j<N;j++){
			E=0;
			REAL(D)[i]=0;
			for (k=0;k<BATCH_SIZE;k++){
				ridx=gsl_rng_uniform(r)*n;
				sample_column=gsl_vector_view_array(REAL(sample)+ridx*m,m);
				for (l=0;l<m;l++) gsl_vector_set(g,l,gsl_ran_gaussian(r,sigma));
				gsl_vector_add(g,&(sample_column.vector));
				gsl_vector_memcpy(diff,&(sample_column.vector));
				gsl_vector_sub(diff,&(xcol.vector)); // diff = x-sample[,j]
				gsl_blas_ddot(diff,diff,&d);  // d = (x-sample[,j])^T * (x-sample[,j])
				d/=gsl_pow_2(sigma);
				E+=exp(-0.5*d)/C;
			}
			f=frac(j,j+1);
			h=E/((double) (j+1)*BATCH_SIZE);
			REAL(D)[i] = f * (REAL(D)[i]) + h;
		}
	}
	gsl_vector_free(diff);
	gsl_vector_free(g);
	gsl_rng_free(r);
	UNPROTECT(1);
	return D;
}
