/*
 *  HarmonicOscillator_gvf.c
 *
 *  GSL C file for the vector field named: HarmonicOscillator
 *
 *  This file was generated by the program VFGEN, version: 2.6.0.dev1
 *  Generated on 19-Oct-2021 at 14:43
 */

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <stdio.h>
#include <stdlib.h>
/*
 *  The vector field.
 */

int HarmonicOscillator_vf(double t, const double y_[], double f_[], void *params)
{
    double v, y;
    double k, c, F;
    double v_flux;
    double *p_;
    int RET=GSL_SUCCESS;
    if (y_ && f_ && params){
      p_ = (double *) params;

      v          = y_[0];
      y          = y_[1];

      k          = p_[0];
      c          = p_[1];
      F          = p_[2];

      v_flux = -k*y+F-v*c;

      f_[0] = v_flux;
      f_[1] = v;
    } else {
      RET=2;
    }
    return RET;
}

/*
 *  The Jacobian.
 */

int HarmonicOscillator_jac(double t, const double y_[], double *jac_, double *dfdt_, void *params)
{
    double v, y;
    double k, c, F;
    double *p_;
    int RET=GSL_SUCCESS;
    
    if (y_ && jac_ && dfdt_ && params){
      p_ = (double *) params;

      v          = y_[0];
      y          = y_[1];

      k          = p_[0];
      c          = p_[1];
      F          = p_[2];

      gsl_matrix_view dfdy_mat = gsl_matrix_view_array(jac_,2,2);
      gsl_matrix *m_ = &dfdy_mat.matrix;

      gsl_matrix_set(m_, 0, 0, -c);
      gsl_matrix_set(m_, 0, 1, -k);
      gsl_matrix_set(m_, 1, 0, 1.0);
      gsl_matrix_set(m_, 1, 1, 0.0);

      dfdt_[0] = 0.0;
      dfdt_[1] = 0.0;
    } else {
      RET=2*2;
    }
    return RET;
}

/*
 *  The Jacobian with respect to the parameters.
 */

int HarmonicOscillator_jacp(double t, const double y_[], double *jacp_, void *params)
{
    double v, y;
    double k, c, F;
    double *p_;
    int RET=GSL_SUCCESS;
    if (y_ && jacp_ && params){
      p_ = (double *) params;

      v          = y_[0];
      y          = y_[1];

      k          = p_[0];
      c          = p_[1];
      F          = p_[2];

      gsl_matrix_view dfdp_mat = gsl_matrix_view_array(jacp_,2,3);
      gsl_matrix *m_ = &dfdp_mat.matrix;

      gsl_matrix_set(m_, 0, 0, -y);
      gsl_matrix_set(m_, 0, 1, -v);
      gsl_matrix_set(m_, 0, 2, 1.0);
      gsl_matrix_set(m_, 1, 0, 0.0);
      gsl_matrix_set(m_, 1, 1, 0.0);
      gsl_matrix_set(m_, 1, 2, 0.0);
    } else {
      RET=2*3;
    }
    return GSL_SUCCESS;
}

/*
 *  User function
 */
int HarmonicOscillator_func(double t, const double y_[], double *f, void *params)
{
	double v, y;
	double k, c, F;
	double v_flux;
	double *p_ = params;
	if (!f || !y_ || !p_) return 1;
	v          = y_[0];
	y          = y_[1];
	k          = p_[0];
	c          = p_[1];
	F          = p_[2];

	v_flux = -k*y+F-v*c;

	f[0] = sqrt(y*y+v*v);
	//fprintf(stderr,"[%s] sqrt((%g)^2 + (%g)^2) = %g\n",__func__,y,v,f[0]);
	//fflush(stderr);
	return GSL_SUCCESS;
}

/* default parameters */
int HarmonicOscillator_default(double t, void *params)
{
	double *p_ = params;
	if (!p_) return 3;
	p_[0]=1.0;
	p_[1]=0.0;
	p_[2]=0.0;
	return GSL_SUCCESS;
}

int HarmonicOscillator_init(double t, double y_[], void *params)
{
	double v, y;
	double *p_ = params;
	if (!y_ || !params) return 2;
	y_[0]=0.0;
	y_[1]=1.0;
	return GSL_SUCCESS;
}
