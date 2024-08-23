#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

/* The error code indicates how to pre-allocate memory
 * for output values such as `f_`. The _vf function returns
 * the number of state variables, if any of the args are `NULL`.
 * evaluation errors can be indicated by negative return values.
 * GSL_SUCCESS (0) is returned when no error occurred.
 */

/* ode vector field: y'=f(t,y;p) */
int HarmonicOscillator_vf(double t, const double y_[], double f_[], void *par)
{
	double *p_=par;
	if (!y_ || !f_) return 2;
	double k=p_[0];
	double c=p_[1];
	double F=p_[2];
	double v=y_[0];
	double y=y_[1];
	f_[0] = F-k*y-c*v;
	f_[1] = v;
	return GSL_SUCCESS;
}
/* ode Jacobian df(t,y;p)/dy */
int HarmonicOscillator_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jac_) return 2*2;
	double k=p_[0];
	double c=p_[1];
	double F=p_[2];
	double v=y_[0];
	double y=y_[1];
/* column 1 (df/dy_0) */
	jac_[0] = -c; /* [0, 0] */
	jac_[2] = 1; /* [1, 0] */
/* column 2 (df/dy_1) */
	jac_[1] = -k; /* [0, 1] */
	jac_[3] = 0; /* [1, 1] */
	return GSL_SUCCESS;
}
/* ode parameter Jacobian df(t,y;p)/dp */
int HarmonicOscillator_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jacp_) return 2*3;
	double k=p_[0];
	double c=p_[1];
	double F=p_[2];
	double v=y_[0];
	double y=y_[1];
/* column 1 (df/dp_0) */
	jacp_[0] = -y; /* [0, 0] */
	jacp_[3] = 0; /* [1, 0] */
/* column 2 (df/dp_1) */
	jacp_[1] = -v; /* [0, 1] */
	jacp_[4] = 0; /* [1, 1] */
/* column 3 (df/dp_2) */
	jacp_[2] = 1; /* [0, 2] */
	jacp_[5] = 0; /* [1, 2] */
	return GSL_SUCCESS;
}
/* ode Functions F(t,y;p) */
int HarmonicOscillator_func(double t, const double y_[], double *func_, void *par)
{
	double *p_=par;
	if (!y_ || !func_) return 1;
	double k=p_[0];
	double c=p_[1];
	double F=p_[2];
	double v=y_[0];
	double y=y_[1];
	func_[0] = sqrt(y*y+v*v); /* norm */
	return GSL_SUCCESS;
}
/* Function Jacobian dF(t,y;p)/dy */
int HarmonicOscillator_funcJac(double t, const double y_[], double *funcJac_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJac_) return 2;
	double k=p_[0];
	double c=p_[1];
	double F=p_[2];
	double v=y_[0];
	double y=y_[1];
/* column 1 (dF/dy_0) */
	funcJac_[0] = v/sqrt(gsl_pow_2(y)+gsl_pow_2(v)); /* [0, 0] */
/* column 2 (dF/dy_1) */
	funcJac_[1] = y/sqrt(gsl_pow_2(y)+gsl_pow_2(v)); /* [0, 1] */
	return GSL_SUCCESS;
}
/* Function parameter Jacobian dF(t,y;p)/dp */
int HarmonicOscillator_funcJacp(double t, const double y_[], double *funcJacp_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJacp_) return 3;
	double k=p_[0];
	double c=p_[1];
	double F=p_[2];
	double v=y_[0];
	double y=y_[1];
/* column 1 (dF/dp_0) */
	funcJacp_[0] = 0; /* [0, 0] */
/* column 2 (dF/dp_1) */
	funcJacp_[1] = 0; /* [0, 1] */
	return GSL_SUCCESS;
}
/* ode default parameters */
int HarmonicOscillator_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 3;
	p_[0] = 0.5;
	p_[1] = 0.0;
	p_[2] = 0.0;
	return GSL_SUCCESS;
}
/* ode initial values */
int HarmonicOscillator_init(double t, double *y_, void *par)
{
	double *p_=par;
	if (!y_) return 2;
	double k=p_[0];
	double c=p_[1];
	double F=p_[2];
	/* the initial value of y may depend on the parameters. */
	y_[0] = 0.0;
	y_[1] = 1.0;
	return GSL_SUCCESS;
}
