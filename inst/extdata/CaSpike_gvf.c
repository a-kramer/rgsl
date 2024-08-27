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
int CaSpike_vf(double t, const double y_[], double f_[], void *par)
{
	double *p_=par;
	if (!y_ || !f_) return 1;
	double tau=p_[0];
	double CaBase=p_[1];
	double dCa=p_[2];
	double Ca=y_[0];
	f_[0] = CaBase/tau - Ca/tau;
	return GSL_SUCCESS;
}
/* Scheduled Event function,
   EventLabel specifies which of the possible transformations to apply,
   dose can specify a scalar intensity for this transformation. */
int CaSpike_event(double t, double y_[], void *par, int EventLabel, double dose)
{
	double *p_=par;
	if (!y_ || !par || EventLabel<0) return 1;
	enum eventLabel { Spike, numEvents }; /* event name indexes */
	enum stateVariable { var_Ca, numStateVar }; /* state variable indexes  */
	enum param { par_tau,par_CaBase,par_dCa, numParam }; /* parameter indexes  */
	double tau=p_[0];
	double CaBase=p_[1];
	double dCa=p_[2];
	double Ca=y_[0];
	switch(EventLabel){
	case Spike:
		y_[var_Ca] = Ca+dCa; /* state variable transformation */
	break;
	}
	return GSL_SUCCESS;
}

/* ode Jacobian df(t,y;p)/dy */
int CaSpike_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jac_) return 1*1;
	double tau=p_[0];
	double CaBase=p_[1];
	double dCa=p_[2];
	double Ca=y_[0];
/* column 1 (df/dy_0) */
	jac_[0] = -(1/tau); /* [0, 0] */
	return GSL_SUCCESS;
}
/* ode parameter Jacobian df(t,y;p)/dp */
int CaSpike_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jacp_) return 1*3;
	double tau=p_[0];
	double CaBase=p_[1];
	double dCa=p_[2];
	double Ca=y_[0];
/* column 1 (df/dp_0) */
	jacp_[0] = Ca/gsl_pow_2(tau)-CaBase/gsl_pow_2(tau); /* [0, 0] */
/* column 2 (df/dp_1) */
	jacp_[1] = 1/tau; /* [0, 1] */
/* column 3 (df/dp_2) */
	jacp_[2] = 0; /* [0, 2] */
	return GSL_SUCCESS;
}
/* ode Functions F(t,y;p) */
int CaSpike_func(double t, const double y_[], double *func_, void *par)
{
	double *p_=par;
	if (!y_ || !func_) return 1;
	double tau=p_[0];
	double CaBase=p_[1];
	double dCa=p_[2];
	double Ca=y_[0];
	func_[0] = Ca; /* Ca */
	return GSL_SUCCESS;
}
/* Function Jacobian dF(t,y;p)/dy */
int CaSpike_funcJac(double t, const double y_[], double *funcJac_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJac_) return 1;
	double tau=p_[0];
	double CaBase=p_[1];
	double dCa=p_[2];
	double Ca=y_[0];
/* column 1 (dF/dy_0) */
	funcJac_[0] = 1; /* [0, 0] */
	return GSL_SUCCESS;
}
/* Function parameter Jacobian dF(t,y;p)/dp */
int CaSpike_funcJacp(double t, const double y_[], double *funcJacp_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJacp_) return 3;
	double tau=p_[0];
	double CaBase=p_[1];
	double dCa=p_[2];
	double Ca=y_[0];
/* column 1 (dF/dp_0) */
	funcJacp_[0] = 0; /* [0, 0] */
	return GSL_SUCCESS;
}
/* ode default parameters */
int CaSpike_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 3;
	p_[0] = 0.3;
	p_[1] = 0.1;
	p_[2] = 1;
	return GSL_SUCCESS;
}
/* ode initial values */
int CaSpike_init(double t, double *y_, void *par)
{
	double *p_=par;
	if (!y_) return 1;
	double tau=p_[0];
	double CaBase=p_[1];
	double dCa=p_[2];
	/* the initial value of y may depend on the parameters. */
	y_[0] = CaBase;
	return GSL_SUCCESS;
}
