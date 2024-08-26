#' Initial Value Problem solution in C
#'
#' This is a wrapper. It uses the .Call function to call the C
#' function r_gsl_odeiv2_outer(). The C program solves a set of ODE intial
#' value problems and returns the trajectory y(t;p) for every real
#' valued parameter vector p. We use the solvers from the GNU
#' Scientific Library module odeiv2.
#'
#' The model will be simulated N×M times, where N is the number of
#' simulation experiments and M the number of columns in the p
#' matrix. The output is an N-sized list of 3 dimensional arrays y,
#' with y[i,j,k] corresponding to state variable i, time t[j], and
#' parameter set p[,k].
#'
#' The model may have output functions that model observable
#' quantities of a system, the c-source can include a function with a
#' name ending in "_func" (e.g. mymodel.so with
#'
#' ```
#' int mymodel_func(double t, double y[], double fout[], void *p)
#' ```
#'
#' where fout will contain the output-function values after the call).
#'
#' @param model_name the name of the ODE model to simulate (a shared library of the same name will be dynamically loaded and needs to be created first)
#' @param experiments a list of N simulation experiments (time, parameters, initial value, events)
#' @param p a matrix of parameters with M columns
#' @param abs.tol absolute tolerance (one real number, defaults to 1e-6)
#' @param rel.tol relative tolerance (one real number, defaults to 1e-5)
#' @param initial.step.size initial value for the step size; the step size will adapt to a value that observes the tolerances; defaults to 1e-3
#' @return the solution trajectories y(t;p) for all experiments, as well as the output functions (if MODEL_func() is present in the .so file)
#' @keywords ODE
#' @useDynLib rgsl, odeiv_outer_f=r_gsl_odeiv2_outer_func
#' @export
#' @examples
#' y0 <- c(0,1)
#' t <- seq(0,1,length.out=100)
#' u <- c(0,0)
#' e <- list(time=t,input=u,initial_value=y0)
#' y <- r_gsl_odeiv2_outer("HarmonicOscillator",t,y0,p=matrix(seq(0,1,length.out=3),ncol=3))
r_gsl_odeiv2_outer <- function(name,experiments,p,abs.tol=1e-6,rel.tol=1e-5,initial.step.size=1e-3){
	if (is.character(comment(name))){
		so <- comment(name)
	} else {
		so <- paste0("./",name,".so")
		comment(name)<-so
	}
	stopifnot(file.exists(so))
	if (!is.matrix(p)) p <- as.matrix(p)
	y <- .Call(odeiv_outer_f,name,experiments,p,abs.tol,rel.tol,initial.step.size)
	return(y)
}

#' Initial Value Problem solution in C
#'
#' This is a wrapper. It uses the .Call function to call the C
#' function r_gsl_odeiv2_outer(). The C program solves a set of ODE
#' intial value problems and returns the trajectory y(t;p) for every
#' real valued parameter vector p. We use the solvers from the GNU
#' Scientific Library module odeiv2. This function is similar to
#' `r_gsl_odeiv2_outer`, but also estimates the sensitivity of the
#' solution: dx/dp (how would the solution change with slightly
#' different parameters).
#'
#' The returned list also contains the items "stateSensitivity" and
#' "funcSensitivity". Calculating sensitivities takes much more time
#' than evaluating the output functions, so this version of the solver
#' treats output functions as mandatory.
#'
#' @param model_name the name of the ODE model to simulate (a shared library of the same name will be dynamically loaded and needs to be created first)
#' @param experiments a list of N simulation experiments (time, parameters, initial value, events)
#' @param p a matrix of parameters with M columns
#' @param abs.tol absolute tolerance (one real number, defaults to 1e-6)
#' @param rel.tol relative tolerance (one real number, defaults to 1e-5)
#' @param initial.step.size initial value for the step size; the step size will adapt to a value that observes the tolerances; defaults to 1e-3
#' @return the solution trajectories y(t;p) for all experiments, as well as the output functions (if MODEL_func() is present in the .so file)
#' @keywords ODE
#' @useDynLib rgsl, odeiv_outer_sens=r_gsl_odeiv2_outer_sens
#' @export
#' @examples
#' y0 <- c(0,1)
#' t <- seq(0,1,length.out=100)
#' u <- c(0,0)
#' e <- list(time=t,input=u,initial_value=y0)
#' y <- r_gsl_odeiv2_outer("HarmonicOscillator",t,y0,p=matrix(seq(0,1,length.out=3),ncol=3))
r_gsl_odeiv2_outer_sens <- function(name,experiments,p,abs.tol=1e-6,rel.tol=1e-5,initial.step.size=1e-3){
	if (is.character(comment(name))){
		so <- comment(name)
	} else {
		so <- paste0("./",name,".so")
		comment(name)<-so
	}
	stopifnot(file.exists(so))
	if (!is.matrix(p)) p <- as.matrix(p)
	y <- .Call(odeiv_outer_sens,name,experiments,p,abs.tol,rel.tol,initial.step.size)
	return(y)
}

#' Initial Value Problem solution in C, no functions
#'
#' This is a wrapper. It uses the .Call function to call the C
#' function r_gsl_odeiv2_outer(). The C program solves a set of ODE
#' intial value problems and returns the trajectory y(t;p) for every
#' real valued parameter vector p. We use the solvers from the GNU
#' Scientific Library module odeiv2. This is similar to
#' `r_gsl_odeiv2_outer` but does not return (or evaluate) the system's
#' output function (model_func()).
#'
#'
#' @param model_name the name of the ODE model to simulate (a shared library of the same name will be dynamically loaded and needs to be created first)
#' @param experiments a list of N simulation experiments (time, parameters, initial value, events)
#' @param p a matrix of parameters with M columns
#' @param abs.tol absolute tolerance (one real number, defaults to 1e-6)
#' @param rel.tol relative tolerance (one real number, defaults to 1e-5)
#' @param initial.step.size initial value for the step size; the step size will adapt to a value that observes the tolerances; defaults to 1e-3
#' @return the solution trajectories y(t;p) for all experiments, as well as the output functions (if MODEL_func() is present in the .so file)
#' @keywords ODE
#' @useDynLib rgsl, odeiv_outer_state=r_gsl_odeiv2_outer
#' @export
#' @examples
#' y0 <- c(0,1)
#' t <- seq(0,1,length.out=100)
#' u <- c(0,0)
#' e <- list(time=t,input=u,initial_value=y0)
#' y <- r_gsl_odeiv2_outer("HarmonicOscillator",t,y0,p=matrix(seq(0,1,length.out=3),ncol=3))
r_gsl_odeiv2_outer_state_only <- function(name,experiments,p,abs.tol=1e-6,rel.tol=1e-5,initial.step.size=1e-3){
	if (is.character(comment(name))){
		so <- comment(name)
	} else {
		so <- paste0("./",name,".so")
		comment(name)<-so
	}
	stopifnot(file.exists(so))
	if (!is.matrix(p)) p <- as.matrix(p)
	y <- .Call(odeiv_outer_state,name,experiments,p,abs.tol,rel.tol,initial.step.size)
	return(y)
}

