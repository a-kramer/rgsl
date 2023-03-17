#' Initial Value Problem solution in C
#'
#' This is a wrapper. It uses the .Call function to call the C
#' function r_gsl_odeiv2(). The C program solves a set of ODE intial
#' value problems and returns the trajectory y(t;p) for every real
#' valued parameter vector p. We use the solvers from the GNU
#' Scientific Library module odeiv2.
#'
#' The model will be simulated N times, where N is the number of
#' columns in p. The output is a 3 dimensional array y, with y[i,j,k]
#' corresponding to state variable i, time t[j], and parameter set
#' p[,k].
#'
#' @param model_name the name of the ODE model to simulate (a shared library of the same name will be dynamically loaded and needs to be created first)
#' @param t a vector of time values for the desired output
#' @param y0 initial value of the state vector, at time t[1]
#' @param p a matrix, each column is a valid parameter set for the model
#' @return the solution trajectories y(t;p) for all p[,k] (3-dim-array)
#' @keywords ODE
#' @useDynLib rgsl, odeiv=r_gsl_odeiv2
#' @export
#' @examples
#' y0 <- c(0,1)
#' t <- seq(0,1,length.out=100)
#' p <- c(1,0,0)
#' y <- r_gsl_odeiv2("HarmonicOscillator",t,y0,p)
r_gsl_odeiv2 <- function(name,t,y0,p,events=NULL){
	if (is.character(comment(name))){
		so <- comment(name)
	} else {
		so <- paste0("./",name,".so")
		comment(name)<-so
	}
	stopifnot(file.exists(so))
	if ("time" %in% names(events)) {
		stopifnot(is.vector(p))
		events <- list(events)
	}
	if (is.vector(y0)){
		dim(y0)<-c(length(y0),1)
	}
	if (is.vector(p)) {
		dim(p)<-c(length(p),1)
	}
	if (is.character(colnames(p)) && is.character(names(events))) {
		stopifnot(all(names(events) %in% colnames(p)))
	}
	y <- .Call(odeiv,name,t,y0,p,events)
	dimnames(y) <- list(rownames(y0),names(t),colnames(p))
	return(y)
}


#' Initial Value Problem solution in C
#'
#' This is a wrapper. It uses the .Call function to call the C
#' function r_gsl_odeiv2(). The C program solves a set of ODE intial
#' value problems and returns the trajectory y(t;p) for every real
#' valued parameter vector p. We use the solvers from the GNU
#' Scientific Library module odeiv2.
#'
#' The model will be simulated N times, where N is the number of
#' columns in p. The output is a 3 dimensional array y, with y[i,j,k]
#' corresponding to state variable i, time t[j], and parameter set
#' p[,k].
#'
#' @param model_name the name of the ODE model to simulate (a shared library of the same name will be dynamically loaded and needs to be created first)
#' @param experiments a list of simulation experiments (time, parameters, initial value, events)
#' @return the solution trajectories y(t;p) for all experiments
#' @keywords ODE
#' @useDynLib rgsl, simulate=r_gsl_odeiv2_simulate
#' @export
#' @examples
#' y0 <- c(0,1)
#' t <- seq(0,1,length.out=100)
#' p <- c(1,0,0)
#' e <- list(time=t,parameters=p,initial_value=y0)
#' y <- r_gsl_odeiv2("HarmonicOscillator",t,y0,p)
r_gsl_odeiv2_sim <- function(name,experiments){
	if (is.character(comment(name))){
		so <- comment(name)
	} else {
		so <- paste0("./",name,".so")
		comment(name)<-so
	}
	stopifnot(file.exists(so))
	y <- .Call(simulate,name,experiments)
	return(y)
}


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
#' @return the solution trajectories y(t;p) for all experiments, as well as the output functions (if MODEL_func() is present in the .so file)
#' @keywords ODE
#' @useDynLib rgsl, odeiv_outer_e=r_gsl_odeiv2_outer
#' @export
#' @examples
#' y0 <- c(0,1)
#' t <- seq(0,1,length.out=100)
#' u <- c(0,0)
#' e <- list(time=t,input=u,initial_value=y0)
#' y <- r_gsl_odeiv2_outer("HarmonicOscillator",t,y0,p=matrix(seq(0,1,length.out=3),ncol=3))
r_gsl_odeiv2_outer <- function(name,experiments,p){
	if (is.character(comment(name))){
		so <- comment(name)
	} else {
		so <- paste0("./",name,".so")
		comment(name)<-so
	}
	stopifnot(file.exists(so))
	stopifnot(is.matrix(p))
	stopifnot(any(c('outputTimes','time') %in% names(experiments[[1]])))
	stopifnot(is.double(experiments[[1]]$outputTimes) || is.double(experiments[[1]]$time))
	y <- .Call(odeiv_outer_e,name,experiments,p)
	return(y)
}

