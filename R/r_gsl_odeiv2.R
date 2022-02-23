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
    so <- paste0(name,".so")
    stopifnot(file.exists(so))
    if ("time" %in% names(events)) {
        stopifnot(is.vector(p))
        events=list(events)
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
    so <- paste0(name,".so")
    stopifnot(file.exists(so))
    y <- .Call(simulate,name,experiments)
    return(y)
}

