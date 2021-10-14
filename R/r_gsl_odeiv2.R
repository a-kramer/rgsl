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
r_gsl_odeiv2 <- function(name,t,y0,p){
    so <- paste0(name,".so")
    stopifnot(file.exists(so))
    #rgsl <- dyn.load("r_gsl_odeiv2.so")
    y <- .Call(odeiv,name,t,y0,as.matrix(p))
    return(y)
}
