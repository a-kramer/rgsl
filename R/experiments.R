#' Compile the ODE model if necessary
#'
#' Given the name of a C file, create a shared object in the current
#' directory, if it doesn't exist yet.
#' @param file.c a file with the model sources, ending in _gvf.c
#' @return the name of the shared object file
#' @export
model.so <- function(file.c){
	LIBS <- "-lgsl -lgslcblas -lm"
	CFLAGS <- "-shared -fPIC -Wfatal-errors -O2"
	so <- paste0("./",sub("_gvf[.]c",".so",basename(file.c)))
	if (!file.exists(so)){
		system2("cc",sprintf("%s -o %s '%s' %s",CFLAGS,so,file.c,LIBS))
	}
	stopifnot(file.exists(so))
	return(so)
}

#' create a repeating linear transformation
#'
#' a linear tranmsformation, with shift: y <- A %*% y + b is
#' characterized by the matrix A, the vector b and time at which the
#' transformation happens. If the same transformation is applied
#' several times, A and b need to be expanded into arrays with three
#' indices.
#' @param lt length of the transformation time vector
#' @param A transformation matrix
#' @param b shift vector
#' @return a list of two items, a transformation array A and a shift
#'     array b, both three dimensional: length(dim([Ab]))==3. A and b
#'     repeat lt times.
#' @export
affine.transform <- function(lt=1,A=1,b=0){
	n <- nrow(as.matrix(A))
	m <- nrow(as.matrix(b))
	A <- array(A,dim=c(n,n,lt))
	b <- array(b,dim=c(m,1,lt));
	return(list(A=A,b=b))
}

#' A scheduled transformation event
#'
#' This function creates a list of all necessary characteristics of a
#' repeating transformation event.
#'
#' A scheduled event happens at a predetermined time t. It affects
#' both model state and parameters.
#'
#' The transformations are affine (linear + shift). If no
#' transformation shall occur, an identity matrix can be used.
#' They are encoded as a series of matrices A and a series of
#' vectors b, both series have the same length as t. The third index
#' of A and b corresponds to the time: t[i], A[,,i], and b[,,i] belong
#' together.
#'
#' @param t the times at which transformations occur
#' @param state.tf the transformation that affects the state
#' @param param.tf the transformation that affects the parameters
#' @return a list with two components: "time" and "tf";
#' tf itself is a list of two components: "state" and "param";
#' both are transformations with the components A and b.
event.tf <- function(t,state.tf,param.tf){
	tf <- list(state=state.tf,param=param.tf)
	return(list(time=t,tf=tf))
}

#' Create a simulation experiment
#'
#' Creates a list of simulation instructions, with names compatible
#' with the solver. The transformations can be made using the function
#' affine.transform().
#'
#' @param t time vector
#' @param y0 initial value of the state variables
#' @param event.t times at which events happen
#' @param state.tf state transformation to be applied at times event.t
#' @param param.tf parameter changes to apply at times event.t
#' @return a list of simulation experiments
make.experiment <- function(t,y0,par,event.t=NULL,state.tf=NULL,param.tf=NULL){
	if (!is.null(event.t)){
		stopifnot(length(event.t) == dim(state.tf$A)[3])
		stopifnot(length(event.t) == dim(state.tf$b)[2])
		stopifnot(length(event.t) == dim(param.tf$A)[3])
		stopifnot(length(event.t) == dim(param.tf$b)[2])
		ev <- event.tf(event.t,state.tf,param.tf)
		return(list(time=t,parameters=par,initial_value=y0,events=ev))
	}
	return(list(time=t,parameters=par,initial_value=y0))
}
