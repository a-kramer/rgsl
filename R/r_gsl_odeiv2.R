#' Selects the integration method, by name
#'
#' These are the methods in the gsl library (documented in the
#' official documentation), but in reverse order, as they are
#' approximately ordered by complexity, with more complex methods
#' usually being better.
#'
#' It is therefore a reasonable approach to try methods from the more
#' complex end of the list first and try the next method if the
#' solutions are too slow. But ensure checking the accuracy/stability of the
#' result. The mapping between method names and keys:
#'
#' ```
#'      msbdf:	 0
#'    msadams:	 1
#'      bsimp:	 2
#'     rk4imp:	 3
#'     rk2imp:	 4
#'     rk1imp:	 5
#'      rk8pd:	 6
#'       rkck:	 7
#'      rkf45:	 8
#'        rk4:	 9
#'        rk2:	10
#' ```
#'
#' The returned value must be passed to the solver function calls, e.g.: r_gsl_odeiv2_outer(...,method = integrationMaethod("msbdf"))
#' @export
#' @param charMethod a scalar character string, defaulting to msbdf
#' @return an integer offset, which will select the named method in the solvers
integrationMethod <- function(charMethod=c("msbdf","msadams","bsimp","rk4imp","rk2imp","rk1imp","rk8pd","rkck","rkf45","rk4","rk2")){
	return(
		switch(charMethod[1],
			msbdf=0,
			msadams=1,
			bsimp=2,
			rk4imp=3,
			rk2imp=4,
			rk1imp=5,
			rk8pd=6,
			rkck=7,
			rkf45=8,
			rk4=9,
			rk2=10,
			0
		)
	)
}

#' Reverse look-up of method name from key
#'
#' These are the methods in the gsl library (documented in the
#' official documentation), but in reverse order, as they are
#' approximately ordered by complexity, with more complex methods
#' usually being better (but slower).
#'
#' It is therefore a reasonable approach to try methods from the more
#' complex end of the list first and try the next method if the
#' solutions are too slow. But we need to check the accuracy/stability of the
#' result. The mapping between method names and keys:
#'
#' ```
#'      msbdf:	 0
#'    msadams:	 1
#'      bsimp:	 2
#'     rk4imp:	 3
#'     rk2imp:	 4
#'     rk1imp:	 5
#'      rk8pd:	 6
#'       rkck:	 7
#'      rkf45:	 8
#'        rk4:	 9
#'        rk2:	10
#' ```
#'
#' The returned value is an integer index.
#' @export
#' @param key an integer from 0 to 10 (this is used as an offset in c, for 11 items)
#' @return a string representation of the integration method.
nameMethod <- function(key){
	charMethod <- c("msbdf","msadams","bsimp","rk4imp","rk2imp","rk1imp","rk8pd","rkck","rkf45","rk4","rk2")
	k <- round(key+1)
	if (all(0 < k) && all(k <= length(charMethod))){
		return(charMethod[k])
	} else {
		warning(sprintf("invalid key %i",key[k==0 | k>length(charMethod)]))
		return("")
	}
}


#' Initial Value Problem solution in C
#'
#' This is a wrapper. It uses the .Call function to call the C
#' function `r_gsl_odeiv2_outer()`. The C program solves a set of ODE intial
#' value problems and returns the trajectory `y(t;p)` for every real
#' valued parameter vector `p`. We use the solvers from the GNU
#' Scientific Library module odeiv2.
#'
#' The model will be simulated N×M times, where N is the number of
#' simulation experiments and M the number of columns in the p
#' matrix. The output is an N-sized list of 3 dimensional arrays y,
#' with y[i,j,k] corresponding to state variable i, time t[j], and
#' parameter set p[,k].
#'
#' The model must have output functions (observable quantities of a
#' system), the C-source can include a function with a name ending in
#' "_func" (e.g. mymodel.so with
#'
#' ```
#' int mymodel_func(double t, double y[], double fout[], void *p)
#' ```
#'
#' where fout will contain the output-function values after the call).
#'
#' @param name the name of the ODE model to simulate (a shared library of the same name will be dynamically loaded and needs to be created first)
#' @param experiments a list of `N` simulation experiments (time, parameters, initial value, events)
#' @param p a matrix of parameters with M columns
#' @param abs.tol absolute tolerance, real scalar
#' @param rel.tol relative tolerance, real scalar
#' @param initial.step.size initial value for the step size; the step size will adapt to a value that observes the tolerances, real scalar
#' @return a list of the solution trajectories `y(t;p)` for all experiments (named like the experiments), as well as the output functions
#' @keywords ODE
#' @useDynLib rgsl, odeiv_outer_f=r_gsl_odeiv2_outer_func
#' @export
r_gsl_odeiv2_outer <- function(name,experiments,p,abs.tol=1e-6,rel.tol=1e-5,initial.step.size=1e-3,method=0){
	if (is.character(comment(name))){
		so <- comment(name)
	} else {
		so <- paste0("./",name,".so")
		comment(name)<-so
	}
	if (!file.exists(so)){
            warning(sprintf("[r_gsl_odeiv2_outer] for model name «%s», in directory «%s» file «%s» not found.",name,getwd(),so))
	}
	if (!is.matrix(p)) p <- as.matrix(p)
	y <- .Call(odeiv_outer_f,name,experiments,p,abs.tol,rel.tol,initial.step.size,method)
	for (i in seq_along(experiments)){
		if ("initialState" %in% names(experiments[[i]])){
			dimnames(y[[i]]$state) <- list(names(experiments[[i]]$initialState),NULL,NULL)
		}
		if ("outputValues" %in% names(experiments[[i]])){
			dimnames(y[[i]]$func) <- list(names(experiments[[i]]$outputValues),NULL,NULL)
		}
	}
	return(y)
}

#' Initial Value Problem solution in C with sensitivity approximation
#'
#' This is a wrapper. It uses the .Call function to call the C
#' function r_gsl_odeiv2_outer(). The C program solves a set of ODE
#' intial value problems and returns the trajectory y(t;p) for every
#' real valued parameter vector p (column of a matrix). We use the
#' solvers from the GNU Scientific Library module odeiv2. This
#' function is similar to `r_gsl_odeiv2_outer`, but also estimates the
#' sensitivity of the solution: dx/dp (i.e.: how would the solution
#' change with slightly different parameters). For all functions in
#' this package, it is important to understand that we aim to solve
#' the model for several parameter vectors in one go. For this reason,
#' `p` is a matrix, even though each simulation only needs a vector
#' `p[,l]` (wich 1<l<=M). The returned value represents the outer
#' product of simulation experiments and columns of `p`: (N * M)
#' solutions.
#'
#' The returned list also contains the items "stateSensitivity" and
#' "funcSensitivity". Calculating sensitivities takes much more time
#' than evaluating the output functions, so this version of the solver
#' treats output functions as mandatory and there is no 'state_only'
#' version of this.
#'
#' A test for the quality of the sensitivity approximation is
#' calculating the relative error $E$ of a small shift in parameters
#' $\delta$:
#'
#' $$\frac{\| y(t;p+\delta) - (y(t;p) + S(t;p) \delta) \|}{\| y(t;p) \|} = E(t;p)\,,$$
#'
#' where $\delta$ is a small change in parameters (e.g. a Gaussian
#' random vector with $\mu=0$). $E$ can be calculated with only one
#' additional simulation; whenever $E(t;p)$ is small, then the
#' sensitivity approximation has worked.
#'
#' @param model_name the name of the ODE model to simulate (a shared
#'     library of the same name will be dynamically loaded and needs
#'     to be created first)
#' @param experiments a list of N simulation experiments (time,
#'     parameters, initial value, events)
#' @param p a matrix of parameters with M columns, possibly 1
#' @param abs.tol absolute tolerance, scalar of type double (numeric)
#' @param rel.tol relative tolerance, scalar of type double (numeric)
#' @param initial.step.size initial value for the step size; the step
#'     size will _adapt_ to a value that observes the tolerances; real
#'     scalar
#' @return a list containing the solution trajectories y(t;p) for all
#'     experiments (named like the experiments), as well as the output
#'     functions, and the sensitivities of both. The sensitivities are
#'     returned as a list of 3d-arrays, where the list index
#'     enumerates the different columns of `p`, the sensitivities have
#'     the structure `stateSensitivity[[l]][i,j,k]`, where i is for
#'     the state-variables/functions, j is for the parameter index in
#'     $dy_i/dx_j$, and k is for time; and finally `l` is for the
#'     different columns of `p`.
#' @keywords ODE
#' @useDynLib rgsl, odeiv_outer_sens=r_gsl_odeiv2_outer_sens
#' @export
#' @examples
#' y0 <- c(0,1)
#' t <- seq(0,1,length.out=100)
#' u <- c(0,0)
#' e <- list(time=t,input=u,initial_value=y0)
#' y <- r_gsl_odeiv2_outer("HarmonicOscillator",t,y0,p=matrix(seq(0,1,length.out=3),ncol=3))
r_gsl_odeiv2_outer_sens <- function(name,experiments,p,abs.tol=1e-6,rel.tol=1e-5,initial.step.size=1e-3,method=0){
	if (is.character(comment(name))){
		so <- comment(name)
	} else {
		so <- paste0("./",name,".so")
		comment(name)<-so
	}
	stopifnot(file.exists(so))
	if (!is.matrix(p)) p <- as.matrix(p)
	y <- .Call(odeiv_outer_sens,name,experiments,p,abs.tol,rel.tol,initial.step.size,method)
	for (i in seq_along(experiments)){
		if ("initialState" %in% names(experiments[[i]])){
			dimnames(y[[i]]$state) <- list(names(experiments[[i]]$initialState),NULL,NULL)
			for (j in seq_along(y[[i]]$stateSensitivity)){
			 dimnames(y[[i]]$stateSensitivity[[j]]) <- list(names(experiments[[i]]$initialState),NULL,NULL)
			}
		}
		if ("outputValues" %in% names(experiments[[i]])){
			dimnames(y[[i]]$func) <- list(names(experiments[[i]]$outputValues),NULL,NULL)
			for (j in seq_along(y[[i]]$funcSensitivity)){
				dimnames(y[[i]]$funcSensitivity[[j]]) <- list(names(experiments[[i]]$outputValues),NULL,NULL)
			}
		}
	}
	return(y)
}

#' Initial Value Problem solution in C, no functions
#'
#' This is a wrapper. It uses the .Call function to call the C
#' implementation. The C program solves a set of ODE intial value
#' problems and returns the trajectory y(t;p) for every real valued
#' parameter vector p. We use the solvers from the GNU Scientific
#' Library module odeiv2. This is similar to `r_gsl_odeiv2_outer` but
#' does not return (or evaluate) the system's output function
#' (MODEL_func()).
#'
#' Details: the c-code can be found in the `src/r_gsl_odeiv2.c`, the C
#' function is called `r_gsl_odeiv2_outer_state_only()`.
#'
#' @param name the name of the ODE model to simulate (a shared library of the same name will be dynamically loaded and needs to be created first)
#' @param experiments a list of N simulation experiments (time, parameters, initial value, events)
#' @param p a matrix of parameters with M columns
#' @param abs.tol absolute tolerance (one real number, defaults to 1e-6)
#' @param rel.tol relative tolerance (one real number, defaults to 1e-5)
#' @param initial.step.size initial value for the step size; the step size will adapt to a value that observes the tolerances; defaults to 1e-3
#' @return a list with solution trajectories y(t;p) for all experiments
#' @keywords ODE
#' @useDynLib rgsl, odeiv_outer_state=r_gsl_odeiv2_outer_state_only
#' @export
r_gsl_odeiv2_outer_state_only <- function(name,experiments,p,abs.tol=1e-6,rel.tol=1e-5,initial.step.size=1e-3,method=0){
	if (is.character(comment(name))){
		so <- comment(name)
	} else {
		so <- paste0("./",name,".so")
		comment(name)<-so
	}
	stopifnot(file.exists(so))
	if (!is.matrix(p)) p <- as.matrix(p)
	y <- .Call(odeiv_outer_state,name,experiments,p,abs.tol,rel.tol,initial.step.size,method)
	for (i in seq_along(experiments)){
		if ("initialState" %in% names(experiments[[i]])){
			dimnames(y[[i]]$state) <- list(names(experiments[[i]]$initialState),NULL,NULL)
		}
	}
	return(y)
}

