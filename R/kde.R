check_arg <- function(x,d){
	if (!is.matrix(x)) {
		z <- matrix(x,nrow=d)
	} else if (nrow(x)!=d && ncol(x)==d) {
		z <- t(x)
	} else {
		z <- x
	}
	if (!(nrow(as.matrix(z)) == d)){
		stop(sprintf("Sizes of argument and sample do not match: sample is «%i» dimensional (number of variables), argument x has «%i» variables.",nrow(Z),nrow(as.matrix(x))))
	}
	return(x)
}


#' kernel density estimator
#'
#' written in C, the prupose of this function is to make very fast
#' kernel density estimates, with no user supplied bandwidth
#' parameters.
#'
#' @useDynLib rgsl, density_c=dkde
#' @param sample a matrix with columns of observations of continuous
#'     variables from a target density.
#' @return A closure, function of x, where x is the matrix of columns
#'     for which to estimate densities
kde.c <- function(sample,normalize=FALSE){
	d <- dim(as.matrix(sample))
	if (d[1]>d[2]){
		Z=t(sample)
	} else {
		Z=sample
	}
	if (normalize){
		C<-chol(cov(t(Z)))
		Z<-solve(C,Z)
	}
	sigma <- estimate.sigma(Z)
	print(sigma)
	if (normalize) {
		pd <- function(x){
			z <- solve(C,x)
			y <- .Call(density_c,Z,z,sigma)
			return(y)
		}
	} else {
		pd <- function(x){
			z <- as.matrix(x)
			y <- .Call(density_c,Z,z,sigma)
			return(y)
		}
	}
	return(pd)
}

#' estimates the sigma parameter of the kernel
#'
#' this function will try to estimate a scalar sigma, assuming that
#' the sample variables are all very similar to one another.
#'
#' @useDynLib rgsl estimate_sigma
#' @param sample a sample matrix from the target density
estimate.sigma <- function(sample){
	d <- dim(as.matrix(sample))
	if (d[1]>d[2]){
		Z=t(sample)
	} else {
		Z=sample
	}
	s <- .Call(estimate_sigma,sample)
	return(s)
}
