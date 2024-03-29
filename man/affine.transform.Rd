% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experiments.R
\name{affine.transform}
\alias{affine.transform}
\title{create a repeating linear transformation}
\usage{
affine.transform(lt = 1, A = 1, b = 0)
}
\arguments{
\item{lt}{length of the transformation time vector}

\item{A}{transformation matrix}

\item{b}{shift vector}
}
\value{
a list of two items, a transformation array A and a shift
array b, both three dimensional: length(dim(\link{Ab}))==3. A and b
repeat lt times.
}
\description{
a linear tranmsformation, with shift: y <- A \%*\% y + b is
characterized by the matrix A, the vector b and time at which the
transformation happens. If the same transformation is applied
several times, A and b need to be expanded into arrays with three
indices.
}
