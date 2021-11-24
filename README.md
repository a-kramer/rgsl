# rgsl - an interface between R and gsl_odeiv2

This R package solves a series of initial value problems given as an
ordinary differential equation, an initial state and a numerical
parameter vector. The C code calls [gsl odeiv2](https://www.gnu.org/software/gsl/doc/html/ode-initval.html)
module functions to solve the problem or problems.

The parameter vector can be replaced by an _n×m_ matrix of several
parameterisations (each column is a distinct parameter vector). The
solver will be called _m_ times. If the parameters _p_ are a matrix
then so should be the initial value _y0_, with identical numbers of
columns.

OpenMP is used to make this set of simulations run in parallel.

The result is returned to R as a 3-dimensional array.

This project was [funded by the EU](./ACKNOWLEDGMENTS.md) as part of the [Human Brain Project](https://www.humanbrainproject.eu/en/) and supported by [EBRAINS](https://ebrains.eu/) infrastructure.

The 2 in `odeiv2` is from the GSL, there is no older R package of this
name.

## Install

Using the `remotes` package:

```R
remotes::install_github("a-kramer/rgsl")
```

Ensure that the [GNU Scientific Library](https://www.gnu.org/software/gsl/doc/html/index.html) (gsl) is installed in your system and pkg-config can find it:

```bash
$ pkg-config --libs gsl
```

The above command should print something similar to:
```
-lgsl -lgslcblas -lm
```

## Usage

The simplest call is:
```R
y<-r_gsl_odeiv2("NameOfModel",t,y0,p)
```

where `NameOfModel` indicates a file `NameOfModel.so` in the current
directory, `t` is the vector of output time-points, with `t[1]`
corresponding to `y0` (the initial conditions). The parameters `p`
must be numeric (not arbitrary content).

Both `y0` and `p` can be matrices, if they are more than one
simulation is performed.

The simulations can be more intricate if the event system is used, see
below.

### Events

Often, a simulation requires an instantaneous intervention at a known
point in time. For this purpose, a list of scheduled events can be fed
into the solver. This does not cover all needs and use cases, but this
concept can be used to mimic some experimental setups without calling
the R functions multiple times.

Each simulation can receive an event of this kind:
```R
te <- 3
event[[i]] <- list(time=te,A,a,B,b)
```
where at time `te` both the state and the parameters will be transformed, like this:
```R
y <- A %*% y + a
p <- B %*% p + b
```
but in C. Events can also be named, with names corresponding to the
column-names of `p`. In this case experiments with no events don't
require a blank entry in the events list. The events list contains
only meaningful (non-empty) events and they are accessed by name.

If `time` is a vector, then all transformation quantities can be 3
dimensional arrays, where the third dimension corresponds to time:

```R
# at time te[j], the transformation is
y <- A[,,j] %*% y + a[,1,j]
p <- B[,,j] %*% p + b[,1,j]
```
(the second dimension of `a` and `b` remains unused).

If all `A` and `B` are diagonal, then the second dimension of these
matrices can be 1 and understood as a list of diagonals:
```R
# if A and B are diagonals, then for all time[j]
y <- A[,1,j] * y + a[,1,j]
p <- B[,1,j] * p + b[,1,j]
```
but in C (again).

If `p` is a vector, then the list of events is _allowed_ to be just
one event, _not_ within a _one element list_ (but that is also ok).

This event system is less powerful than a triggered event system,
where an instantaneous event occurs when a condition is met.

The same effect could be achieved by calling `r_gsl_odeiv2` and
integrating up to the event time, perform the transformation in R and
then continue the integration from the transformed state using a
second call to `r_gsl_odeiv2`.

## Purpose

In this package we consider the initial value problem:

<span class="math display" style="width: 15em; margin: auto; padding: 2em"><em>ẏ</em> = <em>f</em>(<em>y</em>, <em>t</em>; <em>p</em>)</span>,<span style="width=4em"/><span class="math display" style="width: 15em; margin: auto; padding: 2em"><em>y</em>(<em>t</em><sub>0</sub>) = <em>y</em><sub>0</sub></span>

where _p_ is a parameter vector. 


We think that it is helpful to distinguish between _constants_
(parameters that never change), _unknown parameters_, and _input
parameters_.

|Type|Symbol|Description|
|---:|:----:|:----------|
|constants|NONE|never change, are built into the model |
|unknown| k | may be inferred from data |
|known| u | these distinguish different experiments |

The overall vector _p_ contains both the unknown _k_ and the known
_u_, e.g. as `p=c(k,u)` (aka a direct sum). But any ordering is possible (here _c_ means
concatenate, like in R).

To determine the parameters _k_, the model _f_ needs to be simulated
in the provided experimental scenarios and the parameters are found
through optimisation, or sampling (or perhaps another, similar
method). 

Typically, these methods need to solve the model many times. This is
especially true for Bayesian methods (sampling). So, solving a batch
of similar problems in one function call is desirable. In the next two
Sections we discuss two common Scenarios.

### MCMC sampling

Let us assume that we have obtained a data set _D_, that includes
measured observables form different experiments. Each experiment is
distinguished from the others by the values in the input to the model.

The input may well be dynamic, and if it is, the model needs to
contain the possible input dynamics within the model files. But, let
us assume here that we can switch between input dynamics by setting
specific input parameter values: *u₁, u₂, u₃, …*. We shall evaluate
the model for a given suggested parameter vector _k_. Then we can
construct a _p_ for each input: *p₁=c(k,u₁)*. Example R code:

```R
# two experiments, "Stim" and "Control"
u<-load_inputs_from_file("my_file_of_inputs.dat")
k<-c(1,0,2,3) # let's say this makes sense
p<-matrix(c(k,u$Stim,k,u$Control),ncol=2)
t=seq(0,10,length.out=12)
y0=c(1,0,0) # initial state vector

library(rgsl)
y<-r_gsl_odeiv2("NameOfModel",t,y0,p)

# the below plot will correspont to parameters c(k,u$Stim)
plot(t,y[1,,1])
# and this will be the ratio between Stim and Control:
lines(t,y[1,,1]/y[1,,2])
```

The result is a three dimensional array _y_, where `y[i,j,l]`
will be the state variable _i_ at time _j_ for parameter vector _l_.

### Sample Post-Processing

Similarly, if an MCMC sample _K_ for a given model already exists, the
user may want to simulate a prediction based on this sample. To do
this the whole sample may be provided to the C solver used here and
the result will be a batch of trajectories, one for each sampled
parameter vector. The sample matrix _K_ is a set of column-vectors,
each a valid parameterisation of the model. We want to predict the
system's behaviour for a given input of interest: _u_.

We construct a parameter matrix:

```R
K<-load_sample_from_file("sample_file.hdf5")
u<-c(1,2,3)
n<-dim(K)
P<-rbind(K,matrix(u,nrow=3,ncol=n(2)))
```

So, the input is now repeated for each column in _K_. Then we call the
`r_gsl_odeiv2()` interface function to simulate the scenario _u_ for
each parameter vector. The result is a matrix `y[,,l]` where `l` can
be any number in the range `1:n(2)`.


### Models 

This package requires the model to be provided as a shared library:
`${MODEL}.so`. The function names within the file need to correspond to the model name:

```bash
${MODEL}_vf();
${MODEL}_jac();
```

### Example File

See the example file [HarmonicOscillator](./HarmonicOscillator_gvf.c)
to inspect an auto-generated right-hand-side and Jacobian function
compatible with the solvers from the gsl.

The example script [HarmonicOscillator.R](./HarmonicOscillator.R)
contains a `demo()` function.

These C code files (for the gsl module
[odeiv2](https://www.gnu.org/software/gsl/doc/html/ode-initval.html))
can be automatically obtained using from
[VFGEN](https://github.com/WarrenWeckesser/vfgen) or written manually.

## SBtabVFGEN

This project integrates well with
[SBtabVFGEN](https://github.com/a-kramer/SBtabVFGEN) as this package
prints vfgen files from source files written using the SBtab
format. SBtab is a format in systems biology and is suited for
writing/drafting biological models. The workflow could be:


```
  User written:                           generated          Simulation
  +-----------+      +-----------+      +------------+      +----------+
  |           |      |           |      |  (CVODE)   |      |          |
  | SBtab (M) +--+-->+   VFGEN   +----->+  ODE code  +--+-->+   rgsl   |
  |           |  |   |           |      | +jacobian  |  |   |          |
  +-----------+  |   +-----------+      +------------+  |   +----------+
                 |                                      |
            +----+--------------+                 +----------------+
            |                   |                 |                |
            | sbtab_to_vfgen()  |                 | r_gsl_odeiv2() |
            |                   |                 |                |
            +-------------------+                 +----------------+

```
