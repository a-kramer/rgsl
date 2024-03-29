# Purpose

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

## MCMC sampling

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

## Sample Post-Processing

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

