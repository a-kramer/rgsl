# rgsl - an interface between R and gsl_odeiv2

This R package solves a series of initial value problems given as an
[ordinary differential equation](MODELS.md), an initial state and a numerical
parameter vector. 

This project was [funded by the EU](./ACKNOWLEDGMENTS.md) as part of the [Human Brain Project](https://www.humanbrainproject.eu/en/) and supported by [EBRAINS](https://ebrains.eu/) infrastructure.

The 2 in `odeiv2` is from the GSL, there is no older R package of this
name.

To [install](INSTALL.md):

```R
remotes::install_github("a-kramer/rgsl")
```

## Description and Purpose

The C code calls
[gsl_odeiv2](https://www.gnu.org/software/gsl/doc/html/ode-initval.html)
module functions to solve the problem or problems. The goal is to
offload as much work as possible to the C code and keep the overhead
minimal. That is why this package expects to solve a set of problems,
rather than one, for the same model file (varying in e.g.: initial
conditions, or parameters).

The package contains 3 interface functions, they accept different ways
of defining a set of problems,each with their own drawbacks and
advantages. The interface functions are described in the following
Sections.

The ODE has to exist as a [shared library](MODELS.md) (`.so`) file (currently in
the current working directory: `?setwd` and `?getwd`). There are some
assumptions we make about the contents of the shared library file.

Here we assume that the solutions serve some scientific purpose and
the lab experiments come with _observables_, some measureable values
that depend on the system's state (but are not the full state
vector). We call the part of the model that calculates the observables
`${ModelName}_func()` (vfgen also calls them *Functions* of the model). 

### Note

This package cannot use an ODE written in R, which the reader
may consider impractical. For plain R models, there is of course
[deSolve](https://cran.r-project.org/web/packages/deSolve/) that
remains an alternative. Our goal is to auto-generate the shared
library whenever possible and avoid writing the model files by
hand. This package does not do this in any way, but
[vfgen](https://warrenweckesser.github.io/vfgen/) does. But, any
method of auto-generating code will do, as long as it creates somewhat
compatible c files.

## Interface function `r_gsl_odeiv2`

The [purpose](PURPOSE.md) of this interface is very similar to the
others, and is easiest to use for simple problems.

Usage:

```R
library("rgsl")

ny <- 3
np <- 4
nt <- 100
m <- 20

t <- seq(0,1,length.out=nt)
y0 <- matrix(c(1,2,3),nrow=ny,ncol=m)
p <- matrix( ... ,nrow=np,ncol=m)
y <- r_gsl_odeiv2("ModelName",t,y0,p)
```

The function `r_gsl_odeiv2()` will look for the file `ModelName.so` in
the current directory.

The parameters can be an _n×m_ matrix of several parameterisations
(each column is a distinct parameter vector of length _n_). 

The solver will be looped _m_ times (using each column) within the C
function of the same name. If the parameters _p_ are a matrix then
`y0` should be as well (initial values), with identical numbers of
columns.

OpenMP is used to make this set of simulations run in parallel, if
possible.

The result is returned to R as a 3-dimensional array of size _ny×nt×m_.

The simulations can be more intricate if the [event system](EVENTS.md) is used.

### Usage Notes

The time `t` is a vector of output time-points, with `t[1]`
corresponding to `y0` (the initial conditions). The different
simulation scannot differ in the amount or values of the time
points.

The parameters `p` must be numeric (not some arbitrary content passed
to the ODE right-hand-side functions).

Both `y0` and `p` can be matrices; if they are, more than one
simulation is performed, but keep in mind:
```R
stopifnot(ncol(y0)==ncol(p))
```

## List of Experiments: `r_gsl_odeiv2_sim`

This interface assumes that the user has a list of *N* simulation
experiments and uses this call structure:

```R
library("rgsl")
y<-r_gsl_odeiv2_sim(NameOfModel,experiments)
```

Each list item is a list of named properties that define a simulation
run, e.g.:

```R
experiments[[1]][["time"]] <- seq(0,1,length.out=100)
experiments[[1]][["initial_value"]] <- c(0,0,0)
experiments[[1]][["parameters"]] <- c(1,2,3,4,5)
```
The entries are found by their name, so they naming is not arbitrary:

| name | optional | meaning |
|-----:|:-------: |:--------|
| `time` | no | output time |
| `initial_value` | no | a state space vector `y(t=t0)` |
| `parameters` | no | a numeric vector of parameters |
| `events` | yes | an event description for this run |

The simulation experiments will be done in parallel.  This interface
is useful whenever your problem lends itself to this description.

Particularly, if you want to simulate a subset of experiments, you can
remove unwanted items from the list easily (which is not as easy with
plain `r_gsl_odeiv2()`, see above).

## Interface Function `r_gsl_odeiv2_outer`

This interface is a mixture of the first two. We assume that the user
has a set of *M* parameter vectors to test (e.g. for model fitting or
hypothesis testing or whatever), but otherwise there is a specific
list of *N* simulation instructions (or _protocols_, or
_scenarios_). This interface will perform the outer product of *N×M*
simulations:

```R
y<-r_gsl_odeiv2_sim(NameOfModel, experiments, parameters)
```

The `parameters` are a matrix, each column a possible model
parameterization. The list of `experiments` may contain parameters,
but they will be ignored, this matrix will be used.

This interface has the advantage that the number of parameter columns
does not have to correspond to the length of the experiments list. So,
once again the list of experiments can be truncated to any subset of
scenarios.

## SBtabVFGEN

We are working on models in systems biology and define models using
the SBtab format. Therefore, this project integrates somwhat well with
[SBtabVFGEN](https://github.com/a-kramer/SBtabVFGEN): it contains
`sbtab_to_vfgen`, which prints vfgen files (given SBtab). 

SBtab is a format within the realm of systems biology (but less common
than SBML) and is suited for writing/drafting biological models. The
workflow could be:


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
