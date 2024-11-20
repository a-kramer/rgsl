# rgsl - an interface between R and gsl_odeiv2

This R package solves a series of initial value problems given as an
[ordinary differential equation](MODELS.md), an initial state and a numerical
parameter vector.

This project was [funded by the EU](./ACKNOWLEDGMENTS.md) as part of the [Human Brain Project](https://www.humanbrainproject.eu/en/) and supported by [EBRAINS](https://ebrains.eu/) infrastructure.

The 2 in `odeiv2` is from the GSL, there is no older R package of this
name.

The [integration methods](./integrators.md) (stepping functions) are also described in the official documentation.

To [install](INSTALL.md):

```R
remotes::install_github("a-kramer/rgsl")

```

OR

```R
remotes::install_github("icpm-kth/rgsl") # official fork
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
[vfgen](https://warrenweckesser.github.io/vfgen/) does. Any
method of auto-generating code will do, as long as it creates somewhat
compatible c files. We use [icpm-kth/RPN-derivative](icpm-kth/RPN-derivative).

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
t0 <- -1
t <- seq(0,1,length.out=nt)
y0 <- matrix(c(1,2,3),nrow=ny,ncol=m)
p <- matrix( ... ,nrow=np,ncol=m)

experiments <- list(
  a=list(initialState=y0, outputTimes=t, initialTime=t0, input=c(1,2)),
  b=list(initialState=y0, outputTimes=t, initialTime=t0, input=c(0,0))
)
y <- r_gsl_odeiv2_outer("ModelName",experiments,p)
```

This will simulate an outer product of simulation `experiments` with every column of `p`.

The function `r_gsl_odeiv2_outer()` will look for the file `ModelName.so` in
the current directory. Otherwise, a specific file-path can be selected using a `comment` on the model name:

```R
modelName <- "apoptosis"
comment(modelName) <- "../apoptosis.so"
```

The parameters can be an _n×m_ matrix of several parameterisations
(each column is a distinct parameter vector of length _n_). 

The solver will be looped _m_ times (using each column) within the C
code.

The simulations can be more intricate if the [event system](EVENTS.md) is used.

### Returns

A 3-dimensional array `y` of size _ny×nt×m_ (state, time, parameter set).

### Usage Notes

The time `t` is a vector of output time-points.

The parameters `p` must be numeric (not some arbitrary content passed
to the ODE right-hand-side functions).

## List of Experiments

This interface assumes that the user has a list of *N* simulation
experiments and uses this call structure:

```R
library("rgsl")
y<-r_gsl_odeiv2_outer(NameOfModel,experiments)
```

Each list item is a list of named properties that define a simulation
run, e.g.:

```R
experiments[[1]][["time"]] <- seq(0,1,length.out=100)
experiments[[1]][["initial_value"]] <- c(0,0,0)
experiments[[1]][["parameters"]] <- c(1,2,3,4,5)
```

The entries are found by their name, so they naming is not arbitrary
(but there are some alternative spellings):

|                           name | optional | meaning                           |
|-------------------------------:|:--------:|:----------------------------------|
|           `time`,`outputTimes` | no       | output time                       |
| `initial_value`,`initialState` | no       | a state space vector `y(t=t0)`    |
|     `parameters`,`param`,`par` | no       | a numeric vector of parameters    |
|     `events`,`scheduledEvents` | yes      | an event description for this run |

This interface is useful whenever your problem lends itself to this
description.

Particularly, if you want to simulate a subset of experiments, you can
remove unwanted items from the list easily (which is not as easy with
plain `r_gsl_odeiv2()`, see above).

### Returns

a list `y` of the same size as the experiment list, each entry has a
`[["state"]]` component and often a `[["func"]]` (unless specifically
suppressed) component (it is filled in if `_func` exists). Each
`y[[i]]$state` is a 3d-array, where the second index corresponds to
the output time points and the third enumerates the different
parameter vectors (possibly only 1).

The `func` component contains the output functions, as calculated by
the `${MODEL}_func` function. If that function is undefined, use 

```R
rgsl::r_gsl_odeiv2_outer_state_only # ${MODEL}_func is not used
```

### Structure of Simulation Experiments

`experiments`:

- `time` simulation output time
- `initial_value` *y(t=t₀)*
- `input` partial parameter vector
- `events` OPTIONAL
    + `time` event occurence times (nt)
    + `label` transformation label, a 0-based offset (used as `int`) which selects the transformation to apply
    + `dose` a scalar floating point variable (used as `double`)

### Returns

a list `y` of the same size as the experiment list, each entry has a
`[["state"]]` component and a `[["func"]]` component (it is filled in
if `_func` exists). Each `y[[i]][["state"]]` is a 3d array, where the
last dimension corresponds to the columns of `parameters`.

The `func` component contains the output functions, as calculated by
the `${MODEL}_func` function.

## SBtabVFGEN

We are working on models in systems biology and define models using
the SBtab format (loosely). Therefore, this project integrates somwhat well with
[SBtabVFGEN](https://github.com/a-kramer/SBtabVFGEN): it contains
`sbtab_to_vfgen`, which prints vfgen files (given SBtab). 

SBtab is a format within the realm of systems biology (but less common
than SBML) and is suited for writing/drafting biological models. The
workflow could be:


```
  User written:                           generated          Simulation
  +-----------+      +-----------+      +------------+      +----------+
  |           |      |  vfgen    |      |  (GSL)     |      |          |
  | SBtab (M) +--+-->+    OR     +----->+  ODE code  +--+-->+   rgsl   |
  |           |  |   |  ode.sh   |      | +jacobian  |  |   |          |
  +-----------+  |   +-----------+      +------------+  |   +----------+
                 |                                      |
            +----+--------------+                 +----------------+
            |                   |                 |                |
            | sbtab_to_vfgen()  |                 | r_gsl_odeiv2() |
            |                   |                 |                |
            +-------------------+                 +----------------+

```

VFGEN can be difficult to install. As an alternative, one can also
generate the ode code via `ode.sh` in the
[RPN-derivative/sh](github.com/icpm-kth/RPN-derivative) repository
(ode.sh accepts the sam efile format as vfgen, but only generates code
in C or R; it also lacks the time delay features of vfgen).

# Notes

In earlier versions we used OpenMP for parallelization, this has not
improved simulation times (parallelization is hard) and
`parallel::mclapply` turned out much easier to use (as are probably
the other functions in the parallel package).

OpenMP had the advantage that the returned value was shaped exactly as
we wanted it to be, while with `parallel::mclapply`, the result is
shaped by that function and need to be reshaped to recreate the old
values.

