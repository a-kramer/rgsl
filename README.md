# rgsl - an interface between R and gsl_odeiv2

This R package solves a series of initial value problems for ordinary
differential equations formulated as a list of simulation
experiments to be integrated using [C-functions](MODELS.md).

This project was [funded by the EU](./ACKNOWLEDGMENTS.md) as part of
the [Human Brain Project](https://www.humanbrainproject.eu/en/) and
supported by [EBRAINS](https://ebrains.eu/) infrastructure.

The 2 in `odeiv2` is from the GSL, there is no older R package of this
name.

The [integration methods](./integrators.md) (stepping functions) we interface with are
also described in the official [documentation](https://www.gnu.org/software/gsl/doc/html/ode-initval.html).

To [install](INSTALL.md):

```R
remotes::install_github("a-kramer/rgsl")

```

OR

```R
remotes::install_github("icpm-kth/rgsl") # official fork
```

### Note on Operating Systems

This package will work on UNIX-like systems (GNU Linux, BSD, MACOS),
but almost certainly not on windows, because of the shared libraries,
compilers, etc.

This is not just plain R code.

It may work coincidentally, but we have not tried.

## Description and Purpose

The C code calls
[gsl_odeiv2](https://www.gnu.org/software/gsl/doc/html/ode-initval.html)
module functions to solve the problem or problems. The goal is to
offload as much work as possible to the C code and keep the overhead
minimal. That is why this package expects to solve a set of problems,
rather than one, for the same model file (varying in e.g.: initial
conditions, or parameters).

The package contains only high-level interface functions, each
interface function returns a list of solutions. The interface
functions are described in the following Sections.

The ODE has to exist as a [shared library](MODELS.md) (`.so`)
file. There are some assumptions we make about the contents of the
shared library file, and the naming of the functions.

Here we assume that the solutions serve some scientific purpose and
the lab experiments come with a description of _observables_: measureable values
that depend on the system's state. We call the part of the model that calculates the observables
`${ModelName}_func()` (vfgen also calls them *Functions* of the model).

### Note

We aim to find the perfect balance between speed and
feasibility. This package cannot use an ODE written in R, which the
reader may consider impractical, but it is faster. For plain R
models, there is of course
[deSolve](https://cran.r-project.org/web/packages/deSolve/) that
remains an alternative. Our strategy is to auto-generate the shared
library whenever possible and avoid writing the model files by
hand. This package does not itself auto-generate code, but
[vfgen](https://warrenweckesser.github.io/vfgen/) does. Any method of
auto-generating code will do, as long as it creates somewhat
compatible c files. We use
[icpm-kth/RPN-derivative](icpm-kth/RPN-derivative).

The GSL
[odeiv2](https://www.gnu.org/software/gsl/doc/html/ode-initval.html)
module is well documented. In addition, we have decided on this rule:
the C functions return the number of principal return values, when
called with `NULL` pointers.

We have considiered the SUNDIALS solvers like Cvode, but they turned
out too complex for us, the documentation is difficult to read, and
the types, macros, and function names change with major versions. We
consider it not feasible (for us) to maintain an interface to Cvode
and chose the GSL solvers instead.

## List of Experiments `gsl_odeiv2_outer` and `gsl_odeiv2_outer_state_only`

This interface assumes that the user has a list of *N* simulation
experiments and uses this call structure:

```R
library("rgsl")
comment(NameOfModel) <- "path/to/sharedLibrary.so"  # optional
y <- r_gsl_odeiv2_outer(NameOfModel,experiments,p)
```

Where `p` refers to the parameters of the model. This may be a matrix
of _M_ parameter vectors, each column will result in a simulation of
every experiment. The result is the outer product of _M_×_M_
simulations.

Each list item of `experiments` is itself a list of named properties that define a simulation
run, e.g.:

```R
experiments[[1]][["time"]] <- seq(0,1,length.out=100)
experiments[[1]][["initial_value"]] <- c(0,0,0)
experiments[[1]][["input"]] <- c(1,2,3,4,5)
```

The entries are found by name, so their naming is not arbitrary
(but there are some alternative spellings):

|                           name | optional | meaning                           |
|-------------------------------:|:--------:|:----------------------------------|
|           `time`,`outputTimes` | no       | output time                       |
| `initial_value`,`initialState` | no       | a state space vector `y(t=t0)`    |
|     `events`,`scheduledEvents` | yes      | an event description for this run |
|                        `input` | yes      | additional parameters             |

The model knows how many parameters it needs to simulate, it uses the
suuplied parameter matrix `p` to fill the internal parameter buffer
from the front and the input parameters from the `experiments` to fill
the buffer from the back. In the normal case, these parameters complement each other:

```R
# something like this will happen in the C code, but with memcpy()
# p: supplied parameters
# u: input parameters
modelPar <- c(p,u) # effective model parameters
```

The idea is that `p` represents the model's internal parameters and
the values come from optimization, sampling, or various competing
model hypotheses; but the input parameters `u` are related to the
actual experiments. These input parameters are what makes several
experiments different (in addition to initial values and scheduled
events).

It is easy to simulate a subset of experiments by subsetting the experiments list.
If the experiments are named then so is the return value.

### Returns

a list `y` of the same size as the experiment list, each entry has a
`[["state"]]` component and a `[["func"]]` (the function
`r_gsl_odeiv2_outer_state_only` suppresses this) component (it is
filled in if `_func` exists). Each `y[[i]]$state` is a 3d-array, where
the second index corresponds to the output time points and the third
enumerates the _M_ different parameter vectors (possibly only 1).

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

