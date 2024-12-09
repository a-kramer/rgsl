# Models

When calling the solver functions, the supplied model name will
determine the expected name of the C functions:

```R
y <- r_gsl_odeiv2_outer("xyzq",experiments)
```

Implies that there is a file called `xyzq.so` with functions called
`xyzq_vf()` and other underscore-suffixes.

The file name can differ from the model name, for this, the C code will read the attached comment:

```R
modelName <- "xyzq"
comment(modelName) <- "../a.out.so"
y <- r_gsl_odeiv2_outer(modelName,experiments)
```

This way, the C-code will open `a.out.so` in the parent directory and
not guess.

This package requires the model to be provided as a shared library,
there is no alternative to this, For the remainder of our
explanations, we will use posix-shell syntax to talk about the
function names: `${MODEL}.so`.

The function names within the file need to correspond to the given model
name `MODEL`:

```sh
${MODEL}_vf();
${MODEL}_jac();
${MODEL}_jacp();
${MODEL}_func();
```

The functions are found by name (cannot be arbitrary, in prefix_ or _suffix). The
program [vfgen](https://warrenweckesser.github.io/vfgen/) can be
helpful in creating the `_vf` and `_jac` functions (GSL's odeiv is one
of vfgen's output formats), but not `_func` (it's not one of the
things that vfgen makes, it makes scalar output functions).

The last entry: `*_func` calculates the observables. If you think of
the model as an input/output system then `_func` refers to the
output. This is trivial to write by hand, if needed.

Since the functions are loaded from a shared library, it is not easy
to find out the dimensionality of either the *state space* or *parameter
space*.

To make working with these model files easier for us, we use more
specific return status codes.

This is an example header file, where `MODEL=HarmonicOscillator`

```c
int HarmonicOscillator_vf(double, const double [], double [], void *);
int HarmonicOscillator_jac(double, const double [], double *, double *, void *);
int HarmonicOscillator_jacp(double, const double [], double *, void *);
int HarmonicOscillator_func(double, const double [], double*, void *);
```

All of the routines return an integer status value. On success they
must return `GSL_SUCCESS` (0). But, on failure we make them return the size
of the array they expect in the output slot:

```c
/* this function returns ẏ=f(...) in the third slot */
int HarmonicOscillator_vf(
 double t,  /* time t */
 const double y_[], /* state vector */
 double f_[], /* OUTPUT:  dy[i]/dt = f_[i] */
 void *params) /* parameter values */
{
...
```

To find out how much memory needs to be allocated for `f_` we call
this function with `NULL` in the third slot (`y_` may also be `NULL`
then), and it shall return the integer 2 (because this model has
two-dimensional state-vectors).

VFGEN does not make this happen. Modify the code by hand, semi-automated, or use
[RPN-derivative/sh/ode.sh](icpm-kth/RPN-derivative).

This strategy makes the C files a bit more autonomous as sizes are not
entirely taken on faith (from the supplied initial conditions).

```c
/* The vector field. */
int HarmonicOscillator_vf(double t, const double y_[], double f_[], void *params)
{
    double v, y;
    double k, c, F;
    double v_flux;
    double *p_;
    int RET=GSL_SUCCESS;
    if (y_ && f_){
      p_ = (double *) params;

      v          = y_[0];
      y          = y_[1];

      k          = p_[0];
      c          = p_[1];
      F          = p_[2];

      v_flux = -k*y+F-v*c;

      f_[0] = v_flux;
      f_[1] = v;
    } else {
      RET=2;
    }
    return RET;
}
```

Or use an early return, like RPN-derivative does.

## Example File

See the example file [HarmonicOscillator](./inst/extdata/HarmonicOscillator_gvf.c)
to inspect an auto-generated right-hand-side and Jacobian function
compatible with the solvers from the gsl.
