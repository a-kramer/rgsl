# Models 

This package requires the model to be provided as a shared library:
`${MODEL}.so`. The function names within the file need to correspond to the model name:

```sh
${MODEL}_vf();
${MODEL}_jac();
${MODEL}_jacp();
${MODEL}_func();
```

The functions are found by name and cannot have arbitrary names. The
program [vfgen](https://warrenweckesser.github.io/vfgen/) can be
helpful in creating the `_vf` and `_jac` functions (GSL's odeiv is one
of vfgen's output formats), but not `_func` (it's not one of the
things that vfgen makes, it makes scalar output functions).

The last entry: `*_func` calculates the observables. If you think of
the model as an input/output system then `_func` refers to the output.

Since the functions are loaded from a shared library, it is not easy
to find out the dimensionality of either the *state space* or *parameter
space*. 

To make working with these model files easier for us, we use more
specific return status codes.

This is an example header file:
```c
int HarmonicOscillator_vf(double, const double [], double [], void *);
int HarmonicOscillator_jac(double, const double [], double *, double *, void *);
int HarmonicOscillator_jacp(double, const double [], double *, void *);
int HarmonicOscillator_func(double, const double [], double*, void *);
```

All of the routines return an integer status value. On success they
must return `GSL_SUCCESS`. But on failure we make them return the size
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

VFGEN does not make this happen: we modified the function to do this
by hand.

This strategy makes the c files a bit more autonomous as sizes are not
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


## Example File

See the example file [HarmonicOscillator](./HarmonicOscillator_gvf.c)
to inspect an auto-generated right-hand-side and Jacobian function
compatible with the solvers from the gsl.

The example script [HarmonicOscillator.R](./HarmonicOscillator.R)
contains a `test.something()` functions (where `something` is one of:
plain, events, experiments, outer).

These C code files (for the gsl module
[odeiv2](https://www.gnu.org/software/gsl/doc/html/ode-initval.html))
can be automatically obtained using from
[VFGEN](https://github.com/WarrenWeckesser/vfgen) or written manually.
