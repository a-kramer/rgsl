# Excerpt from the official documentation

GSL solvers for Ordinary Differential Equations (ODEs) are documented
on
[www.gnu.org](https://www.gnu.org/software/gsl/doc/html/ode-initval.html).

This package has two utility functions to select the integration method:

```R
rgsl::integrationMethod("msbdf") == 0
rgsl::nameMethod(0) == "msbdf"
```

The value of `integrationMethod` can be used here:

```R
rgsl::r_gsl_odeiv2_outer(...,method=integrationMethod("msbdf"))
rgsl::r_gsl_odeiv2_outer_sens(...,method=integrationMethod("rkf45"))
rgsl::r_gsl_odeiv2_outer_state_only(...,method=integrationMethod("msadams"))
```

## rk2

```
gsl_odeiv2_step_type *gsl_odeiv2_step_rk2
```

Explicit embedded Runge-Kutta (2, 3) method.

## rk4

```c
gsl_odeiv2_step_type *gsl_odeiv2_step_rk4
```

Explicit 4th order (classical) Runge-Kutta. Error estimation is
carried out by the step doubling method. For more efficient estimate
of the error, use the embedded methods described below.

## rkf45

```c
gsl_odeiv2_step_type *gsl_odeiv2_step_rkf45
```

Explicit embedded Runge-Kutta-Fehlberg (4, 5) method. This method is a
good general-purpose integrator.

## rkck

```c
gsl_odeiv2_step_type *gsl_odeiv2_step_rkck
```

Explicit embedded Runge-Kutta Cash-Karp (4, 5) method.

## rk8pd

```c
gsl_odeiv2_step_type *gsl_odeiv2_step_rk8pd
```

Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method.

## rk1imp

```c
gsl_odeiv2_step_type *gsl_odeiv2_step_rk1imp
```

Implicit Gaussian first order Runge-Kutta. Also known as implicit
Euler or backward Euler method. Error estimation is carried out by the
step doubling method. This algorithm requires the Jacobian and access
to the driver object via `gsl_odeiv2_step_set_driver()`.

## rk2imp

```c
gsl_odeiv2_step_type *gsl_odeiv2_step_rk2imp
```

Implicit Gaussian second order Runge-Kutta. Also known as implicit
mid-point rule. Error estimation is carried out by the step doubling
method. This stepper requires the Jacobian and access to the driver
object via `gsl_odeiv2_step_set_driver()`.

## rk4imp

```c
gsl_odeiv2_step_type *gsl_odeiv2_step_rk4imp
```

Implicit Gaussian 4th order Runge-Kutta. Error estimation is carried
out by the step doubling method. This algorithm requires the Jacobian
and access to the driver object via `gsl_odeiv2_step_set_driver()`.

## bsimp

```c
gsl_odeiv2_step_type *gsl_odeiv2_step_bsimp
```

Implicit Bulirsch-Stoer method of Bader and Deuflhard. The method is
generally suitable for stiff problems. This stepper requires the
Jacobian.

## msadams

```c
gsl_odeiv2_step_type *gsl_odeiv2_step_msadams
```

A variable-coefficient linear multistep Adams method in Nordsieck
form. This stepper uses explicit Adams-Bashforth (predictor) and
implicit Adams-Moulton (corrector) methods in $P(EC)^m$ functional
iteration mode. Method order varies dynamically between 1 and 12. This
stepper requires the access to the driver object via
`gsl_odeiv2_step_set_driver()`.

## msbdf

```c
gsl_odeiv2_step_type *gsl_odeiv2_step_msbdf
```

A variable-coefficient linear multistep backward differentiation
formula (BDF) method in Nordsieck form. This stepper uses the explicit
BDF formula as predictor and implicit BDF formula as corrector. A
modified Newton iteration method is used to solve the system of
non-linear equations. Method order varies dynamically between 1
and 5. The method is generally suitable for stiff problems. This
stepper requires the Jacobian and the access to the driver object via
`gsl_odeiv2_step_set_driver()`.

# List of all integrators in reverse order

Taken from this file, using coreutils:

```sh
printf "{" ; egrep "^gsl_" integrators.md | sed -r 's/^gsl_[^ ]+ [*]//' | tac | sed 's/$/, /' | tr -d '\n'; printf "NULL}\n"
```

```
{gsl_odeiv2_step_msbdf, gsl_odeiv2_step_msadams, gsl_odeiv2_step_bsimp, gsl_odeiv2_step_rk4imp, gsl_odeiv2_step_rk2imp, gsl_odeiv2_step_rk1imp, gsl_odeiv2_step_rk8pd, gsl_odeiv2_step_rkck, gsl_odeiv2_step_rkf45, gsl_odeiv2_step_rk4, gsl_odeiv2_step_rk2, NULL}
```

The list is `NULL` terminated.
