# Events

Often, a simulation requires an instantaneous intervention at a known
point in time. We don't want to interupt simulations and return to R
for simple enough cases.

So, for this purpose, a list of scheduled events can be fed into the
solver. This does not cover all possible needs and use cases, but this
concept can be used to mimic some experimental setups well enough.

Each simulation can receive an event of this kind:
```R
# r_gsl_odeiv2 with m parameter columns
event[[i]] <- list(
                time=event.times,
                tf=list(
				    state=list(A=...,b=...),
				    param=list(A=...,b=...)
					)
				)
```

where at time `event.times[j]` both the state (and the parameters)
will be transformed, like this:

```R
y <- A %*% y + b
```

but in C (using
[BLAS](https://www.gnu.org/software/gsl/doc/html/blas.html)).

## Events within the plain interface `r_gsl_odeiv2`

When using the plain interface `r_gsl_odeiv2(...)`, `events` can be
supplied as a list with as many entries as there are parameter
columns, or a shorter *named* list, with names corresponding to the
column-names of `p`. In this case experiments with no events don't
require a blank entry in the events list (if they do, then it's a
bug). The events list contains only meaningful (non-empty) events and
they are accessed by name.

If `time` is a vector, then all transformation quantities can be 3
dimensional arrays, where the third dimension corresponds to time:

```R
# at time point j, the transformation is
y <- A[,,j] %*% y + b[,1,j]
# similarly for p, but with p's A and b
```

(Note: the second dimension of `b` remains unused).

If for all j the matrices `A[,,j]` are diagonal, then the second dimension of these
matrices can be 1 and understood as a list of diagonals:

```R
# if A and B are diagonals, then for all time[j]
y <- A[,1,j] * y + a[,1,j]
p <- B[,1,j] * p + b[,1,j]
```

but in C (again).

If `p` is a vector (rather than a matrix), then the list of events is _allowed_ to be just
one event, not _necessarily_ within a _one element list_ (but that is also ok).

This event system is less powerful than a triggered event system,
where an instantaneous event occurs when a condition is met.

The same effect could be achieved by calling `r_gsl_odeiv2` and
integrating up to the event time, perform the transformation in R and
then continue the integration from the transformed state using a
second call to `r_gsl_odeiv2`.

## Events for the other interfaces (*sim* and *outer*)

For `r_gsl_odeiv2_sim` and `r_gsl_odeiv2_outer` the events structure
is always part of the experiments list:

```R
ev.t <- c(-1,0,1)
nt<-length(ev.t)

experiments[[i]][["events"]][["time"]] <- ev.t
experiments[[i]][["events"]][["tf"]][["param"]][["A"]] <- array(0,dim=c(np,np,nt))
experiments[[i]][["events"]][["tf"]][["param"]][["b"]] <- array(0.34,dim=c(np,1,nt))
```

sets the `i`th experiment's parameter transformation matrices `A` to
zero for all 3 event times (-1.0, 0.0, and 1.0) and the offset `b` to
0.34 for all parameter values and all event times (a vector of length
`nt`).

The `events` list-structure and names:

- `time` time schedule for sudden events
- `tf` transformations 
    + `state` transforms y
	    * `A` a matrix (multiplicative)
		* `b` a vector (additive)
	+ `param` transforms p
	    * `A` a matrix (multiplicative)
		* `b` a vector (additive)

