# Investigation of Nested Parallelization

These are timing results when OpenMP is used together with `parallel::mclapply` or each individually.

## Just OpenMP 

Here we test how efficiently the cores are used, the maximum number of threads started by OpenMP can be controlled through the environment variable `OMP_NUM_THREADS`, like this:

```sh
$ export OMP_NUM_THREADS=2
$ R
```

The timing is done for the *Harmonic Oscillator* model, 
supplied alongside the code, using the `test.outer2()` function 
which is OpenMP-parallel in the columns of the supplied 
parallelizations to the model.

The time is measured as `system.time(y<-test.outer2(5000))`, where 5000 is the amount of work (parameter columns to simulate with).

The machine has 8 cores of this type:

```
Intel(R) Core(TM) i7-7700HQ CPU @ 2.80GHz
```
Nevertheless, a value of 16 was still permitted and just lead to all 8 cores being busy (a higher number is not illegal it seems). 

### test.outer2

Workload: `5000`

| `OMP_NUM_THREADS`/cores | elapsed time in s |
| -----------------:|:-----------------:|
|16/8|9.299|
|8/8|9.297|
|6/8|8.130|
|4/8|5.679|
|2/8|5.691|
|1/8|3.015|

This is quite surprising, the number of threads adversely influences the run-time. This could mean that the overhead is quite substantial.



### test.plain

Workload: `5000`

| `OMP_NUM_THREADS`/cores | elapsed time in s |
|-----------------:|:-----------------:|
|8/8|8.723|
|6/8|6.397|
|4/8|4.859|
|2/8|4.676|
|1/8|2.630|

In this case, parallelization has an adverse effect.

## Just mclapply

The command: 
```R
n<-parallel::detectCores()
x<-rep(5000/n,n)
system.time(y<-mclapply(x,FUN=test.plain,mc.cores=n))
```

|`mc.cores`|time in s|
|---------:|:-------:|
|8|0.665|
|5|0.740|
|2|1.248|
|1|2.245|

R's parallel package parallizes well, this is what we would expect.


## Both methods together:

The command:

```sh
$ export OMP_NUM_THREADS=2
$ R
```
and then
```R
source("HarmonicOscillator.R")
n<-parallel::detectCores()/2
x<-rep(5000/n,n)
system.time(y<-mclapply(x,FUN=test.plain,mc.cores=n))
```

|`OMP_NUM_THREADS`|`mc.cores`|time in s|
|----------------:|:--------:|:-------:|
|1|8|0.677|
|2|4|1.172|
|4|2|3.570|
|8|8|14.293|

OpenMP makes it worse.

Using `test.outer2` instead of `test.plain` leads to similar results.

## Turning OpenMP off by removing the pragma

OpenMP needs to be activated in the C code by statements such as

```C
/* on */
#pragma omp parallel for

/* off */
//#pragma omp parallel for
```

Is there an improvement if this is turned off (commented out)? 

```R
n<-parallel::detectCores()
print(n)
  
x<-rep(5000/n,n)
system.time(y<-mclapply(x,FUN=test.outer2,mc.cores=n))
[1] 8
   user  system elapsed 
  5.371   0.103   0.816 
```

no improvement over the previous test.

```R  
x<-rep(5000/n,n)
system.time(y<-mclapply(x,FUN=test.plain,mc.cores=n))
[1] 8
   user  system elapsed 
  4.559   0.114   0.696 
```


# Irreproducible effects

On another machine (Xeon Workstation), the call to mclapply never returns if OpenMP is turned on. 

This needs to be investigated further.



