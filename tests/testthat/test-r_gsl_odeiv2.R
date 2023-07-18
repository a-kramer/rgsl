test_that("plain solver interface works",{
	file.c <- rgsl.example()
	so <- model.so(file.c)
	model.name <- sub("_gvf.c","",basename(file.c))
	comment(model.name)<- so
	expect_true(file.exists(so))
	t0  <- 0.0
	t <- seq(0.0,13,length.out=120)
	ny <- 2
	np <- 3
	I3 <- diag(1,3,3)
	I2 <- diag(1,2,2)

	N <- 12 # arbitrary number of simulations
	y0 <- matrix(c(0,1),nrow=ny,ncol=N)
	p <- matrix(c(1,0,0),nrow=np,ncol=N)
	p[2,]=seq(0,N-1,length.out=N) + rnorm(N,mean=0,sd=0.05);

	y <- r_gsl_odeiv2(model.name,t=t,t0=t0,y0=y0,p=p)
	expect_type(y,"double")
	expect_equal(dim(y),c(ny,length(t),N))
	expect_true(all(is.finite(y)))
	plot(t,y[2,,1],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
	for (l in 1:N) {
		lines(t,y[2,,l],lty=l)
	}
})

test_that("experiment solver interface works",{
	file.c <- rgsl.example()
	so <- model.so(file.c)
	model.name <- sub("_gvf.c","",basename(file.c))
	comment(model.name) <- so
	expect_true(file.exists(so))

	t <- seq(0.0,13,length.out=120)

	I3 <- diag(1,3,3)
	I2 <- diag(1,2,2)

	event.t <- c(1.0,2.0,3.0) # schedule of events

	state.tf <- affine.transform(length(event.t),I2,c(1,0))
	## experiment 1
	param.tf <- affine.transform(
			length(event.t),I3,c(0,0,0)
	)
	no.friction <- list(time=t,t0=0.0,parameters=c(1,0,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))
	## experiment 2
	param.tf <- affine.transform(
			length(event.t),I3,c(0,0.1,0)
	)
	low.friction <- list(time=t,t0=0.0,parameters=c(1,1,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))
	## experiment 3
	param.tf <- affine.transform(
			length(event.t),I3,c(0,0.2,0)
	)
	medium.friction <- list(time=t,t0=-10.0,parameters=c(1,2,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))
	## experiment 4
	medium.friction.no.events <- list(time=t,parameters=c(1,2,0),initial_value=c(0,1))

	experiments=list(a=no.friction,b=low.friction,c=medium.friction,d=medium.friction.no.events)

	N<-length(experiments)

	y <- r_gsl_odeiv2_sim(model.name,experiments)

	expect_type(y,"list")
	expect_length(y,N)
	expect_named(y,c("a","b","c","d"))
	expect_named(y[['a']],c("state","func"))
	expect_true(all(is.finite(y$a$state)))
	expect_true(all(is.finite(y$a$func)))
	expect_true(all(is.finite(unlist(y))))
	expect_gt(max(unlist(y)),0.0)
	expect_lt(min(unlist(y)),0.0)

	plot(t,y[[1]][["state"]][2,],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
	for (l in 1:N) {
		lines(t,y[[l]][["state"]][2,],lty=l)
	}
	pty=1:N;
	pty[2:N]<-NA;
	legend("bottomright",legend=names(experiments),lty=1:N,pch=pty)
})

test_that("outer product interface works",{
	file.c <- rgsl.example()
	so <- model.so(file.c)
	model.name <- sub("_gvf.c","",basename(file.c))
	comment(model.name) <- so
	expect_true(file.exists(so))

	t <- seq(0.0,13,length.out=120)
	## 2-dim and 3-dim identity matrix
	I3 <- diag(1,3,3)
	I2 <- diag(1,2,2)
	M <- 40 # arbitrary number of random parameter vectors
	event.t <- c(1.0)
	paramG <- matrix(rnorm(M,mean=1,sd=0.2),nrow=1,ncol=M)
	state.tf <- affine.transform(length(event.t),I2,c(1,0))

	param.tf <- affine.transform(length(event.t),I3,c(0,0,0))
	no.friction <- list(time=t,t0=0.0,input=c(0,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))

	experiments=list(a=no.friction,b=no.friction)
	N <- length(experiments)

	Y <- r_gsl_odeiv2_outer(model.name,experiments,paramG)
	expect_type(Y,"list")
	expect_length(Y,2)
	expect_named(Y,c("a","b"))
	expect_true(all(is.finite(unlist(Y))))
	for (y in Y){
		dev.new()
		plot(t,y$state[2,,1],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
		for (l in 1:N) {
			for (k in 1:M) {
				lines(t,y$state[2,,k],lty=l)
			}
		}
	pty=1:M;
	pty[2:M]<-NA;
	}
})

test_that("failed simulation is different from a successful simulation",{
	file.c <- rgsl.example()
	so <- model.so(file.c)
	model.name <- sub("_gvf.c","",basename(file.c))
	comment(model.name) <- so
	expect_true(file.exists(so))

	t <- seq(0.0,13,length.out=120)
	## 2-dim and 3-dim identity matrix
	I3 <- diag(1,3,3)
	I2 <- diag(1,2,2)
	M <- 40 # arbitrary number of random parameter vectors
	event.t <- c(1.0)
	paramG <- matrix(rnorm(M,mean=1,sd=0.2),nrow=1,ncol=M)
	ok.state.tf <- affine.transform(length(event.t),I2,c(1,0))
	bad.state.tf <- affine.transform(length(event.t),I2,c(Inf,1.0))
	param.tf <- affine.transform(length(event.t),I3,c(0,0,0))
	probably.ok <- list(time=t,t0=0.0,input=c(0,0),initial_value=c(0,1),events=event.tf(event.t,ok.state.tf,param.tf))
	probably.fails <- list(time=t,t0=0.0,input=c(0.0,-100.0),initial_value=c(0,1),events=event.tf(event.t,bad.state.tf,param.tf))

	experiments=list(a=probably.ok,b=probably.fails)
	N <- length(experiments)

	Y <- r_gsl_odeiv2_outer(model.name,experiments,paramG)
	expect_type(Y,"list")
	expect_length(Y,2)
	expect_named(Y,c("a","b"))
	expect_false(all(is.finite(unlist(Y[[2]]))))
	expect_true(all(is.finite(unlist(Y[[1]]))))
	for (y in Y){
		dev.new()
		plot(t,y$state[2,,1],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
		for (l in 1:N) {
			for (k in 1:M) {
				lines(t,y$state[2,,k],lty=l)
			}
		}
	pty=1:M;
	pty[2:M]<-NA;
	}
})
