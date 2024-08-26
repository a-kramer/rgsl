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
	expect_true(all(is.finite(Y[[1]]$state)))
	expect_true(all(is.finite(Y[[2]]$state)))
	for (y in Y){
		dev.new()
		plot(t,y$state[2,,1],main="HO: Outer Product function",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
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
	expect_false(all(is.finite(Y[[2]]$state)))
	expect_true(all(is.finite(Y[[1]]$state)))
	for (y in Y){
		dev.new()
		if (all(is.finite(y$state))){
			mt <- "HO: did not fail"
		} else {
			mt <- "HO: fails at t=1"
		}
		plot(t,y$state[2,,1],main=mt,sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
		for (l in 1:N) {
			for (k in 1:M) {
				lines(t,y$state[2,,k],lty=l)
			}
		}
	pty=1:M;
	pty[2:M]<-NA;
	}
})
