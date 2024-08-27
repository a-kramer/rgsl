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


test_that("event function works",{
	modelName <- "CaSpike"
	model.c <- rgsl.example(modelName)
	comment(modelName) <- model.so(model.c)
	y0 <- 0.7
	p <- c(tau=0.01)
	t_out <- seq(0,1,length.out=100)
	ev.t <- seq(0,1,by=0.05)
	event <- list(label=integer(length(ev.t)),time=ev.t,dose=integer(length(ev.t)))
	a <- list(input=c(dCa=0.5),event=event,outputTimes=t_out,t0=0)
	b <- list(input=c(dCa=1.5),event=event,outputTimes=t_out,t0=0)
	c <- list(input=c(dCa=3.5),event=event,outputTimes=t_out,t0=0)
	experiments <- list(a,b,c)
	y <- r_gsl_odeiv2_outer(modelName,experiments,as.matrix(p))
	expect_type(y,"list")
	expect_length(y,length(experiments))
	expect_named(y,names(experiments))
	expect_true(all(is.finite(y[[1]]$state)))
	expect_true(all(is.finite(y[[2]]$state)))
	dev.new()
	plot(t_out,y[[1]]$state[1,,1],main=names(experiments)[1],sub="y(t) = [t>t.ev] * dCa*exp(-t/tau) + CaBase",xlab="time",ylab="state y[1]",type="p",ylim=c(0,max(y[[3]]$state)))
	for (i in seq_along(y)){
		Ca <- y[[i]]$state[1,,1]
		expect_length(Ca,length(t_out))
		lines(t_out,Ca)
	}
	legend("topright",legend=c("",names(experiments)),lty=c(NA,1,2,3),pch=c(1,NA,NA,NA))
})
