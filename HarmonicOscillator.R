
test.plain <-function(N=3){
	name <- "HarmonicOscillator"
	comment(name)  <- "./HarmonicOscillator.so"
	t <- seq(0,13,length.out=120)
	ny <- 2
	np <- 3
	I3 <- diag(1,3,3)
	I2 <- diag(1,2,2)

	y0 <- matrix(c(0,1),nrow=ny,ncol=N)
	p <- matrix(c(1,0,0),nrow=np,ncol=N)

	p[2,]=seq(0,N-1,length.out=N) + rnorm(N,mean=0,sd=0.05);

	if (require("rgsl")){
		y <- r_gsl_odeiv2(name,t=t,y0=y0,p=p)
		plot(t,y[2,,1],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
		for (l in 1:N) {
			lines(t,y[2,,l],lty=l)
		}
		pty=1:N;
		pty[2:N]<-NA;
##		legend("bottomright",legend=colnames(p),lty=1:N,pch=pty)
	}
	return(y)
}


test.events <-function(){
	name <- "HarmonicOscillator"
	t <- seq(0,13,length.out=120)
	N <- 3
	ny <- 2
	np <- 3
	I3 <- diag(1,3,3)
	I2 <- diag(1,2,2)

	y0 <- matrix(c(0,1),nrow=ny,ncol=N)
	p <- matrix(c(1,0,0),nrow=np,ncol=N)

	for (i in 1:N) {
		p[2,i]=(i-1); #/N
	}
	colnames(p) <- c("no friction","tiny amount of friction","lots of friction")
	## events:
	all.events=list()
	event.t <- c(1.0,2.0,3.0)
	state.tf <- affine.transform(length(event.t),I2,c(1,0))

	param.tf <- affine.transform(length(event.t),I3,c(0,0,0))
	all.events[[1]] <- event.tf(event.t,state.tf,param.tf)

	param.tf <- affine.transform(length(event.t),I3,c(0,0.1,0))
	all.events[[2]] <- event.tf(event.t,state.tf,param.tf)

	param.tf <- affine.transform(length(event.t),I3,c(0,0.2,0))
	all.events[[3]] <- event.tf(event.t,state.tf,param.tf)

	names(all.events) <- colnames(p);

	if (require("rgsl")){
		y <- r_gsl_odeiv2(name,t=t,y0=y0,p=p,all.events)
		plot(t,y[2,,1],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
		for (l in 1:N) {
			lines(t,y[2,,l],lty=l)
		}
		pty=1:N;
		pty[2:N]<-NA;
		legend("bottomright",legend=colnames(p),lty=1:N,pch=pty)
	}
	return(y)
}



test.experiments <-function(){
	name <- "HarmonicOscillator"
	t <- seq(0,13,length.out=120)
	## 2-dim and 3-dim identity matrix
	I3 <- diag(1,3,3)
	I2 <- diag(1,2,2)
	N <- 4
	event.t <- c(1.0,2.0,3.0)

	state.tf <- affine.transform(length(event.t),I2,c(1,0))

	param.tf <- affine.transform(length(event.t),I3,c(0,0,0))
	no.friction <- list(time=t,parameters=c(1,0,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))

	param.tf <- affine.transform(length(event.t),I3,c(0,0.1,0))
	low.friction <- list(time=t,parameters=c(1,1,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))

	param.tf <- affine.transform(length(event.t),I3,c(0,0.2,0))
	medium.friction <- list(time=t,parameters=c(1,2,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))

	medium.friction.no.events <- list(time=t,parameters=c(1,2,0),initial_value=c(0,1))

	experiments=list(a=no.friction,b=low.friction,c=medium.friction,d=medium.friction.no.events)
	if (require("rgsl")){
		y <- r_gsl_odeiv2_sim(name,experiments)
		plot(t,y[[1]][["state"]][2,],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
		for (l in 1:N) {
			lines(t,y[[l]][["state"]][2,],lty=l)
		}
		pty=1:N;
		pty[2:N]<-NA;
		legend("bottomright",legend=names(experiments),lty=1:N,pch=pty)
	}
	return(y)
}

test.experiments2 <-function(){
	name <- "HarmonicOscillator"
	t <- as.double(seq(0,13,length.out=120))
	## 2-dim and 3-dim identity matrix
	I3 <- diag(1,3,3)
	I2 <- diag(1,2,2)
	N <- 3
	event.t <- c(1.0,2.0,3.0)

	state.tf <- affine.transform(length(event.t),I2,c(1,0))

	param.tf <- affine.transform(length(event.t),I3,c(0,0,0))
	no.friction <- list(outputTimes=t,parameters=c(1,0,0),initialState=c(0,1),scheduledEvents=event.tf(event.t,state.tf,param.tf))

	param.tf <- affine.transform(length(event.t),I3,c(0,0.1,0))
	low.friction <- list(outputTimes=t,parameters=c(1,1,0),initialState=c(0,1),events=event.tf(event.t,state.tf,param.tf))

	param.tf <- affine.transform(length(event.t),I3,c(0,0.2,0))
	medium.friction <- list(outputTimes=t,parameters=c(1,2,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))
	experiments=list(a=no.friction,b=low.friction,c=medium.friction)
	if (require("rgsl")){
		y <- r_gsl_odeiv2_sim(name,experiments)
		plot(t,y[[1]][["state"]][2,],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
		for (l in 1:N) {
			lines(t,y[[l]][["state"]][2,],lty=l)
		}
		pty=1:N;
		pty[2:N]<-NA;
		legend("bottomright",legend=names(experiments),lty=1:N,pch=pty)
	}
	return(y)
}

test.outer <-function(M=4){
	name <- "HarmonicOscillator"
	t <- seq(0,13,length.out=120)
	## 2-dim and 3-dim identity matrix
	I3 <- diag(1,3,3)
	I2 <- diag(1,2,2)
	N <- 1
	event.t <- c(1.0)
	paramG <- matrix(rnorm(M,mean=1,sd=0.2),nrow=1,ncol=M)
	paramU <- matrix(runif(M,min=0.5,max=1.5),nrow=1,ncol=M)
	state.tf <- affine.transform(length(event.t),I2,c(1,0))

	param.tf <- affine.transform(length(event.t),I3,c(0,0,0))
	no.friction <- list(time=t,input=c(0,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))

	experiments=list(a=no.friction)

	if (require("rgsl")){
		param<-list(Uniform=paramU,Gaussian=paramG)
		Y <- parallel::mclapply(X=param, function(p) r_gsl_odeiv2_outer(name,experiments,p))
		for (y in Y){
		dev.new()
		plot(t,y[[1]][["state"]][2,,1],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
		for (l in 1:N) {
		for (k in 1:M){
			lines(t,y[[l]][["state"]][2,,k],lty=l)
		}
		}
		pty=1:M;
		pty[2:M]<-NA;
		}
	}
	return(y)
}

test.outer.parallel <- function(M=4){
	require(parallel)
	name <- "HarmonicOscillator"
	comment(name) <- "./HarmonicOscillator.so"
	t <- seq(0,13,length.out=120)
	## 2-dim and 3-dim identity matrix
	I3 <- diag(1,3,3)
	I2 <- diag(1,2,2)
	## 1 free parameter, 2 input parameters
	paramG <- matrix(rnorm(1*M,mean=1,sd=0.2),nrow=1,ncol=M)

	event.t <- c(1.0,2.0,3.0)
	state.tf <- affine.transform(length(event.t),I2,c(1,0))

	param.tf <- affine.transform(length(event.t),I3,c(0,0,0))
	no.friction <- list(time=t,initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf),input=c(0,0))

	param.tf <- affine.transform(length(event.t),I3,c(0,0.1,0))
	low.friction <- list(time=t,initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf),input=c(0,1))

	param.tf <- affine.transform(length(event.t),I3,c(0,0.2,0))
	medium.friction <- list(time=t,initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf),input=c(0,2))

	medium.friction.no.events <- list(time=t,initial_value=c(0,1),input=c(0,2))

	experiments <- rep(list(a=no.friction,b=low.friction,c=medium.friction,d=medium.friction.no.events),4)

	N <- length(experiments)
	if (require("rgsl")){
		cat("not at all parallel\n")
		print(system.time(Y <- r_gsl_odeiv2_outer(name,experiments,paramG)))
		cat("mclapply parallel\n")
		print(system.time(Y <- unlist(mclapply(1:N, function(i) r_gsl_odeiv2_outer(name,experiments[i],paramG)),recursive=FALSE)))
		if (M<100){
		for (y in Y){
		dev.new()
		plot(t,y[["state"]][2,,1],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
		for (k in 1:M) {
			lines(t,y[["state"]][2,,k],lty=1)
		}
		pty=1:M;
		pty[2:M]<-NA;
		}
		}
	}
	return(Y)
}

test.outer.wo.events<-function(){
	name <- "HarmonicOscillator"
	t <- seq(0,13,length.out=120)
	N <- 1
	M <- 20
	param <- matrix(rnorm(M,mean=1,sd=0.2),nrow=1,ncol=M)
	no.friction <- list(time=t,input=c(0,0),initial_value=c(0,1))
	experiments<-list(a=no.friction)
	if (require("rgsl")){
		y <- r_gsl_odeiv2_outer(name,experiments,param)
		plot(t,y[[1]][["state"]][2,,1],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
		for (l in 1:N) {
		for (k in 1:M){
			lines(t,y[[l]][["state"]][2,,k],lty=l)
		}
		}
		pty=1:M;
		pty[2:M]<-NA;
		legend("bottomright",legend=sprintf("%4.3g",param),lty=1:N,pch=pty)

	}
	return(y)

}
