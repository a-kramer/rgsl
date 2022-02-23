model.so <- function(name){
    LIBS <- "-lgsl -lgslcblas -lm"
    CFLAGS <- "-shared -fPIC -Wall -O2"
    so <- sprintf("%s.so",name)
    if (!file.exists(so)){
        system2("gcc",sprintf("%s -o %s %s_gvf.c %s",CFLAGS,so,name,LIBS))
    }
    return(so)
}


test.events <-function(){
    name <- "HarmonicOscillator"
    t <- seq(0,13,length.out=120)
    N <- 3
    ny <- 2
    np <- 3
    y0 <- matrix(c(0,1),nrow=ny,ncol=N)
    p <- matrix(c(1,0,0),nrow=np,ncol=N)
    for (i in 1:N) {
        p[2,i]=(i-1); #/N
    }
    colnames(p) <- c("no firction","tiny amount of friction","lots of friction")
    ## events:
    t.event <- c(1.0,2.0,3.0)
    lt <- length(t.event)
    ## 2-dim and 3-dim identity matrix
    I3 <- diag(1,3,3)
    I2 <- diag(1,2,2)
    ## the state is reset to initial
    M <- array(I2,dim=c(ny,ny,lt))
    K <- array(I3,dim=c(np,np,lt))
    print(dim(M))
    ## y[i] <- v[i]*y[i] + a[i]
    ## p <- B %*% p + b
    all.events=list()
    for (i in 1:N){
        a <- array(c(1,0),dim=c(ny,1,lt));
        b <- array(c(0,i*0.1,0),dim=c(np,1,lt));
        state.tf <- list(A=M,b=a)
        param.tf <- list(A=K,b=b)
        tf <- list(state=state.tf,param=param.tf)
        ev <- list(time=t.event,tf=tf)
        print(ev)
        all.events[[i]]  <- ev
    }
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

transform <- function(lt=1,A=1,b=0){
    n <- nrow(as.matrix(A))
    m <- nrow(as.matrix(b))
    A <- array(A,dim=c(n,n,lt))
    b <- array(b,dim=c(m,1,lt));
    return(list(A=A,b=b))
}

event.tf <- function(t,state.tf,param.tf){
    tf <- list(state=state.tf,param=param.tf)
    return(list(time=t,tf=tf))
}


test.experiments <-function(){
    name <- "HarmonicOscillator"
    t <- seq(0,13,length.out=120)
    ## 2-dim and 3-dim identity matrix
    I3 <- diag(1,3,3)
    I2 <- diag(1,2,2)
    N <- 3    
    event.t <- c(1.0,2.0,3.0)

    state.tf <- transform(length(event.t),I2,c(1,0))

    param.tf <- transform(length(event.t),I3,c(0,0,0))
    no.friction <- list(time=t,parameters=c(1,0,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))

    param.tf <- transform(length(event.t),I3,c(0,0.1,0))
    low.friction <- list(time=t,parameters=c(1,1,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))

    param.tf <- transform(length(event.t),I3,c(0,0.2,0))
    medium.friction <- list(time=t,parameters=c(1,2,0),initial_value=c(0,1),events=event.tf(event.t,state.tf,param.tf))
    
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
