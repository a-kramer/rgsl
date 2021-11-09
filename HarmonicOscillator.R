demo<-function(){
    name <- "HarmonicOscillator"
    LIBS <- "-lgsl -lgslcblas -lm"
    CFLAGS <- "-shared -fPIC -Wall -O2"
    so <- sprintf("%s.so",name)
    if (!file.exists(so)){
        system2("gcc",sprintf("%s -o %s %s_gvf.c %s",CFLAGS,so,name,LIBS))
    }
    t <- seq(0,13,length.out=120)
    N <- 6
    ny <- 2
    np <- 3
    y0 <- matrix(c(0,1),nrow=ny,ncol=N)
    p <- matrix(c(1,0,0),nrow=np,ncol=N)
    for (i in 1:N) {
        p[2,i]=(i-1)/N
    }
    ## events:
    t.event <- c(1.0,2.0,3.0)
    lt <- length(t.event)
    ## the state is reset to initial
    M <- array(0,dim=c(ny,ny,lt))
    ## friction is increased
    I3 <- diag(1,3,3)
    K <- array(c(I3,I3,I3),dim=c(np,np,lt))
    print(dim(M))
    ev <- list(time=t.event,A=M,B=K,a=y0,b=c(0,0.1,0))

    all.events=list()
    for (i in 1:N){
        all.events[[i]]  <- ev
    }

    
    if (require("rgsl")){
        y <- r_gsl_odeiv2(name,t=t,y0=y0,p=p,all.events)
        plot(t,y[2,,1],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
        for (l in 1:N) {
            lines(t,y[2,,l],lty=l)
        }
        pty=1:N;
        pty[2:N]<-NA;
        legend("bottomright",legend=sprintf("c=%g",p[2,]),lty=1:N,pch=pty)
    }
}
