demo<-function(){
    name <- "HarmonicOscillator"
    LIBS <- "-lgsl -lgslcblas -lm"
    CFLAGS <- "-shared -fPIC -Wall -O2"
    so <- sprintf("%s.so",name)
    if (!file.exists(so)){
        system2("gcc",sprintf("%s -o %s %s_gvf.c %s",CFLAGS,so,name,LIBS))
    }
    t <- seq(0,13,length.out=120)
    N <- 10

    y0 <- matrix(c(0,1),nrow=2,ncol=N)
    p <- matrix(c(1,0,0.1),nrow=3,ncol=N)
    for (i in 1:N) {
        p[2,i]=(i-1)/N
    }
    ## events:
    t.event <- c(1,2,3)
    lt <- length(t)
    M <- array(c(diag(2),diag(2),diag(2)),dim=c(2,2,3))
    K <- array(c(diag(3),diag(3),diag(3)),dim=c(3,3,3))
    print(dim(M))
    
    if (require("rgsl")){
        y <- r_gsl_odeiv2(name,t=t,y0=y0,p=p,M=M,K=K,t.event=t.event)
        plot(t,y[2,,1],main="Damped Harmonic Oscillator",sub="y'' = -ky -cy' with varying damping c (dy/dt=y')",xlab="time",ylab="state y(t;c)")
        for (l in 1:N) {
            lines(t,y[2,,l],lty=l)
        }
        pty=1:N;
        pty[2:N]<-NA;
        legend("bottomright",legend=sprintf("c=%g",p[2,]),lty=1:N,pch=pty)
    }
}
