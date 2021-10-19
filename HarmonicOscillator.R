demo<-function(){
    name <- "HarmonicOscillator"
    LIBS <- "-lgsl -lgslcblas -lm"
    CFLAGS <- "-shared -fPIC"
    so <- sprintf("%s.so",name)
    if (!file.exists(so)){
        system2("gcc",sprintf("%s -o %s %s_gvf.c %s",CFLAGS,so,name,LIBS))
    }
    t <- seq(0,130,length.out=1200)
    y0 <- c(0,1)
    N <- 64
    p <- matrix(c(1,0,0.1),nrow=3,ncol=N)
    for (i in 1:N) {
        p[2,i]=(i-1)/N
    }
    if (require("rgsl")){
        y <- r_gsl_odeiv2(name,t,y0,p)
        plot(t,y[2,,1],main="Harmonic Oscillator",sub="with varying damping c",xlab="time",ylab="y(t) amplitude")
        for (l in 1:N) {
            lines(t,y[2,,l],lty=l)
        }
        legend("bottomright",legend=sprintf("c=%g",p[2,]),lty=1:N)
    }
}
