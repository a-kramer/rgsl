test_that("sensitivity is correct",{
	file.c <- rgsl.example()
	so <- model.so(file.c)
	model.name <- sub("_gvf.c","",basename(file.c))
	comment(model.name) <- so
	expect_true(file.exists(so))
	nt <- 120
	t_ <- seq(0.0,2.0,length.out=nt)
	## 2-dim and 3-dim identity matrix
	I3 <- diag(1,3,3)
	I2 <- diag(1,2,2)
	M <- 2 # arbitrary number of random parameter vectors
	l <- 2 # number of parameters
	ny <- 2
	paramG <- matrix(rnorm(l*M,mean=c(3.0,0.5),sd=c(0.3,0.03)),nrow=l,ncol=M)
	e1 <- list(time=t_,t0=0.0,input=0.0,initial_value=c(0.0,1.0))

	experiments=list(a=e1)
	N <- length(experiments)

	Y <- r_gsl_odeiv2_outer_sens(model.name,experiments,paramG)
	y <- Y[[1]]
	expect_type(Y,"list")
	expect_true(all(is.finite(Y[[1]]$state)))
	expect_length(Y,N)
	expect_length(y,5)
	expect_equal(names(y),c('state','func','stateSensitivity','funcSensitivity','cpuSeconds'))
	sState <- y$stateSensitivity
	sFunc <- y$funcSensitivity
	expect_true(is.list(sState))
	expect_true(is.list(sFunc))
	expect_length(sState,M)
	expect_length(sFunc,M)
	Delta <- (paramG[,2] - paramG[,1])
	expect_gt(norm(Delta,'2'),0.0)
	expect_lt(norm(Delta,'2'),1.0)
	predict_y_p2 <- y$state[,,1]
	true_y_p2 <- y$state[,,2]
	for (j in seq(nt)){
	 predict_y_p2[,j] <- y$state[,j,1] + sState[[1]][,seq(l),j] %*% Delta
	}
	expect_lt(norm(true_y_p2 - predict_y_p2,"2")/nt,5e-1)

	predict_f_p2 <- matrix(y$func[,,1],1,nt)
	predict_y_p2 <- matrix(y$state[,,1],ny,nt)
	true_f_p2 <- matrix(y$func[,,2],1,nt)
	for (j in seq(nt)){
		sy <- matrix(sState[[1]][,seq(l),j],ny,l)
		sf <- matrix(sFunc[[1]][,seq(l),j],1,l)
		predict_f_p2[,j] <- y$func[,j,1] + sf %*% Delta
		predict_y_p2[,j] <- y$state[,j,1] + sy %*% Delta
	}
	expect_lt(norm(true_f_p2 - predict_f_p2,"2")/nt,5e-1)
	dev.new()
	plot(t_,true_f_p2,main="prediction of dots via funcSensitivity",xlab='t',ylab='func',pch=1)
	lines(t_,predict_f_p2,lty=1)
	lines(t_,y$func[,,1],lty=2)
	legend('bottomright',c('true values','predicted','basis of prediction'),lty=c(NA,1,2),pch=c(1,NA,NA))
	dev.new()
	plot(t_,true_y_p2[1,],main="prediction of dots via stateSensitivity",xlab='t',ylab='state 1',pch=1)
	lines(t_,predict_y_p2[1,],lty=1)
	lines(t_,y$state[1,,1],lty=2)
	legend('bottomright',c('true values','predicted','basis of prediction'),lty=c(NA,1,2),pch=c(1,NA,NA))
})
