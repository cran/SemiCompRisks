
	simSurv	<- function(x, beta.true, alpha.true, kappa.true, cens){
		
		n <- dim(x)[1]
		p <- dim(x)[2]		
		
		LP	<- as.vector(beta.true %*% t(x)) 

		T	<- rweibull(n, shape = alpha.true, scale = exp(-(log(kappa.true) + LP)/alpha.true))

		delta <- rep(NA, n)
		y		<- T
		Cen		<- runif(n, cens[1], cens[2])
		
		ind1	<- which(T < Cen)
		y[ind1] <- T[ind1]
		delta[ind1] <- 1
		
		ind0	<- which(T >= Cen)
		y[ind0] <- Cen[ind0]
		delta[ind0] <- 0
		
		ret <- list(y = y, delta = delta)
		class(ret) <- "simSurv"
		return(ret)
				
		}	
	
	
	