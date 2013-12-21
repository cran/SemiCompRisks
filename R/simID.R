
	simID	<- function(x1, x2, x3, beta1.true, beta2.true, beta3.true,
						alpha1.true, alpha2.true, alpha3.true, 
						kappa1.true, kappa2.true, kappa3.true, gamma.true, cens){
		
		n <- dim(x1)[1]
		p1 <- dim(x1)[2]
		p2 <- dim(x2)[2]
		p3 <- dim(x3)[2]						
		
		LP1	<- as.vector(beta1.true %*% t(x1)) 
		LP2	<- as.vector(beta2.true %*% t(x2))
		LP3	<- as.vector(beta3.true %*% t(x3))
		
		##
		
		Rind	<- NULL
		
		# at time = 0,
		
		R	<- rweibull(n, shape = alpha1.true, scale = exp(-(log(kappa1.true) + LP1 + log(gamma.true))/alpha1.true))
		D	<- rweibull(n, shape = alpha2.true, scale = exp(-(log(kappa2.true) + LP2 + log(gamma.true))/alpha2.true))
		
		yesR	<- R < D & R
		D[yesR]	<- R[yesR] + rweibull(sum(yesR), shape = alpha3.true, scale = exp(-(log(kappa3.true) + LP3[yesR] + log(gamma.true[yesR]))/alpha3.true))
			
		delta1 <- rep(NA, n)
		delta2 <- rep(NA, n)
		y1		<- R
		y2		<- D
		Cen		<- runif(n, cens[1], cens[2])
		
		ind01	<- which(D < R & D < Cen)
		y1[ind01] <- D[ind01]
		delta1[ind01] <- 0
		delta2[ind01] <- 1
		
		ind10 <- which(R < D & R < Cen & D >= Cen)
		y2[ind10] <- Cen[ind10]
		delta1[ind10] <- 1
		delta2[ind10] <- 0	
		
		ind00 <- which(R >= Cen & D >= Cen)	
		y1[ind00] <- Cen[ind00]
		y2[ind00] <- Cen[ind00]		
		delta1[ind00] <- 0
		delta2[ind00] <- 0	
		
		ind11 <- which(R < Cen & D < Cen & R < D)
		delta1[ind11] <- 1
		delta2[ind11] <- 1					
		
		ret <- list(y1 = y1, y2 = y2, delta1 = delta1, delta2 = delta2)
		class(ret) <- "simID"
		return(ret)				
		
		}	
	
	
	