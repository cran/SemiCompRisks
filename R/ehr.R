	

ehr <- function(Cov2 = NULL, Cov3 = NULL, x, xlim = c(0, 10)){
	
if(class(x) == "BayesID"){
	
	if(x$setup$model == "Markov"){
		
		nChain = x$setup$nChain
			
		if(x$setup$type == "semi-parametric"){
				
			nStore <- length(x$chain1$J1.p)
			numSpl <- nStore * nChain			

			s2.p <- x$chain1$time_lambda2
			s3.p <- x$chain1$time_lambda3
	
			s2.n	<- length(s2.p)
			s3.n	<- length(s3.p)	
	
			time <- sort(unique(c(s2.p, s3.p)))
			time <- time[time <= min(max(s2.p), max(s3.p))]
	
			J 	 <- length(time) - 1
	
			lambda2.fin	<- matrix(NA, nStore*nChain, J+1)
			lambda3.fin	<- matrix(NA, nStore*nChain, J+1)
	
			ind2.lam	<- rep(NA, J+1)
			ind3.lam	<- rep(NA, J+1)
	
			ehrSubj		<- matrix(NA, numSpl, J+1)
			
			j = 1
			for(i in 1:length(time)){
				while(time[i] > s2.p[j]){
					j = j+1
				}
				ind2.lam[i] <-  j
			}
			j = 1
			for(i in 1:length(time)){
				while(time[i] > s3.p[j]){
					j = j+1
				}
				ind3.lam[i] <-  j
			}
			

			if(!is.null(Cov2)){
				beta2.p <- x$chain1$beta2.p
				}
			if(!is.null(Cov3)){
				beta3.p <- x$chain1$beta3.p
				}							
			lambda2.fin	<- x$chain1$lambda2.fin
			lambda3.fin	<- x$chain1$lambda3.fin				
			
			if(nChain > 1){
				for(i in 2:nChain){
					nam <- paste("chain", i, sep="")
				if(!is.null(Cov2)){					
					beta2.p <- rbind(beta2.p, x[[nam]]$beta2.p)
					}
				if(!is.null(Cov3)){					
					beta3.p <- rbind(beta3.p, x[[nam]]$beta3.p)
					}
					lambda2.fin <- rbind(lambda2.fin, x[[nam]]$lambda2.fin)
					lambda3.fin <- rbind(lambda3.fin, x[[nam]]$lambda3.fin)											
				}							
			}

	
			for(i in 1:numSpl){
				if(!is.null(Cov2)){
					LP2	<- as.vector(Cov2 %*% beta2.p[i, ])	
					}
				if(is.null(Cov2)){
					LP2	<- 0
					}
				if(!is.null(Cov3)){
					LP3	<- as.vector(Cov3 %*% beta3.p[i, ])	
					}
				if(is.null(Cov3)){
					LP3	<- 0
					}					
		
				bhr	<- exp(lambda3.fin[i, ind2.lam] - lambda2.fin[i, ind2.lam]) 
				ehrSubj[i,]	<- bhr * exp(LP3 - LP2)		
				}	
			}
		else{
				
			nStore <- length(x$chain1$alpha1.p)
			numSpl <- nStore * nChain				
			
			time <- seq(xlim[1], xlim[2], diff(xlim)/100)
							
			ehrSubj		<- matrix(NA, numSpl, length(time))
			
			
			if(!is.null(Cov2)){
				beta2.p <- x$chain1$beta2.p
				}
			if(!is.null(Cov3)){
				beta3.p <- x$chain1$beta3.p
				}							
			alpha2.p	<- x$chain1$alpha2.p
			alpha3.p	<- x$chain1$alpha3.p
			kappa2.p	<- x$chain1$kappa2.p
			kappa3.p	<- x$chain1$kappa3.p						
			
			if(nChain > 1){
				for(i in 2:nChain){
					nam <- paste("chain", i, sep="")
				if(!is.null(Cov2)){					
					beta2.p <- rbind(beta2.p, x[[nam]]$beta2.p)
					}
				if(!is.null(Cov3)){					
					beta3.p <- rbind(beta3.p, x[[nam]]$beta3.p)
					}
					alpha2.p	<- c(alpha2.p, x[[nam]]$alpha2.p)
					alpha3.p	<- c(alpha3.p, x[[nam]]$alpha3.p)
					kappa2.p	<- c(kappa2.p, x[[nam]]$kappa2.p)
					kappa3.p	<- c(kappa3.p, x[[nam]]$kappa3.p)										
				}							
			}			
			
	
			for(i in 1:numSpl){
				bhr	<- (alpha3.p[i]*kappa3.p[i]*time^(alpha3.p[i] - 1))/(alpha2.p[i]*kappa2.p[i]*time^(alpha2.p[i] - 1))
				if(!is.null(Cov2)){
					LP2	<- as.vector(Cov2 %*% beta2.p[i, ])	
					}
				if(is.null(Cov2)){
					LP2	<- 0
					}
				if(!is.null(Cov3)){
					LP3	<- as.vector(Cov3 %*% beta3.p[i, ])	
					}
				if(is.null(Cov3)){
					LP3	<- 0
					}						
				
				ehrSubj[i,]	<- bhr * exp(LP3 - LP2)		
				if(xlim[1] == 0){
					ehrSubj[i,1] <- 1	
					}				
				}
				
		}
		ret <- list(ehr = ehrSubj, time = time, type = x$setup$type)	
		class(ret) <- "ehr"
		return(ret)					
	}
	else print("Error: the object is not returned by 'Markov model'.")
	}
else print("Error: the class of the object is not 'BayesID'.")
}
	
		
		
		
		
		
		
	

