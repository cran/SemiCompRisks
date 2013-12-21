








summary.BayesID <- function(object, ...)
{
	x <- object
	
	nChain = x$setup$nChain

	cat("\n Type of analysis: ", x$setup$type,"\n")
	cat("\n Type of model: ", x$setup$model,"\n")
	
	##
	
if(length(x$chain1$beta1.p) != 0){
	
	cat("\n Regression parameter (beta1): posterior median (PM)\n")	
	
	p1	= dim(x$chain1$beta1.p)[2]	
	beta.p <- x$chain1$beta1.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			beta.p <- rbind(beta.p, x[[nam]]$beta1.p)											
			}							
		}	
	
		
	beta.pMed <- apply(beta.p, 2, median)
	beta.pUb <- apply(beta.p, 2, quantile, prob = 0.975)		
	beta.pLb <- apply(beta.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, p1, 4)
	rownames(tbl) <- x$chain1$covNames
	colnames(tbl) <- c("PM", "exp(PM)", "lower .975", "upper .975")	
			
	tbl[,1]	<- beta.pMed
	tbl[,2]	<- exp(beta.pMed)
	tbl[,3]	<- exp(beta.pLb)	
	tbl[,4]	<- exp(beta.pUb)
		
	print(round(tbl, 2))	
}

if(length(x$chain1$beta2.p) != 0){
	
	cat("\n Regression parameter (beta2): posterior median (PM)\n")	
	
	p2	= dim(x$chain1$beta2.p)[2]	
	beta.p <- x$chain1$beta2.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			beta.p <- rbind(beta.p, x[[nam]]$beta2.p)											
			}							
		}	
	
		
	beta.pMed <- apply(beta.p, 2, median)
	beta.pUb <- apply(beta.p, 2, quantile, prob = 0.975)		
	beta.pLb <- apply(beta.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, p2, 4)
	rownames(tbl) <- x$chain1$covNames
	colnames(tbl) <- c("PM", "exp(PM)", "lower .975", "upper .975")	
			
	tbl[,1]	<- beta.pMed
	tbl[,2]	<- exp(beta.pMed)
	tbl[,3]	<- exp(beta.pLb)	
	tbl[,4]	<- exp(beta.pUb)
		
	print(round(tbl, 2))	
}


if(length(x$chain1$beta3.p) != 0){
	
	cat("\n Regression parameter (beta3): posterior median (PM)\n")	
	
	p3	= dim(x$chain1$beta3.p)[2]	
	beta.p <- x$chain1$beta3.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			beta.p <- rbind(beta.p, x[[nam]]$beta3.p)											
			}							
		}	
	
		
	beta.pMed <- apply(beta.p, 2, median)
	beta.pUb <- apply(beta.p, 2, quantile, prob = 0.975)		
	beta.pLb <- apply(beta.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, p3, 4)
	rownames(tbl) <- x$chain1$covNames
	colnames(tbl) <- c("PM", "exp(PM)", "lower .975", "upper .975")	
			
	tbl[,1]	<- beta.pMed
	tbl[,2]	<- exp(beta.pMed)
	tbl[,3]	<- exp(beta.pLb)	
	tbl[,4]	<- exp(beta.pUb)
		
	print(round(tbl, 2))	
}

	cat("\n theta: posterior median (PM)\n")

	theta.p <- x$chain1$theta.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			theta.p <- rbind(theta.p, x[[nam]]$theta.p)											
			}							
		}	
	
		
	theta.pMed <- apply(theta.p, 2, median)
	theta.pUb <- apply(theta.p, 2, quantile, prob = 0.975)		
	theta.pLb <- apply(theta.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- theta.pMed
	tbl[,2]	<- theta.pLb	
	tbl[,3]	<- theta.pUb
	
	print(round(tbl, 2))	
	

if(x$setup$type == "parametric"){

	cat("\n alpha1: posterior median (PM)\n")

	alpha.p <- x$chain1$alpha1.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			alpha.p <- rbind(alpha.p, x[[nam]]$alpha1.p)											
			}							
		}	
	
		
	alpha.pMed <- apply(alpha.p, 2, median)
	alpha.pUb <- apply(alpha.p, 2, quantile, prob = 0.975)		
	alpha.pLb <- apply(alpha.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- alpha.pMed
	tbl[,2]	<- alpha.pLb	
	tbl[,3]	<- alpha.pUb
	
	print(round(tbl, 2))
	
	cat("\n alpha2: posterior median (PM)\n")	
	
	alpha.p <- x$chain1$alpha2.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			alpha.p <- rbind(alpha.p, x[[nam]]$alpha2.p)											
			}							
		}	
	
		
	alpha.pMed <- apply(alpha.p, 2, median)
	alpha.pUb <- apply(alpha.p, 2, quantile, prob = 0.975)		
	alpha.pLb <- apply(alpha.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- alpha.pMed
	tbl[,2]	<- alpha.pLb	
	tbl[,3]	<- alpha.pUb
	
	print(round(tbl, 2))
	

	cat("\n alpha3: posterior median (PM)\n")
	
	alpha.p <- x$chain1$alpha3.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			alpha.p <- rbind(alpha.p, x[[nam]]$alpha3.p)											
			}							
		}	
	
		
	alpha.pMed <- apply(alpha.p, 2, median)
	alpha.pUb <- apply(alpha.p, 2, quantile, prob = 0.975)		
	alpha.pLb <- apply(alpha.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- alpha.pMed
	tbl[,2]	<- alpha.pLb	
	tbl[,3]	<- alpha.pUb
	
	print(round(tbl, 2))	
				
				
				

	cat("\n kappa1: posterior median (PM)\n")

	kappa.p <- x$chain1$kappa1.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			kappa.p <- rbind(kappa.p, x[[nam]]$kappa1.p)											
			}							
		}	
		
	kappa.pMed <- apply(kappa.p, 2, median)
	kappa.pUb <- apply(kappa.p, 2, quantile, prob = 0.975)		
	kappa.pLb <- apply(kappa.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- kappa.pMed
	tbl[,2]	<- kappa.pLb	
	tbl[,3]	<- kappa.pUb
	
	print(round(tbl, 2))	
	
	

	cat("\n kappa2: posterior median (PM)\n")

	kappa.p <- x$chain1$kappa2.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			kappa.p <- rbind(kappa.p, x[[nam]]$kappa2.p)											
			}							
		}	
		
	kappa.pMed <- apply(kappa.p, 2, median)
	kappa.pUb <- apply(kappa.p, 2, quantile, prob = 0.975)		
	kappa.pLb <- apply(kappa.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- kappa.pMed
	tbl[,2]	<- kappa.pLb	
	tbl[,3]	<- kappa.pUb
	
	print(round(tbl, 2))	
	
	

	cat("\n kappa3: posterior median (PM)\n")

	kappa.p <- x$chain1$kappa3.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			kappa.p <- rbind(kappa.p, x[[nam]]$kappa3.p)											
			}							
		}	
		
	kappa.pMed <- apply(kappa.p, 2, median)
	kappa.pUb <- apply(kappa.p, 2, quantile, prob = 0.975)		
	kappa.pLb <- apply(kappa.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- kappa.pMed
	tbl[,2]	<- kappa.pLb	
	tbl[,3]	<- kappa.pUb
	
	print(round(tbl, 2))	

}


if(x$setup$type == "semi-parametric"){
	
	cat("\n mu_lam1: posterior median (PM)\n")

	mu_lam.p <- x$chain1$mu_lam1.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			mu_lam.p <- rbind(mu_lam.p, x[[nam]]$mu_lam1.p)											
			}							
		}	
	
		
	mu_lam.pMed <- apply(mu_lam.p, 2, median)
	mu_lam.pUb <- apply(mu_lam.p, 2, quantile, prob = 0.975)		
	mu_lam.pLb <- apply(mu_lam.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- mu_lam.pMed
	tbl[,2]	<- mu_lam.pLb	
	tbl[,3]	<- mu_lam.pUb
	
	print(round(tbl, 2))	
	
	
	cat("\n mu_lam2: posterior median (PM)\n")

	mu_lam.p <- x$chain1$mu_lam2.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			mu_lam.p <- rbind(mu_lam.p, x[[nam]]$mu_lam2.p)											
			}							
		}	
	
		
	mu_lam.pMed <- apply(mu_lam.p, 2, median)
	mu_lam.pUb <- apply(mu_lam.p, 2, quantile, prob = 0.975)		
	mu_lam.pLb <- apply(mu_lam.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- mu_lam.pMed
	tbl[,2]	<- mu_lam.pLb	
	tbl[,3]	<- mu_lam.pUb
	
	print(round(tbl, 2))	
	
	
	cat("\n mu_lam3: posterior median (PM)\n")

	mu_lam.p <- x$chain1$mu_lam3.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			mu_lam.p <- rbind(mu_lam.p, x[[nam]]$mu_lam3.p)											
			}							
		}	
	
		
	mu_lam.pMed <- apply(mu_lam.p, 2, median)
	mu_lam.pUb <- apply(mu_lam.p, 2, quantile, prob = 0.975)		
	mu_lam.pLb <- apply(mu_lam.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- mu_lam.pMed
	tbl[,2]	<- mu_lam.pLb	
	tbl[,3]	<- mu_lam.pUb
	
	print(round(tbl, 2))	
	
			
	
	cat("\n sigSq_lam1: posterior median (PM)\n")

	sigSq_lam.p <- x$chain1$sigSq_lam1.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			sigSq_lam.p <- rbind(sigSq_lam.p, x[[nam]]$sigSq_lam1.p)											
			}							
		}	
	
		
	sigSq_lam.pMed <- apply(sigSq_lam.p, 2, median)
	sigSq_lam.pUb <- apply(sigSq_lam.p, 2, quantile, prob = 0.975)		
	sigSq_lam.pLb <- apply(sigSq_lam.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- sigSq_lam.pMed
	tbl[,2]	<- sigSq_lam.pLb	
	tbl[,3]	<- sigSq_lam.pUb
	
	print(round(tbl, 2))
	
	cat("\n sigSq_lam2: posterior median (PM)\n")

	sigSq_lam.p <- x$chain1$sigSq_lam2.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			sigSq_lam.p <- rbind(sigSq_lam.p, x[[nam]]$sigSq_lam2.p)											
			}							
		}	
	
		
	sigSq_lam.pMed <- apply(sigSq_lam.p, 2, median)
	sigSq_lam.pUb <- apply(sigSq_lam.p, 2, quantile, prob = 0.975)		
	sigSq_lam.pLb <- apply(sigSq_lam.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- sigSq_lam.pMed
	tbl[,2]	<- sigSq_lam.pLb	
	tbl[,3]	<- sigSq_lam.pUb
	
	print(round(tbl, 2))	
	
	cat("\n sigSq_lam3: posterior median (PM)\n")

	sigSq_lam.p <- x$chain1$sigSq_lam3.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			sigSq_lam.p <- rbind(sigSq_lam.p, x[[nam]]$sigSq_lam3.p)											
			}							
		}	
	
		
	sigSq_lam.pMed <- apply(sigSq_lam.p, 2, median)
	sigSq_lam.pUb <- apply(sigSq_lam.p, 2, quantile, prob = 0.975)		
	sigSq_lam.pLb <- apply(sigSq_lam.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- sigSq_lam.pMed
	tbl[,2]	<- sigSq_lam.pLb	
	tbl[,3]	<- sigSq_lam.pUb
	
	print(round(tbl, 2))	

	cat("\n J1: posterior median (PM)\n")

	J.p <- x$chain1$J1.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			J.p <- rbind(J.p, x[[nam]]$J1.p)											
			}							
		}	
	
		
	J.pMed <- apply(J.p, 2, median)
	J.pUb <- apply(J.p, 2, quantile, prob = 0.975)		
	J.pLb <- apply(J.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- J.pMed
	tbl[,2]	<- J.pLb	
	tbl[,3]	<- J.pUb
	
	print(round(tbl, 2))
	
	cat("\n J2: posterior median (PM)\n")

	J.p <- x$chain1$J2.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			J.p <- rbind(J.p, x[[nam]]$J2.p)											
			}							
		}	
	
		
	J.pMed <- apply(J.p, 2, median)
	J.pUb <- apply(J.p, 2, quantile, prob = 0.975)		
	J.pLb <- apply(J.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- J.pMed
	tbl[,2]	<- J.pLb	
	tbl[,3]	<- J.pUb
	
	print(round(tbl, 2))	
	
	cat("\n J3: posterior median (PM)\n")

	J.p <- x$chain1$J3.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			J.p <- rbind(J.p, x[[nam]]$J3.p)											
			}							
		}	
	
		
	J.pMed <- apply(J.p, 2, median)
	J.pUb <- apply(J.p, 2, quantile, prob = 0.975)		
	J.pLb <- apply(J.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- J.pMed
	tbl[,2]	<- J.pLb	
	tbl[,3]	<- J.pUb
	
	print(round(tbl, 2))	
}

}






summary.BayesSurv <- function(object, ...)
{
	x <- object
	
	nChain = x$setup$nChain

	cat("\n Type of analysis: ", x$setup$type,"\n")

if(length(x$chain1$beta.p) != 0){
	
	cat("\n Regression parameter (beta): posterior median (PM)\n")	
	
	p	= dim(x$chain1$beta.p)[2]	
	beta.p <- x$chain1$beta.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			beta.p <- rbind(beta.p, x[[nam]]$beta.p)											
			}							
		}	
	
		
	beta.pMed <- apply(beta.p, 2, median)
	beta.pUb <- apply(beta.p, 2, quantile, prob = 0.975)		
	beta.pLb <- apply(beta.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, p, 4)
	rownames(tbl) <- x$chain1$covNames
	colnames(tbl) <- c("PM", "exp(PM)", "lower .975", "upper .975")	
			
	tbl[,1]	<- beta.pMed
	tbl[,2]	<- exp(beta.pMed)
	tbl[,3]	<- exp(beta.pLb)	
	tbl[,4]	<- exp(beta.pUb)
		
	print(round(tbl, 2))	
}

if(x$setup$type == "parametric"){

	cat("\n alpha: posterior median (PM)\n")

	alpha.p <- x$chain1$alpha.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			alpha.p <- rbind(alpha.p, x[[nam]]$alpha.p)											
			}							
		}	
	
		
	alpha.pMed <- apply(alpha.p, 2, median)
	alpha.pUb <- apply(alpha.p, 2, quantile, prob = 0.975)		
	alpha.pLb <- apply(alpha.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- alpha.pMed
	tbl[,2]	<- alpha.pLb	
	tbl[,3]	<- alpha.pUb
	
	print(round(tbl, 2))	

	cat("\n kappa: posterior median (PM)\n")

	kappa.p <- x$chain1$kappa.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			kappa.p <- rbind(kappa.p, x[[nam]]$kappa.p)											
			}							
		}	
		
	kappa.pMed <- apply(kappa.p, 2, median)
	kappa.pUb <- apply(kappa.p, 2, quantile, prob = 0.975)		
	kappa.pLb <- apply(kappa.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- kappa.pMed
	tbl[,2]	<- kappa.pLb	
	tbl[,3]	<- kappa.pUb
	
	print(round(tbl, 2))	

}


if(x$setup$type == "semi-parametric"){
	
	cat("\n mu_lam: posterior median (PM)\n")

	mu_lam.p <- x$chain1$mu_lam.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			mu_lam.p <- rbind(mu_lam.p, x[[nam]]$mu_lam.p)											
			}							
		}	
	
		
	mu_lam.pMed <- apply(mu_lam.p, 2, median)
	mu_lam.pUb <- apply(mu_lam.p, 2, quantile, prob = 0.975)		
	mu_lam.pLb <- apply(mu_lam.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- mu_lam.pMed
	tbl[,2]	<- mu_lam.pLb	
	tbl[,3]	<- mu_lam.pUb
	
	print(round(tbl, 2))	
	
	cat("\n sigSq_lam: posterior median (PM)\n")

	sigSq_lam.p <- x$chain1$sigSq_lam.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			sigSq_lam.p <- rbind(sigSq_lam.p, x[[nam]]$sigSq_lam.p)											
			}							
		}	
	
		
	sigSq_lam.pMed <- apply(sigSq_lam.p, 2, median)
	sigSq_lam.pUb <- apply(sigSq_lam.p, 2, quantile, prob = 0.975)		
	sigSq_lam.pLb <- apply(sigSq_lam.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- sigSq_lam.pMed
	tbl[,2]	<- sigSq_lam.pLb	
	tbl[,3]	<- sigSq_lam.pUb
	
	print(round(tbl, 2))

	cat("\n J: posterior median (PM)\n")

	J.p <- x$chain1$J.p

	if(nChain > 1){
		for(i in 2:nChain){
			nam <- paste("chain", i, sep="")
			J.p <- rbind(J.p, x[[nam]]$J.p)											
			}							
		}	
	
		
	J.pMed <- apply(J.p, 2, median)
	J.pUb <- apply(J.p, 2, quantile, prob = 0.975)		
	J.pLb <- apply(J.p, 2, quantile, prob = 0.025)
		
	tbl <- matrix(NA, 1, 3)
	colnames(tbl) <- c("PM", "lower .975", "upper .975")	
			
	tbl[,1]	<- J.pMed
	tbl[,2]	<- J.pLb	
	tbl[,3]	<- J.pUb
	
	print(round(tbl, 2))
}

}









print.BayesID <- function(x, ...)
{
	nChain = x$setup$nChain
	
  ##
	cat("\n Type of analysis: ", x$setup$type,"\n")
	
	if(x$setup$type == "semi-parametric"){
		cat("\n Type of model: ", x$setup$model,"\n")
	}
	
  ##
	cat("\n Number of scans: ", x$setup$numReps,"\n")
  ##
	cat("\n Number of chains: ", nChain,"\n")
  ##
	cat("\n Thinning: ",x$setup$thin,"\n")
  ##
	cat("\n Percentage of burnin: ", x$setup$burninPerc*100, "%\n")

# convergence diagnostics

if(nChain > 1){

	cat("\n ###################################")	
	cat("\n Potential scale reduction fator\n")		

	
	if(length(x$chain1$beta1.p) != 0){
	
	#beta1
		
	p1	= dim(x$chain1$beta1.p)[2]	
	
	psrfBeta <- rep(NA, p1)	
	for(j in 1:p1){
		
		#namPara = paste("beta_", j, sep = "")
		
		beta1 <- x$chain1$beta1[,j]
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			beta1 <- cbind(beta1, x[[nam]]$beta1[,j])
		}
		psrfBeta[j] <- calcPSR(beta1)
	}
		cat("\n beta1: p1 =", p1, " elements \n")
		print(round(psrfBeta, 2))		
	}	
	
	if(length(x$chain1$beta2.p) != 0){
	
	#beta2
		
	p2	= dim(x$chain1$beta2.p)[2]	
	
	psrfBeta <- rep(NA, p2)	
	for(j in 1:p2){
		
		#namPara = paste("beta_", j, sep = "")
		
		beta2 <- x$chain1$beta2[,j]
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			beta2 <- cbind(beta2, x[[nam]]$beta2[,j])
		}
		psrfBeta[j] <- calcPSR(beta2)
	}
		cat("\n beta2: p2 =", p2, " elements \n")
		print(round(psrfBeta, 2))		
	}
	
	if(length(x$chain1$beta3.p) != 0){
	
	#beta3
		
	p3	= dim(x$chain1$beta3.p)[2]	
	
	psrfBeta <- rep(NA, p3)	
	for(j in 1:p3){
		
		#namPara = paste("beta_", j, sep = "")
		
		beta3 <- x$chain1$beta3[,j]
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			beta3 <- cbind(beta3, x[[nam]]$beta3[,j])
		}
		psrfBeta[j] <- calcPSR(beta3)
	}
		cat("\n beta3: p3 =", p3, " elements \n")
		print(round(psrfBeta, 2))		
	}	
	
		theta <- x$chain1$theta.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			theta <- cbind(theta, x[[nam]]$theta.p)
		}
		psrftheta <- calcPSR(theta)
		
		cat("\n theta: \n")
		print(round(psrftheta, 2))		
		
		
if(x$setup$type == "semi-parametric"){
	
	ntime1  = length(x$chain1$time_lambda1)
	ntime2  = length(x$chain1$time_lambda3)
	ntime3  = length(x$chain1$time_lambda3)			
		
	# lambda's	
	
	psrfLam <- rep(NA, ntime1)
	
	for(j in 1:ntime1){
		
		lambda1 <- x$chain1$lambda1.fin[,j]
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			lambda1 <- cbind(lambda1, x[[nam]]$lambda1.fin[,j])
		}
		psrfLam[j] <- calcPSR(lambda1)
	}
		
		cat("\n lambda1: summary statistics", "\n")
		print(round(summary(psrfLam), 2))

	psrfLam <- rep(NA, ntime1)
	
	for(j in 1:ntime2){
		
		lambda2 <- x$chain1$lambda2.fin[,j]
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			lambda2 <- cbind(lambda2, x[[nam]]$lambda2.fin[,j])
		}
		psrfLam[j] <- calcPSR(lambda2)
	}
		
		cat("\n lambda2: summary statistics", "\n")
		print(round(summary(psrfLam), 2))

	psrfLam <- rep(NA, ntime1)
	
	for(j in 1:ntime3){
		
		lambda3 <- x$chain1$lambda3.fin[,j]
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			lambda3 <- cbind(lambda3, x[[nam]]$lambda3.fin[,j])
		}
		psrfLam[j] <- calcPSR(lambda3)
	}
		
		cat("\n lambda3: summary statistics", "\n")
		print(round(summary(psrfLam), 2))


	# mu_lam
		
		mu <- x$chain1$mu_lam1.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			mu <- cbind(mu, x[[nam]]$mu_lam1.p)
		}
		psrfMu <- calcPSR(mu)
		
		cat("\n mu_lam1: \n")
		print(round(psrfMu, 2))		
		
		mu <- x$chain1$mu_lam2.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			mu <- cbind(mu, x[[nam]]$mu_lam2.p)
		}
		psrfMu <- calcPSR(mu)
		
		cat("\n mu_lam2: \n")
		print(round(psrfMu, 2))	
		
		mu <- x$chain1$mu_lam3.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			mu <- cbind(mu, x[[nam]]$mu_lam3.p)
		}
		psrfMu <- calcPSR(mu)
		
		cat("\n mu_lam3: \n")
		print(round(psrfMu, 2))					
		
	# sigSq_lam
		
		sig <- x$chain1$sigSq_lam1.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			sig <- cbind(sig, x[[nam]]$sigSq_lam1.p)
		}
		psrfSig <- calcPSR(sig)
		
		cat("\n SigmaSq_lam1: \n")
		print(round(psrfSig, 2))		
		
		sig <- x$chain1$sigSq_lam2.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			sig <- cbind(sig, x[[nam]]$sigSq_lam2.p)
		}
		psrfSig <- calcPSR(sig)
		
		cat("\n SigmaSq_lam2: \n")
		print(round(psrfSig, 2))	
		
		sig <- x$chain1$sigSq_lam3.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			sig <- cbind(sig, x[[nam]]$sigSq_lam3.p)
		}
		psrfSig <- calcPSR(sig)
		
		cat("\n SigmaSq_lam3: \n")
		print(round(psrfSig, 2))					
		
	# J
		
		J <- x$chain1$J1.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			J <- cbind(J, x[[nam]]$J1.p)
		}
		psrfJ <- calcPSR(J)
		
		cat("\n J1: \n")
		print(round(psrfJ, 2))	
		
		J <- x$chain1$J2.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			J <- cbind(J, x[[nam]]$J2.p)
		}
		psrfJ <- calcPSR(J)
		
		cat("\n J2: \n")
		print(round(psrfJ, 2))		
		
		J <- x$chain1$J3.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			J <- cbind(J, x[[nam]]$J3.p)
		}
		psrfJ <- calcPSR(J)
		
		cat("\n J3: \n")
		print(round(psrfJ, 2))		
		
} # type = semi-parametric	

if(x$setup$type == "parametric"){
	
	# alpha
		
		alpha <- x$chain1$alpha1.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			alpha <- cbind(alpha, x[[nam]]$alpha1.p)
		}
		psrfAlpha <- calcPSR(alpha)
		
		cat("\n alpha1: \n")
		print(round(psrfAlpha, 2))	
		
		alpha <- x$chain1$alpha2.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			alpha <- cbind(alpha, x[[nam]]$alpha2.p)
		}
		psrfAlpha <- calcPSR(alpha)
		
		cat("\n alpha2: \n")
		print(round(psrfAlpha, 2))		
		
		alpha <- x$chain1$alpha3.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			alpha <- cbind(alpha, x[[nam]]$alpha3.p)
		}
		psrfAlpha <- calcPSR(alpha)
		
		cat("\n alpha3: \n")
		print(round(psrfAlpha, 2))					

	# kappa
		
		kappa <- x$chain1$kappa1.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			kappa <- cbind(kappa, x[[nam]]$kappa1.p)
		}
		psrfKappa <- calcPSR(kappa)
		
		cat("\n kappa1: \n")
		print(round(psrfKappa, 2))
		
		kappa <- x$chain1$kappa2.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			kappa <- cbind(kappa, x[[nam]]$kappa2.p)
		}
		psrfKappa <- calcPSR(kappa)
		
		cat("\n kappa2: \n")
		print(round(psrfKappa, 2))	
		
		kappa <- x$chain1$kappa3.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			kappa <- cbind(kappa, x[[nam]]$kappa3.p)
		}
		psrfKappa <- calcPSR(kappa)
		
		cat("\n kappa3: \n")
		print(round(psrfKappa, 2))			
		
} # type = parametric			

	cat("\n ###################################")

}
else if(nChain == 1){
	cat("The potential scale reduction fator cannot be printed. \n")	
	cat("The number of chains must be larger than 1. \n")	
}

	invisible()
}













print.BayesSurv <- function(x, ...)
{
	nChain = x$setup$nChain
	
  ##
	cat("\n Type of analysis: ", x$setup$type,"\n")
  ##
	cat("\n Number of scans: ", x$setup$numReps,"\n")
  ##
	cat("\n Number of chains: ", nChain,"\n")
  ##
	cat("\n Thinning: ",x$setup$thin,"\n")
  ##
	cat("\n Percentage of burnin: ", x$setup$burninPerc*100, "%\n")

# convergence diagnostics

if(nChain > 1){

	cat("\n ###################################")	
	cat("\n Potential scale reduction fator\n")		

	
	if(length(x$chain1$beta.p) != 0){
	
	#beta	
		
	p	= dim(x$chain1$beta.p)[2]	
	
	psrfBeta <- rep(NA, p)	
	for(j in 1:p){
		
		#namPara = paste("beta_", j, sep = "")
		
		beta <- x$chain1$beta[,j]
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			beta <- cbind(beta, x[[nam]]$beta[,j])
		}
		psrfBeta[j] <- calcPSR(beta)
	}
		cat("\n beta: p =",p, " elements \n")
		print(round(psrfBeta, 2))		
	}	
		
if(x$setup$type == "semi-parametric"){
	
	ntime  = length(x$chain1$time_lambda)	
		
	# lambda's	
	
	psrfLam <- rep(NA, ntime)
	
	for(j in 1:ntime){
		
		namPara = paste("beta_", j, sep = "")
		
		lambda <- x$chain1$lambda.fin[,j]
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			lambda <- cbind(lambda, x[[nam]]$lambda.fin[,j])
		}
		psrfLam[j] <- calcPSR(lambda)
	}
		
		cat("\n lambda: summary statistics", "\n")
		print(round(summary(psrfLam), 2))


	# mu_lam
		
		mu <- x$chain1$mu_lam.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			mu <- cbind(mu, x[[nam]]$mu_lam.p)
		}
		psrfMu <- calcPSR(mu)
		
		cat("\n mu_lam: \n")
		print(round(psrfMu, 2))		
		
	# sigSq_lam
		
		sig <- x$chain1$sigSq_lam.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			sig <- cbind(sig, x[[nam]]$sigSq_lam.p)
		}
		psrfSig <- calcPSR(sig)
		
		cat("\n SigmaSq_lam: \n")
		print(round(psrfSig, 2))		
		
	# J
		
		J <- x$chain1$J.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			J <- cbind(J, x[[nam]]$J.p)
		}
		psrfJ <- calcPSR(J)
		
		cat("\n J: \n")
		print(round(psrfJ, 2))	
		
} # type = semi-parametric	

if(x$setup$type == "parametric"){
	
	# alpha
		
		alpha <- x$chain1$alpha.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			alpha <- cbind(alpha, x[[nam]]$alpha.p)
		}
		psrfAlpha <- calcPSR(alpha)
		
		cat("\n alpha: \n")
		print(round(psrfAlpha, 2))			

	# kappa
		
		kappa <- x$chain1$kappa.p
		for(i in 2:nChain){
			nam <- paste("chain", i, sep = "")
			kappa <- cbind(kappa, x[[nam]]$kappa.p)
		}
		psrfKappa <- calcPSR(kappa)
		
		cat("\n kappa: \n")
		print(round(psrfKappa, 2))
		
} # type = parametric			

	cat("\n ###################################")

}
else if(nChain == 1){
	cat("The potential scale reduction fator cannot be printed. \n")	
	cat("The number of chains must be larger than 1. \n")		
}

	invisible()
}






plot.BayesSurv <- function(x, tlim = c(0, 10), ...){
	
nChain = x$setup$nChain	
	
if(x$setup$type == "semi-parametric"){
	
	lambda.fin	<- x$chain1$lambda.fin	
	
	if(nChain > 1){
		for(i in 2:nChain){		
			nam <- paste("chain", i, sep="")			
			lambda.fin <- rbind(lambda.fin, x[[nam]]$lambda.fin)													
			}							
		}

	time <- x$chain1$time_lambda
		
	basehazMed <- exp(apply(lambda.fin, 2, median))
	basehazUb <- exp(apply(lambda.fin, 2, quantile, prob = 0.975))
	basehazLb <- exp(apply(lambda.fin, 2, quantile, prob = 0.025))
	
	plot(time, basehazMed, ylim = c(min(basehazLb), max(basehazUb)), 
		type = "l", xlim = c(0, max(time)),
		main = expression(paste("Estimates (posterior median) of ", h[0](t), "")),
		ylab = "Baseline hazard", xlab = "Time")
	lines(time, basehazUb, lty = 2)
	lines(time, basehazLb, lty = 2)	
}	


if(x$setup$type == "parametric"){
	
	time <- seq(tlim[1], tlim[2], diff(tlim)/100)
	nStore <- length(x$chain1$alpha.p)
	numSpl <- nStore * nChain				
	
	basehaz <- matrix(NA, numSpl, length(time))
	
	alpha.p	<- x$chain1$alpha.p
	kappa.p	<- x$chain1$kappa.p	
	
	if(nChain > 1){
		for(i in 2:nChain){		
			nam <- paste("chain", i, sep="")			
			alpha.p <- c(alpha.p, x[[nam]]$alpha.p)													
			kappa.p <- c(alpha.p, x[[nam]]$kappa.p)			
			}							
		}		
	for(i in 1:numSpl){
		basehaz[i, ] <- alpha.p[i] * kappa.p[i] * time^(alpha.p[i] - 1)
	}
	if(tlim[1] == 0){
		basehaz[,1] <- 0	
	}
	
	basehazMed <- apply(basehaz, 2, median)
	basehazUb <- apply(basehaz, 2, quantile, prob = 0.975)	
	basehazLb <- apply(basehaz, 2, quantile, prob = 0.025)
	
	plot(time, basehazMed, ylim = c(min(basehazLb), max(basehazUb)), 
		type = "l",
		main = expression(paste("Estimates (posterior median) of ", h[0](t), "")),
		ylab = "Baseline hazard", xlab = "Time")
	lines(time, basehazUb, lty = 2)
	lines(time, basehazLb, lty = 2)
}		
		
}










plot.BayesID <- function(x, tlim = c(0, 10), ...){	
	
nChain = x$setup$nChain	
		
if(x$setup$type == "semi-parametric"){
	
	time1 <- x$chain1$time_lambda1
	time2 <- x$chain1$time_lambda2
	time3 <- x$chain1$time_lambda3	
	
	lambda1.fin	<- x$chain1$lambda1.fin
	lambda2.fin	<- x$chain1$lambda2.fin
	lambda3.fin	<- x$chain1$lambda3.fin				
	
	if(nChain > 1){
		for(i in 2:nChain){		
			nam <- paste("chain", i, sep="")			
			lambda1.fin <- rbind(lambda1.fin, x[[nam]]$lambda1.fin)
			lambda2.fin <- rbind(lambda2.fin, x[[nam]]$lambda2.fin)
			lambda3.fin <- rbind(lambda3.fin, x[[nam]]$lambda3.fin)
			}							
		}	

	basehaz1Med <- exp(apply(lambda1.fin, 2, median))
	basehaz1Ub <- exp(apply(lambda1.fin, 2, quantile, prob = 0.975))
	basehaz1Lb <- exp(apply(lambda1.fin, 2, quantile, prob = 0.025))		

	basehaz2Med <- exp(apply(lambda2.fin, 2, median))
	basehaz2Ub <- exp(apply(lambda2.fin, 2, quantile, prob = 0.975))
	basehaz2Lb <- exp(apply(lambda2.fin, 2, quantile, prob = 0.025))		
	
	basehaz3Med <- exp(apply(lambda3.fin, 2, median))
	basehaz3Ub <- exp(apply(lambda3.fin, 2, quantile, prob = 0.975))
	basehaz3Lb <- exp(apply(lambda3.fin, 2, quantile, prob = 0.025))
				
	ylim = c(min(basehaz1Lb, basehaz2Lb, basehaz3Lb), max(basehaz1Ub, basehaz2Ub, basehaz3Ub))
	
	par(mfrow = c(1,3))
	plot(time1, basehaz1Med, ylim = ylim, 
		type = "l",
		main = expression(paste("Estimates (posterior median) of ", h[0][1](t), "")),
		ylab = "Baseline hazard", xlab = "Time")
	lines(time1, basehaz1Ub, lty = 2)
	lines(time1, basehaz1Lb, lty = 2)	
		
	plot(time2, basehaz2Med, ylim = ylim, 
		type = "l",
		main = expression(paste("Estimates (posterior median) of ", h[0][2](t), "")),
		ylab = "Baseline hazard", xlab = "Time")
	lines(time2, basehaz2Ub, lty = 2)
	lines(time2, basehaz2Lb, lty = 2)				
	
if(x$setup$model == "Markov"){
			
	plot(time3, basehaz3Med, ylim = ylim, 
		type = "l",
		main = expression(paste("Estimates (posterior median) of ", h[0][3](t), "")),
		ylab = "Baseline hazard", xlab = "Time")
	lines(time3, basehaz3Ub, lty = 2)
	lines(time3, basehaz3Lb, lty = 2)	
}	

if(x$setup$model == "semi-Markov"){
	
	plot(time3, basehaz3Med, ylim = ylim, 
		type = "l",
		main = expression(paste("Estimates (posterior median) of ", h[0][3](t), "")),
		ylab = "Baseline hazard", xlab = "Time since non-terminal event")
	lines(time3, basehaz3Ub, lty = 2)
	lines(time3, basehaz3Lb, lty = 2)	
}		
}	


if(x$setup$type == "parametric"){
	time <- seq(tlim[1], tlim[2], diff(tlim)/100)
	nStore <- length(x$chain1$alpha1.p)
	numSpl <- nStore * nChain
		
	basehaz1 <- matrix(NA, numSpl, length(time))
	basehaz2 <- matrix(NA, numSpl, length(time))
	basehaz3 <- matrix(NA, numSpl, length(time))
	
	alpha1.p	<- x$chain1$alpha1.p
	alpha2.p	<- x$chain1$alpha2.p
	alpha3.p	<- x$chain1$alpha3.p		
	kappa1.p	<- x$chain1$kappa1.p
	kappa2.p	<- x$chain1$kappa2.p
	kappa3.p	<- x$chain1$kappa3.p
	
	if(nChain > 1){
		for(i in 2:nChain){		
			nam <- paste("chain", i, sep="")			
			alpha1.p <- c(alpha1.p, x[[nam]]$alpha1.p)
			alpha2.p <- c(alpha2.p, x[[nam]]$alpha2.p)
			alpha3.p <- c(alpha3.p, x[[nam]]$alpha3.p)
																						
			kappa1.p <- c(kappa1.p, x[[nam]]$kappa1.p)
			kappa2.p <- c(kappa2.p, x[[nam]]$kappa2.p)
			kappa3.p <- c(kappa3.p, x[[nam]]$kappa3.p)									
			}							
		}						
	
	for(i in 1:numSpl){
		basehaz1[i, ] <- alpha1.p[i] * kappa1.p[i] * time^(alpha1.p[i] - 1)
		basehaz2[i, ] <- alpha2.p[i] * kappa2.p[i] * time^(alpha2.p[i] - 1)
		basehaz3[i, ] <- alpha3.p[i] * kappa3.p[i] * time^(alpha3.p[i] - 1)				
	}
	
	if(tlim[1] == 0){
		basehaz1[,1] <- 0	
		basehaz2[,1] <- 0	
		basehaz3[,1] <- 0					
	}
	
	basehaz1Med <- apply(basehaz1, 2, median)
	basehaz1Ub <- apply(basehaz1, 2, quantile, prob = 0.975)	
	basehaz1Lb <- apply(basehaz1, 2, quantile, prob = 0.025)
	basehaz2Med <- apply(basehaz2, 2, median)
	basehaz2Ub <- apply(basehaz2, 2, quantile, prob = 0.975)	
	basehaz2Lb <- apply(basehaz2, 2, quantile, prob = 0.025)
	basehaz3Med <- apply(basehaz3, 2, median)
	basehaz3Ub <- apply(basehaz3, 2, quantile, prob = 0.975)	
	basehaz3Lb <- apply(basehaz3, 2, quantile, prob = 0.025)		
	
	ylim = c(min(basehaz1Lb, basehaz2Lb, basehaz3Lb), max(basehaz1Ub, basehaz2Ub, basehaz3Ub))
	
	par(mfrow = c(1,3))
	plot(time, basehaz1Med, ylim = ylim, 
		type = "l",
		main = expression(paste("Estimates (posterior median) of ", h[0][1](t), "")),
		ylab = "Baseline hazard", xlab = "Time")
	lines(time, basehaz1Ub, lty = 2)
	lines(time, basehaz1Lb, lty = 2)	
		
	plot(time, basehaz2Med, ylim = ylim, 
		type = "l",
		main = expression(paste("Estimates (posterior median) of ", h[0][2](t), "")),
		ylab = "Baseline hazard", xlab = "Time")
	lines(time, basehaz2Ub, lty = 2)
	lines(time, basehaz2Lb, lty = 2)
	
	plot(time, basehaz3Med, ylim = ylim, 
		type = "l",
		main = expression(paste("Estimates (posterior median) of ", h[0][3](t), "")),
		ylab = "Baseline hazard", xlab = "Time")
	lines(time, basehaz3Ub, lty = 2)
	lines(time, basehaz3Lb, lty = 2)	
}		
		
}



plot.ehr <- function(x, tlim = c(0, 10), ...){
	if(x$type == "semi-parametric"){
		time 	<- x$time
	
		ehrMed <- apply(x$ehr, 2, median)
		ehrUb <- apply(x$ehr, 2, quantile, prob = 0.975)	
		ehrLb <- apply(x$ehr, 2, quantile, prob = 0.025)		
	
		plot(time, ehrMed, ylim = c(0, max(ehrUb)), 
			type = "l",
			main = "Estimates (posterior median) of EHR",
			ylab = "Explanatory hazard ratio", xlab = "Time")
		lines(time, ehrUb, lty = 3)
		lines(time, ehrLb, lty = 3)
		abline(h = 1, col = "red", lty = 4)			
	}
	else{
		time 	<- x$time
	
		ehrMed <- apply(x$ehr, 2, median)
		ehrUb <- apply(x$ehr, 2, quantile, prob = 0.975)	
		ehrLb <- apply(x$ehr, 2, quantile, prob = 0.025)		
	
		plot(time, ehrMed, ylim = c(0, max(ehrUb)), 
			type = "l",
			main = "Estimates (posterior median) of EHR",
			ylab = "Explanatory hazard ratio", xlab = "Time")
		lines(time, ehrUb, lty = 3)
		lines(time, ehrLb, lty = 3)
		abline(h = 1, col = "red", lty = 4)				
	}
}








