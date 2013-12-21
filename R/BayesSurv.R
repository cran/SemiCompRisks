

### type = "semi-parametric" or "parametric"

BayesSurv <- function(survData, 
						hyperParams,
						startValues,								
						mcmcParams,
						numReps,
						thin,
						path,
						burninPerc=0.5,
						type = "semi-parametric",
						nChain = 1)
{
	
if(class(startValues) == "list" & length(startValues) == nChain){	
	
	chain = 1
	ret <- list()
	
	while(chain <= nChain){
		
	cat("chain: ", chain, "\n")		
	
	startV <- as.vector(unlist(startValues[chain]))	
	nam = paste("chain", chain, sep="")	
	
	# type = "semi-parametric"
	
	if(type == "semi-parametric"){

		###
		n	<- dim(survData)[1]
		p	<- dim(survData)[2] - 2
	
		J_max = mcmcParams[4]

		###	
		nStore <- numReps/thin * (1 - burninPerc)

		nTime_lambda = mcmcParams[6]
		time_lambda	= mcmcParams[(length(mcmcParams)-nTime_lambda+1):length(mcmcParams)]
	
		mcmc <- .C("BpeSurvmcmc",
					survData            = as.double(as.matrix(survData)),
					n                   = as.integer(n),
					p                   = as.integer(p),
					hyperParams         = as.double(hyperParams),
					startValues         = as.double(startV),
					mcmcParams          = as.double(mcmcParams),
					numReps             = as.integer(numReps),
					thin                = as.integer(thin),
					burninPerc          = as.double(burninPerc),
					samples_beta        = as.double(rep(0, nStore*p)),
					samples_mu_lam      = as.double(rep(0, nStore*1)),
					samples_sigSq_lam	= as.double(rep(0, nStore*1)),
					samples_J           = as.double(rep(0, nStore*1)),
					samples_s           = as.double(rep(0, nStore*(J_max + 1))),
					samples_misc        = as.double(rep(0, p + 2)),
					lambda_fin			= as.double(rep(0, nStore*nTime_lambda)))

		if(p > 0){
     		beta.p 			<- matrix(mcmc$samples_beta, nrow = nStore, byrow = TRUE)
       		}
		if(p == 0){
        	beta.p 			<- NULL
        	}  
		lambda.fin 		<- matrix(mcmc$lambda_fin, nrow = nStore, byrow = TRUE)			
		mu_lam.p 		<- matrix(mcmc$samples_mu_lam, nrow = nStore, byrow = TRUE)
		sigSq_lam.p 	<- matrix(mcmc$samples_sigSq_lam, nrow = nStore, byrow = TRUE)	
		J.p 			<- matrix(mcmc$samples_J, nrow = nStore, byrow = TRUE)
		s.p 			<- matrix(mcmc$samples_s, nrow = nStore, byrow = TRUE)			
    	if(p > 0){
       	 	accept.beta 	<- as.vector(mcmc$samples_misc[1:p])
        	}
    	if(p == 0){
       	 	accept.beta 	<- NULL
        	}
		accept.BI		<- as.vector(mcmc$samples_misc[(p+1)])
		accept.DI		<- as.vector(mcmc$samples_misc[(p+2)])	
	
    	if(p > 0){
			covNames = colnames(survData)[-c(1,2)]
        	}	
    	if(p == 0){
			covNames = NULL
        	}	
		
		ret[[nam]] <- list(beta.p = beta.p, lambda.fin = lambda.fin, mu_lam.p = mu_lam.p, sigSq_lam.p = sigSq_lam.p, 
					J.p = J.p, s.p = s.p, accept.beta = accept.beta, accept.BI = accept.BI, accept.DI = accept.DI,
					covNames = covNames, time_lambda = time_lambda, type = type)

		}


	# type = "parametric"
	
	if(type == "parametric"){

		###
		n	<- dim(survData)[1]
		p	<- dim(survData)[2] - 2

		###
	
		nStore <- numReps/thin * (1 - burninPerc)

		mcmc <- .C("BweibSurvmcmc",
						survData 		= as.double(as.matrix(survData)),
						n				= as.integer(n),
						p				= as.integer(p),
						hyperParams 	= as.double(hyperParams),
						mcmcParams		= as.double(mcmcParams),
						startValues 	= as.double(startV),
						numReps			= as.integer(numReps),
						thin			= as.integer(thin),
						burninPerc      = as.double(burninPerc),
						samples_beta 	= as.double(rep(0, nStore*p)),
						samples_alpha 	= as.double(rep(0, nStore*1)),
						samples_kappa 	= as.double(rep(0, nStore*1)),
						samples_misc	= as.double(rep(0, p + 1)))
    
		if(p > 0){
        	beta.p 			<- matrix(mcmc$samples_beta, nrow = nStore, byrow = TRUE)
        	}
		if(p == 0){
   	     	beta.p 			<- NULL
   	     	}    
		alpha.p 		<- matrix(mcmc$samples_alpha, nrow = nStore, byrow = TRUE)
		kappa.p 		<- matrix(mcmc$samples_kappa, nrow = nStore, byrow = TRUE)
   		if(p > 0){
        	accept.beta 	<- as.vector(mcmc$samples_misc[1:p])
        	}
    	if(p == 0){
        	accept.beta 	<- NULL
        	}
		accept.alpha	<- as.vector(mcmc$samples_misc[(p+1)])
    	if(p > 0){
			covNames = colnames(survData)[-c(1,2)]
        	}	
    	if(p == 0){
			covNames = NULL
        	}
			
		ret[[nam]] <- list(beta.p = beta.p, alpha.p = alpha.p, kappa.p = kappa.p, accept.beta = accept.beta, 
					accept.alpha = accept.alpha, covNames = covNames, type = type)

		
		}
				
	chain = chain + 1	
	}		
	
	ret[["setup"]]	<- list(hyperParams = hyperParams, startValues = startValues, mcmcParams = mcmcParams, numReps = numReps, thin = thin, path = path, burninPerc = burninPerc, type = type, nChain = nChain)	

	class(ret) <- "BayesSurv"
	return(ret)

}
else{
	print("The 'startValues' should be the list of length equal to 'nChain'.")
}

}





























