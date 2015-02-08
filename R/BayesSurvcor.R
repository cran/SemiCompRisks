

### hz.type = "Weibull" or "PEM"
### re.type = "Normal" or "DPM"



BayesSurvcor <- function(survData, 
						hyperParams,
						startValues,								
						mcmcParams,
						numReps,
						thin,
						path = "results/",
						burninPerc=0.5,
                        hz.type = "Weibull",
                        re.type = "Normal",
                        storeV = TRUE,
						nChain = 1)
{
	
if(class(startValues) == "list" & length(startValues) == nChain){
    
	dir.create(paste(path))
	
	chain = 1
	ret <- list()
	
	while(chain <= nChain){
		
	cat("chain: ", chain, "\n")		
	
	startV <- as.vector(unlist(startValues[chain]))	
	nam = paste("chain", chain, sep="")	
	
    # hz.type = "Weibull"
	
	if(hz.type == "Weibull"){
        
        # re.type = "Normal"
        
        #######################################
        ############ Weibull-Normal ###########
        #######################################
        
        if(re.type == "Normal"){
            
            ###
            n	<- dim(survData)[1]
            p	<- dim(survData)[2] - 3
            
            J	<- length(unique(survData[,3]))
            
            nj	<- rep(NA, J)
            
            for(i in 1:J){
                nj[i]	<- length(which(survData[,3] == i))
            }
            
            ###
            
            nStore <- numReps/thin * (1 - burninPerc)
            
            mcmc <- .C("BweibCorSurvmcmc",
            survData 		= as.double(as.matrix(survData)),
            n				= as.integer(n),
            p				= as.integer(p),
            J				= as.integer(J),
            nj				= as.double(nj),
            hyperParams 	= as.double(hyperParams),
            mcmcParams		= as.double(mcmcParams),
            startValues 	= as.double(startV),
            numReps			= as.integer(numReps),
            thin			= as.integer(thin),
            burninPerc      = as.double(burninPerc),
            samples_beta 	= as.double(rep(0, nStore*p)),
            samples_alpha 	= as.double(rep(0, nStore*1)),
            samples_kappa 	= as.double(rep(0, nStore*1)),
            samples_V		= as.double(rep(0, nStore*J)),
            samples_zeta	= as.double(rep(0, nStore*1)),
            samples_misc	= as.double(rep(0, p+1+J)))

            
            
            if(p > 0){
                beta.p 		<- matrix(mcmc$samples_beta, nrow = nStore, byrow = TRUE)
            }
            if(p == 0){
                beta.p 		<- NULL
            }
            
            
            alpha.p 		<- matrix(mcmc$samples_alpha, nrow = nStore, byrow = TRUE)
            kappa.p 		<- matrix(mcmc$samples_kappa, nrow = nStore, byrow = TRUE)
            V.p             <- matrix(mcmc$samples_V, nrow = nStore, byrow = TRUE)
            zeta.p 			<- matrix(mcmc$samples_zeta, nrow = nStore, byrow = TRUE)
            
            if(p > 0){
                accept.beta 	<- as.vector(mcmc$samples_misc[1:p])
            }
            if(p == 0){
                accept.beta 	<- NULL
            }
            
            accept.alpha	<- as.vector(mcmc$samples_misc[p+1])
            accept.V 	<- as.vector(mcmc$samples_misc[(p+1+1):(p+1+J)])
            
            Vsummary <- as.matrix(apply(V.p, 2, summary))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.975))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.025))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, sd))
            rownames(Vsummary)[7:9] <- c("0.975", "0.025", "sd")	
            
            if(storeV == TRUE)
            {
                save(V.p, file = paste(path, "VPch", chain, ".RData", sep = ""))
            }
            
            
            if(p > 0){
                covNames = colnames(survData)[c(4:(3+p))]
            }	
            if(p == 0){
                covNames = NULL
            }
            
            ret[[nam]] <- list(beta.p = beta.p, alpha.p = alpha.p, kappa.p = kappa.p, zeta.p = zeta.p, accept.beta = accept.beta, accept.alpha = accept.alpha, accept.V = accept.V, Vsum = Vsummary, covNames = covNames)

		
		} ## end: if Weibull-Normal
        
        # re.type = "DPM"
        
        #################################################
        ############ Weibull-DPM (univariate) ###########
        #################################################
        
        if(re.type == "DPM"){
            
            ###
            n	<- dim(survData)[1]
            p	<- dim(survData)[2] - 3
            
            J	<- length(unique(survData[,3]))
            nj	<- rep(NA, J)
            
            for(i in 1:J){
                nj[i]	<- length(which(survData[,3] == i))
            }
            
            ###
            
            nStore <- numReps/thin * (1 - burninPerc)
            
            mcmc <- .C("BweibDpCorSurvmcmc",
            survData 		= as.double(as.matrix(survData)),
            n				= as.integer(n),
            p				= as.integer(p),
            J				= as.integer(J),
            nj				= as.double(nj),
            hyperParams 	= as.double(hyperParams),
            mcmcParams		= as.double(mcmcParams),
            startValues 	= as.double(startV),
            numReps			= as.integer(numReps),
            thin			= as.integer(thin),
            burninPerc      = as.double(burninPerc),
            samples_beta 	= as.double(rep(0, nStore*p)),
            samples_alpha 	= as.double(rep(0, nStore*1)),
            samples_kappa 	= as.double(rep(0, nStore*1)),
            samples_V		= as.double(rep(0, nStore*J)),
            samples_c		= as.double(rep(0, nStore*J)),
            samples_mu		= as.double(rep(0, nStore*J)),
            samples_zeta	= as.double(rep(0, nStore*J)),
            samples_tau     = as.double(rep(0, nStore*1)),
            samples_misc	= as.double(rep(0, p+1+J)))
            
            if(p > 0){
                beta.p 		<- matrix(mcmc$samples_beta, nrow = nStore, byrow = TRUE)
            }
            if(p == 0){
                beta.p 		<- NULL
            }
            
            
            alpha.p 		<- matrix(mcmc$samples_alpha, nrow = nStore, byrow = TRUE)
            kappa.p 		<- matrix(mcmc$samples_kappa, nrow = nStore, byrow = TRUE)
            V.p             <- matrix(mcmc$samples_V, nrow = nStore, byrow = TRUE)
            c.p				<- matrix(mcmc$samples_c, nrow = nStore, byrow = TRUE)
            mu.p            <- matrix(mcmc$samples_mu, nrow = nStore, byrow = TRUE)
            zeta.p          <- matrix(mcmc$samples_zeta, nrow = nStore, byrow = TRUE)
            tau.p           <- matrix(mcmc$samples_tau, nrow = nStore, byrow = TRUE)
            
            if(p > 0){
                accept.beta 	<- as.vector(mcmc$samples_misc[1:p])
            }
            if(p == 0){
                accept.beta 	<- NULL
            }
            
            accept.alpha	<- as.vector(mcmc$samples_misc[p+1])
            accept.V 	<- as.vector(mcmc$samples_misc[(p+1+1):(p+1+J)])
            
            
            Vsummary <- as.matrix(apply(V.p, 2, summary))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.975))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.025))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, sd))
            rownames(Vsummary)[7:9] <- c("0.975", "0.025", "sd")
            
            if(storeV == TRUE)
            {
                save(V.p, file = paste(path, "VPch", chain, ".RData", sep = ""))
            }
            
            if(p > 0){
                covNames = colnames(survData)[c(4:(3+p))]
            }	
            if(p == 0){
                covNames = NULL
            }
            
            ret[[nam]] <- list(beta.p = beta.p, alpha.p = alpha.p, kappa.p = kappa.p, c.p = c.p, mu.p = mu.p, zeta.p = zeta.p, tau.p = tau.p, accept.beta = accept.beta, accept.alpha = accept.alpha, accept.V = accept.V, covNames = covNames, Vsum = Vsummary)
            
        } ## end: if Weibull-DPM
        
    } ## end: if Weibull
        
        
    # hz.type = "PEM"
        
    if(hz.type == "PEM"){
            
        # re.type = "Normal"
            
        ###################################
        ############ PEM-Normal ###########
        ###################################
            
        if(re.type == "Normal"){
            
            ###
            n	<- dim(survData)[1]
            p	<- dim(survData)[2] - 3
            
            K_max = mcmcParams[4]
            
            nTime_lambda = mcmcParams[6]
            
            time_lambda	= mcmcParams[(length(mcmcParams)-nTime_lambda):(length(mcmcParams)-1)]
            
            J	<- length(unique(survData[,3]))
            
            nj	<- rep(NA, J)
            
            for(i in 1:J){
                nj[i]	<- length(which(survData[,3] == i))
            }
                       
            ###
            
            nStore <- numReps/thin * (1 - burninPerc)
            
            mcmc <- .C("BpeMvnCorSurvmcmc",
            survData 		= as.double(as.matrix(survData)),
            n				= as.integer(n),
            p				= as.integer(p),
            J				= as.integer(J),
            nj				= as.double(nj),
            hyperParams 	= as.double(hyperParams),
            mcmcParams		= as.double(mcmcParams),
            startValues 	= as.double(startV),
            numReps			= as.integer(numReps),
            thin			= as.integer(thin),
            burninPerc      = as.double(burninPerc),
            samples_beta 	= as.double(rep(0, nStore*p)),
            samples_mu_lam  = as.double(rep(0, nStore*1)),
            samples_sigSq_lam  = as.double(rep(0, nStore*1)),
            samples_K          = as.double(rep(0, nStore*1)),
            samples_s          = as.double(rep(0, nStore*(K_max + 1))),
            samples_V		= as.double(rep(0, nStore*J)),
            samples_zeta	= as.double(rep(0, nStore*1)),
            samples_misc	= as.double(rep(0, p+2+J)),
            lambda_fin     = as.double(rep(0, nStore*nTime_lambda)),
            dev     		= as.double(rep(0, nStore*1)))

            
            
            if(p > 0){
                beta.p 		<- matrix(mcmc$samples_beta, nrow = nStore, byrow = TRUE)
            }
            if(p == 0){
                beta.p 		<- NULL
            }
            
            lambda.fin 	<- matrix(mcmc$lambda_fin, nrow = nStore, byrow = TRUE)
            
            mu_lam.p 		<- matrix(mcmc$samples_mu_lam, nrow = nStore, byrow = TRUE)
            sigSq_lam.p 	<- matrix(mcmc$samples_sigSq_lam, nrow = nStore, byrow = TRUE)
            K.p 			<- matrix(mcmc$samples_K, nrow = nStore, byrow = TRUE)
            s.p 			<- matrix(mcmc$samples_s, nrow = nStore, byrow = TRUE)
            V.p             <- matrix(mcmc$samples_V, nrow = nStore, byrow = TRUE)
            zeta.p 			<- matrix(mcmc$samples_zeta, nrow = nStore, byrow = TRUE)
            
            if(p > 0){
                accept.beta 	<- as.vector(mcmc$samples_misc[1:p])
            }
            if(p == 0){
                accept.beta 	<- NULL
            }
            
            accept.BI		<- as.vector(mcmc$samples_misc[(p)+1])
            accept.DI		<- as.vector(mcmc$samples_misc[(p)+2])
            accept.V 	<- as.vector(mcmc$samples_misc[(p+3):(p+3+J)])
            
            dev.p 			<- matrix(mcmc$dev, nrow = nStore, byrow = TRUE)
            
            Vsummary <- as.matrix(apply(V.p, 2, summary))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.975))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.025))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, sd))
            rownames(Vsummary)[7:9] <- c("0.975", "0.025", "sd")
            
            if(storeV == TRUE)
            {
                save(V.p, file = paste(path, "VPch", chain, ".RData", sep = ""))
            }
            
            
            if(p > 0){
                covNames = colnames(survData)[c(4:(3+p))]
            }	
            if(p == 0){
                covNames = NULL
            }
            
            ret[[nam]] <- list(beta.p = beta.p, lambda.fin = lambda.fin, mu_lam.p = mu_lam.p, sigSq_lam.p = sigSq_lam.p, K.p = K.p, s.p = s.p, zeta.p = zeta.p, accept.beta = accept.beta, accept.BI = accept.BI, accept.DI = accept.DI, time_lambda = time_lambda, accept.V = accept.V, covNames = covNames, dev.p = dev.p, Vsum = Vsummary)
            
        }   ## end: if PEM-Normal
        
        # re.type = "DPM"
        
        ###################################
        ############ PEM-DPM ###########
        ###################################
        
        if(re.type == "DPM"){
            
            ###
            n	<- dim(survData)[1]
            p	<- dim(survData)[2] - 3
            
            K_max = mcmcParams[4]
            
            nTime_lambda = mcmcParams[6]
            
            time_lambda	= mcmcParams[(length(mcmcParams)-nTime_lambda):(length(mcmcParams)-1)]
            
            J	<- length(unique(survData[,3]))

            nj	<- rep(NA, J)
            
            for(i in 1:J){
                nj[i]	<- length(which(survData[,3] == i))
            }
            
            ###
            
            nStore <- numReps/thin * (1 - burninPerc)
            
            mcmc <- .C("BpeDpCorSurvmcmc",
            survData 		= as.double(as.matrix(survData)),
            n				= as.integer(n),
            p				= as.integer(p),
            J				= as.integer(J),
            nj				= as.double(nj),
            hyperParams 	= as.double(hyperParams),
            mcmcParams		= as.double(mcmcParams),
            startValues 	= as.double(startV),
            numReps			= as.integer(numReps),
            thin			= as.integer(thin),
            burninPerc      = as.double(burninPerc),
            samples_beta 	= as.double(rep(0, nStore*p)),
            samples_mu_lam  = as.double(rep(0, nStore*1)),
            samples_sigSq_lam  = as.double(rep(0, nStore*1)),
            samples_K          = as.double(rep(0, nStore*1)),
            samples_s          = as.double(rep(0, nStore*(K_max + 1))),
            samples_V		= as.double(rep(0, nStore*J)),
            samples_c		= as.double(rep(0, nStore*J)),
            samples_mu		= as.double(rep(0, nStore*J)),
            samples_zeta	= as.double(rep(0, nStore*J)),
            samples_tau     = as.double(rep(0, nStore*1)),
            samples_misc	= as.double(rep(0, p+2+J)),
            lambda_fin     = as.double(rep(0, nStore*nTime_lambda)),
            dev     		= as.double(rep(0, nStore*1)))
            
            
            
            if(p > 0){
                beta.p 		<- matrix(mcmc$samples_beta, nrow = nStore, byrow = TRUE)
            }
            if(p == 0){
                beta.p 		<- NULL
            }
            
            lambda.fin 	<- matrix(mcmc$lambda_fin, nrow = nStore, byrow = TRUE)
            
            mu_lam.p 		<- matrix(mcmc$samples_mu_lam, nrow = nStore, byrow = TRUE)
            sigSq_lam.p 	<- matrix(mcmc$samples_sigSq_lam, nrow = nStore, byrow = TRUE)
            K.p 			<- matrix(mcmc$samples_K, nrow = nStore, byrow = TRUE)
            s.p 			<- matrix(mcmc$samples_s, nrow = nStore, byrow = TRUE)
            V.p             <- matrix(mcmc$samples_V, nrow = nStore, byrow = TRUE)
            c.p				<- matrix(mcmc$samples_c, nrow = nStore, byrow = TRUE)
            mu.p            <- matrix(mcmc$samples_mu, nrow = nStore, byrow = TRUE)
            zeta.p          <- matrix(mcmc$samples_zeta, nrow = nStore, byrow = TRUE)
            tau.p           <- matrix(mcmc$samples_tau, nrow = nStore, byrow = TRUE)
            
            if(p > 0){
                accept.beta 	<- as.vector(mcmc$samples_misc[1:p])
            }
            if(p == 0){
                accept.beta 	<- NULL
            }
            
            accept.BI		<- as.vector(mcmc$samples_misc[(p)+1])
            accept.DI		<- as.vector(mcmc$samples_misc[(p)+2])
            accept.V 	<- as.vector(mcmc$samples_misc[(p+3):(p+2+J)])
            
            dev.p 			<- matrix(mcmc$dev, nrow = nStore, byrow = TRUE)
            
            Vsummary <- as.matrix(apply(V.p, 2, summary))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.975))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.025))
            Vsummary <- rbind(Vsummary, apply(V.p, 2, sd))
            rownames(Vsummary)[7:9] <- c("0.975", "0.025", "sd")
            
            if(storeV == TRUE)
            {
                save(V.p, file = paste(path, "VPch", chain, ".RData", sep = ""))
            }
            
            if(p > 0){
                covNames = colnames(survData)[c(4:(3+p))]
            }	
            if(p == 0){
                covNames = NULL
            }
            
            ret[[nam]] <- list(beta.p = beta.p, lambda.fin = lambda.fin, mu_lam.p = mu_lam.p, sigSq_lam.p = sigSq_lam.p, K.p = K.p, s.p = s.p, c.p = c.p, mu.p = mu.p, zeta.p = zeta.p, tau.p = tau.p, accept.beta = accept.beta, accept.BI = accept.BI, accept.DI = accept.DI, time_lambda = time_lambda, accept.V = accept.V, covNames = covNames, dev.p = dev.p, Vsum = Vsummary)
            
            
        }   ## end: if PEM-DPM
        
        
    } ## end: if PEM
        
        
				
	chain = chain + 1	
	
    }## end: while(chain <= nChain)
	
	
    
    ret[["setup"]]	<- list(hyperParams = hyperParams, startValues = startValues, mcmcParams = mcmcParams, numReps = numReps, thin = thin, path = path, burninPerc = burninPerc, hz.type = hz.type, re.type = re.type, nChain = nChain)

	class(ret) <- "BayesSurvcor"
	return(ret)

}  ## end: if(class(startValues) == "list" & length(startValues) == nChain)

else{
	print("The 'startValues' should be the list of length equal to 'nChain'.")
}

} # end of function "BayesSurvcor"





























