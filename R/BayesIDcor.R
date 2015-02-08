

### hz.type = "Weibull" or "PEM"
### re.type = "MVN" or "DPM"
### model = "Markov" or "semi-Markov"

BayesIDcor <- function(survData, 
					nCov,					
					hyperParams,
					startValues,								
					mcmcParams,
					nGam_save,					
					numReps,
					thin,
					path = "results/",
					burninPerc=0.5,
					hz.type = "Weibull",
                    re.type = "MVN",
					model = "Markov",
                    storeV=rep(TRUE, 3),
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
    
        # re.type = "MVN"
	
        if(re.type == "MVN"){    
    
            # model = "Markov"

            #######################################
            ############ Weibull-MVN-M ############
            #######################################
		
            if(model == "Markov"){
	
			###
            n	<- dim(survData)[1]
            p1	<- nCov[1]
            p2	<- nCov[2]
            p3	<- nCov[3]
            
            J	<- length(unique(survData[,5]))
            
            nj	<- rep(NA, J)
            
            for(i in 1:J){
                nj[i]	<- length(which(survData[,5] == i))
            }            
	
			nStore <- numReps/thin * (1 - burninPerc)
                
                startValues1 <- startV[1:(p1+p2+p3+n+7)]
                startValues2 <- startV[(p1+p2+p3+n+8):length(startV)]
                startV <- c(startValues1, 1, 1, startValues2)
                
                mcmcParams1 <- mcmcParams[1:4]
                mcmcParams2 <- mcmcParams[5:7]
                
                mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
	

			mcmc <- .C("BweibMvnCorScrmcmc",
						survData 		= as.double(as.matrix(survData)),
						n				= as.integer(n),
						p1				= as.integer(p1),
						p2				= as.integer(p2),
						p3				= as.integer(p3),
						J				= as.integer(J),
						nj				= as.double(nj),
                        hyperParams 	= as.double(hyperParams),
						mcmcParams		= as.double(mcmcParams),
						startValues 	= as.double(startV),
						numReps			= as.integer(numReps),
						thin			= as.integer(thin),
						burninPerc      = as.double(burninPerc),
						nGam_save		= as.integer(nGam_save),
						samples_beta1 	= as.double(rep(0, nStore*p1)),
						samples_beta2 	= as.double(rep(0, nStore*p2)),
						samples_beta3 	= as.double(rep(0, nStore*p3)),
						samples_alpha1 	= as.double(rep(0, nStore*1)),
						samples_alpha2 	= as.double(rep(0, nStore*1)),
						samples_alpha3 	= as.double(rep(0, nStore*1)),
						samples_kappa1 	= as.double(rep(0, nStore*1)),
						samples_kappa2 	= as.double(rep(0, nStore*1)),
						samples_kappa3 	= as.double(rep(0, nStore*1)),
						samples_nu2 	= as.double(rep(0, nStore*1)),
						samples_nu3 	= as.double(rep(0, nStore*1)),
						samples_theta 	= as.double(rep(0, nStore*1)),
						samples_V1		= as.double(rep(0, nStore*J)),
						samples_V2		= as.double(rep(0, nStore*J)),
						samples_V3		= as.double(rep(0, nStore*J)),	
						samples_Sigma_V	= as.double(rep(0, nStore*3*3)),
						samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
						samples_misc	= as.double(rep(0, (p1+p2+p3+6+n+1+n+J+J+J))),
                        gammaP			= as.double(rep(0, n)),
                        dev     			= as.double(rep(0, nStore*1)),
                        invLH = as.double(rep(0, n)),
                        logLH_fin = as.double(0),
                        lpml     			= as.double(rep(0, nStore*1)),
                        lpml2     			= as.double(rep(0, nStore*1)))

            if(p1 > 0){
                beta1.p 		<- matrix(mcmc$samples_beta1, nrow = nStore, byrow = TRUE)		
            }					
            if(p1 == 0){
                beta1.p 		<- NULL
            }
            if(p2 > 0){
                beta2.p 		<- matrix(mcmc$samples_beta2, nrow = nStore, byrow = TRUE)		
            }					
            if(p2 == 0){
                beta2.p 		<- NULL
            }
            if(p3 > 0){
                beta3.p 		<- matrix(mcmc$samples_beta3, nrow = nStore, byrow = TRUE)		
            }					
            if(p3 == 0){
                beta3.p 		<- NULL
            }		

            alpha1.p 		<- matrix(mcmc$samples_alpha1, nrow = nStore, byrow = TRUE)		
            alpha2.p 		<- matrix(mcmc$samples_alpha2, nrow = nStore, byrow = TRUE)
            alpha3.p 		<- matrix(mcmc$samples_alpha3, nrow = nStore, byrow = TRUE)				
            kappa1.p 		<- matrix(mcmc$samples_kappa1, nrow = nStore, byrow = TRUE)
            kappa2.p 		<- matrix(mcmc$samples_kappa2, nrow = nStore, byrow = TRUE)
            kappa3.p 		<- matrix(mcmc$samples_kappa3, nrow = nStore, byrow = TRUE)
            theta.p 		<- matrix(mcmc$samples_theta, nrow = nStore, byrow = TRUE)
            gamma.p 		<- matrix(mcmc$samples_gamma, nrow = nStore, byrow = TRUE)

            V1.p            <- matrix(mcmc$samples_V1, nrow = nStore, byrow = TRUE)
            V2.p            <- matrix(mcmc$samples_V2, nrow = nStore, byrow = TRUE)
            V3.p            <- matrix(mcmc$samples_V3, nrow = nStore, byrow = TRUE)
            Sigma_V.p		<- array(as.vector(mcmc$samples_Sigma_V), c(3, 3, nStore))	
		
            if(p1 > 0){
                accept.beta1 	<- as.vector(mcmc$samples_misc[1:(p1)])		
            }
            if(p1 == 0){
                accept.beta1 	<- NULL
            }
            if(p2 > 0){
                accept.beta2 	<- as.vector(mcmc$samples_misc[(p1+1):(p1+p2)])
            }
            if(p2 == 0){
                accept.beta2 	<- NULL
            }
            if(p3 > 0){
                accept.beta3 	<- as.vector(mcmc$samples_misc[(p1+p2+1):(p1+p2+p3)])		
            }
            if(p3 == 0){
                accept.beta3 	<- NULL
            }		

            accept.alpha1	<- as.vector(mcmc$samples_misc[(p1+p2+p3+1)])
            accept.alpha2	<- as.vector(mcmc$samples_misc[(p1+p2+p3+2)])
            accept.alpha3	<- as.vector(mcmc$samples_misc[(p1+p2+p3+3)])
            accept.theta	<- as.vector(mcmc$samples_misc[(p1+p2+p3+4)])
            accept.V1       <- as.vector(mcmc$samples_misc[(p1+p2+p3+7+n+n+1):(p1+p2+p3+7+n+n+J)])
            accept.V2       <- as.vector(mcmc$samples_misc[(p1+p2+p3+7+n+n+J+1):(p1+p2+p3+7+n+n+J+J)])
            accept.V3       <- as.vector(mcmc$samples_misc[(p1+p2+p3+7+n+n+J+J+1):(p1+p2+p3+7+n+n+J+J+J)])
	
            V1summary <- as.matrix(apply(V1.p, 2, summary))
            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
	
            V2summary <- as.matrix(apply(V2.p, 2, summary))
            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
	
            V3summary <- as.matrix(apply(V3.p, 2, summary))
            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")
        
            if(storeV[1] == TRUE)
            {
                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
            }
            if(storeV[2] == TRUE)
            {
                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
            }
            if(storeV[3] == TRUE)
            {
                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
            }
    
            save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))

	
            if(p1 > 0){
                covNames1 = colnames(survData)[c(5:(4+p1))]
                }	
            if(p1 == 0){
                covNames1 = NULL
                }
        
            if(p2 > 0){
                covNames2 = colnames(survData)[c((4+p1+1):(4+p1+p2))]
                }	
            if(p2 == 0){
                covNames2 = NULL
                }     
            if(p3 > 0){
                covNames3 = colnames(survData)[c((4+p1+p2+1):(4+p1+p2+p3))]
                }	
            if(p3 == 0){
                covNames3 = NULL
                }
    
            ### posterior predictive checks ###
    
            ## 1. log pseudo-marginal likelihood
    
            invLH.p <- matrix(mcmc$invLH, nrow = n, byrow = TRUE)
    
            cpo     <- 1/invLH.p
            
            LPML <- sum(log(cpo))
    
            # or
    
            lpml.p <- matrix(mcmc$lpml, nrow = nStore, byrow = TRUE)
    
            lpml2.p <- matrix(mcmc$lpml2, nrow = nStore, byrow = TRUE)
    
            ## 2. deviance information criterion
    
            gamma_mean <- matrix(mcmc$gammaP, nrow = n, byrow = TRUE)
            dev.p 	   <- matrix(mcmc$dev, nrow = nStore, byrow = TRUE)
            Dbar        <- mean(dev.p)
            pD          <- Dbar - (-2*mcmc$logLH_fin)
    
            DIC <- pD + Dbar    	
	
			
			ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, theta.p = theta.p, Sigma_V.p = Sigma_V.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, accept.V1 = accept.V1, accept.V2 = accept.V2, accept.V3 = accept.V3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, dev.p = dev.p, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, cpo = cpo, LPML = LPML, gamma_mean = gamma_mean, DIC = DIC, val = mcmc$logLH_fin, lpml.p = lpml.p, lpml2.p = lpml2.p)

			} ## end: if Weibull-MVN-M
            
            
            # model = "semi-Markov"
            
            #######################################
            ############ 2 Weibull-MVN-SM #########
            #######################################	
		
            if(model == "semi-Markov"){
            
            ###
            n	<- dim(survData)[1]
            p1	<- nCov[1]
            p2	<- nCov[2]
            p3	<- nCov[3]
	
            J	<- length(unique(survData[,5]))
	
            nj	<- rep(NA, J)
	
            for(i in 1:J){
                nj[i]	<- length(which(survData[,5] == i))
            }

	
            nStore <- numReps/thin * (1 - burninPerc)
                
                startValues1 <- startV[1:(p1+p2+p3+n+7)]
                startValues2 <- startV[(p1+p2+p3+n+8):length(startV)]
                startV <- c(startValues1, 1, 1, startValues2)
                
                mcmcParams1 <- mcmcParams[1:4]
                mcmcParams2 <- mcmcParams[5:7]
                
                mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
	
            mcmc <- .C("BweibMvnCorScrSMmcmc",
						survData 		= as.double(as.matrix(survData)),
						n				= as.integer(n),
						p1				= as.integer(p1),
						p2				= as.integer(p2),
						p3				= as.integer(p3),
						J				= as.integer(J),
						nj				= as.double(nj),
                        hyperParams 	= as.double(hyperParams),
						mcmcParams		= as.double(mcmcParams),
						startValues 	= as.double(startV),
						numReps			= as.integer(numReps),
						thin			= as.integer(thin),
						burninPerc      = as.double(burninPerc),
						nGam_save		= as.integer(nGam_save),
						samples_beta1 	= as.double(rep(0, nStore*p1)),
						samples_beta2 	= as.double(rep(0, nStore*p2)),
						samples_beta3 	= as.double(rep(0, nStore*p3)),
						samples_alpha1 	= as.double(rep(0, nStore*1)),
						samples_alpha2 	= as.double(rep(0, nStore*1)),
						samples_alpha3 	= as.double(rep(0, nStore*1)),
						samples_kappa1 	= as.double(rep(0, nStore*1)),
						samples_kappa2 	= as.double(rep(0, nStore*1)),
						samples_kappa3 	= as.double(rep(0, nStore*1)),
						samples_nu2 	= as.double(rep(0, nStore*1)),
						samples_nu3 	= as.double(rep(0, nStore*1)),
						samples_theta 	= as.double(rep(0, nStore*1)),
						samples_V1		= as.double(rep(0, nStore*J)),
						samples_V2		= as.double(rep(0, nStore*J)),
						samples_V3		= as.double(rep(0, nStore*J)),	
						samples_Sigma_V	= as.double(rep(0, nStore*3*3)),
						samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                        samples_misc	= as.double(rep(0, (p1+p2+p3+6+n+1+n+J+J+J))),
                        gammaP			= as.double(rep(0, n)),
                        dev     			= as.double(rep(0, nStore*1)),
                        invLH = as.double(rep(0, n)),
                        logLH_fin = as.double(0),
                        lpml     			= as.double(rep(0, nStore*1)),
                        lpml2     			= as.double(rep(0, nStore*1)))
				
				
				
            if(p1 > 0){
                beta1.p 		<- matrix(mcmc$samples_beta1, nrow = nStore, byrow = TRUE)
            }					
            if(p1 == 0){
                beta1.p 		<- NULL
            }
            if(p2 > 0){
                beta2.p 		<- matrix(mcmc$samples_beta2, nrow = nStore, byrow = TRUE)		
            }					
            if(p2 == 0){
                beta2.p 		<- NULL
            }
            if(p3 > 0){
                beta3.p 		<- matrix(mcmc$samples_beta3, nrow = nStore, byrow = TRUE)		
            }					
            if(p3 == 0){
                beta3.p 		<- NULL
            }		

            alpha1.p 		<- matrix(mcmc$samples_alpha1, nrow = nStore, byrow = TRUE)		
            alpha2.p 		<- matrix(mcmc$samples_alpha2, nrow = nStore, byrow = TRUE)		
            alpha3.p 		<- matrix(mcmc$samples_alpha3, nrow = nStore, byrow = TRUE)				
            kappa1.p 		<- matrix(mcmc$samples_kappa1, nrow = nStore, byrow = TRUE)
            kappa2.p 		<- matrix(mcmc$samples_kappa2, nrow = nStore, byrow = TRUE)
            kappa3.p 		<- matrix(mcmc$samples_kappa3, nrow = nStore, byrow = TRUE)
            theta.p 		<- matrix(mcmc$samples_theta, nrow = nStore, byrow = TRUE)
            gamma.p 		<- matrix(mcmc$samples_gamma, nrow = nStore, byrow = TRUE)
            V1.p            <- matrix(mcmc$samples_V1, nrow = nStore, byrow = TRUE)
            V2.p            <- matrix(mcmc$samples_V2, nrow = nStore, byrow = TRUE)
            V3.p            <- matrix(mcmc$samples_V3, nrow = nStore, byrow = TRUE)
            Sigma_V.p		<- array(as.vector(mcmc$samples_Sigma_V), c(3, 3, nStore))	
		
            if(p1 > 0){
                accept.beta1 	<- as.vector(mcmc$samples_misc[1:(p1)])		
            }
            if(p1 == 0){
                accept.beta1 	<- NULL
            }
            if(p2 > 0){
                accept.beta2 	<- as.vector(mcmc$samples_misc[(p1+1):(p1+p2)])
            }
            if(p2 == 0){
                accept.beta2 	<- NULL
            }
            if(p3 > 0){
                accept.beta3 	<- as.vector(mcmc$samples_misc[(p1+p2+1):(p1+p2+p3)])		
            }
            if(p3 == 0){
                accept.beta3 	<- NULL
            }		

            accept.alpha1	<- as.vector(mcmc$samples_misc[(p1+p2+p3+1)])
            accept.alpha2	<- as.vector(mcmc$samples_misc[(p1+p2+p3+2)])
            accept.alpha3	<- as.vector(mcmc$samples_misc[(p1+p2+p3+3)])
            accept.theta	<- as.vector(mcmc$samples_misc[(p1+p2+p3+4)])
            accept.V1       <- as.vector(mcmc$samples_misc[(p1+p2+p3+7+n+n+1):(p1+p2+p3+7+n+n+J)])
            accept.V2       <- as.vector(mcmc$samples_misc[(p1+p2+p3+7+n+n+J+1):(p1+p2+p3+7+n+n+J+J)])
            accept.V3       <- as.vector(mcmc$samples_misc[(p1+p2+p3+7+n+n+J+J+1):(p1+p2+p3+7+n+n+J+J+J)])
	
            V1summary <- as.matrix(apply(V1.p, 2, summary))
            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
	
            V2summary <- as.matrix(apply(V2.p, 2, summary))
            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
	
            V3summary <- as.matrix(apply(V3.p, 2, summary))
            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")
	
            if(storeV[1] == TRUE)
            {
                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
            }
            if(storeV[2] == TRUE)
            {
                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
            }
            if(storeV[3] == TRUE)
            {
                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
            }
    
            save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
    
	
	
            if(p1 > 0){
                covNames1 = colnames(survData)[c(6:(5+p1))]
                }	
            if(p1 == 0){
                covNames1 = NULL
                }
        
            if(p2 > 0){
                covNames2 = colnames(survData)[c((5+p1+1):(5+p1+p2))]
                }	
            if(p2 == 0){
                covNames2 = NULL
                }     
            if(p3 > 0){
                covNames3 = colnames(survData)[c((5+p1+p2+1):(5+p1+p2+p3))]
                }	
            if(p3 == 0){
                covNames3 = NULL
                }
    
    
    
            ### posterior predictive checks ###

            ## 1. log pseudo-marginal likelihood
    
            invLH.p <- matrix(mcmc$invLH, nrow = n, byrow = TRUE)
    
            cpo     <- 1/invLH.p
    
            LPML <- sum(log(cpo))
    
            # or
    
            lpml.p <- matrix(mcmc$lpml, nrow = nStore, byrow = TRUE)
    
            lpml2.p <- matrix(mcmc$lpml2, nrow = nStore, byrow = TRUE)

            ## 2. deviance information criterion
    
            gamma_mean <- matrix(mcmc$gammaP, nrow = n, byrow = TRUE)
            dev.p 	   <- matrix(mcmc$dev, nrow = nStore, byrow = TRUE)
            Dbar        <- mean(dev.p)
            pD          <- Dbar - (-2*mcmc$logLH_fin)
    
            DIC <- pD + Dbar
          		
            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, theta.p = theta.p, Sigma_V.p = Sigma_V.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, accept.V1 = accept.V1, accept.V2 = accept.V2, accept.V3 = accept.V3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, dev.p = dev.p, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, cpo = cpo, LPML = LPML, gamma_mean = gamma_mean, DIC = DIC, val = mcmc$logLH_fin, lpml.p = lpml.p, lpml2.p = lpml2.p)

                        
            }  ## end: if Weibull-MVN-SM          

        } ## end: if Weibull-MVN
        
        
        # re.type = "DPM"
	
        if(re.type == "DPM"){    
    
            # model = "Markov"

            #######################################
            ############ Weibull-DPM-M ############
            #######################################
		
            if(model == "Markov"){  
            
            ###
            n	<- dim(survData)[1]
            p1	<- nCov[1]
            p2	<- nCov[2]
            p3	<- nCov[3]
	
            J	<- length(unique(survData[,5]))
            nj	<- rep(NA, J)
	
            for(i in 1:J){
                nj[i]	<- length(which(survData[,5] == i))
            }
	
            ### p_1+p_2+p_3+n+4J+8
	
            nStore <- numReps/thin * (1 - burninPerc)
                
                startValues1 <- startV[1:(p1+p2+p3+n+7)]
                startValues2 <- startV[(p1+p2+p3+n+8):length(startV)]
                startV <- c(startValues1, 1, 1, startValues2)
                
                mcmcParams1 <- mcmcParams[1:4]
                mcmcParams2 <- mcmcParams[5:7]
                
                mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
                
	
            mcmc <- .C("BweibDpCorScrmcmc",
						survData 		= as.double(as.matrix(survData)),
						n				= as.integer(n),
						p1				= as.integer(p1),
						p2				= as.integer(p2),
						p3				= as.integer(p3),
						J				= as.integer(J),
						nj				= as.double(nj),
                        hyperParams 	= as.double(hyperParams),
						mcmcParams		= as.double(mcmcParams),
						startValues 	= as.double(startV),
						numReps			= as.integer(numReps),
						thin			= as.integer(thin),
						burninPerc      = as.double(burninPerc),
						nGam_save		= as.integer(nGam_save),
						samples_beta1 	= as.double(rep(0, nStore*p1)),
						samples_beta2 	= as.double(rep(0, nStore*p2)),
						samples_beta3 	= as.double(rep(0, nStore*p3)),
						samples_alpha1 	= as.double(rep(0, nStore*1)),
						samples_alpha2 	= as.double(rep(0, nStore*1)),
						samples_alpha3 	= as.double(rep(0, nStore*1)),
						samples_kappa1 	= as.double(rep(0, nStore*1)),
						samples_kappa2 	= as.double(rep(0, nStore*1)),
						samples_kappa3 	= as.double(rep(0, nStore*1)),
						samples_nu2 	= as.double(rep(0, nStore*1)),
						samples_nu3 	= as.double(rep(0, nStore*1)),
						samples_theta 	= as.double(rep(0, nStore*1)),
						samples_V1		= as.double(rep(0, nStore*J)),
						samples_V2		= as.double(rep(0, nStore*J)),
						samples_V3		= as.double(rep(0, nStore*J)),	
						samples_c		= as.double(rep(0, nStore*J)),
						samples_mu		= as.double(rep(0, nStore*3*J)),
						samples_Sigma	= as.double(rep(0, nStore*3*3*J)),
                        samples_tau     = as.double(rep(0, nStore*1)),
						samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                        samples_misc	= as.double(rep(0, (p1+p2+p3+6+n+1+n+J))),
                        gammaP			= as.double(rep(0, n)),
                        dev     			= as.double(rep(0, nStore*1)),
                        invLH = as.double(rep(0, n)),
                        logLH_fin = as.double(0),
                        lpml     			= as.double(rep(0, nStore*1)),
                        lpml2     			= as.double(rep(0, nStore*1)))
				
				
				
            if(p1 > 0){
                beta1.p 		<- matrix(mcmc$samples_beta1, nrow = nStore, byrow = TRUE)
            }					
            if(p1 == 0){
                beta1.p 		<- NULL
            }
            if(p2 > 0){
                beta2.p 		<- matrix(mcmc$samples_beta2, nrow = nStore, byrow = TRUE)		
            }
            if(p2 == 0){
                beta2.p 		<- NULL
            }
            if(p3 > 0){
                beta3.p 		<- matrix(mcmc$samples_beta3, nrow = nStore, byrow = TRUE)		
            }					
            if(p3 == 0){
                beta3.p 		<- NULL
            }		

            alpha1.p 		<- matrix(mcmc$samples_alpha1, nrow = nStore, byrow = TRUE)		
            alpha2.p 		<- matrix(mcmc$samples_alpha2, nrow = nStore, byrow = TRUE)		
            alpha3.p 		<- matrix(mcmc$samples_alpha3, nrow = nStore, byrow = TRUE)				
            kappa1.p 		<- matrix(mcmc$samples_kappa1, nrow = nStore, byrow = TRUE)
            kappa2.p 		<- matrix(mcmc$samples_kappa2, nrow = nStore, byrow = TRUE)
            kappa3.p 		<- matrix(mcmc$samples_kappa3, nrow = nStore, byrow = TRUE)
            nu2.p           <- matrix(mcmc$samples_nu2, nrow = nStore, byrow = TRUE)
            nu3.p           <- matrix(mcmc$samples_nu3, nrow = nStore, byrow = TRUE)
            theta.p 		<- matrix(mcmc$samples_theta, nrow = nStore, byrow = TRUE)
            gamma.p 		<- matrix(mcmc$samples_gamma, nrow = nStore, byrow = TRUE)

            V1.p            <- matrix(mcmc$samples_V1, nrow = nStore, byrow = TRUE)
            V2.p            <- matrix(mcmc$samples_V2, nrow = nStore, byrow = TRUE)
            V3.p            <- matrix(mcmc$samples_V3, nrow = nStore, byrow = TRUE)
            c.p            <- matrix(mcmc$samples_c, nrow = nStore, byrow = TRUE)
            mu.p            <- matrix(mcmc$samples_mu, nrow = nStore, byrow = TRUE)
            Sigma.p         <- array(as.vector(mcmc$samples_Sigma), c(3, 3 * J, nStore))
            tau.p           <- matrix(mcmc$samples_tau, nrow = nStore, byrow = TRUE)		
		
            if(p1 > 0){
                accept.beta1 	<- as.vector(mcmc$samples_misc[1:(p1)])		
            }
            if(p1 == 0){
                accept.beta1 	<- NULL
            }
            if(p2 > 0){
                accept.beta2 	<- as.vector(mcmc$samples_misc[(p1+1):(p1+p2)])
            }
            if(p2 == 0){
                accept.beta2 	<- NULL
            }
            if(p3 > 0){
                accept.beta3 	<- as.vector(mcmc$samples_misc[(p1+p2+1):(p1+p2+p3)])		
            }
            if(p3 == 0){
                accept.beta3 	<- NULL
            }		

            accept.alpha1	<- as.vector(mcmc$samples_misc[(p1+p2+p3+1)])
            accept.alpha2	<- as.vector(mcmc$samples_misc[(p1+p2+p3+2)])
            accept.alpha3	<- as.vector(mcmc$samples_misc[(p1+p2+p3+3)])
            accept.theta	<- as.vector(mcmc$samples_misc[(p1+p2+p3+4)])
            accept.nu2      <- as.vector(mcmc$samples_misc[(p1+p2+p3+5)])
            accept.nu3      <- as.vector(mcmc$samples_misc[(p1+p2+p3+6)])
            accept.gamma    <- as.vector(mcmc$samples_misc[(p1+p2+p3+6+1):(p1+p2+p3+6+n)])
            lastChgProp     <- as.vector(mcmc$samples_misc[(p1+p2+p3+6+n+1)])
            mhGam_chk       <- as.vector(mcmc$samples_misc[(p1+p2+p3+6+n+1+1):(p1+p2+p3+6+n+1+n)])
            accept.V        <- as.vector(mcmc$samples_misc[(p1+p2+p3+7+n+n+1):(p1+p2+p3+7+n+n+J)])


            V1summary <- as.matrix(apply(V1.p, 2, summary))
            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")	
	
            V2summary <- as.matrix(apply(V2.p, 2, summary))
            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")	

            V3summary <- as.matrix(apply(V3.p, 2, summary))
            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")
    		
	
            if(storeV[1] == TRUE)
            {
                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
            }
            if(storeV[2] == TRUE)
            {
                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
            }
            if(storeV[3] == TRUE)
            {
                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
            }
    
            save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
    	
	
            if(p1 > 0){
                covNames1 = colnames(survData)[c(5:(4+p1))]
                }	
            if(p1 == 0){
                covNames1 = NULL
                }
        
            if(p2 > 0){
                covNames2 = colnames(survData)[c((4+p1+1):(4+p1+p2))]
                }	
            if(p2 == 0){
                covNames2 = NULL
                }     
            if(p3 > 0){
                covNames3 = colnames(survData)[c((4+p1+p2+1):(4+p1+p2+p3))]
                }	
            if(p3 == 0){
                covNames3 = NULL
                }
    
    
            ### posterior predictive checks ###
    
            ## 1. log pseudo-marginal likelihood
    
            invLH.p <- matrix(mcmc$invLH, nrow = n, byrow = TRUE)

            cpo     <- 1/invLH.p
    
            LPML <- sum(log(cpo))
    
            # or
    
            lpml.p <- matrix(mcmc$lpml, nrow = nStore, byrow = TRUE)
            
            lpml2.p <- matrix(mcmc$lpml2, nrow = nStore, byrow = TRUE)
    
            ## 2. deviance information criterion
    
            gamma_mean <- matrix(mcmc$gammaP, nrow = n, byrow = TRUE)
            dev.p 	   <- matrix(mcmc$dev, nrow = nStore, byrow = TRUE)
            Dbar        <- mean(dev.p)
            pD          <- Dbar - (-2*mcmc$logLH_fin)
    
            DIC <- pD + Dbar
    
    
    
          		
            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, nu2.p = nu2.p, nu3.p = nu3.p, theta.p = theta.p, V1.p = V1.p, V2.p = V2.p, V3.p = V3.p, c.p = c.p, mu.p = mu.p, Sigma.p = Sigma.p, tau.p = tau.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, accept.nu2 = accept.nu2, accept.nu3 = accept.nu3, accept.gamma = accept.gamma, lastChgProp = lastChgProp, mhGam_chk = mhGam_chk, accept.V = accept.V, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, dev.p = dev.p, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, cpo = cpo, LPML = LPML, gamma_mean = gamma_mean, DIC = DIC, val = mcmc$logLH_fin, lpml.p = lpml.p, lpml2.p = lpml2.p)
            
            
            }  ## end: if Weibull-DPM-M      
            
                
            # model = "semi-Markov"

            #######################################
            ############ Weibull-DPM-SM ############
            #######################################
		
            if(model == "semi-Markov"){  
            
            ###
            n	<- dim(survData)[1]
            p1	<- nCov[1]
            p2	<- nCov[2]
            p3	<- nCov[3]
	
            J	<- length(unique(survData[,5]))
		
            nj	<- rep(NA, J)
	
            for(i in 1:J){
                nj[i]	<- length(which(survData[,5] == i))
            }

            ###
	
            nStore <- numReps/thin * (1 - burninPerc)
                
                startValues1 <- startV[1:(p1+p2+p3+n+7)]
                startValues2 <- startV[(p1+p2+p3+n+8):length(startV)]
                startV <- c(startValues1, 1, 1, startValues2)
                
                mcmcParams1 <- mcmcParams[1:4]
                mcmcParams2 <- mcmcParams[5:7]
                
                mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
	
            mcmc <- .C("BweibDpCorScrSMmcmc",
						survData 		= as.double(as.matrix(survData)),
						n				= as.integer(n),
						p1				= as.integer(p1),
						p2				= as.integer(p2),
						p3				= as.integer(p3),
						J				= as.integer(J),
						nj				= as.double(nj),
                        hyperParams 	= as.double(hyperParams),
						mcmcParams		= as.double(mcmcParams),
						startValues 	= as.double(startV),
						numReps			= as.integer(numReps),
						thin			= as.integer(thin),
						burninPerc      = as.double(burninPerc),
						nGam_save		= as.integer(nGam_save),
						samples_beta1 	= as.double(rep(0, nStore*p1)),
						samples_beta2 	= as.double(rep(0, nStore*p2)),
						samples_beta3 	= as.double(rep(0, nStore*p3)),
						samples_alpha1 	= as.double(rep(0, nStore*1)),
						samples_alpha2 	= as.double(rep(0, nStore*1)),
						samples_alpha3 	= as.double(rep(0, nStore*1)),
						samples_kappa1 	= as.double(rep(0, nStore*1)),
						samples_kappa2 	= as.double(rep(0, nStore*1)),
						samples_kappa3 	= as.double(rep(0, nStore*1)),
						samples_nu2 	= as.double(rep(0, nStore*1)),
						samples_nu3 	= as.double(rep(0, nStore*1)),
						samples_theta 	= as.double(rep(0, nStore*1)),
						samples_V1		= as.double(rep(0, nStore*J)),
						samples_V2		= as.double(rep(0, nStore*J)),
						samples_V3		= as.double(rep(0, nStore*J)),	
						samples_c		= as.double(rep(0, nStore*J)),
						samples_mu		= as.double(rep(0, nStore*3*J)),
						samples_Sigma	= as.double(rep(0, nStore*3*3*J)),
                        samples_tau     = as.double(rep(0, nStore*1)),
						samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
						samples_misc	= as.double(rep(0, (p1+p2+p3+6+n+1+n+J))),
                        gammaP			= as.double(rep(0, n)),
                        dev     			= as.double(rep(0, nStore*1)),
                        invLH = as.double(rep(0, n)),
                        logLH_fin = as.double(0),
                        lpml     			= as.double(rep(0, nStore*1)),
                        lpml2     			= as.double(rep(0, nStore*1)))
            
				
            if(p1 > 0){
                beta1.p 		<- matrix(mcmc$samples_beta1, nrow = nStore, byrow = TRUE)
            }					
            if(p1 == 0){
                beta1.p 		<- NULL
            }
            if(p2 > 0){
                beta2.p 		<- matrix(mcmc$samples_beta2, nrow = nStore, byrow = TRUE)		
            }					
            if(p2 == 0){
                beta2.p 		<- NULL
            }
            if(p3 > 0){
                beta3.p 		<- matrix(mcmc$samples_beta3, nrow = nStore, byrow = TRUE)
            }					
            if(p3 == 0){
                beta3.p 		<- NULL
            }		

            alpha1.p 		<- matrix(mcmc$samples_alpha1, nrow = nStore, byrow = TRUE)		
            alpha2.p 		<- matrix(mcmc$samples_alpha2, nrow = nStore, byrow = TRUE)		
            alpha3.p 		<- matrix(mcmc$samples_alpha3, nrow = nStore, byrow = TRUE)
            kappa1.p 		<- matrix(mcmc$samples_kappa1, nrow = nStore, byrow = TRUE)
            kappa2.p 		<- matrix(mcmc$samples_kappa2, nrow = nStore, byrow = TRUE)
            kappa3.p 		<- matrix(mcmc$samples_kappa3, nrow = nStore, byrow = TRUE)
            nu2.p           <- matrix(mcmc$samples_nu2, nrow = nStore, byrow = TRUE)
            nu3.p           <- matrix(mcmc$samples_nu3, nrow = nStore, byrow = TRUE)
            theta.p 		<- matrix(mcmc$samples_theta, nrow = nStore, byrow = TRUE)
            gamma.p 		<- matrix(mcmc$samples_gamma, nrow = nStore, byrow = TRUE)

            V1.p            <- matrix(mcmc$samples_V1, nrow = nStore, byrow = TRUE)
            V2.p            <- matrix(mcmc$samples_V2, nrow = nStore, byrow = TRUE)
            V3.p            <- matrix(mcmc$samples_V3, nrow = nStore, byrow = TRUE)
            c.p            <- matrix(mcmc$samples_c, nrow = nStore, byrow = TRUE)
            mu.p            <- matrix(mcmc$samples_mu, nrow = nStore, byrow = TRUE)
            Sigma.p         <- array(as.vector(mcmc$samples_Sigma), c(3, 3 * J, nStore))
            tau.p           <- matrix(mcmc$samples_tau, nrow = nStore, byrow = TRUE)		
		
            if(p1 > 0){
                accept.beta1 	<- as.vector(mcmc$samples_misc[1:(p1)])		
            }
            if(p1 == 0){
                accept.beta1 	<- NULL
            }
            if(p2 > 0){
                accept.beta2 	<- as.vector(mcmc$samples_misc[(p1+1):(p1+p2)])
            }
            if(p2 == 0){
                accept.beta2 	<- NULL
            }
            if(p3 > 0){
                accept.beta3 	<- as.vector(mcmc$samples_misc[(p1+p2+1):(p1+p2+p3)])		
            }
            if(p3 == 0){
                accept.beta3 	<- NULL
            }		

            accept.alpha1	<- as.vector(mcmc$samples_misc[(p1+p2+p3+1)])
            accept.alpha2	<- as.vector(mcmc$samples_misc[(p1+p2+p3+2)])
            accept.alpha3	<- as.vector(mcmc$samples_misc[(p1+p2+p3+3)])
            accept.theta	<- as.vector(mcmc$samples_misc[(p1+p2+p3+4)])
            accept.nu2      <- as.vector(mcmc$samples_misc[(p1+p2+p3+5)])
            accept.nu3      <- as.vector(mcmc$samples_misc[(p1+p2+p3+6)])
            accept.gamma    <- as.vector(mcmc$samples_misc[(p1+p2+p3+6+1):(p1+p2+p3+6+n)])
            lastChgProp     <- as.vector(mcmc$samples_misc[(p1+p2+p3+6+n+1)])
            mhGam_chk       <- as.vector(mcmc$samples_misc[(p1+p2+p3+6+n+1+1):(p1+p2+p3+6+n+1+n)])
            accept.V        <- as.vector(mcmc$samples_misc[(p1+p2+p3+7+n+n+1):(p1+p2+p3+7+n+n+J)])
    
            V1summary <- as.matrix(apply(V1.p, 2, summary))
            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")	
	
            V2summary <- as.matrix(apply(V2.p, 2, summary))
            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")	
	
            V3summary <- as.matrix(apply(V3.p, 2, summary))
            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")	
    		
	
            if(storeV[1] == TRUE)
            {
                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
            }
            if(storeV[2] == TRUE)
            {
                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
            }
            if(storeV[3] == TRUE)
            {
                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
            }
    
            save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
    	
	
            if(p1 > 0){
                covNames1 = colnames(survData)[c(6:(5+p1))]
                }	
            if(p1 == 0){
                covNames1 = NULL
                }
        
            if(p2 > 0){
                covNames2 = colnames(survData)[c((5+p1+1):(5+p1+p2))]
                }	
            if(p2 == 0){
                covNames2 = NULL
                }     
            if(p3 > 0){
                covNames3 = colnames(survData)[c((5+p1+p2+1):(5+p1+p2+p3))]
                }	
            if(p3 == 0){
                covNames3 = NULL
                }         
    
    
            ### posterior predictive checks ###
    
            ## 1. log pseudo-marginal likelihood
    
            invLH.p <- matrix(mcmc$invLH, nrow = n, byrow = TRUE)
    
            cpo     <- 1/invLH.p
    
            LPML <- sum(log(cpo))
    
            # or
    
            lpml.p <- matrix(mcmc$lpml, nrow = nStore, byrow = TRUE)
    
            lpml2.p <- matrix(mcmc$lpml2, nrow = nStore, byrow = TRUE)
    
            ## 2. deviance information criterion
    
            gamma_mean <- matrix(mcmc$gammaP, nrow = n, byrow = TRUE)
            dev.p 	   <- matrix(mcmc$dev, nrow = nStore, byrow = TRUE)
            Dbar        <- mean(dev.p)
            pD          <- Dbar - (-2*mcmc$logLH_fin)
    
            DIC <- pD + Dbar
    
    

    
    
            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, nu2.p = nu2.p, nu3.p = nu3.p, theta.p = theta.p, V1.p = V1.p, V2.p = V2.p, V3.p = V3.p, c.p = c.p, mu.p = mu.p, Sigma.p = Sigma.p, tau.p = tau.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, accept.nu2 = accept.nu2, accept.nu3 = accept.nu3, accept.gamma = accept.gamma, lastChgProp = lastChgProp, mhGam_chk = mhGam_chk, accept.V = accept.V, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, dev.p = dev.p, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, cpo = cpo, LPML = LPML, gamma_mean = gamma_mean, DIC = DIC, val = mcmc$logLH_fin, lpml.p = lpml.p, lpml2.p = lpml2.p)            
            
            
            }  ## end: if Weibull-DPM-SM
            
        } ## end: if Weibull-DPM                                                            

    } ## end: if Weibull
			

        
        
        
        
        

        
        
    # hz.type = "PEM"
        
    if(hz.type == "PEM"){
            
        # re.type = "MVN"
            
        if(re.type == "MVN"){
                
            # model = "Markov"
                
            #######################################
            ############ PEM-MVN-M ############
            #######################################
                
            if(model == "Markov"){
                
                ###
                n	<- dim(survData)[1]
                p1	<- nCov[1]
                p2	<- nCov[2]
                p3	<- nCov[3]
                
                K1_max = mcmcParams[10]
                K2_max = mcmcParams[11]
                K3_max = mcmcParams[12]
                
                nTime_lambda1 = mcmcParams[16]
                nTime_lambda2 = mcmcParams[17]
                nTime_lambda3 = mcmcParams[18]
                
                time_lambda3	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3+1):(length(mcmcParams)-4)]
                time_lambda2	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2+1):(length(mcmcParams)-4-nTime_lambda3)]
                time_lambda1	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2-nTime_lambda1+1):(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2)]
                
                
                
                J	<- length(unique(survData[,5]))
                
                nj	<- rep(NA, J)
                
                for(i in 1:J){
                    nj[i]	<- length(which(survData[,5] == i))
                }
                
                ###
                
                nStore <- numReps/thin * (1 - burninPerc)
                
                K_1 <- startV[p1+p2+p3+1]
                K_2 <- startV[p1+p2+p3+2]
                K_3 <- startV[p1+p2+p3+3]
                
                num_s_propBI1 <- mcmcParams[7]
                num_s_propBI2 <- mcmcParams[8]
                num_s_propBI3 <- mcmcParams[9]
                
                
                startValues1 <- startV[1:(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1))]
                startValues2 <- startV[(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1)+1):length(startV)]
                startV <- c(startValues1, 1, 1, startValues2)
                
                mcmcParams1 <- mcmcParams[1:(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+1)]
                mcmcParams2 <- mcmcParams[(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+2):(length(mcmcParams)-1)]
                
                mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
                
                
                
                mcmc <- .C("BpeMvnCorScrmcmc",
                survData 		= as.double(as.matrix(survData)),
                n				= as.integer(n),
                p1				= as.integer(p1),
                p2				= as.integer(p2),
                p3				= as.integer(p3),
                J				= as.integer(J),
                nj				= as.double(nj),
                hyperParams 	= as.double(hyperParams),
                startValues 	= as.double(startV),
                mcmcParams		= as.double(mcmcParams),
                numReps			= as.integer(numReps),
                thin			= as.integer(thin),
                burninPerc      = as.double(burninPerc),
                nGam_save		= as.integer(nGam_save),
                samples_beta1 	= as.double(rep(0, nStore*p1)),
                samples_beta2 	= as.double(rep(0, nStore*p2)),
                samples_beta3 	= as.double(rep(0, nStore*p3)),
                samples_mu_lam1     = as.double(rep(0, nStore*1)),
                samples_mu_lam2     = as.double(rep(0, nStore*1)),
                samples_mu_lam3     = as.double(rep(0, nStore*1)),
                samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
                samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
                samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
                samples_K1          = as.double(rep(0, nStore*1)),
                samples_K2          = as.double(rep(0, nStore*1)),
                samples_K3          = as.double(rep(0, nStore*1)),
                samples_s1          = as.double(rep(0, nStore*(K1_max + 1))),
                samples_s2          = as.double(rep(0, nStore*(K2_max + 1))),
                samples_s3          = as.double(rep(0, nStore*(K3_max + 1))),
                samples_nu2 	= as.double(rep(0, nStore*1)),
                samples_nu3 	= as.double(rep(0, nStore*1)),
                samples_theta 	= as.double(rep(0, nStore*1)),
                samples_V1		= as.double(rep(0, nStore*J)),
                samples_V2		= as.double(rep(0, nStore*J)),
                samples_V3		= as.double(rep(0, nStore*J)),
                samples_Sigma_V	= as.double(rep(0, nStore*3*3)),
                samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                samples_gamma_last = as.double(rep(0, n)),
                samples_misc	= as.double(rep(0, (p1+p2+p3+9+n+1+n+J+J+J))),
                lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
                lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
                lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)),
                gammaP			= as.double(rep(0, n)),
                dev     			= as.double(rep(0, nStore*1)),
                invLH = as.double(rep(0, n)),
                logLH_fin = as.double(0),
                lpml     			= as.double(rep(0, nStore*1)),
                lpml2     			= as.double(rep(0, nStore*1)))
                
                
                
                if(p1 > 0){
                    beta1.p 		<- matrix(mcmc$samples_beta1, nrow = nStore, byrow = TRUE)
                }
                if(p1 == 0){
                    beta1.p 		<- NULL
                }
                if(p2 > 0){
                    beta2.p 		<- matrix(mcmc$samples_beta2, nrow = nStore, byrow = TRUE)
                }
                if(p2 == 0){
                    beta2.p 		<- NULL
                }
                if(p3 > 0){
                    beta3.p 		<- matrix(mcmc$samples_beta3, nrow = nStore, byrow = TRUE)
                }
                if(p3 == 0){
                    beta3.p 		<- NULL
                }
                
                lambda1.fin 	<- matrix(mcmc$lambda1_fin, nrow = nStore, byrow = TRUE)
                lambda2.fin 	<- matrix(mcmc$lambda2_fin, nrow = nStore, byrow = TRUE)
                lambda3.fin 	<- matrix(mcmc$lambda3_fin, nrow = nStore, byrow = TRUE)
                
                mu_lam1.p 		<- matrix(mcmc$samples_mu_lam1, nrow = nStore, byrow = TRUE)
                mu_lam2.p 		<- matrix(mcmc$samples_mu_lam2, nrow = nStore, byrow = TRUE)
                mu_lam3.p 		<- matrix(mcmc$samples_mu_lam3, nrow = nStore, byrow = TRUE)
                sigSq_lam1.p 	<- matrix(mcmc$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)
                sigSq_lam2.p 	<- matrix(mcmc$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)
                sigSq_lam3.p 	<- matrix(mcmc$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)
                
                K1.p 			<- matrix(mcmc$samples_K1, nrow = nStore, byrow = TRUE)
                K2.p 			<- matrix(mcmc$samples_K2, nrow = nStore, byrow = TRUE)
                K3.p 			<- matrix(mcmc$samples_K3, nrow = nStore, byrow = TRUE)
                s1.p 			<- matrix(mcmc$samples_s1, nrow = nStore, byrow = TRUE)
                s2.p 			<- matrix(mcmc$samples_s2, nrow = nStore, byrow = TRUE)
                s3.p 			<- matrix(mcmc$samples_s3, nrow = nStore, byrow = TRUE)
                
                theta.p 		<- matrix(mcmc$samples_theta, nrow = nStore, byrow = TRUE)
                gamma.p 		<- matrix(mcmc$samples_gamma, nrow = nStore, byrow = TRUE)
                
                V1.p            <- matrix(mcmc$samples_V1, nrow = nStore, byrow = TRUE)
                V2.p            <- matrix(mcmc$samples_V2, nrow = nStore, byrow = TRUE)
                V3.p            <- matrix(mcmc$samples_V3, nrow = nStore, byrow = TRUE)
                Sigma_V.p		<- array(as.vector(mcmc$samples_Sigma_V), c(3, 3, nStore))
                
                if(p1 > 0){
                    accept.beta1 	<- as.vector(mcmc$samples_misc[1:(p1)])
                }
                if(p1 == 0){
                    accept.beta1 	<- NULL
                }
                if(p2 > 0){
                    accept.beta2 	<- as.vector(mcmc$samples_misc[(p1+1):(p1+p2)])
                }
                if(p2 == 0){
                    accept.beta2 	<- NULL
                }
                if(p3 > 0){
                    accept.beta3 	<- as.vector(mcmc$samples_misc[(p1+p2+1):(p1+p2+p3)])
                }
                if(p3 == 0){
                    accept.beta3 	<- NULL
                }
                
                accept.BI1		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+1])
                accept.DI1		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+2])
                accept.BI2		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+3])
                accept.DI2		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+4])
                accept.BI3		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+5])
                accept.DI3		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+6])
                accept.theta	<- as.vector(mcmc$samples_misc[(p1+p2+p3+7)])
                accept.V1       <- as.vector(mcmc$samples_misc[(p1+p2+p3+10+n+n+1):(p1+p2+p3+10+n+n+J)])
                accept.V2       <- as.vector(mcmc$samples_misc[(p1+p2+p3+10+n+n+J+1):(p1+p2+p3+10+n+n+J+J)])
                accept.V3       <- as.vector(mcmc$samples_misc[(p1+p2+p3+10+n+n+J+J+1):(p1+p2+p3+10+n+n+J+J+J)])


                
                V1summary <- as.matrix(apply(V1.p, 2, summary))
                V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                
                V2summary <- as.matrix(apply(V2.p, 2, summary))
                V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                
                V3summary <- as.matrix(apply(V3.p, 2, summary))
                V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")			
                
                if(storeV[1] == TRUE)
                {
                    save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                }
                if(storeV[2] == TRUE)
                {
                    save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                }
                if(storeV[3] == TRUE)
                {
                    save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                }
                
                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                
                
                if(p1 > 0){
                    covNames1 = colnames(survData)[c(6:(5+p1))]
                }	
                if(p1 == 0){
                    covNames1 = NULL
                }
                
                if(p2 > 0){
                    covNames2 = colnames(survData)[c((5+p1+1):(5+p1+p2))]
                }	
                if(p2 == 0){
                    covNames2 = NULL
                }     
                if(p3 > 0){
                    covNames3 = colnames(survData)[c((5+p1+p2+1):(5+p1+p2+p3))]
                }	
                if(p3 == 0){
                    covNames3 = NULL
                }       
                
                
                
                ### posterior predictive checks ###
                
                ## 1. log pseudo-marginal likelihood
                
                invLH.p <- matrix(mcmc$invLH, nrow = n, byrow = TRUE)		
                
                cpo     <- 1/invLH.p
                
                LPML <- sum(log(cpo))
                
                # or
                
                lpml.p <- matrix(mcmc$lpml, nrow = nStore, byrow = TRUE)
                
                lpml2.p <- matrix(mcmc$lpml2, nrow = nStore, byrow = TRUE)	
                
                ## 2. deviance information criterion
                
                gamma_mean <- matrix(mcmc$gammaP, nrow = n, byrow = TRUE)
                dev.p 	   <- matrix(mcmc$dev, nrow = nStore, byrow = TRUE)	
                Dbar        <- mean(dev.p)
                pD          <- Dbar - (-2*mcmc$logLH_fin)
                
                DIC <- pD + Dbar
                
          		
                ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, K1.p = K1.p, K2.p = K2.p, K3.p = K3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, theta.p = theta.p, Sigma_V.p = Sigma_V.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, accept.V1 = accept.V1, accept.V2 = accept.V2, accept.V3 = accept.V3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, dev.p = dev.p, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, cpo = cpo, LPML = LPML, gamma_mean = gamma_mean, DIC = DIC, val = mcmc$logLH_fin, lpml.p = lpml.p, lpml2.p = lpml2.p)
     
                
            }   ## end: if PEM-MVN-M
            
            # model = "semi-Markov"
            
            #######################################
            ############ PEM-MVN-SM ############
            #######################################
            
            if(model == "semi-Markov"){
                
                ###
                n	<- dim(survData)[1]
                p1	<- nCov[1]
                p2	<- nCov[2]
                p3	<- nCov[3]
                
                K1_max = mcmcParams[10]
                K2_max = mcmcParams[11]
                K3_max = mcmcParams[12]
                
                nTime_lambda1 = mcmcParams[16]
                nTime_lambda2 = mcmcParams[17]
                nTime_lambda3 = mcmcParams[18]
                
                time_lambda3	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3+1):(length(mcmcParams)-4)]
                time_lambda2	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2+1):(length(mcmcParams)-4-nTime_lambda3)]
                time_lambda1	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2-nTime_lambda1+1):(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2)]
                
                
                J	<- length(unique(survData[,5]))
                
                nj	<- rep(NA, J)
                
                for(i in 1:J){
                    nj[i]	<- length(which(survData[,5] == i))
                }
                
                ###
                
                nStore <- numReps/thin * (1 - burninPerc)
                
                K_1 <- startV[p1+p2+p3+1]
                K_2 <- startV[p1+p2+p3+2]
                K_3 <- startV[p1+p2+p3+3]
                
                num_s_propBI1 <- mcmcParams[7]
                num_s_propBI2 <- mcmcParams[8]
                num_s_propBI3 <- mcmcParams[9]
                
                
                startValues1 <- startV[1:(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1))]
                startValues2 <- startV[(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1)+1):length(startV)]
                startV <- c(startValues1, 1, 1, startValues2)
                
                mcmcParams1 <- mcmcParams[1:(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+1)]
                mcmcParams2 <- mcmcParams[(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+2):(length(mcmcParams)-1)]
                
                mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
                
                
                
                mcmc <- .C("BpeMvnCorScrSMmcmc",
                survData 		= as.double(as.matrix(survData)),
                n				= as.integer(n),
                p1				= as.integer(p1),
                p2				= as.integer(p2),
                p3				= as.integer(p3),
                J				= as.integer(J),
                nj				= as.double(nj),
                hyperParams 	= as.double(hyperParams),
                startValues 	= as.double(startV),
                mcmcParams		= as.double(mcmcParams),
                numReps			= as.integer(numReps),
                thin			= as.integer(thin),
                burninPerc      = as.double(burninPerc),
                nGam_save		= as.integer(nGam_save),
                samples_beta1 	= as.double(rep(0, nStore*p1)),
                samples_beta2 	= as.double(rep(0, nStore*p2)),
                samples_beta3 	= as.double(rep(0, nStore*p3)),
                samples_mu_lam1     = as.double(rep(0, nStore*1)),
                samples_mu_lam2     = as.double(rep(0, nStore*1)),
                samples_mu_lam3     = as.double(rep(0, nStore*1)),
                samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
                samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
                samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
                samples_K1          = as.double(rep(0, nStore*1)),
                samples_K2          = as.double(rep(0, nStore*1)),
                samples_K3          = as.double(rep(0, nStore*1)),
                samples_s1          = as.double(rep(0, nStore*(K1_max + 1))),
                samples_s2          = as.double(rep(0, nStore*(K2_max + 1))),
                samples_s3          = as.double(rep(0, nStore*(K3_max + 1))),
                samples_nu2 	= as.double(rep(0, nStore*1)),
                samples_nu3 	= as.double(rep(0, nStore*1)),
                samples_theta 	= as.double(rep(0, nStore*1)),
                samples_V1		= as.double(rep(0, nStore*J)),
                samples_V2		= as.double(rep(0, nStore*J)),
                samples_V3		= as.double(rep(0, nStore*J)),
                samples_Sigma_V	= as.double(rep(0, nStore*3*3)),
                samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                samples_gamma_last = as.double(rep(0, n)),
                samples_misc	= as.double(rep(0, (p1+p2+p3+9+n+1+n+J+J+J))),
                lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
                lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
                lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)),
                gammaP			= as.double(rep(0, n)),
                dev     			= as.double(rep(0, nStore*1)),
                invLH = as.double(rep(0, n)),
                logLH_fin = as.double(0),
                lpml     			= as.double(rep(0, nStore*1)))

                
                
                if(p1 > 0){
                    beta1.p 		<- matrix(mcmc$samples_beta1, nrow = nStore, byrow = TRUE)
                }
                if(p1 == 0){
                    beta1.p 		<- NULL
                }
                if(p2 > 0){
                    beta2.p 		<- matrix(mcmc$samples_beta2, nrow = nStore, byrow = TRUE)
                }
                if(p2 == 0){
                    beta2.p 		<- NULL
                }
                if(p3 > 0){
                    beta3.p 		<- matrix(mcmc$samples_beta3, nrow = nStore, byrow = TRUE)
                }
                if(p3 == 0){
                    beta3.p 		<- NULL
                }
                
                lambda1.fin 	<- matrix(mcmc$lambda1_fin, nrow = nStore, byrow = TRUE)
                lambda2.fin 	<- matrix(mcmc$lambda2_fin, nrow = nStore, byrow = TRUE)
                lambda3.fin 	<- matrix(mcmc$lambda3_fin, nrow = nStore, byrow = TRUE)
                
                mu_lam1.p 		<- matrix(mcmc$samples_mu_lam1, nrow = nStore, byrow = TRUE)
                mu_lam2.p 		<- matrix(mcmc$samples_mu_lam2, nrow = nStore, byrow = TRUE)
                mu_lam3.p 		<- matrix(mcmc$samples_mu_lam3, nrow = nStore, byrow = TRUE)
                sigSq_lam1.p 	<- matrix(mcmc$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)
                sigSq_lam2.p 	<- matrix(mcmc$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)
                sigSq_lam3.p 	<- matrix(mcmc$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)
                
                K1.p 			<- matrix(mcmc$samples_K1, nrow = nStore, byrow = TRUE)
                K2.p 			<- matrix(mcmc$samples_K2, nrow = nStore, byrow = TRUE)
                K3.p 			<- matrix(mcmc$samples_K3, nrow = nStore, byrow = TRUE)
                s1.p 			<- matrix(mcmc$samples_s1, nrow = nStore, byrow = TRUE)
                s2.p 			<- matrix(mcmc$samples_s2, nrow = nStore, byrow = TRUE)
                s3.p 			<- matrix(mcmc$samples_s3, nrow = nStore, byrow = TRUE)
                
                theta.p 		<- matrix(mcmc$samples_theta, nrow = nStore, byrow = TRUE)
                gamma.p 		<- matrix(mcmc$samples_gamma, nrow = nStore, byrow = TRUE)
                
                V1.p            <- matrix(mcmc$samples_V1, nrow = nStore, byrow = TRUE)
                V2.p            <- matrix(mcmc$samples_V2, nrow = nStore, byrow = TRUE)
                V3.p            <- matrix(mcmc$samples_V3, nrow = nStore, byrow = TRUE)
                Sigma_V.p		<- array(as.vector(mcmc$samples_Sigma_V), c(3, 3, nStore))
                
                if(p1 > 0){
                    accept.beta1 	<- as.vector(mcmc$samples_misc[1:(p1)])
                }
                if(p1 == 0){
                    accept.beta1 	<- NULL
                }
                if(p2 > 0){
                    accept.beta2 	<- as.vector(mcmc$samples_misc[(p1+1):(p1+p2)])
                }
                if(p2 == 0){
                    accept.beta2 	<- NULL
                }
                if(p3 > 0){
                    accept.beta3 	<- as.vector(mcmc$samples_misc[(p1+p2+1):(p1+p2+p3)])
                }
                if(p3 == 0){
                    accept.beta3 	<- NULL
                }
                
                accept.BI1		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+1])
                accept.DI1		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+2])
                accept.BI2		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+3])
                accept.DI2		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+4])
                accept.BI3		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+5])
                accept.DI3		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+6])
                accept.theta	<- as.vector(mcmc$samples_misc[(p1+p2+p3+7)])
                accept.V1       <- as.vector(mcmc$samples_misc[(p1+p2+p3+10+n+n+1):(p1+p2+p3+10+n+n+J)])
                accept.V2       <- as.vector(mcmc$samples_misc[(p1+p2+p3+10+n+n+J+1):(p1+p2+p3+10+n+n+J+J)])
                accept.V3       <- as.vector(mcmc$samples_misc[(p1+p2+p3+10+n+n+J+J+1):(p1+p2+p3+10+n+n+J+J+J)])
               
                
                V1summary <- as.matrix(apply(V1.p, 2, summary))
                V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                
                V2summary <- as.matrix(apply(V2.p, 2, summary))
                V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                
                V3summary <- as.matrix(apply(V3.p, 2, summary))
                V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")			
                
                if(storeV[1] == TRUE)
                {
                    save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                }
                if(storeV[2] == TRUE)
                {
                    save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                }
                if(storeV[3] == TRUE)
                {
                    save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                }
                
                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                
                
                if(p1 > 0){
                    covNames1 = colnames(survData)[c(6:(5+p1))]
                }	
                if(p1 == 0){
                    covNames1 = NULL
                }
                
                if(p2 > 0){
                    covNames2 = colnames(survData)[c((5+p1+1):(5+p1+p2))]
                }	
                if(p2 == 0){
                    covNames2 = NULL
                }     
                if(p3 > 0){
                    covNames3 = colnames(survData)[c((5+p1+p2+1):(5+p1+p2+p3))]
                }	
                if(p3 == 0){
                    covNames3 = NULL
                }
                
                
                
                ### posterior predictive checks ###
                
                ## 1. log pseudo-marginal likelihood
                
                invLH.p <- matrix(mcmc$invLH, nrow = n, byrow = TRUE)		
                
                cpo     <- 1/invLH.p
                
                LPML <- sum(log(cpo))
                
                # or
                
                lpml.p <- matrix(mcmc$lpml, nrow = nStore, byrow = TRUE)
                
                
                ## 2. deviance information criterion
                
                gamma_mean <- matrix(mcmc$gammaP, nrow = n, byrow = TRUE)
                dev.p 	   <- matrix(mcmc$dev, nrow = nStore, byrow = TRUE)	
                Dbar        <- mean(dev.p)
                pD          <- Dbar - (-2*mcmc$logLH_fin)
                
                DIC <- pD + Dbar              
                                
                ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, K1.p = K1.p, K2.p = K2.p, K3.p = K3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, theta.p = theta.p, Sigma_V.p = Sigma_V.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, accept.V1 = accept.V1, accept.V2 = accept.V2, accept.V3 = accept.V3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, dev.p = dev.p, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, cpo = cpo, LPML = LPML, gamma_mean = gamma_mean, DIC = DIC, val = mcmc$logLH_fin, lpml.p = lpml.p)
                
                
            }   ## end: if PEM-MVN-SM
            
        }   ## end: if PEM-MVN
        
        
        # re.type = "DPM"
        
        if(re.type == "DPM"){
            
            # model = "Markov"
            
            #######################################
            ############ PEM-DPM-M ############
            #######################################
            
            if(model == "Markov"){
                
                ###
                n	<- dim(survData)[1]
                p1	<- nCov[1]
                p2	<- nCov[2]
                p3	<- nCov[3]
                
                K1_max = mcmcParams[10]
                K2_max = mcmcParams[11]
                K3_max = mcmcParams[12]
                
                nTime_lambda1 = mcmcParams[16]
                nTime_lambda2 = mcmcParams[17]
                nTime_lambda3 = mcmcParams[18]
                
                time_lambda3	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3+1):(length(mcmcParams)-4)]
                time_lambda2	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2+1):(length(mcmcParams)-4-nTime_lambda3)]
                time_lambda1	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2-nTime_lambda1+1):(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2)]
                
                
                
                J	<- length(unique(survData[,5]))
                
                nj	<- rep(NA, J)
                
                for(i in 1:J){
                    nj[i]	<- length(which(survData[,5] == i))
                }
                
                ###
                
                nStore <- numReps/thin * (1 - burninPerc)
                
                
                K_1 <- startV[p1+p2+p3+1]
                K_2 <- startV[p1+p2+p3+2]
                K_3 <- startV[p1+p2+p3+3]
                
                num_s_propBI1 <- mcmcParams[7]
                num_s_propBI2 <- mcmcParams[8]
                num_s_propBI3 <- mcmcParams[9]
                
                
                startValues1 <- startV[1:(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1))]
                startValues2 <- startV[(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1)+1):length(startV)]
                startV <- c(startValues1, 1, 1, startValues2)
                
                mcmcParams1 <- mcmcParams[1:(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+1)]
                mcmcParams2 <- mcmcParams[(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+2):(length(mcmcParams)-1)]
                
                mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
                
                
                mcmc <- .C("BpeDpCorScrmcmc",
                survData 		= as.double(as.matrix(survData)),
                n				= as.integer(n),
                p1				= as.integer(p1),
                p2				= as.integer(p2),
                p3				= as.integer(p3),
                J				= as.integer(J),
                nj				= as.double(nj),
                hyperParams 	= as.double(hyperParams),
                startValues 	= as.double(startV),
                mcmcParams		= as.double(mcmcParams),
                numReps			= as.integer(numReps),
                thin			= as.integer(thin),
                burninPerc      = as.double(burninPerc),
                nGam_save		= as.integer(nGam_save),
                samples_beta1 	= as.double(rep(0, nStore*p1)),
                samples_beta2 	= as.double(rep(0, nStore*p2)),
                samples_beta3 	= as.double(rep(0, nStore*p3)),
                samples_mu_lam1     = as.double(rep(0, nStore*1)),
                samples_mu_lam2     = as.double(rep(0, nStore*1)),
                samples_mu_lam3     = as.double(rep(0, nStore*1)),
                samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
                samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
                samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
                samples_K1          = as.double(rep(0, nStore*1)),
                samples_K2          = as.double(rep(0, nStore*1)),
                samples_K3          = as.double(rep(0, nStore*1)),
                samples_s1          = as.double(rep(0, nStore*(K1_max + 1))),
                samples_s2          = as.double(rep(0, nStore*(K2_max + 1))),
                samples_s3          = as.double(rep(0, nStore*(K3_max + 1))),
                samples_nu2 	= as.double(rep(0, nStore*1)),
                samples_nu3 	= as.double(rep(0, nStore*1)),
                samples_theta 	= as.double(rep(0, nStore*1)),
                samples_V1		= as.double(rep(0, nStore*J)),
                samples_V2		= as.double(rep(0, nStore*J)),
                samples_V3		= as.double(rep(0, nStore*J)),
                samples_c		= as.double(rep(0, nStore*J)),
                samples_mu		= as.double(rep(0, nStore*3*J)),
                samples_Sigma	= as.double(rep(0, nStore*3*3*J)),
                samples_tau     = as.double(rep(0, nStore*1)),
                samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                samples_gamma_last = as.double(rep(0, n)),
                samples_misc	= as.double(rep(0, (p1+p2+p3+9+n+1+n+J))),
                lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
                lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
                lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)),
                gammaP			= as.double(rep(0, n)),
                dev     			= as.double(rep(0, nStore*1)),
                invLH = as.double(rep(0, n)),
                logLH_fin = as.double(0),
                lpml     			= as.double(rep(0, nStore*1)))
				
                
                
                if(p1 > 0){
                    beta1.p 		<- matrix(mcmc$samples_beta1, nrow = nStore, byrow = TRUE)
                }
                if(p1 == 0){
                    beta1.p 		<- NULL
                }
                if(p2 > 0){
                    beta2.p 		<- matrix(mcmc$samples_beta2, nrow = nStore, byrow = TRUE)
                }
                if(p2 == 0){
                    beta2.p 		<- NULL
                }
                if(p3 > 0){
                    beta3.p 		<- matrix(mcmc$samples_beta3, nrow = nStore, byrow = TRUE)
                }
                if(p3 == 0){
                    beta3.p 		<- NULL
                }
                
                lambda1.fin 	<- matrix(mcmc$lambda1_fin, nrow = nStore, byrow = TRUE)
                lambda2.fin 	<- matrix(mcmc$lambda2_fin, nrow = nStore, byrow = TRUE)
                lambda3.fin 	<- matrix(mcmc$lambda3_fin, nrow = nStore, byrow = TRUE)
                
                mu_lam1.p 		<- matrix(mcmc$samples_mu_lam1, nrow = nStore, byrow = TRUE)
                mu_lam2.p 		<- matrix(mcmc$samples_mu_lam2, nrow = nStore, byrow = TRUE)
                mu_lam3.p 		<- matrix(mcmc$samples_mu_lam3, nrow = nStore, byrow = TRUE)
                sigSq_lam1.p 	<- matrix(mcmc$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)
                sigSq_lam2.p 	<- matrix(mcmc$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)
                sigSq_lam3.p 	<- matrix(mcmc$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)
                
                K1.p 			<- matrix(mcmc$samples_K1, nrow = nStore, byrow = TRUE)
                K2.p 			<- matrix(mcmc$samples_K2, nrow = nStore, byrow = TRUE)
                K3.p 			<- matrix(mcmc$samples_K3, nrow = nStore, byrow = TRUE)
                s1.p 			<- matrix(mcmc$samples_s1, nrow = nStore, byrow = TRUE)
                s2.p 			<- matrix(mcmc$samples_s2, nrow = nStore, byrow = TRUE)
                s3.p 			<- matrix(mcmc$samples_s3, nrow = nStore, byrow = TRUE)
                
                theta.p 		<- matrix(mcmc$samples_theta, nrow = nStore, byrow = TRUE)
                gamma.p 		<- matrix(mcmc$samples_gamma, nrow = nStore, byrow = TRUE)
                
                V1.p            <- matrix(mcmc$samples_V1, nrow = nStore, byrow = TRUE)
                V2.p            <- matrix(mcmc$samples_V2, nrow = nStore, byrow = TRUE)
                V3.p            <- matrix(mcmc$samples_V3, nrow = nStore, byrow = TRUE)
                c.p            <- matrix(mcmc$samples_c, nrow = nStore, byrow = TRUE)
                mu.p            <- matrix(mcmc$samples_mu, nrow = nStore, byrow = TRUE)
                Sigma.p         <- array(as.vector(mcmc$samples_Sigma), c(3, 3 * J, nStore))
                tau.p           <- matrix(mcmc$samples_tau, nrow = nStore, byrow = TRUE)
                
                if(p1 > 0){
                    accept.beta1 	<- as.vector(mcmc$samples_misc[1:(p1)])
                }
                if(p1 == 0){
                    accept.beta1 	<- NULL
                }
                if(p2 > 0){
                    accept.beta2 	<- as.vector(mcmc$samples_misc[(p1+1):(p1+p2)])
                }
                if(p2 == 0){
                    accept.beta2 	<- NULL
                }
                if(p3 > 0){
                    accept.beta3 	<- as.vector(mcmc$samples_misc[(p1+p2+1):(p1+p2+p3)])
                }
                if(p3 == 0){
                    accept.beta3 	<- NULL
                }
                
                accept.BI1		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+1])
                accept.DI1		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+2])
                accept.BI2		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+3])
                accept.DI2		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+4])
                accept.BI3		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+5])
                accept.DI3		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+6])
                accept.theta	<- as.vector(mcmc$samples_misc[(p1+p2+p3+7)])
                accept.V       <- as.vector(mcmc$samples_misc[(p1+p2+p3+10+n+n+1):(p1+p2+p3+10+n+n+J)])
                                
                V1summary <- as.matrix(apply(V1.p, 2, summary))
                V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                
                V2summary <- as.matrix(apply(V2.p, 2, summary))
                V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                
                V3summary <- as.matrix(apply(V3.p, 2, summary))
                V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")			
                
                if(storeV[1] == TRUE)
                {
                    save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                }
                if(storeV[2] == TRUE)
                {
                    save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                }
                if(storeV[3] == TRUE)
                {
                    save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                }
                
                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                
                
                if(p1 > 0){
                    covNames1 = colnames(survData)[c(6:(5+p1))]
                }	
                if(p1 == 0){
                    covNames1 = NULL
                }
                
                if(p2 > 0){
                    covNames2 = colnames(survData)[c((5+p1+1):(5+p1+p2))]
                }	
                if(p2 == 0){
                    covNames2 = NULL
                }     
                if(p3 > 0){
                    covNames3 = colnames(survData)[c((5+p1+p2+1):(5+p1+p2+p3))]
                }	
                if(p3 == 0){
                    covNames3 = NULL
                }
                
                
                
                ### posterior predictive checks ###
                
                ## 1. log pseudo-marginal likelihood
                
                invLH.p <- matrix(mcmc$invLH, nrow = n, byrow = TRUE)		
                
                cpo     <- 1/invLH.p
                
                LPML <- sum(log(cpo))
                
                # or
                
                lpml.p <- matrix(mcmc$lpml, nrow = nStore, byrow = TRUE)	
                
                
                ## 2. deviance information criterion
                
                gamma_mean <- matrix(mcmc$gammaP, nrow = n, byrow = TRUE)
                dev.p 	   <- matrix(mcmc$dev, nrow = nStore, byrow = TRUE)	
                Dbar        <- mean(dev.p)
                pD          <- Dbar - (-2*mcmc$logLH_fin)
                
                DIC <- pD + Dbar        
                
                
          		
                ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, K1.p = K1.p, K2.p = K2.p, K3.p = K3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, theta.p = theta.p, c.p = c.p, mu.p = mu.p, Sigma.p = Sigma.p, tau.p = tau.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, accept.V = accept.V, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, dev.p = dev.p, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, cpo = cpo, LPML = LPML, gamma_mean = gamma_mean, DIC = DIC, val = mcmc$logLH_fin, lpml.p = lpml.p)
           
            
            }   ## end: if PEM-DPM-M
            
            
            
            # model = "semi-Markov"
            
            #######################################
            ############ PEM-DPM-SM ############
            #######################################
            
            if(model == "semi-Markov"){
                
                ###
                n	<- dim(survData)[1]
                p1	<- nCov[1]
                p2	<- nCov[2]
                p3	<- nCov[3]
                
                K1_max = mcmcParams[10]
                K2_max = mcmcParams[11]
                K3_max = mcmcParams[12]
                
                nTime_lambda1 = mcmcParams[16]
                nTime_lambda2 = mcmcParams[17]
                nTime_lambda3 = mcmcParams[18]
                
                time_lambda3	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3+1):(length(mcmcParams)-4)]
                time_lambda2	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2+1):(length(mcmcParams)-4-nTime_lambda3)]
                time_lambda1	= mcmcParams[(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2-nTime_lambda1+1):(length(mcmcParams)-4-nTime_lambda3-nTime_lambda2)]
                
                
                J	<- length(unique(survData[,5]))
                
                nj	<- rep(NA, J)
                
                for(i in 1:J){
                    nj[i]	<- length(which(survData[,5] == i))
                }

                ###
                
                nStore <- numReps/thin * (1 - burninPerc)
                
                K_1 <- startV[p1+p2+p3+1]
                K_2 <- startV[p1+p2+p3+2]
                K_3 <- startV[p1+p2+p3+3]
                
                num_s_propBI1 <- mcmcParams[7]
                num_s_propBI2 <- mcmcParams[8]
                num_s_propBI3 <- mcmcParams[9]
                
                
                startValues1 <- startV[1:(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1))]
                startValues2 <- startV[(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1)+1):length(startV)]
                startV <- c(startValues1, 1, 1, startValues2)
                
                mcmcParams1 <- mcmcParams[1:(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+1)]
                mcmcParams2 <- mcmcParams[(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+2):(length(mcmcParams)-1)]
                
                mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
                
                mcmc <- .C("BpeDpCorScrSMmcmc",
                survData 		= as.double(as.matrix(survData)),
                n				= as.integer(n),
                p1				= as.integer(p1),
                p2				= as.integer(p2),
                p3				= as.integer(p3),
                J				= as.integer(J),
                nj				= as.double(nj),
                hyperParams 	= as.double(hyperParams),
                startValues 	= as.double(startV),
                mcmcParams		= as.double(mcmcParams),
                numReps			= as.integer(numReps),
                thin			= as.integer(thin),
                burninPerc      = as.double(burninPerc),
                nGam_save		= as.integer(nGam_save),
                samples_beta1 	= as.double(rep(0, nStore*p1)),
                samples_beta2 	= as.double(rep(0, nStore*p2)),
                samples_beta3 	= as.double(rep(0, nStore*p3)),
                samples_mu_lam1     = as.double(rep(0, nStore*1)),
                samples_mu_lam2     = as.double(rep(0, nStore*1)),
                samples_mu_lam3     = as.double(rep(0, nStore*1)),
                samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
                samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
                samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
                samples_K1          = as.double(rep(0, nStore*1)),
                samples_K2          = as.double(rep(0, nStore*1)),
                samples_K3          = as.double(rep(0, nStore*1)),
                samples_s1          = as.double(rep(0, nStore*(K1_max + 1))),
                samples_s2          = as.double(rep(0, nStore*(K2_max + 1))),
                samples_s3          = as.double(rep(0, nStore*(K3_max + 1))),
                samples_nu2 	= as.double(rep(0, nStore*1)),
                samples_nu3 	= as.double(rep(0, nStore*1)),
                samples_theta 	= as.double(rep(0, nStore*1)),
                samples_V1		= as.double(rep(0, nStore*J)),
                samples_V2		= as.double(rep(0, nStore*J)),
                samples_V3		= as.double(rep(0, nStore*J)),
                samples_c		= as.double(rep(0, nStore*J)),
                samples_mu		= as.double(rep(0, nStore*3*J)),
                samples_Sigma	= as.double(rep(0, nStore*3*3*J)),
                samples_tau     = as.double(rep(0, nStore*1)),
                samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                samples_gamma_last = as.double(rep(0, n)),
                samples_misc	= as.double(rep(0, (p1+p2+p3+9+n+1+n+J))),
                lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
                lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
                lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)),
                gammaP			= as.double(rep(0, n)),
                dev     			= as.double(rep(0, nStore*1)),
                invLH = as.double(rep(0, n)),
                logLH_fin = as.double(0),
                lpml     			= as.double(rep(0, nStore*1)))
				                
                
                if(p1 > 0){
                    beta1.p 		<- matrix(mcmc$samples_beta1, nrow = nStore, byrow = TRUE)
                }
                if(p1 == 0){
                    beta1.p 		<- NULL
                }
                if(p2 > 0){
                    beta2.p 		<- matrix(mcmc$samples_beta2, nrow = nStore, byrow = TRUE)
                }
                if(p2 == 0){
                    beta2.p 		<- NULL
                }
                if(p3 > 0){
                    beta3.p 		<- matrix(mcmc$samples_beta3, nrow = nStore, byrow = TRUE)
                }
                if(p3 == 0){
                    beta3.p 		<- NULL
                }
                
                lambda1.fin 	<- matrix(mcmc$lambda1_fin, nrow = nStore, byrow = TRUE)
                lambda2.fin 	<- matrix(mcmc$lambda2_fin, nrow = nStore, byrow = TRUE)
                lambda3.fin 	<- matrix(mcmc$lambda3_fin, nrow = nStore, byrow = TRUE)
                
                mu_lam1.p 		<- matrix(mcmc$samples_mu_lam1, nrow = nStore, byrow = TRUE)
                mu_lam2.p 		<- matrix(mcmc$samples_mu_lam2, nrow = nStore, byrow = TRUE)
                mu_lam3.p 		<- matrix(mcmc$samples_mu_lam3, nrow = nStore, byrow = TRUE)
                sigSq_lam1.p 	<- matrix(mcmc$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)
                sigSq_lam2.p 	<- matrix(mcmc$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)
                sigSq_lam3.p 	<- matrix(mcmc$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)
                
                K1.p 			<- matrix(mcmc$samples_K1, nrow = nStore, byrow = TRUE)
                K2.p 			<- matrix(mcmc$samples_K2, nrow = nStore, byrow = TRUE)
                K3.p 			<- matrix(mcmc$samples_K3, nrow = nStore, byrow = TRUE)
                s1.p 			<- matrix(mcmc$samples_s1, nrow = nStore, byrow = TRUE)
                s2.p 			<- matrix(mcmc$samples_s2, nrow = nStore, byrow = TRUE)
                s3.p 			<- matrix(mcmc$samples_s3, nrow = nStore, byrow = TRUE)
                
                theta.p 		<- matrix(mcmc$samples_theta, nrow = nStore, byrow = TRUE)
                gamma.p 		<- matrix(mcmc$samples_gamma, nrow = nStore, byrow = TRUE)
                
                V1.p            <- matrix(mcmc$samples_V1, nrow = nStore, byrow = TRUE)
                V2.p            <- matrix(mcmc$samples_V2, nrow = nStore, byrow = TRUE)
                V3.p            <- matrix(mcmc$samples_V3, nrow = nStore, byrow = TRUE)
                c.p            <- matrix(mcmc$samples_c, nrow = nStore, byrow = TRUE)
                mu.p            <- matrix(mcmc$samples_mu, nrow = nStore, byrow = TRUE)
                Sigma.p         <- array(as.vector(mcmc$samples_Sigma), c(3, 3 * J, nStore))
                tau.p           <- matrix(mcmc$samples_tau, nrow = nStore, byrow = TRUE)
                
                if(p1 > 0){
                    accept.beta1 	<- as.vector(mcmc$samples_misc[1:(p1)])
                }
                if(p1 == 0){
                    accept.beta1 	<- NULL
                }
                if(p2 > 0){
                    accept.beta2 	<- as.vector(mcmc$samples_misc[(p1+1):(p1+p2)])
                }
                if(p2 == 0){
                    accept.beta2 	<- NULL
                }
                if(p3 > 0){
                    accept.beta3 	<- as.vector(mcmc$samples_misc[(p1+p2+1):(p1+p2+p3)])
                }
                if(p3 == 0){
                    accept.beta3 	<- NULL
                }
                
                accept.BI1		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+1])
                accept.DI1		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+2])
                accept.BI2		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+3])
                accept.DI2		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+4])
                accept.BI3		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+5])
                accept.DI3		<- as.vector(mcmc$samples_misc[(p1+p2+p3)+6])
                accept.theta	<- as.vector(mcmc$samples_misc[(p1+p2+p3+7)])
                accept.V       <- as.vector(mcmc$samples_misc[(p1+p2+p3+10+n+n+1):(p1+p2+p3+10+n+n+J)])
                
                
                V1summary <- as.matrix(apply(V1.p, 2, summary))
                V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                
                V2summary <- as.matrix(apply(V2.p, 2, summary))
                V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                
                V3summary <- as.matrix(apply(V3.p, 2, summary))
                V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")			
                
                if(storeV[1] == TRUE)
                {
                    save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                }
                if(storeV[2] == TRUE)
                {
                    save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                }
                if(storeV[3] == TRUE)
                {
                    save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                }
                
                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                
                
                if(p1 > 0){
                    covNames1 = colnames(survData)[c(6:(5+p1))]
                }	
                if(p1 == 0){
                    covNames1 = NULL
                }
                
                if(p2 > 0){
                    covNames2 = colnames(survData)[c((5+p1+1):(5+p1+p2))]
                }	
                if(p2 == 0){
                    covNames2 = NULL
                }     
                if(p3 > 0){
                    covNames3 = colnames(survData)[c((5+p1+p2+1):(5+p1+p2+p3))]
                }	
                if(p3 == 0){
                    covNames3 = NULL
                }
                
                
                
                
                ### posterior predictive checks ###
                
                ## 1. log pseudo-marginal likelihood
                
                invLH.p <- matrix(mcmc$invLH, nrow = n, byrow = TRUE)		
                
                cpo     <- 1/invLH.p
                
                LPML <- sum(log(cpo))
                
                # or
                
                lpml.p <- matrix(mcmc$lpml, nrow = nStore, byrow = TRUE)    
                
                
                
                ## 2. deviance information criterion
                
                gamma_mean <- matrix(mcmc$gammaP, nrow = n, byrow = TRUE)
                dev.p 	   <- matrix(mcmc$dev, nrow = nStore, byrow = TRUE)	
                Dbar        <- mean(dev.p)
                pD          <- Dbar - (-2*mcmc$logLH_fin)
                
                DIC <- pD + Dbar              
                
                
                
          		
                ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, K1.p = K1.p, K2.p = K2.p, K3.p = K3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, theta.p = theta.p, c.p = c.p, mu.p = mu.p, Sigma.p = Sigma.p, tau.p = tau.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, accept.V = accept.V, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, dev.p = dev.p, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, cpo = cpo, LPML = LPML, gamma_mean = gamma_mean, DIC = DIC, val = mcmc$logLH_fin, lpml.p = lpml.p)
                
                
            }   ## end: if PEM-DPM-SM
            
            
        }   ## end: if PEM-DPM
        
        
    }   ## end: if PEM
        

		chain = chain + 1	
        

	}	## end: while(chain <= nChain)
    
    
		ret[["setup"]]	<- list(nCov = nCov, hyperParams = hyperParams, startValues = startValues, mcmcParams = mcmcParams, nGam_save = nGam_save, numReps = numReps, thin = thin, path = path, burninPerc = burninPerc, hz.type = hz.type, re.type = re.type, model = model, nChain = nChain)	   
			   
		class(ret) <- "BayesIDcor"
		return(ret)
		
        
        
}  ## end: if(class(startValues) == "list" & length(startValues) == nChain)	



else{
	print("The 'startValues' should be the list of length equal to 'nChain'.")
}







} # end of function "BayesIDcor"





























