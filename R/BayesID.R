

### hz.type = "PEM" (semi-parametric) or "Weibull" (parametric)
### model = "Markov" or "semi-Markov"

BayesID <- function(survData, 
					nCov,					
					hyperParams,
					startValues,								
					mcmcParams,
					nGam_save,					
					numReps,
					thin,
					path,
					burninPerc=0.5,
					hz.type = "Weibull",
					model = "Markov",
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

	# hz.type = "PEM"
	
	if(hz.type == "PEM"){
		
		# model = "Markov"		
		
		if(model == "Markov"){
	
			###
			n	<- dim(survData)[1]
			p1	<- nCov[1]
			p2	<- nCov[2]
			p3	<- nCov[3]
	
			J1_max = mcmcParams[10]
			J2_max = mcmcParams[11]
			J3_max = mcmcParams[12]		
	
			nStore <- numReps/thin * (1 - burninPerc)
	
			nTime_lambda1 = mcmcParams[16]
			nTime_lambda2 = mcmcParams[17]
			nTime_lambda3 = mcmcParams[18]		
	
			time_lambda3	= mcmcParams[(length(mcmcParams)-1-nTime_lambda3+1):(length(mcmcParams)-1)]	
			time_lambda2	= mcmcParams[(length(mcmcParams)-1-nTime_lambda3-nTime_lambda2+1):(length(mcmcParams)-1-nTime_lambda3)]
			time_lambda1	= mcmcParams[(length(mcmcParams)-1-nTime_lambda3-nTime_lambda2-nTime_lambda1+1):(length(mcmcParams)-1-nTime_lambda3-nTime_lambda2)]			
	

			mcmc <- .C("BpeScrmcmc",
						survData            = as.double(as.matrix(survData)),
						n                   = as.integer(n),
						p1					= as.integer(p1),
						p2					= as.integer(p2),
						p3					= as.integer(p3),
						hyperParams         = as.double(hyperParams),
						startValues         = as.double(startV),
						mcmcParams          = as.double(mcmcParams),
						numReps             = as.integer(numReps),
						thin                = as.integer(thin),
						burninPerc          = as.double(burninPerc),
						nGam_save			= as.integer(nGam_save),
						samples_beta1 		= as.double(rep(0, nStore*p1)),
						samples_beta2 		= as.double(rep(0, nStore*p2)),
						samples_beta3 		= as.double(rep(0, nStore*p3)),
						samples_mu_lam1     = as.double(rep(0, nStore*1)),
						samples_mu_lam2     = as.double(rep(0, nStore*1)),
						samples_mu_lam3     = as.double(rep(0, nStore*1)),
						samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
						samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
						samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
						samples_J1          = as.double(rep(0, nStore*1)),
						samples_J2          = as.double(rep(0, nStore*1)),
						samples_J3          = as.double(rep(0, nStore*1)),
						samples_s1          = as.double(rep(0, nStore*(J1_max + 1))),
						samples_s2          = as.double(rep(0, nStore*(J2_max + 1))),
						samples_s3          = as.double(rep(0, nStore*(J3_max + 1))),
						samples_theta 		= as.double(rep(0, nStore*1)),
						samples_gamma 		= as.double(rep(0, nStore*nGam_save)),
						samples_misc        = as.double(rep(0, p1 + p2 + p3 + 7)),
						lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
						lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
						lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)))

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
			gamma.p 		<- matrix(mcmc$samples_gamma, nrow = nStore, byrow = TRUE)	
			theta.p 		<- matrix(mcmc$samples_theta, nrow = nStore, byrow = TRUE)						
			mu_lam1.p 		<- matrix(mcmc$samples_mu_lam1, nrow = nStore, byrow = TRUE)
			mu_lam2.p 		<- matrix(mcmc$samples_mu_lam2, nrow = nStore, byrow = TRUE)
			mu_lam3.p 		<- matrix(mcmc$samples_mu_lam3, nrow = nStore, byrow = TRUE)		
			sigSq_lam1.p 	<- matrix(mcmc$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)	
			sigSq_lam2.p 	<- matrix(mcmc$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)	
			sigSq_lam3.p 	<- matrix(mcmc$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)			
			J1.p 			<- matrix(mcmc$samples_J1, nrow = nStore, byrow = TRUE)
			J2.p 			<- matrix(mcmc$samples_J2, nrow = nStore, byrow = TRUE)
			J3.p 			<- matrix(mcmc$samples_J3, nrow = nStore, byrow = TRUE)		
			s1.p 			<- matrix(mcmc$samples_s1, nrow = nStore, byrow = TRUE)
			s2.p 			<- matrix(mcmc$samples_s2, nrow = nStore, byrow = TRUE)
			s3.p 			<- matrix(mcmc$samples_s3, nrow = nStore, byrow = TRUE)	
							
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
			accept.theta	<- as.vector(mcmc$samples_misc[(p1+p2+p3)+7])	
		
			write.table(gamma.p, file = paste(path,"gammaPch", chain, ".txt", sep = ""))
	
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
	
			
			ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, theta.p = theta.p, J1.p = J1.p, J2.p = J2.p, J3.p = J3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, hz.type = hz.type, model = model)

        } #if(model == "Markov")
			
		# model = "semi-Markov"			
		
		if(model == "semi-Markov"){
			
			###
			n	<- dim(survData)[1]
			p1	<- nCov[1]
			p2	<- nCov[2]
			p3	<- nCov[3]
	
			J1_max = mcmcParams[10]
			J2_max = mcmcParams[11]
			J3_max = mcmcParams[12]		
	
			nStore <- numReps/thin * (1 - burninPerc)
	
			nTime_lambda1 = mcmcParams[16]
			nTime_lambda2 = mcmcParams[17]
			nTime_lambda3 = mcmcParams[18]		
	
			time_lambda3	= mcmcParams[(length(mcmcParams)-1-nTime_lambda3+1):(length(mcmcParams)-1)]	
			time_lambda2	= mcmcParams[(length(mcmcParams)-1-nTime_lambda3-nTime_lambda2+1):(length(mcmcParams)-1-nTime_lambda3)]
			time_lambda1	= mcmcParams[(length(mcmcParams)-1-nTime_lambda3-nTime_lambda2-nTime_lambda1+1):(length(mcmcParams)-1-nTime_lambda3-nTime_lambda2)]			
	

			mcmc <- .C("BpeScrSMmcmc",
						survData            = as.double(as.matrix(survData)),
						n                   = as.integer(n),
						p1					= as.integer(p1),
						p2					= as.integer(p2),
						p3					= as.integer(p3),
						hyperParams         = as.double(hyperParams),
						startValues         = as.double(startV),
						mcmcParams          = as.double(mcmcParams),
						numReps             = as.integer(numReps),
						thin                = as.integer(thin),
						burninPerc          = as.double(burninPerc),
						nGam_save			= as.integer(nGam_save),						
						samples_beta1 		= as.double(rep(0, nStore*p1)),
						samples_beta2 		= as.double(rep(0, nStore*p2)),
						samples_beta3 		= as.double(rep(0, nStore*p3)),
						samples_mu_lam1     = as.double(rep(0, nStore*1)),
						samples_mu_lam2     = as.double(rep(0, nStore*1)),
						samples_mu_lam3     = as.double(rep(0, nStore*1)),
						samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
						samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
						samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
						samples_J1          = as.double(rep(0, nStore*1)),
						samples_J2          = as.double(rep(0, nStore*1)),
						samples_J3          = as.double(rep(0, nStore*1)),
						samples_s1          = as.double(rep(0, nStore*(J1_max + 1))),
						samples_s2          = as.double(rep(0, nStore*(J2_max + 1))),
						samples_s3          = as.double(rep(0, nStore*(J3_max + 1))),
						samples_theta 		= as.double(rep(0, nStore*1)),
						samples_gamma 		= as.double(rep(0, nStore*nGam_save)),
						samples_misc        = as.double(rep(0, p1 + p2 + p3 + 7)),
						lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
						lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
						lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)))

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
			gamma.p 		<- matrix(mcmc$samples_gamma, nrow = nStore, byrow = TRUE)	
			theta.p 		<- matrix(mcmc$samples_theta, nrow = nStore, byrow = TRUE)						
			mu_lam1.p 		<- matrix(mcmc$samples_mu_lam1, nrow = nStore, byrow = TRUE)
			mu_lam2.p 		<- matrix(mcmc$samples_mu_lam2, nrow = nStore, byrow = TRUE)
			mu_lam3.p 		<- matrix(mcmc$samples_mu_lam3, nrow = nStore, byrow = TRUE)		
			sigSq_lam1.p 	<- matrix(mcmc$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)	
			sigSq_lam2.p 	<- matrix(mcmc$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)	
			sigSq_lam3.p 	<- matrix(mcmc$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)			
			J1.p 			<- matrix(mcmc$samples_J1, nrow = nStore, byrow = TRUE)
			J2.p 			<- matrix(mcmc$samples_J2, nrow = nStore, byrow = TRUE)
			J3.p 			<- matrix(mcmc$samples_J3, nrow = nStore, byrow = TRUE)		
			s1.p 			<- matrix(mcmc$samples_s1, nrow = nStore, byrow = TRUE)
			s2.p 			<- matrix(mcmc$samples_s2, nrow = nStore, byrow = TRUE)
			s3.p 			<- matrix(mcmc$samples_s3, nrow = nStore, byrow = TRUE)
								
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
			accept.theta	<- as.vector(mcmc$samples_misc[(p1+p2+p3)+7])	
	
	
			write.table(gamma.p, file = paste(path,"gammaPch", chain, ".txt", sep = ""))
	
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
	
			
			ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, theta.p = theta.p, J1.p = J1.p, J2.p = J2.p, J3.p = J3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, hz.type = hz.type, model = model)

        }	# if(model == "semi-Markov")

    } # if(hz.type == "PEM")


	# hz.type = "Weibull"
	
	if(hz.type == "Weibull"){
        
		# model = "Markov"
        
		if(model == "Markov"){
            n	<- dim(survData)[1]
            p1	<- nCov[1]
            p2	<- nCov[2]
            p3	<- nCov[3]
            
            ###
            
            nStore <- numReps/thin * (1 - burninPerc)
            
            mcmc <- .C("BweibScrmcmc",
            survData 		= as.double(as.matrix(survData)),
            n				= as.integer(n),
            p1				= as.integer(p1),
            p2				= as.integer(p2),
            p3				= as.integer(p3),
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
            samples_theta 	= as.double(rep(0, nStore*1)),
            samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
            samples_misc	= as.double(rep(0, (p1+p2+p3+4))))
            
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
            
            write.table(gamma.p, file = paste(path,"gammaPch", chain, ".txt", sep = ""))
            
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
            
            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, theta.p = theta.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, hz.type = hz.type, model = model)
        } # if(model == "Markov")

		# model = "semi-Markov"        
        
        if(model == "semi-Markov"){
            n	<- dim(survData)[1]
            p1	<- nCov[1]
            p2	<- nCov[2]
            p3	<- nCov[3]
            
            ###
            
            nStore <- numReps/thin * (1 - burninPerc)
            
            mcmc <- .C("BweibScrSMmcmc",
            survData 		= as.double(as.matrix(survData)),
            n				= as.integer(n),
            p1				= as.integer(p1),
            p2				= as.integer(p2),
            p3				= as.integer(p3),
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
            samples_theta 	= as.double(rep(0, nStore*1)),
            samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
            samples_misc	= as.double(rep(0, (p1+p2+p3+4))))
            
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
            
            write.table(gamma.p, file = paste(path,"gammaPch", chain, ".txt", sep = ""))
            
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
            
            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, theta.p = theta.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, hz.type = hz.type, model = model)
        } # if(model == "semi-Markov")
	
		
    } # if(hz.type == "Weibull")
		
		chain = chain + 1	
	}	# while(chain <= nChain)
    
		ret[["setup"]]	<- list(nCov = nCov, hyperParams = hyperParams, startValues = startValues, mcmcParams = mcmcParams, nGam_save = nGam_save, numReps = numReps, thin = thin, path = path, burninPerc = burninPerc, hz.type = hz.type, model = model, nChain = nChain)
			   
		class(ret) <- "BayesID"
		return(ret)
		
}	
else{
	print("The 'startValues' should be the list of length equal to 'nChain'.")
}

}





























