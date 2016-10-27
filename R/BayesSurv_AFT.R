

BayesSurv_AFT <- function(Y,
lin.pred,
data,
model = "LN",
hyperParams,
startValues,
mcmcParams,
path = NULL)
{
    mcmcList    <- mcmcParams
    
    if((mcmcList$run$numReps / mcmcList$run$thin * mcmcList$run$burninPerc) %% 1 == 0)
    {
        nChain <- length(startValues)
        
        hz.type 	<- model[1]
        
        Xmat <- model.frame(lin.pred, data=data)
        
        p <- ncol(Xmat)
        
        ##
        
        if(p == 0){
            survData <- Y
        }
        
        if(p > 0){
            survData <- cbind(Y, Xmat)
        }
        
        n	<- dim(survData)[1]
        
        Y[,1] <- log(Y[,1])
        Y[,2] <- log(Y[,2])
        Y[,3] <- log(Y[,3])
        
        yLInf <- rep(0, n)
        for(i in 1:n) if(Y[i,1] == -Inf)
        {
            Y[i,1] <- -9.9e10
            yLInf[i] <- 1
        }
        
        yUInf <- rep(0, n)
        for(i in 1:n) if(Y[i,2] == Inf)
        {
            Y[i,2] <- 9.9e10
            yUInf[i] <- 1
        }
        
        c0Inf <- rep(0, n)
        for(i in 1:n) if(Y[i,3] == -Inf)
        {
            Y[i,3] <- -9.9e10
            c0Inf[i] <- 1
        }
        
        if(!is.null(path)){
            dir.create(paste(path), recursive = TRUE, showWarnings = FALSE)
        }
        
        
        ### setting hyperparameters
        
        if(hz.type == "DPM")
        {
            hyperP  <- as.vector(c(hyperParams$DPM$DPM.ab, hyperParams$DPM$Tau.ab, hyperParams$DPM$DPM.mu, hyperParams$DPM$DPM.sigSq))
        }
        
        if(hz.type == "LN")
        {
            hyperP  <- as.vector(c(hyperParams$LN$LN.ab))
        }
        
        ### mcmc setting
        
        mcmcP   <- as.vector(c(mcmcParams$tuning$beta.prop.var, mcmcParams$tuning$mu.prop.var, mcmcParams$tuning$zeta.prop.var))
        
        
        chain = 1
        ret <- list()
        
        while(chain <= nChain){
            
            cat("chain: ", chain, "\n")
            nam = paste("chain", chain, sep="")
            
            temp <- startValues[[chain]]
            
            ### setting starting values
            
            if(hz.type == "DPM")
            {
                startV <- as.vector(c(y=temp$common$y, beta=temp$common$beta, r=temp$DPM$DPM.class, tau=temp$DPM$DPM.tau, mu=temp$DPM$DPM.mu, zeta=temp$DPM$DPM.zeta))
            }
            
            if(hz.type == "LN")
            {
                startV <- as.vector(c(y=temp$common$y, beta=temp$common$beta, mu=temp$LN$LN.mu, sigSq=temp$LN$LN.sigSq))
            }
            
            # hz.type = "LN"
            
            if(hz.type == "LN"){
                
                numReps     <- mcmcParams$run$numReps
                thin        <- mcmcParams$run$thin
                burninPerc  <- mcmcParams$run$burninPerc
                nStore <- round(numReps/thin*(1-burninPerc))
                
                mcmcRet     <- .C("BAFTunimcmc",
                Ymat            = as.double(as.matrix(Y)),
                yUInf			= as.double(yUInf),
                c0Inf			= as.double(c0Inf),
                Xmat           	= as.double(as.matrix(Xmat)),
                n				= as.integer(n),
                p				= as.integer(p),
                hyperP          = as.double(hyperP),
                mcmcP           = as.double(mcmcP),
                startValues 	= as.double(startV),
                numReps			= as.integer(numReps),
                thin			= as.integer(thin),
                burninPerc      = as.double(burninPerc),
                samples_y       = as.double(rep(0, nStore*n)),
                samples_beta    = as.double(rep(0, nStore*p)),
                samples_beta0   = as.double(rep(0, nStore*1)),
                samples_sigSq   = as.double(rep(0, nStore*1)),
                samples_misc    = as.double(rep(0, p+1+1)))
                
                y.p <- matrix(as.vector(mcmcRet$samples_y), nrow=nStore, byrow=T)
                if(p >0)
                {
                    beta.p <- matrix(as.vector(mcmcRet$samples_beta), nrow=nStore, byrow=T)
                }else
                {
                    beta.p <- NULL
                }
                
                mu.p <- matrix(as.vector(mcmcRet$samples_beta0), nrow=nStore, byrow=T)
                sigSq.p <- matrix(as.vector(mcmcRet$samples_sigSq), nrow=nStore, byrow=T)
                
                if(p >0)
                {
                    accept.beta	 <- as.vector(mcmcRet$samples_misc[1:p])
                }else
                {
                    accept.beta <- NULL
                }
                
                accept.mu	 <- as.vector(mcmcRet$samples_misc[p+1])
                accept.sigSq	 <- as.vector(mcmcRet$samples_misc[p+2])
                
                if(p > 0){
                    covNames = colnames(Xmat)
                }
                if(p == 0){
                    covNames = NULL
                }
                
                ret[[nam]] <- list(y.p = y.p, beta.p = beta.p, mu.p=mu.p, sigSq.p = sigSq.p, accept.beta = accept.beta, accept.mu = accept.mu, accept.sigSq = accept.sigSq, covNames = covNames, model = hz.type)
                
            }   # if(hz.type == "LN")
            
            
            
            # hz.type = "DPM"
            
            if(hz.type == "DPM"){
                
                numReps     <- mcmcParams$run$numReps
                thin        <- mcmcParams$run$thin
                burninPerc  <- mcmcParams$run$burninPerc
                nStore <- round(numReps/thin*(1-burninPerc))
                
                mcmcRet     <- .C("BAFT_DPunimcmc",
                Ymat            = as.double(as.matrix(Y)),
                yLInf			= as.double(yLInf),
                yUInf			= as.double(yUInf),
                c0Inf			= as.double(c0Inf),
                Xmat           	= as.double(as.matrix(Xmat)),
                n				= as.integer(n),
                p				= as.integer(p),
                hyperP          = as.double(hyperP),
                mcmcP           = as.double(mcmcP),
                startValues 	= as.double(startV),
                numReps			= as.integer(numReps),
                thin			= as.integer(thin),
                burninPerc      = as.double(burninPerc),
                samples_y       = as.double(rep(0, nStore*n)),
                samples_beta    = as.double(rep(0, nStore*p)),
                samples_r		= as.double(rep(0, nStore*n)),
                samples_mu      = as.double(rep(0, nStore*n)),
                samples_sigSq   = as.double(rep(0, nStore*n)),
                samples_tau     = as.double(rep(0, nStore*1)),
                samples_misc    = as.double(rep(0, p+2)))
                
                y.p <- matrix(as.vector(mcmcRet$samples_y), nrow=nStore, byrow=T)
                if(p >0)
                {
                    beta.p <- matrix(as.vector(mcmcRet$samples_beta), nrow=nStore, byrow=T)
                }else
                {
                    beta.p <- NULL
                }
                r.p     <- matrix(mcmcRet$samples_r, nrow = nStore, byrow = T)
                mu.p    <- matrix(mcmcRet$samples_mu, nrow = nStore, byrow = T)
                sigSq.p <- matrix(as.vector(mcmcRet$samples_sigSq), nrow=nStore, byrow=T)
                tau.p   <- matrix(mcmcRet$samples_tau, nrow = nStore, byrow = T)
                
                if(p >0)
                {
                    accept.beta	 <- as.vector(mcmcRet$samples_misc[1:p])
                }else
                {
                    accept.beta <- NULL
                }
                accept.mu	 <- as.vector(mcmcRet$samples_misc[p+1])
                accept.sigSq	 <- as.vector(mcmcRet$samples_misc[p+2])
                
                if(p > 0){
                    covNames = colnames(Xmat)
                }
                if(p == 0){
                    covNames = NULL
                }
                
                ret[[nam]] <- list(y.p = y.p, beta.p = beta.p, r.p = r.p, mu.p = mu.p, sigSq.p = sigSq.p, tau.p = tau.p, accept.beta = accept.beta, accept.mu = accept.mu, accept.sigSq=accept.sigSq, covNames = covNames, model = hz.type)
            }
            
            chain = chain + 1
        }
        
        ret[["setup"]]	<- list(hyperParams = hyperParams, startValues = startValues, mcmcParams = mcmcParams, numReps = numReps, thin = thin, path = path, burninPerc = burninPerc, model = hz.type, nChain = nChain)
        
        if(hz.type == "LN")
        {
            class(ret) <- c("Bayes_AFT", "Surv", "Ind", "LN")
        }
        if(hz.type == "DPM")
        {
            class(ret) <- c("Bayes_AFT", "Surv", "Ind", "DPM")
        }
        
        return(ret)
        
    }
    else{
        print(" (numReps * burninPerc) must be divisible by (thin)")
    }
    
}





























