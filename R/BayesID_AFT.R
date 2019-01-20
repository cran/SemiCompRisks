BayesID_AFT <- function(Formula,
data,
model = "LN",
hyperParams,
startValues,
mcmcParams,
na.action = "na.fail",
subset=NULL,
path=NULL)
{
    mcmcList    <- mcmcParams
    
    if((mcmcList$run$numReps / mcmcList$run$thin * mcmcList$run$burninPerc) %% 1 == 0)
    {
        nChain <- length(startValues)
        
        hz.type     <- model[1]
        
        if(na.action != "na.fail" & na.action != "na.omit")
        {
            stop("na.action should be either na.fail or na.omit")
        }
        
        form2 <- as.Formula(paste(Formula[2], Formula[1], Formula[3], sep = ""))
        
        if(hz.type == "DPM")
        {
            for(i in 1:nChain)
            {
                nam1 <- paste("DPM.class1ch", i, sep = "")
                data[[nam1]] <- startValues[[i]]$DPM$DPM.class1
                nam2 <- paste("DPM.class2ch", i, sep = "")
                data[[nam2]] <- startValues[[i]]$DPM$DPM.class2
                nam3 <- paste("DPM.class3ch", i, sep = "")
                data[[nam3]] <- startValues[[i]]$DPM$DPM.class3
                
                nam4 <- paste("DPM.mu1ch", i, sep = "")
                data[[nam4]] <- startValues[[i]]$DPM$DPM.mu1
                nam5 <- paste("DPM.mu2ch", i, sep = "")
                data[[nam5]] <- startValues[[i]]$DPM$DPM.mu2
                nam6 <- paste("DPM.mu3ch", i, sep = "")
                data[[nam6]] <- startValues[[i]]$DPM$DPM.mu3
                
                nam7 <- paste("DPM.zeta1ch", i, sep = "")
                data[[nam7]] <- startValues[[i]]$DPM$DPM.zeta1
                nam8 <- paste("DPM.zeta2ch", i, sep = "")
                data[[nam8]] <- startValues[[i]]$DPM$DPM.zeta2
                nam9 <- paste("DPM.zeta3ch", i, sep = "")
                data[[nam9]] <- startValues[[i]]$DPM$DPM.zeta3
                
                form2 <- as.Formula(paste(form2[2], form2[1], form2[3], "| ", nam1, "| ", nam2, "| ", nam3, "| ", nam4, "| ", nam5, "| ", nam6, "| ", nam7, "| ", nam8, "| ", nam9, sep = ""))
            }
        }
        
        for(i in 1:nChain)
        {
            nam1 <- paste("y1ch", i, sep = "")
            data[[nam1]] <- startValues[[i]]$common$y1
            nam2 <- paste("y2ch", i, sep = "")
            data[[nam2]] <- startValues[[i]]$common$y2
            nam3<- paste("gamch", i, sep = "")
            data[[nam3]] <- startValues[[i]]$common$gamma
            
            form2 <- as.Formula(paste(form2[2], form2[1], form2[3], "| ", nam1, "| ", nam2, "| ", nam3, sep = ""))
        }
    
        data <- model.frame(form2, data=data, na.action = na.action, subset = subset)
        
        if(hz.type == "DPM")
        {
            for(i in 1:nChain)
            {
                nam1 <- paste("DPM.class1ch", i, sep = "")
                startValues[[i]]$DPM$DPM.class1 <- data[[nam1]]
                nam2 <- paste("DPM.class2ch", i, sep = "")
                startValues[[i]]$DPM$DPM.class2 <- data[[nam2]]
                nam3 <- paste("DPM.class3ch", i, sep = "")
                startValues[[i]]$DPM$DPM.class3 <- data[[nam3]]
                
                nam4 <- paste("DPM.mu1ch", i, sep = "")
                startValues[[i]]$DPM$DPM.mu1 <- data[[nam4]]
                nam5 <- paste("DPM.mu2ch", i, sep = "")
                startValues[[i]]$DPM$DPM.mu2 <- data[[nam5]]
                nam6 <- paste("DPM.mu3ch", i, sep = "")
                startValues[[i]]$DPM$DPM.mu3 <- data[[nam6]]
                
                nam7 <- paste("DPM.zeta1ch", i, sep = "")
                startValues[[i]]$DPM$DPM.zeta1 <- data[[nam7]]
                nam8 <- paste("DPM.zeta2ch", i, sep = "")
                startValues[[i]]$DPM$DPM.zeta2 <- data[[nam8]]
                nam9 <- paste("DPM.zeta3ch", i, sep = "")
                startValues[[i]]$DPM$DPM.zeta3 <- data[[nam9]]
            }
        }

        for(i in 1:nChain)
        {
            nam1 <- paste("y1ch", i, sep = "")
            startValues[[i]]$common$y1 <- data[[nam1]]
            nam2 <- paste("y2ch", i, sep = "")
            startValues[[i]]$common$y2 <- data[[nam2]]
            nam3<- paste("gamch", i, sep = "")
            startValues[[i]]$common$gamma <- data[[nam3]]
        }
        
        LT <- model.part(Formula, data=data, lhs=1)
        y1.mat <- model.part(Formula, data=data, lhs=2)
        y2.mat <- model.part(Formula, data=data, lhs=3)
        
        Y <- cbind(y1.mat, y2.mat, LT)
        
        Xmat1 <- model.frame(formula(Formula, lhs=0, rhs=1), data=data)
        Xmat2 <- model.frame(formula(Formula, lhs=0, rhs=2), data=data)
        Xmat3 <- model.frame(formula(Formula, lhs=0, rhs=3), data=data)
        
        p1 <- ncol(Xmat1)
        p2 <- ncol(Xmat2)
        p3 <- ncol(Xmat3)
        
        nCov <- c(p1, p2, p3)
        
        ##
        
        if(p1 == 0 & p2 == 0 & p3 == 0){
            survData <- Y
        }
        if(p1 == 0 & p2 > 0 & p3 > 0){
            survData <- cbind(Y, Xmat2, Xmat3)
        }
        if(p1 > 0 & p2 == 0 & p3 > 0){
            survData <- cbind(Y, Xmat1, Xmat3)
        }
        if(p1 > 0 & p2 > 0 & p3 == 0){
            survData <- cbind(Y, Xmat1, Xmat2)
        }
        if(p1 > 0 & p2 == 0 & p3 == 0){
            survData <- cbind(Y, Xmat1)
        }
        if(p1 == 0 & p2 > 0 & p3 == 0){
            survData <- cbind(Y, Xmat2)
        }
        if(p1 == 0 & p2 == 0 & p3 > 0){
            survData <- cbind(Y, Xmat3)
        }
        if(p1 > 0 & p2 > 0 & p3 > 0){
            survData <- cbind(Y, Xmat1, Xmat2, Xmat3)
        }
        
        n	<- dim(survData)[1]
        
        Y[,1] <- log(Y[,1])
        Y[,2] <- log(Y[,2])
        Y[,3] <- log(Y[,3])
        Y[,4] <- log(Y[,4])
        Y[,5] <- log(Y[,5])
        
        y1LInf <- rep(0, n)
        for(i in 1:n) if(Y[i,1] == -Inf)
        {
            Y[i,1] <- -9.9e10
            y1LInf[i] <- 1
        }
        y2LInf <- rep(0, n)
        for(i in 1:n) if(Y[i,3] == -Inf)
        {
            Y[i,3] <- -9.9e10
            y2LInf[i] <- 1
        }
        
        y1UInf <- rep(0, n)
        for(i in 1:n) if(Y[i,2] == Inf)
        {
            Y[i,2] <- 9.9e10
            y1UInf[i] <- 1
        }
        
        y2UInf <- rep(0, n)
        for(i in 1:n) if(Y[i,4] == Inf)
        {
            Y[i,4] <- 9.9e10
            y2UInf[i] <- 1
        }
        
        
        c0Inf <- rep(0, n)
        for(i in 1:n) if(Y[i,5] == -Inf)
        {
            Y[i,5] <- -9.9e10
            c0Inf[i] <- 1
        }
        
        if(!is.null(path)){
            dir.create(paste(path), recursive = TRUE, showWarnings = FALSE)
        }
        
        ### setting hyperparameters
        
        if(hz.type == "DPM")
        {
            hyperP  <- as.vector(c(hyperParams$theta, hyperParams$DPM$DPM.ab1, hyperParams$DPM$DPM.ab2, hyperParams$DPM$DPM.ab3, hyperParams$DPM$Tau.ab1, hyperParams$DPM$Tau.ab2, hyperParams$DPM$Tau.ab3, hyperParams$DPM$DPM.mu1, hyperParams$DPM$DPM.mu2, hyperParams$DPM$DPM.mu3, hyperParams$DPM$DPM.sigSq1, hyperParams$DPM$DPM.sigSq2, hyperParams$DPM$DPM.sigSq3))
        }
        
        if(hz.type == "LN")
        {
            hyperP  <- as.vector(c(hyperParams$theta, hyperParams$LN$LN.ab1,hyperParams$LN$LN.ab2, hyperParams$LN$LN.ab3))
        }
        
        ### mcmc setting
        
        mcmcP   <- as.vector(c(mcmcParams$tuning$betag.prop.var, mcmcParams$tuning$mug.prop.var, mcmcParams$tuning$zetag.prop.var, mcmcParams$tuning$gamma.prop.var))
        
        chain = 1
        ret <- list()
        
        while(chain <= nChain){
            
            cat("chain: ", chain, "\n")
            nam = paste("chain", chain, sep="")
            
            temp <- startValues[[chain]]
            
            ### setting starting values
            
            if(hz.type == "DPM")
            {
                startV <- as.vector(c(y1=temp$common$y1, y2=temp$common$y2, beta1=temp$common$beta1, beta2=temp$common$beta2, beta3=temp$common$beta3, r1=temp$DPM$DPM.class1, r2=temp$DPM$DPM.class2, r3=temp$DPM$DPM.class3, tau=temp$DPM$DPM.tau, gamma=temp$common$gamma, theta=temp$common$theta, mu1=temp$DPM$DPM.mu1, mu2=temp$DPM$DPM.mu2, mu3=temp$DPM$DPM.mu3, zeta1=temp$DPM$DPM.zeta1, zeta2=temp$DPM$DPM.zeta2, zeta3=temp$DPM$DPM.zeta3))
            }
            
            if(hz.type == "LN")
            {
                startV <- as.vector(c(y1=temp$common$y1, y2=temp$common$y2, beta1=temp$common$beta1, beta2=temp$common$beta2, beta3=temp$common$beta3, mu=temp$LN$LN.mu, sigSq=temp$LN$LN.sigSq, gamma=temp$common$gamma, theta=temp$common$theta))
            }
            
            
            # hz.type = "LN"
            
            if(hz.type == "LN")
            {
                nGam_save   <- mcmcParams$storage$nGam_save
                nY1_save = mcmcParams$storage$nY1_save
                nY2_save = mcmcParams$storage$nY2_save
                nY1_NA_save = mcmcParams$storage$nY1.NA_save
                
                numReps     <- mcmcParams$run$numReps
                thin        <- mcmcParams$run$thin
                burninPerc  <- mcmcParams$run$burninPerc
                
                nStore <- round(numReps/thin*(1-burninPerc))
                
                mcmcRet     <- .C("BAFTscrmcmc",
                Ymat            = as.double(as.matrix(Y)),
                y1LInf			= as.double(y1LInf),
                y1UInf			= as.double(y1UInf),
                y2LInf			= as.double(y2LInf),
                y2UInf			= as.double(y2UInf),
                c0Inf			= as.double(c0Inf),
                X1mat          	= as.double(as.matrix(Xmat1)),
                X2mat          	= as.double(as.matrix(Xmat2)),
                X3mat          	= as.double(as.matrix(Xmat3)),
                n				= as.integer(n),
                p1				= as.integer(p1),
                p2				= as.integer(p2),
                p3				= as.integer(p3),
                hyperP          = as.double(hyperP),
                mcmcP           = as.double(mcmcP),
                startValues 		= as.double(startV),
                numReps			= as.integer(numReps),
                thin				= as.integer(thin),
                burninPerc      = as.double(burninPerc),
                nGam_save		= as.integer(nGam_save),
                samples_y1       = as.double(rep(0, nStore*n)),
                samples_y2       = as.double(rep(0, nStore*n)),
                samples_y1_NA    = as.double(rep(0, nStore*n)),
                samples_beta1    = as.double(rep(0, nStore*p1)),
                samples_beta2    = as.double(rep(0, nStore*p2)),
                samples_beta3    = as.double(rep(0, nStore*p3)),
                samples_beta01   = as.double(rep(0, nStore*1)),
                samples_beta02   = as.double(rep(0, nStore*1)),
                samples_beta03   = as.double(rep(0, nStore*1)),
                samples_sigSq1   = as.double(rep(0, nStore*1)),
                samples_sigSq2   = as.double(rep(0, nStore*1)),
                samples_sigSq3   = as.double(rep(0, nStore*1)),
                samples_theta 	= as.double(rep(0, nStore*1)),
                samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                samples_misc    = as.double(rep(0, p1+p2+p3+3+3+n)),
                samples_logLH   = as.double(rep(0, nStore*1)),
                LH_i_mean       = as.double(rep(0, n*1)),
                invLH_i_mean     = as.double(rep(0, n*1)))
                
                y1.p <- matrix(as.vector(mcmcRet$samples_y1), nrow=nStore, byrow=T)
                y2.p <- matrix(as.vector(mcmcRet$samples_y2), nrow=nStore, byrow=T)
                y1.NA <- matrix(as.vector(mcmcRet$samples_y1_NA), nrow=nStore, byrow=T)
                
                if(p1 >0)
                {
                    beta1.p <- matrix(as.vector(mcmcRet$samples_beta1), nrow=nStore, byrow=T)
                }else
                {
                    beta1.p <- NULL
                }
                if(p2 >0)
                {
                    beta2.p <- matrix(as.vector(mcmcRet$samples_beta2), nrow=nStore, byrow=T)
                }else
                {
                    beta2.p <- NULL
                }
                if(p3 >0)
                {
                    beta3.p <- matrix(as.vector(mcmcRet$samples_beta3), nrow=nStore, byrow=T)
                }else
                {
                    beta3.p <- NULL
                }
                
                mu1.p <- matrix(as.vector(mcmcRet$samples_beta01), nrow=nStore, byrow=T)
                mu2.p <- matrix(as.vector(mcmcRet$samples_beta02), nrow=nStore, byrow=T)
                mu3.p <- matrix(as.vector(mcmcRet$samples_beta03), nrow=nStore, byrow=T)
                sigSq1.p <- matrix(as.vector(mcmcRet$samples_sigSq1), nrow=nStore, byrow=T)
                sigSq2.p <- matrix(as.vector(mcmcRet$samples_sigSq2), nrow=nStore, byrow=T)
                sigSq3.p <- matrix(as.vector(mcmcRet$samples_sigSq3), nrow=nStore, byrow=T)
                gamma.p <- matrix(as.vector(mcmcRet$samples_gamma), nrow = nStore, byrow = T)
                theta.p <- matrix(as.vector(mcmcRet$samples_theta), nrow = nStore, byrow = T)
                
                if(p1 >0)
                {
                    accept.beta1	<- as.vector(mcmcRet$samples_misc[1:p1])
                }else
                {
                    accept.beta1 <- NULL
                }
                if(p2 >0)
                {
                    accept.beta2	 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])
                }else
                {
                    accept.beta2 <- NULL
                }
                if(p3 >0)
                {
                    accept.beta3		<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])
                }else
                {
                    accept.beta3 <- NULL
                }
                
                accept.mu1	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+1])
                accept.mu2	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+2])
                accept.mu3	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+3])
                accept.sigSq1	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+3+1])
                accept.sigSq2	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+3+2])
                accept.sigSq3	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+3+3])
                accept.gamma	 <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+3+3+1):(p1+p2+p3+3+3+n)])
                
                logLH <- matrix(as.vector(mcmcRet$samples_logLH), nrow = nStore, byrow = T)
                
                if(p1 > 0){
                    covNames1 = colnames(Xmat1)
                }
                if(p1 == 0){
                    covNames1 = NULL
                }
                if(p2 > 0){
                    covNames2 = colnames(Xmat2)
                }
                if(p2 == 0){
                    covNames2 = NULL
                }
                if(p3 > 0){
                    covNames3 = colnames(Xmat3)
                }
                if(p3 == 0){
                    covNames3 = NULL
                }
                
                ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, mu1.p=mu1.p, mu2.p=mu2.p, mu3.p=mu3.p, sigSq1.p = sigSq1.p, sigSq2.p = sigSq2.p, sigSq3.p = sigSq3.p, gamma.p=gamma.p, theta.p=theta.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.mu1 = accept.mu1, accept.mu2 = accept.mu2, accept.mu3 = accept.mu3, accept.sigSq1 = accept.sigSq1, accept.sigSq2 = accept.sigSq2, accept.sigSq3 = accept.sigSq3, accept.gamma=accept.gamma, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, logLH_mean=mean(logLH), LH_i_mean=mcmcRet$LH_i_mean, invLH_i_mean=mcmcRet$invLH_i_mean)
                
            }# if(hz.type == "LN")
            
            
            # hz.type = "DPM"
            
            if(hz.type == "DPM"){
                
                nGam_save   <- mcmcParams$storage$nGam_save
                nY1_save = mcmcParams$storage$nY1_save
                nY2_save = mcmcParams$storage$nY2_save
                nY1_NA_save = mcmcParams$storage$nY1.NA_save
                
                numReps     <- mcmcParams$run$numReps
                thin        <- mcmcParams$run$thin
                burninPerc  <- mcmcParams$run$burninPerc
                
                nStore <- round(numReps/thin*(1-burninPerc))
                
                mcmcRet     <- .C("BAFT_DPscrmcmc",
                Ymat            = as.double(as.matrix(Y)),
                y1LInf			= as.double(y1LInf),
                y1UInf			= as.double(y1UInf),
                y2LInf			= as.double(y2LInf),
                y2UInf			= as.double(y2UInf),
                c0Inf			= as.double(c0Inf),
                X1mat          	= as.double(as.matrix(Xmat1)),
                X2mat          	= as.double(as.matrix(Xmat2)),
                X3mat          	= as.double(as.matrix(Xmat3)),
                n				= as.integer(n),
                p1				= as.integer(p1),
                p2				= as.integer(p2),
                p3				= as.integer(p3),
                hyperP          = as.double(hyperP),
                mcmcP           = as.double(mcmcP),
                startValues 	= as.double(startV),
                numReps			= as.integer(numReps),
                thin			= as.integer(thin),
                burninPerc      = as.double(burninPerc),
                nGam_save		= as.integer(nGam_save),
                nY1_save		= as.integer(nY1_save),
                nY2_save		= as.integer(nY2_save),
                nY1_NA_save		= as.integer(nY1_NA_save),
                samples_y1       = as.double(rep(0, nStore*nY1_save)),
                samples_y2       = as.double(rep(0, nStore*nY2_save)),
                samples_y1_NA    = as.double(rep(0, nStore*nY1_NA_save)),
                samples_beta1    = as.double(rep(0, nStore*p1)),
                samples_beta2    = as.double(rep(0, nStore*p2)),
                samples_beta3    = as.double(rep(0, nStore*p3)),
                samples_r1		= as.double(rep(0, nStore*n)),
                samples_r2		= as.double(rep(0, nStore*n)),
                samples_r3		= as.double(rep(0, nStore*n)),
                samples_beta01	= as.double(rep(0, nStore*n)),
                samples_beta02	= as.double(rep(0, nStore*n)),
                samples_beta03	= as.double(rep(0, nStore*n)),
                samples_sigSq1   = as.double(rep(0, nStore*n)),
                samples_sigSq2   = as.double(rep(0, nStore*n)),
                samples_sigSq3   = as.double(rep(0, nStore*n)),
                samples_tau1     = as.double(rep(0, nStore*1)),
                samples_tau2     = as.double(rep(0, nStore*1)),
                samples_tau3     = as.double(rep(0, nStore*1)),
                samples_theta 	= as.double(rep(0, nStore*1)),
                samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                samples_misc    = as.double(rep(0, p1+p2+p3+n+3+3)),
                samples_logLH   = as.double(rep(0, nStore*1)),
                LH_i_mean       = as.double(rep(0, n*1)),
                invLH_i_mean     = as.double(rep(0, n*1)))
                
                
                y1.p <- matrix(as.vector(mcmcRet$samples_y1), nrow=nStore, byrow=T)
                y2.p <- matrix(as.vector(mcmcRet$samples_y2), nrow=nStore, byrow=T)
                y1.NA <- matrix(as.vector(mcmcRet$samples_y1_NA), nrow=nStore, byrow=T)
                
                if(p1 >0)
                {
                    beta1.p <- matrix(as.vector(mcmcRet$samples_beta1), nrow=nStore, byrow=T)
                }else
                {
                    beta1.p <- NULL
                }
                if(p2 >0)
                {
                    beta2.p <- matrix(as.vector(mcmcRet$samples_beta2), nrow=nStore, byrow=T)
                }else
                {
                    beta2.p <- NULL
                }
                if(p3 >0)
                {
                    beta3.p <- matrix(as.vector(mcmcRet$samples_beta3), nrow=nStore, byrow=T)
                }else
                {
                    beta3.p <- NULL
                }
                
                mu1.p    <- matrix(mcmcRet$samples_beta01, nrow = nStore, byrow = TRUE)
                mu2.p    <- matrix(mcmcRet$samples_beta02, nrow = nStore, byrow = TRUE)
                mu3.p    <- matrix(mcmcRet$samples_beta03, nrow = nStore, byrow = TRUE)
                sigSq1.p    <- matrix(as.vector(mcmcRet$samples_sigSq1), nrow=nStore, byrow=TRUE)
                sigSq2.p    <- matrix(as.vector(mcmcRet$samples_sigSq2), nrow=nStore, byrow=TRUE)
                sigSq3.p    <- matrix(as.vector(mcmcRet$samples_sigSq3), nrow=nStore, byrow=TRUE)
                
                r1.p     <- matrix(mcmcRet$samples_r1, nrow = nStore, byrow = TRUE)
                r2.p     <- matrix(mcmcRet$samples_r2, nrow = nStore, byrow = TRUE)
                r3.p     <- matrix(mcmcRet$samples_r3, nrow = nStore, byrow = TRUE)
                tau1.p   <- matrix(mcmcRet$samples_tau1, nrow = nStore, byrow = TRUE)
                tau2.p   <- matrix(mcmcRet$samples_tau2, nrow = nStore, byrow = TRUE)
                tau3.p   <- matrix(mcmcRet$samples_tau3, nrow = nStore, byrow = TRUE)
                
                gamma.p <- matrix(as.vector(mcmcRet$samples_gamma), nrow = nStore, byrow = T)
                y1.p <- matrix(as.vector(mcmcRet$samples_y1), nrow = nStore, byrow = T)
                y2.p <- matrix(as.vector(mcmcRet$samples_y2), nrow = nStore, byrow = T)
                y1_NA.p <- matrix(as.vector(mcmcRet$samples_y1_NA), nrow = nStore, byrow = T)
                
                
                theta.p <- matrix(as.vector(mcmcRet$samples_theta), nrow = nStore, byrow = T)
                
                if(p1 >0)
                {
                    accept.beta1	<- as.vector(mcmcRet$samples_misc[1:p1])
                }else
                {
                    accept.beta1 <- NULL
                }
                if(p2 >0)
                {
                    accept.beta2	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])
                }else
                {
                    accept.beta2 <- NULL
                }
                if(p3 >0)
                {
                    accept.beta3	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])
                }else
                {
                    accept.beta3 <- NULL
                }
                
                accept.gamma	 <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+1):(p1+p2+p3+n)])
                
                accept.mu1	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+n+1])
                accept.mu2	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+n+2])
                accept.mu3	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+n+3])
                accept.sigSq1	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+n+4])
                accept.sigSq2	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+n+5])
                accept.sigSq3	 <- as.vector(mcmcRet$samples_misc[p1+p2+p3+n+6])
                
                if(nGam_save > 0 & !is.null(path)){
                    save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                }
                
                if(nY1_save > 0 & !is.null(path))
                {
                    save(y1.p, file = paste(path, "/y1ch", chain, ".RData", sep = ""))
                }
                if(nY2_save > 0 & !is.null(path))
                {
                    save(y2.p, file = paste(path, "/y2ch", chain, ".RData", sep = ""))
                }
                if(nY1_NA_save > 0 & !is.null(path))
                {
                    save(y1_NA.p, file = paste(path, "/y1NAch", chain, ".RData", sep = ""))
                }
                
                
                logLH <- matrix(as.vector(mcmcRet$samples_logLH), nrow = nStore, byrow = T)
                
                if(p1 > 0){
                    covNames1 = colnames(Xmat1)
                }
                if(p1 == 0){
                    covNames1 = NULL
                }
                if(p2 > 0){
                    covNames2 = colnames(Xmat2)
                }
                if(p2 == 0){
                    covNames2 = NULL
                }
                if(p3 > 0){
                    covNames3 = colnames(Xmat3)
                }
                if(p3 == 0){
                    covNames3 = NULL
                }
                
                ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, mu1.p=mu1.p, mu2.p=mu2.p, mu3.p=mu3.p, sigSq1.p=sigSq1.p, sigSq2.p=sigSq2.p, sigSq3.p=sigSq3.p, r1.p = r1.p, r2.p = r2.p, r3.p = r3.p, tau1.p = tau1.p, tau2.p = tau2.p, tau3.p = tau3.p, theta.p=theta.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.gamma=accept.gamma, accept.mu1=accept.mu1, accept.mu2=accept.mu2, accept.mu3=accept.mu3, accept.sigSq1=accept.sigSq1, accept.sigSq2=accept.sigSq2, accept.sigSq3=accept.sigSq3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, logLH_mean=mean(logLH), LH_i_mean=mcmcRet$LH_i_mean, invLH_i_mean=mcmcRet$invLH_i_mean)
                
            } # if(hz.type == "DPM")
            
            chain = chain + 1
        }   # while(chain <= nChain)
        
        ret[["setup"]]	<- list(nCov = nCov, hyperParams = hyperParams, startValues = startValues, mcmcParams = mcmcParams, nGam_save = nGam_save, numReps = numReps, thin = thin, path=path, burninPerc = burninPerc, model = hz.type, nChain = nChain)
        
        class(ret) <- "Bayes_AFT"
        
        if(hz.type == "LN")
        {
            ret$class <- c("Bayes_AFT", "ID", "Ind", "LN")
        }
        if(hz.type == "DPM")
        {
            ret$class <- c("Bayes_AFT", "ID", "Ind", "DPM")
        }
        
        return(ret)
        
    }else{
        warning(" (numReps * burninPerc) must be divisible by (thin)")
    }
    
    
    
    
}
























