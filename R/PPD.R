
PPD <- function(fit, x1, x2, x3, t1, t2)
{
    if(class(fit) != "Bayes_HReg" | fit$class[2] != "ID" | fit$class[4] != "PEM")
    {
        warning("Currently, PPD() can only be used for PEM illness-death models")
    }else if(max(fit$chain1$time_lambda1) < max(t1, t2))
    {
        warning("max(t1, t2) should be less than or equal to max(time_lambda1)")
    }else if(max(fit$chain1$time_lambda2) < max(t1, t2))
    {
        warning("max(t1, t2) should be less than or equal to max(time_lambda2)")
    }else
    {
        MAXt1t2 <- max(t1, t2)
        T1grid = seq(from=MAXt1t2/100, to=MAXt1t2, by=MAXt1t2/100)
        T2grid = seq(from=MAXt1t2/100, to=MAXt1t2, by=MAXt1t2/100)
        T1grid = sort(unique(c(t1, t2, T1grid)))
        T2grid = sort(unique(c(t2, t2, T2grid)))
        
        nChain = fit$setup$nChain
        
        t1ind <- which(T1grid == t1)
        t2ind <- which(T2grid == t2)
        
        time1 = fit$chain1$time_lambda1
        time2 = fit$chain1$time_lambda2
        time3 = fit$chain1$time_lambda3
        
        ntime1  = length(time1)
        ntime2  = length(time2)
        ntime3  = length(time3)
        
        nT1grid <- length(T1grid)
        nT2grid <- length(T2grid)
        
        diffT1 <- diff(c(0, T1grid))
        diffT2 <- diff(c(0, T2grid))
        
        DIFFgrid <- matrix(NA, nT1grid, nT2grid)
        for(i in 1:nT1grid)
        {
            for(j in 1:nT2grid)
            {
                DIFFgrid[i,j] <- T2grid[j] - T1grid[i]
            }
        }
        
        lam1 <- fit$chain1$lambda1.fin
        lam2 <- fit$chain1$lambda2.fin
        lam3 <- fit$chain1$lambda3.fin
        
        beta1 <- fit$chain1$beta1.p
        beta2 <- fit$chain1$beta2.p
        beta3 <- fit$chain1$beta3.p
        
        theta <- fit$chain1$theta.p
        
        if(nChain > 1)
        {
            for(i in 2:nChain){
                nam <- paste("chain", i, sep="")
                lam1 <- rbind(lam1, fit[[nam]]$lambda1.fin)
                lam2 <- rbind(lam2, fit[[nam]]$lambda2.fin)
                lam3 <- rbind(lam3, fit[[nam]]$lambda3.fin)
                
                beta1 <- rbind(beta1, fit[[nam]]$beta1.p)
                beta2 <- rbind(beta2, fit[[nam]]$beta2.p)
                beta3 <- rbind(beta3, fit[[nam]]$beta3.p)
                
                theta <- rbind(theta, fit[[nam]]$theta.p)
            }
        }
        nSample <- dim(theta)[1]
        
        lambda1 <- matrix(NA, nSample, nT1grid)
        lambda2 <- matrix(NA, nSample, nT2grid)
        
        ind1 <- rep(NA, nT1grid)
        for(i in 1:nT1grid)
        {
            ind1[i] <- min(which(T1grid[i] - time1 <= 0))
        }
        ind2 <- rep(NA, nT2grid)
        for(i in 1:nT2grid)
        {
            ind2[i] <- min(which(T2grid[i] - time2 <= 0))
        }
        
        for(i in 1:nSample)
        {
            lambda1[i,] <- lam1[i, ind1]
            lambda2[i,] <- lam2[i, ind2]
        }
        
        ind3 <- matrix(NA, nT1grid, nT2grid)
        for(i in 1:nT1grid)
        {
            for(j in 1:nT2grid)
            {
                ind3[i,j] <- min(which(T2grid[j] - T1grid[i] - time3 <= 0))
            }
        }
        
        f_u <- array(0, c(nSample, nT1grid, nT2grid))
        f_l <- matrix(0, nSample, nT2grid)
        F_u.inc <- array(0, c(nSample, nT1grid, nT2grid))
        F_l.inc <- matrix(0, nSample, nT2grid)
        F_u <- array(0, c(nSample, nT1grid, nT2grid))
        F_l <- matrix(0, nSample, nT2grid)
        
        for(M in 1:nSample)
        {
            if(M %% 1000 == 0)
            {
                cat(paste("Calculating posterior predictive density: ", M, " out of ", nSample, " samples", sep = ""), fill = T)
            }
            
            expXbeta1 <- exp(as.vector(x1 %*% (beta1[M, ])))
            expXbeta2 <- exp(as.vector(x2 %*% (beta2[M, ])))
            expXbeta3 <- exp(as.vector(x3 %*% (beta3[M, ])))
            
            expLam1 <- exp(lambda1[M,])
            expLam2 <- exp(lambda2[M,])
            expLam3 <- exp(lam3[M,])
            
            cumHaz1 <- cumsum(expLam1 * diffT1)
            cumHaz2 <- cumsum(expLam2 * diffT2)
            
            for(j in 1:nT2grid)
            {
                if(j <= t2ind)
                {
                    f_l[M, j] <- expLam2[j]*expXbeta2*(1+theta[M,]*(cumHaz1[j]*expXbeta1 + cumHaz2[j]*expXbeta2))^(-1/theta[M,]-1)
                    F_l.inc[M, j] <- diffT2[j] * f_l[M, j]
                    F_l[M, j] <- sum(F_l.inc[M, 1:j])
                }
            }
            
            if(fit$setup$model == "semi-Markov")
            {
                for(i in 1:nT1grid)
                {
                    
                    for(j in 1:nT2grid)
                    {
                        if(T1grid[i] < T2grid[j] & i <= t1ind & j <= t2ind)
                        {
                            Tstar = T2grid[j] - T1grid[i]
                            ind3 <- min(which(Tstar - time3 <= 0))
                            time3.new <- time3[1:ind3]
                            time3.new[ind3] <- Tstar
                            diff.time3 <- diff(c(0, time3.new))
                            expLam3.new <- expLam3[1:ind3]
                            cumHaz3.diff <- cumsum(expLam3.new * diff.time3)
                            
                            f_u[M, i, j] <- (1+theta[M,])*expLam1[i]*expLam3.new[ind3]*expXbeta1*expXbeta3*
                            (1+theta[M,]*(cumHaz1[i]*expXbeta1 + cumHaz2[i]*expXbeta2 + (cumHaz3.diff[ind3])*expXbeta3))^(-1/theta[M,]-2)
                            F_u.inc[M, i, j] <- diffT1[i] * diffT2[j] * f_u[M, i, j]
                            F_u[M, i, j] <- sum(F_u.inc[M, 1:i, 1:j])
                        }
                    }
                }
            }else
            {
                for(i in 1:nT1grid)
                {
                    for(j in 1:nT2grid)
                    {
                        if(T1grid[i] < T2grid[j] & i <= t1ind & j <= t2ind)
                        {
                            ind3.T1 <- min(which(T1grid[i] - time3 <= 0))
                            ind3.T2 <- min(which(T2grid[j] - time3 <= 0))
                            time3.new1 <- time3[1:ind3.T1]
                            time3.new2 <- time3[1:ind3.T2]
                            time3.new1[ind3.T1] <- T1grid[i]
                            time3.new2[ind3.T2] <- T2grid[j]
                            diff.time3.new1 <- diff(c(0, time3.new1))
                            diff.time3.new2 <- diff(c(0, time3.new2))
                            expLam3.new1 <- expLam3[1:ind3.T1]
                            expLam3.new2 <- expLam3[1:ind3.T2]
                            cumHaz3.T1 <- cumsum(expLam3.new1 * diff.time3.new1)
                            cumHaz3.T2 <- cumsum(expLam3.new2 * diff.time3.new2)
                            
                            f_u[M, i, j] <- (1+theta[M,])*expLam1[i]*expLam3.new2[ind3.T2]*expXbeta1*expXbeta3*
                            (1+theta[M,]*(cumHaz1[i]*expXbeta1 + cumHaz2[i]*expXbeta2 + (cumHaz3.T2[ind3.T2]-cumHaz3.T1[ind3.T1])*expXbeta3))^(-1/theta[M,]-2)
                            F_u.inc[M, i, j] <- diffT1[i] * diffT2[j] * f_u[M, i, j]
                            F_u[M, i, j] <- sum(F_u.inc[M, 1:i, 1:j])
                        }
                    }
                }
            }
            
        }
        #list(F_u = mean(F_u[, t1ind, t2ind]), F_l = mean(F_l[,t2ind]), f_u = f_u, f_l = f_l)
        list(F_u = mean(F_u[, t1ind, t2ind]), F_l = mean(F_l[,t2ind]))
    }
}



