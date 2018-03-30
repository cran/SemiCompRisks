initiate.startValues_AFT <- function(Formula, data, model, nChain=1,
beta1=NULL, beta2=NULL, beta3=NULL, beta=NULL,
gamma=NULL,
theta=NULL,
y1=NULL, y2=NULL, y=NULL,
LN.mu=NULL, LN.sigSq=NULL,
DPM.class1=NULL, DPM.class2=NULL, DPM.class3=NULL, DPM.class=NULL,
DPM.mu1=NULL, DPM.mu2=NULL, DPM.mu3=NULL, DPM.mu=NULL,
DPM.zeta1=NULL, DPM.zeta2=NULL, DPM.zeta3=NULL, DPM.zeta=NULL,
DPM.tau=NULL)
{
    ret <- vector("list", nChain)
    chain = 1
    
    while(chain <= nChain){
    
    ## BayesSurv_AFT
    
    if(length(Formula)[2]==1)
    {
        cat(paste("Start values are initiated for univariate ", model, " model...", sep = ""), cat("\n"))
        
        ##
        LT <- model.part(Formula, data=data, lhs=1)
        y.mat <- model.part(Formula, data=data, lhs=2)
        
        Y <- cbind(y.mat, LT)
        
        Xmat <- model.frame(formula(Formula, lhs=0, rhs=1), data=data)
        
        p <- ncol(Xmat)
        n <- nrow(Xmat)
        
        ##
        if(!is.null(beta1)) stop(paste("'beta1' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(beta2)) stop(paste("'beta2' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(beta3)) stop(paste("'beta3' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(gamma)) stop(paste("'gamma' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(theta)) stop(paste("'theta' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(y1)) stop(paste("'y1' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(y2)) stop(paste("'y2' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.class1)) stop(paste("'DPM.class1' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.class2)) stop(paste("'DPM.class2' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.class3)) stop(paste("'DPM.class3' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.mu1)) stop(paste("'DPM.mu1' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.mu2)) stop(paste("'DPM.mu2' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.mu3)) stop(paste("'DPM.mu3' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.zeta1)) stop(paste("'DPM.zeta1' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.zeta2)) stop(paste("'DPM.zeta2' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.zeta3)) stop(paste("'DPM.zeta3' is for semi-competing risks models so it must be specified as NULL for univariate ", model, " model.", sep = ""))
        
        
        ##
        if(!is.null(beta)){
            if(length(beta) != p) stop(paste("Length of starting value for beta must be", p))
        }
        if(is.null(beta)) beta <- runif(p, -0.1, 0.1)
        
        ##
        if(!is.null(y)){
            if(length(y) != n) stop("Length of starting value for y must be n")
        }
        if(is.null(y))
        {
            y		<- log(Y[,1])
            y[y == -Inf] <- 0
        }
        
        
        ### for LN model
        if(model == "LN")
        {
            ##
            if(!is.null(LN.mu)){
                if(length(LN.mu) != 1) stop("Length of starting value for LN.mu must be 1 for univariate survival analysis")
            }
            if(is.null(LN.mu)) LN.mu <- runif(1, -0.1, 0.1)
            
            ##
            if(!is.null(LN.sigSq)){
                if(length(LN.sigSq) != 1) stop("Length of starting value for LN.sigSq must be 1 for univariate survival analysis")
            }
            if(is.null(LN.sigSq)) LN.sigSq <- runif(1, 0.5, 1.5)
            
        }
        
        
        ### for DPM model
        
        if(model == "DPM")
        {
            ##
            if(!is.null(DPM.class)){
                if(length(DPM.class) != n) stop(paste("Length of starting value for DPM.class must be", n))
            }
            if(is.null(DPM.class)) DPM.class <- sample(1:2, size=n, replace=TRUE)
            
            ##
            if(!is.null(DPM.tau)){
                if(length(DPM.tau) != 1) stop("Length of starting value for DPM.tau must be 1 for univariate survival analysis")
            }
            if(is.null(DPM.tau)) DPM.tau <- c(0.5)
            
            ##
            if(!is.null(DPM.zeta)){
                if(length(DPM.zeta) != n) stop(paste("Length of starting value for DPM.zeta must be", n))
            }
            if(is.null(DPM.zeta)) DPM.zeta <- rep(1/0.01, n)
            
            ##
            if(!is.null(DPM.mu)){
                if(length(DPM.mu) != n) stop(paste("Length of starting value for DPM.mu must be", n))
            }
            if(is.null(DPM.mu)) DPM.mu <- rep(1, n)
        }
        
        ##
        start.common <- list(beta=beta, y=y)
        start.LN     <- list(LN.mu=LN.mu, LN.sigSq=LN.sigSq)
        start.DPM   <- list(DPM.class=DPM.class, DPM.mu=DPM.mu, DPM.zeta=DPM.zeta, DPM.tau=DPM.tau)
        
        ##
        if(model == "LN"){
            value <- list(common=start.common, LN=start.LN)
        }
        if(model == "DPM"){
            value <- list(common=start.common, DPM=start.DPM)
        }
        
    }
    
    ## BayesID_AFT
    
    if(length(Formula)[2]==3)
    {
        cat(paste("Start values are initiated for semi-competing risks ", model[1], " model...", sep = ""), cat("\n"))
        
        ##
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
        n <- nrow(Xmat1)
        
        ##
        if(!is.null(beta)) stop(paste("'beta' is for univariate models so it must be specified as NULL for semi-competing risks ", model, " model.", sep = ""))
        if(!is.null(y)) stop(paste("'y' is for univariate models so it must be specified as NULL for semi-competing risks ", model, " model.", sep = ""))
        if(!is.null(DPM.class)) stop(paste("'DPM.class' is for univariate models so it must be specified as NULL for semi-competing risks ", model, " model.", sep = ""))
        if(!is.null(DPM.mu)) stop(paste("'DPM.mu' is for univariate models so it must be specified as NULL for semi-competing risks ", model, " model.", sep = ""))
        if(!is.null(DPM.zeta)) stop(paste("'DPM.zeta' is for univariate models so it must be specified as NULL for semi-competing risks ", model, " model.", sep = ""))
        
        
        ##
        if(!is.null(beta1)){
            if(length(beta1) != p1) stop(paste("Length of starting value for beta1 must be", p1))
        }
        if(is.null(beta1)) beta1 <- runif(p1, -0.1, 0.1)
        ##
        if(!is.null(beta2)){
            if(length(beta2) != p2) stop(paste("Length of starting value for beta2 must be", p2))
        }
        if(is.null(beta2)) beta2 <- runif(p2, -0.1, 0.1)
        ##
        if(!is.null(beta3)){
            if(length(beta3) != p3) stop(paste("Length of starting value for beta3 must be", p3))
        }
        if(is.null(beta3)) beta3 <- runif(p3, -0.1, 0.1)
        
        ##
        if(!is.null(theta)){
            if(length(theta) != 1) stop("Length of starting value for theta must be 1")
        }
        if(is.null(theta)) theta <- runif(1, 0.1, 1.1)
        
        if(!is.null(gamma)){
            if(length(gamma) != n) stop("Length of starting value for gamma must be n")
        }
        if(is.null(gamma)) gamma <- rnorm(n, 0, sqrt(theta))
        
        if(!is.null(y1)){
            if(length(y1) != n) stop("Length of starting value for y1 must be n")
        }
        if(is.null(y1))
        {
            y1		<- log(Y[,1])
            y1[y1 == -Inf] <- 0
        }
        if(!is.null(y2)){
            if(length(y2) != n) stop("Length of starting value for y2 must be n")
        }
        if(is.null(y2))
        {
            y2		<- log(Y[,3])
            y2[y2 == -Inf] <- 0
        }
        
        ### for LN model
        if(model == "LN")
        {
            ##
            if(!is.null(LN.mu)){
                if(length(LN.mu) != 3) stop("Length of starting value for LN.mu must be 3 for semi-competing risks analysis")
            }
            if(is.null(LN.mu)) LN.mu <- runif(3, -0.1, 0.1)
            
            ##
            if(!is.null(LN.sigSq)){
                if(length(LN.sigSq) != 3) stop("Length of starting value for LN.sigSq must be 3 for semi-competing risks analysis")
            }
            if(is.null(LN.sigSq)) LN.sigSq <- runif(3, 0.5, 1.5)
            
        }
        
        
        ### for DPM model
        
        if(model == "DPM")
        {
            ##
            if(!is.null(DPM.class1)){
                if(length(DPM.class1) != n) stop(paste("Length of starting value for DPM.class1 must be", n))
            }
            if(is.null(DPM.class1)) DPM.class1 <- sample(1:2, size=n, replace=TRUE)
            
            ##
            if(!is.null(DPM.class2)){
                if(length(DPM.class2) != n) stop(paste("Length of starting value for DPM.class2 must be", n))
            }
            if(is.null(DPM.class2)) DPM.class2 <- sample(1:2, size=n, replace=TRUE)
            
            ##
            if(!is.null(DPM.class3)){
                if(length(DPM.class3) != n) stop(paste("Length of starting value for DPM.class3 must be", n))
            }
            if(is.null(DPM.class3)) DPM.class3 <- sample(1:2, size=n, replace=TRUE)
            
            ##
            if(!is.null(DPM.tau)){
                if(length(DPM.tau) != 3) stop("Length of starting value for DPM.tau must be 3 for semi-competing risks analysis")
            }
            if(is.null(DPM.tau)) DPM.tau <- c(0.5, 0.5, 0.5)
            
            ##
            if(!is.null(DPM.zeta1)){
                if(length(DPM.zeta1) != n) stop(paste("Length of starting value for DPM.zeta1 must be", n))
            }
            if(is.null(DPM.zeta1)) DPM.zeta1 <- rep(1/0.01, n)
            
            ##
            if(!is.null(DPM.zeta2)){
                if(length(DPM.zeta2) != n) stop(paste("Length of starting value for DPM.zeta2 must be", n))
            }
            if(is.null(DPM.zeta2)) DPM.zeta2 <- rep(1/0.01, n)
            
            ##
            if(!is.null(DPM.zeta3)){
                if(length(DPM.zeta3) != n) stop(paste("Length of starting value for DPM.zeta3 must be", n))
            }
            if(is.null(DPM.zeta3)) DPM.zeta3 <- rep(1/0.01, n)
            
            ##
            if(!is.null(DPM.mu1)){
                if(length(DPM.mu1) != n) stop(paste("Length of starting value for DPM.mu1 must be", n))
            }
            if(is.null(DPM.mu1)) DPM.mu1 <- rep(1, n)
            
            ##
            if(!is.null(DPM.mu2)){
                if(length(DPM.mu2) != n) stop(paste("Length of starting value for DPM.mu2 must be", n))
            }
            if(is.null(DPM.mu2)) DPM.mu2 <- rep(1, n)
            
            ##
            if(!is.null(DPM.mu3)){
                if(length(DPM.mu3) != n) stop(paste("Length of starting value for DPM.mu3 must be", n))
            }
            if(is.null(DPM.mu3)) DPM.mu3 <- rep(1, n)
            
            
        }
        
        
        
        ##
        start.common <- list(beta1=beta1, beta2=beta2, beta3=beta3, gamma=gamma, theta=theta, y1=y1, y2=y2)
        start.LN     <- list(LN.mu=LN.mu, LN.sigSq=LN.sigSq)
        start.DPM   <- list(DPM.class1=DPM.class1, DPM.class2=DPM.class2, DPM.class3=DPM.class3, DPM.mu1=DPM.mu1, DPM.mu2=DPM.mu2, DPM.mu3=DPM.mu3, DPM.zeta1=DPM.zeta1, DPM.zeta2=DPM.zeta2, DPM.zeta3=DPM.zeta3, DPM.tau=DPM.tau)
        
        ##
        if(model == "LN"){
            value <- list(common=start.common, LN=start.LN)
        }
        if(model == "DPM"){
            value <- list(common=start.common, DPM=start.DPM)
        }
        
    }
    
    ret[[chain]] <- value
    chain = chain + 1
    
    } # while(chain <= nChain)
    
    
    ##
    return(ret)
}






