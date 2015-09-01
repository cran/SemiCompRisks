initiate.startValues <- function(Y, lin.pred, data, model, cluster=NULL,
                                 beta1=NULL, beta2=NULL, beta3=NULL, beta=NULL,
                                 gamma.ji=NULL,
                                 theta=NULL,
                                 V.j1=NULL, V.j2=NULL, V.j3=NULL, V.j = NULL,
                                 WB.alpha=NULL, WB.kappa=NULL,
                                 MVN.SigmaV=NULL, Normal.zeta=NULL,
                                 DPM.class=NULL, DPM.tau=NULL)
{
    
    
    ## BayesSurvcor
    
    if(length(model)==2 & class(lin.pred)=="formula")
    {
      print(paste("Start values are initiated for univariate ", model[1],"-", model[2], " model...", sep = ""), cat("\n"))
      ##
      if(is.null(cluster)) stop(paste("'cluster' must be given for ", model[1],"-", model[2], " model.", sep = ""))
      
        ##
        n <- nrow(Y)
        J <- length(unique(cluster))
        p <- ncol(model.frame(lin.pred, data=data))
        
        ##
        if(!is.null(beta1)) stop(paste("'beta1' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(beta2)) stop(paste("'beta2' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(beta3)) stop(paste("'beta3' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(gamma.ji)) stop(paste("'gamma.ji' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(theta)) stop(paste("'theta' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(V.j1)) stop(paste("'V.j1' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(V.j2)) stop(paste("'V.j2' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(V.j3)) stop(paste("'V.j3' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(MVN.SigmaV)) stop(paste("'MVN.SigmaV' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        
        
        ##
        if(!is.null(beta)){
            if(length(beta) != p) stop(paste("Length of starting value for beta must be", p))
        }
        if(is.null(beta)) beta <- runif(p, -0.1, 0.1)
        
        ##
        if(!is.null(V.j)){
            if(length(V.j) != J) stop(paste("Length of starting values for V.j must be", J))
        }
        if(is.null(V.j)) V.j <- runif(J, -0.1, 0.1)
        
        ##
        if(!is.null(WB.alpha)){
            if(length(WB.alpha) != 1) stop("Length of starting value for WB.alpha must be 1")
        }
        if(is.null(WB.alpha)) WB.alpha <- 1
        ##
        if(!is.null(WB.kappa)){
            if(length(WB.kappa) != 1) stop("Length of starting value for WB.kappa should be 1")
        }
        if(is.null(WB.kappa)) WB.kappa <- 0.01
        
        ##
        if(!is.null(Normal.zeta)){
            if(length(Normal.zeta) != 1) stop("Length of starting value for Normal.zeta should be 1")
        }
        if(is.null(Normal.zeta)) Normal.zeta <- 1
        
        ##
        if(!is.null(DPM.class)){
            if(length(DPM.class) != J) stop(paste("Length of starting value for DPM.class must be", J))
        }
        if(is.null(DPM.class)) DPM.class <- sample(1:3, size=J, replace=TRUE)
        ##
        if(!is.null(DPM.tau)){
            if(length(DPM.tau) != 1) stop("Length of starting value for DPM.tau must be 1")
        }
        if(is.null(DPM.tau)) DPM.tau <- 0.5
        
        ##
        start.common <- list(beta=beta, V.j=V.j)
        start.WB     <- list(WB.alpha=WB.alpha, WB.kappa=WB.kappa)
        start.Normal    <- list(Normal.zeta=Normal.zeta)
        start.DPM    <- list(DPM.class=DPM.class, DPM.tau=DPM.tau)
        
        ##
        if(model[1] == "Weibull"){
            if(model[2] == "Normal") value <- list(common=start.common, WB=start.WB, Normal=start.Normal)
            if(model[2] == "DPM") value <- list(common=start.common, WB=start.WB, DPM=start.DPM)
        }
        if(model[1] == "PEM"){
            if(model[2] == "Normal") value <- list(common=start.common, Normal=start.Normal)
            if(model[2] == "DPM") value <- list(common=start.common, DPM=start.DPM)
        }
    }
    
    
    ## BayesID
    
    if(length(model)==2 & class(lin.pred)=="list")
    {
        print(paste("Start values are initiated for semi-competing risks ", model[2], " model...", sep = ""), cat("\n"))
        ##
        if(!is.null(cluster)) print(paste("Warning: 'cluster' is not required for ", model[2], " model so it is ignored", sep = ""))
        
        ##     
        n <- nrow(Y)
        p1 <- ncol(model.matrix(lin.pred[[1]], data=data)) - 1
        p2 <- ncol(model.matrix(lin.pred[[2]], data=data)) - 1
        p3 <- ncol(model.matrix(lin.pred[[3]], data=data)) - 1
        
        ##
        if(!is.null(beta)) stop(paste("'beta' is for univariate models so it must be specified as NULL for semi-competing risks ", model[2], " model.", sep = ""))
        if(!is.null(V.j)) stop(paste("'V.j' is for univariate models so it must be specified as NULL for semi-competing risks ", model[2], " model.", sep = ""))
        if(!is.null(Normal.zeta)) stop(paste("'Normal.zeta' is for univariate models so it must be specified as NULL for semi-competing risks ", model[2], " model.", sep = ""))
        ##
        if(!is.null(V.j1)) stop(paste("'V.j1' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j2)) stop(paste("'V.j2' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j3)) stop(paste("'V.j3' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j)) stop(paste("'V.j' is for univariate models so it must be specified as NULL for semi-competing risks ", model, " model.", sep = ""))
        ##
        if(!is.null(MVN.SigmaV)) stop(paste("'MVN.SigmaV' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(Normal.zeta)) stop(paste("'Normal.zeta' is for univariate models so it must be specified as NULL for semi-competing risks ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.class)) stop(paste("'DPM.class' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.tau)) stop(paste("'DPM.tau' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        
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
        
        if(!is.null(gamma.ji)){
            if(length(gamma.ji) != n) stop("Length of starting value for gamma.ji must be n")
        }
        if(is.null(gamma.ji)) gamma.ji <- rgamma(n, 1/theta, 1/theta)
        
        ##
        if(!is.null(WB.alpha)){
            if(length(WB.alpha) != 3) stop("Length of starting value for WB.alpha must be 3")
        }
        if(is.null(WB.alpha)) WB.alpha <- c(1, 1, 1)
        ##
        if(!is.null(WB.kappa)){
            if(length(WB.kappa) != 3) stop("Length of starting value for WB.kappa should be 3")
        }
        if(is.null(WB.kappa)) WB.kappa <- c(0.01, 0.01, 0.01)
        
        
        ##
        start.common <- list(beta1=beta1, beta2=beta2, beta3=beta3, gamma.ji=gamma.ji, theta=theta)
        start.WB     <- list(WB.alpha=WB.alpha, WB.kappa=WB.kappa)
        
        ##
        if(model[2] == "Weibull"){
            value <- list(common=start.common, WB=start.WB)
        }
        if(model[2] == "PEM"){
            value <- list(common=start.common)
        }

    }
    
    
    
    ## BayesSurv
    
    if(length(model)==1)
    {
        print(paste("Start values are initiated for univariate ", model, " model...", sep = ""), cat("\n"))
        ##
        if(!is.null(cluster)) print(paste("Warning: 'cluster' is not required for ", model, " model so it is ignored", sep = ""))
        
        ##
        n <- nrow(Y)
        p <- ncol(model.frame(lin.pred, data=data))
        
        
        ##
        if(!is.null(beta1)) stop(paste("'beta1' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(beta2)) stop(paste("'beta2' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(beta3)) stop(paste("'beta3' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(gamma.ji)) stop(paste("'gamma.ji' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(theta)) stop(paste("'theta' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j1)) stop(paste("'V.j1' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j2)) stop(paste("'V.j2' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j3)) stop(paste("'V.j3' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j)) stop(paste("'V.j' is for univariate models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(MVN.SigmaV)) stop(paste("'MVN.SigmaV' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(Normal.zeta)) stop(paste("'Normal.zeta' is for univariate models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.class)) stop(paste("'DPM.class' is for univariate models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.tau)) stop(paste("'DPM.tau' is for univariate models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        
        ##
        if(!is.null(beta)){
            if(length(beta) != p) stop(paste("Length of starting value for beta must be", p))
        }
        if(is.null(beta)) beta <- runif(p, -0.1, 0.1)
        
        
        ##
        if(!is.null(WB.alpha)){
            if(length(WB.alpha) != 1) stop("Length of starting value for WB.alpha must be 1")
        }
        if(is.null(WB.alpha)) WB.alpha <- 1
        ##
        if(!is.null(WB.kappa)){
            if(length(WB.kappa) != 1) stop("Length of starting value for WB.kappa should be 1")
        }
        if(is.null(WB.kappa)) WB.kappa <- 0.01
        
        ##
        start.common <- list(beta=beta)
        start.WB     <- list(WB.alpha=WB.alpha, WB.kappa=WB.kappa)
        
        ##
        if(model == "Weibull"){
            value <- list(common=start.common, WB=start.WB)
        }
        if(model == "PEM"){
            value <- list(common=start.common)
        }
        
    }
    
    
    
    
    ## BayesIDcor
    
  	if(length(model)==3)
  	{
      print(paste("Start values are initiated for semi-competing risks ", model[2],"-", model[3], " model...", sep = ""), cat("\n"))
      ##
      if(is.null(cluster)) stop(paste("'cluster' must be given for ", model[2],"-", model[3], " model.", sep = ""))
      ##
      n <- nrow(Y)
      J <- length(unique(cluster))
      p1 <- ncol(model.matrix(lin.pred[[1]], data=data)) - 1
      p2 <- ncol(model.matrix(lin.pred[[2]], data=data)) - 1
      p3 <- ncol(model.matrix(lin.pred[[3]], data=data)) - 1

        ##
        if(!is.null(beta)) stop(paste("'beta' is for univariate models so it must be specified as NULL for ", model[2],"-", model[3], " model.", sep = ""))
    
        if(!is.null(V.j)) stop(paste("'V.j' is for univariate models so it must be specified as NULL for ", model[2],"-", model[3], " model.", sep = ""))
        if(!is.null(Normal.zeta)) stop(paste("'Normal.zeta' is for univariate models so it must be specified as NULL for ", model[2],"-", model[3], " model.", sep = ""))

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
      if(!is.null(V.j1)){
          if(length(V.j1) != J) stop(paste("Length of starting values for V.j1 must be", J))
      }
      if(is.null(V.j1)) V.j1 <- runif(J, -0.1, 0.1)
      ##
      if(!is.null(V.j2)){
          if(length(V.j2) != J) stop(paste("Length of starting values for V.j2 must be", J))
      }
      if(is.null(V.j2)) V.j2 <- runif(J, -0.1, 0.1)
      ##
      if(!is.null(V.j3)){
          if(length(V.j3) != J) stop(paste("Length of starting values for V.j3 must be", J))
      }
      if(is.null(V.j3)) V.j3 <- runif(J, -0.1, 0.1)
      
      ##
      if(!is.null(theta)){
          if(length(theta) != 1) stop("Length of starting value for theta must be 1")
      }
      if(is.null(theta)) theta <- runif(1, 0.1, 1.1)
      
      if(!is.null(gamma.ji)){
          if(length(gamma.ji) != n) stop("Length of starting value for gamma.ji must be n")
      }
      if(is.null(gamma.ji)) gamma.ji <- rgamma(n, 1/theta, 1/theta)
      
      ##
      if(!is.null(WB.alpha)){
          if(length(WB.alpha) != 3) stop("Length of starting value for WB.alpha must be 3")
      }
      if(is.null(WB.alpha)) WB.alpha <- c(1, 1, 1)
      ##
      if(!is.null(WB.kappa)){
          if(length(WB.kappa) != 3) stop("Length of starting value for WB.kappa should be 3")
      }
      if(is.null(WB.kappa)) WB.kappa <- c(0.01, 0.01, 0.01)
      
      ##
      if(!is.null(MVN.SigmaV))
      {
          if(is.matrix(MVN.SigmaV) == FALSE) stop("Starting value for MVN.SigmaV must be a 3x3 matrix")
          if(is.matrix(MVN.SigmaV) == TRUE){
              if(nrow(MVN.SigmaV) != 3 | ncol(MVN.SigmaV) != 3) stop("Starting value for MVN.SigmaV must be a 3x3 matrix")
          }
      }
      if(is.null(MVN.SigmaV)) MVN.SigmaV <- diag(0.1, 3)
      
      ##
      if(!is.null(DPM.class)){
          if(length(DPM.class) != J) stop(paste("Length of starting value for DPM.class must be", J))
      }
      if(is.null(DPM.class)) DPM.class <- sample(1:3, size=J, replace=TRUE)
      ##
      if(!is.null(DPM.tau)){
          if(length(DPM.tau) != 1) stop("Length of starting value for DPM.tau must be 1")
      }
      if(is.null(DPM.tau)) DPM.tau <- 0.5
      
      ##
      start.common <- list(beta1=beta1, beta2=beta2, beta3=beta3, gamma.ji=gamma.ji, V.j1=V.j1, V.j2=V.j2, V.j3=V.j3, theta=theta)
      start.WB     <- list(WB.alpha=WB.alpha, WB.kappa=WB.kappa)
      start.MVN    <- list(MVN.SigmaV=MVN.SigmaV)
      start.DPM    <- list(DPM.class=DPM.class, DPM.tau=DPM.tau)
      
      ##
      if(model[2] == "Weibull"){
          if(model[3] == "MVN") value <- list(common=start.common, WB=start.WB, MVN=start.MVN)
          if(model[3] == "DPM") value <- list(common=start.common, WB=start.WB, DPM=start.DPM)
      }
      if(model[2] == "PEM"){
          if(model[3] == "MVN") value <- list(common=start.common, MVN=start.MVN)
          if(model[3] == "DPM") value <- list(common=start.common, DPM=start.DPM)
      }
  }
  
  


  ##
  return(value)
}






