FreqSurv_HReg <- function(Formula, data, na.action = "na.fail", subset=NULL)
{
    if(na.action != "na.fail" & na.action != "na.omit")
    {
        stop("na.action should be either na.fail or na.omit")
    }
    form2 <- as.Formula(paste(Formula[2], Formula[1], Formula[3], sep = ""))
    data <- model.frame(form2, data=data, na.action = na.action, subset = subset)
    
    ##
    time1 <- model.part(Formula, data=data, lhs=1)    
    Y <- cbind(time1[1], time1[2])
    y     <- as.vector(Y[,1])
    delta <- as.vector(Y[,2])
    Xmat <- model.frame(formula(Formula, lhs=0, rhs=1), data=data)
    
    ##
    fit.survreg <- survreg(as.formula(paste("Surv(y, delta) ", as.character(formula(Formula, lhs=0, rhs=1))[1], as.character(formula(Formula, lhs=0, rhs=1))[2])), dist="weibull", data=data)
    alpha    <- 1 / fit.survreg$scale
    
    ## log(kappa), log(alpha), log(beta)
    startVals <- c(-alpha*coef(fit.survreg)[1], log(alpha), -coef(fit.survreg)[-1]*alpha)
    ##
    fit0 <- suppressWarnings(nlm(logLike.weibull.Uni, p=startVals * runif(length(startVals), 0.9, 1.1),
    y=y, delta=delta, Xmat=as.matrix(Xmat),
    iterlim=1000, hessian=TRUE))
    ##
    if(fit0$code == 1 | fit0$code == 2)
    {
        myLabels <- c("log(kappa)", "log(alpha)", colnames(Xmat))
        value <- list(estimate=fit0$estimate, Finv=solve(fit0$hessian), logLike=-fit0$minimum, myLabels=myLabels)
        value$class <- c("Freq_HReg", "Surv", "Ind", "WB")
        
        class(value) <- "Freq_HReg"
        return(value)
    }
    
    ##
    invisible()
}
