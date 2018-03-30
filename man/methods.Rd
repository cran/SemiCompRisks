\name{methods}
\alias{print.Bayes_HReg}
\alias{print.Bayes_AFT}
\alias{print.summ.Bayes_HReg}
\alias{print.summ.Bayes_AFT}
\alias{print.Freq_HReg}
\alias{print.summ.Freq_HReg}
\alias{summary.Bayes_HReg}
\alias{summary.Bayes_AFT}
\alias{summary.Freq_HReg}
\alias{coef.Bayes_HReg}
\alias{coef.Bayes_AFT}
\alias{coef.Freq_HReg}
\alias{predict.Bayes_HReg}
\alias{plot.pred.Bayes_HReg}
\alias{predict.Bayes_AFT}
\alias{plot.pred.Bayes_AFT}
\alias{predict.Freq_HReg}
\alias{plot.pred.Freq_HReg}
\alias{vcov.Freq_HReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Methods for objects of classes, \code{Bayes_HReg}/\code{Bayes_AFT}/\code{Freq_HReg}.
}
\description{
The \code{Bayes_HReg} class represents results from Bayesian analysis of semi-competing risks or univariate time-to-event data in the context of hazard regression models.\cr
The \code{Bayes_AFT} class represents results from Bayesian analysis of semi-competing risks or univariate time-to-event data in the context of AFT models.\cr
The \code{Freq_HReg} class represents results from Frequentist analysis of semi-competing risks or univariate time-to-event data in the context of hazard regression models.\cr
}
\usage{
\method{print}{Bayes_HReg}(x, digits=3, ...)
\method{print}{Bayes_AFT}(x, digits=3, ...)
\method{print}{summ.Bayes_HReg}(x, digits=3, ...)
\method{print}{summ.Bayes_AFT}(x, digits=3, ...)
\method{print}{Freq_HReg}(x, digits=3, ...)
\method{print}{summ.Freq_HReg}(x, digits=3, ...)
\method{summary}{Bayes_HReg}(object, digits=3, ...)
\method{summary}{Bayes_AFT}(object, digits=3, ...)
\method{summary}{Freq_HReg}(object, digits=3, ...)
\method{coef}{Bayes_HReg}(object, ...)
\method{coef}{Bayes_AFT}(object, ...)
\method{coef}{Freq_HReg}(object, ...)
\method{predict}{Bayes_HReg}(object, xnew=NULL, x1new=NULL, x2new=NULL,
x3new=NULL, tseq=c(0, 5, 10), ...)
\method{plot}{pred.Bayes_HReg}(x, plot.est="Haz", xlab=NULL, ylab=NULL, ...)
\method{predict}{Bayes_AFT}(object, xnew=NULL, x1new=NULL, x2new=NULL,
x3new=NULL, time, tseq=c(0, 5, 10), ...)
\method{plot}{pred.Bayes_AFT}(x, plot.est="Haz", xlab=NULL, ylab=NULL, ...)
\method{predict}{Freq_HReg}(object, xnew=NULL, x1new=NULL, x2new=NULL,
x3new=NULL, tseq=c(0, 5, 10), ...)
\method{plot}{pred.Freq_HReg}(x, plot.est="Haz", xlab=NULL, ylab=NULL, ...)
\method{vcov}{Freq_HReg}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    an object of class \code{Bayes_HReg} or \code{Bayes_AFT} or \code{Freq_HReg}.
}
  \item{digits}{
    a numeric value indicating the number of digits to display.
}
  \item{object}{
    an object of class \code{Bayes_HReg} or \code{Bayes_AFT} or\code{Freq_HReg}.
}
\item{time}{
    the points at which the baseline survival/hazard functions are evaluated.
}
  \item{tseq}{
	the points at which tick-marks are to be drawn. 
	Required only if the object \code{x} is returned by parametric Weibull-HReg/log-Normal-AFT/DPM-AFT models.
}
  \item{plot.est}{
  	used only if \code{plot} is \code{TRUE}. 
  	If \code{Surv} (the default) then estimated survival functions are produced.
  	If \code{Haz} then estimated hazard functions are produced.
}
  \item{xlab}{
	a title for the x axis.
}
  \item{ylab}{
	a title for the y axis.
}
    \item{xnew}{
    a vector of covariate values with which to predict for which to predict for \eqn{h}.
}
    \item{x1new}{
    a vector of covariate values with which to predict for which to predict for \eqn{h_1}.
}
    \item{x2new}{
    a vector of covariate values with which to predict for which to predict for \eqn{h_2}.
}
    \item{x3new}{
    a vector of covariate values with which to predict for which to predict for \eqn{h_3}.
}
	\item{...}{
	additional arguments.
}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{BayesID_HReg}}, \code{\link{BayesID_AFT}}, \code{\link{BayesSurv_HReg}}, \code{\link{BayesSurv_AFT}}, \code{\link{FreqID_HReg}}, \code{\link{FreqSurv_HReg}}.
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ frequentist framework }
