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
\alias{plot.Bayes_HReg}
\alias{plot.Bayes_AFT}
\alias{plot.Freq_HReg}
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
\method{plot}{Freq_HReg}(x, tseq=c(0, 5, 10), plot=TRUE, plot.est="BS", xlab=NULL, ylab=NULL, ...)
\method{plot}{Bayes_HReg}(x, tseq=c(0, 5, 10), plot=TRUE, plot.est="BS", xlab=NULL, ylab=NULL, ...)
\method{plot}{Bayes_AFT}(x, time, tseq=c(0, 5, 10), plot=TRUE, plot.est="BS", xlab=NULL, ylab=NULL, ...)
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
  \item{plot}{
	If \code{TRUE} (the default) then either estimated baseline hazard functions or estimated baseline survival functions are produced depending on the value of \code{plot.est}. 
	If \code{FALSE}, the summaries which the plots are based on are returned.
}
  \item{plot.est}{
  	used only if \code{plot} is \code{TRUE}. 
  	If \code{BS} (the default) then estimated baseline survival functions are produced. 
  	If \code{BH} then estimated baseline hazard functions are produced.
}
  \item{xlab}{
	a title for the x axis.
}
  \item{ylab}{
	a title for the y axis.
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
