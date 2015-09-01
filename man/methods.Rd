\name{methods}
\alias{print.Bayes}
\alias{print.summ.Bayes}
\alias{print.Freq}
\alias{print.summ.Freq}
\alias{summary.Bayes}
\alias{summary.Freq}
\alias{plot.Bayes}
\alias{plot.Freq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Methods for objects of classes, \code{Bayes}/\code{Freq}.
}
\description{
The \code{Bayes} class represents results from Bayesian analysis of semi-competing risks or univariate time-to-event data.\cr
The \code{Freq} class represents results from Frequentist analysis of semi-competing risks or univariate time-to-event data.\cr
}
\usage{
\method{print}{Bayes}(x, digits=3, ...)
\method{print}{summ.Bayes}(x, digits=3, ...)
\method{print}{Freq}(x, digits=3, ...)
\method{print}{summ.Freq}(x, digits=3, ...)
\method{summary}{Bayes}(object, digits=3, ...)
\method{summary}{Freq}(object, digits=3, ...)
\method{plot}{Freq}(x, tseq=c(0, 5, 10), plot=TRUE, plot.est="BS", xlab=NULL, ylab=NULL, ...)
\method{plot}{Bayes}(x, tseq=c(0, 5, 10), plot=TRUE, plot.est="BS", xlab=NULL, ylab=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    an object of class \code{Bayes} or \code{Freq}.
}
  \item{digits}{
    a numeric value indicating the number of digits to display.
}
  \item{object}{
    an object of class \code{Bayes} or \code{Freq}.
}
  \item{tseq}{
	the points at which tick-marks are to be drawn. 
	Required only if the object \code{x} is returned by parametric Weibull models. 
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
	\code{\link{BayesID}}, \code{\link{BayesSurv}}, \code{\link{FreqID}}, \code{\link{FreqSurv}}.
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ frequentist framework }
