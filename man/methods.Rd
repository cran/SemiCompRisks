\name{methods}
\alias{print.BayesSurv}
\alias{print.BayesSurvcor}
\alias{print.BayesID}
\alias{print.BayesIDcor}
\alias{summary.BayesSurv}
\alias{summary.BayesSurvcor}
\alias{summary.BayesID}
\alias{summary.BayesIDcor}
\alias{plot.BayesSurv}
\alias{plot.BayesSurvcor}
\alias{plot.BayesID}
\alias{plot.BayesIDcor}
\alias{plot.ehr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Methods for objects of classes, \cr
    \code{BayesSurv}/\code{BayesSurvcor}/\code{BayesID}/\code{BayesIDcor}
}
\description{
The \code{BayesSurv} (\code{BayesSurvcor}) class represents results from analysis of (cluster-correlated) univariate time-to-event data.\cr
The \code{BayesID} (\code{BayesIDcor}) class represents results from analysis of (cluster-correlated) semi-competing risks data.
}
\usage{
\method{print}{BayesSurv}(x, ...)
\method{print}{BayesSurvcor}(x, ...)
\method{print}{BayesID}(x, ...)
\method{print}{BayesIDcor}(x, ...)
\method{summary}{BayesSurv}(object, ...)
\method{summary}{BayesSurvcor}(object, ...)
\method{summary}{BayesID}(object, ...)
\method{summary}{BayesIDcor}(object, ...)
\method{plot}{BayesSurv}(x, tlim = c(0, 10), ...)
\method{plot}{BayesSurvcor}(x, tlim = c(0, 10), ...)
\method{plot}{BayesID}(x, tlim = c(0, 10), ...)
\method{plot}{BayesIDcor}(x, tlim = c(0, 10), ...)
\method{plot}{ehr}(x, tlim = c(0, 10), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class; \code{BayesSurv}/\code{BayesID}/\code{ehr}.
}
  \item{object}{
    An object of class; \code{BayesSurv} or \code{BayesID}.
}
  \item{tlim}{
	The time (x-axis) limits (x1, x2) of the plot. Only required for an object returned by parametric (Weibull) models.
}
	\item{...}{
	Additional arguments.
}


}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{BayesSurv}}, \code{\link{BayesSurvcor}}, \code{\link{BayesID}}, \code{\link{BayesIDcor}}, \code{\link{ehr}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
