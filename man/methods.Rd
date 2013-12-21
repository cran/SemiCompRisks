\name{methods}
\alias{print.BayesSurv}
\alias{print.BayesID}
\alias{summary.BayesSurv}
\alias{summary.BayesID}
\alias{plot.BayesSurv}
\alias{plot.BayesID}
\alias{plot.ehr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for objects of classes, "BayesSurv"/"BayesID"
}
\description{
The \code{BayesSurv} class represents results from analysis of univariate time-to-event data.\cr
The \code{BayesID} class represents results from analysis of semi-competing risks data.
}
\usage{
\method{print}{BayesSurv}(x, ...)
\method{print}{BayesID}(x, ...)
\method{summary}{BayesSurv}(object, ...)
\method{summary}{BayesID}(object, ...)
\method{plot}{BayesSurv}(x, tlim = c(0, 10), ...)
\method{plot}{BayesID}(x, tlim = c(0, 10), ...)
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
	The time (x-axis) limits (x1, x2) of the plot. Only required for an object returned by parametric models.
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
	\code{\link{BayesSurv}}, \code{\link{BayesID}}, \code{\link{ehr}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
