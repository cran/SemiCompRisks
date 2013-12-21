\name{simSurv}
\alias{simSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to simulate right censored survival data under Weibull model
}
\description{
The function to simulate right censored survival data under Weibull model
}
\usage{
simSurv(x, beta.true, alpha.true, kappa.true, cens)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Covariate matrix, \code{n} observations by \code{p} variables
}
  \item{beta.true}{
True value for \eqn{\beta}
}
  \item{alpha.true}{
True value for \eqn{\alpha}
}
  \item{kappa.true}{
True value for \eqn{\kappa}
}
  \item{cens}{
A vector with two numeric elements. The right censoring times are generated from Uniform(\eqn{cens[1]}, \eqn{cens[2]})
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{simSurv} returns an object of class \code{simSurv} \cr
\item{y}{a vector of \code{n} times to the event}
\item{delta}{a vector of \code{n} censoring indicators for the event time (1=event occurred, 0=censored)}
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2013). 
Bayesian Semi-parametric Analysis of Semi-competing Risks Data: 
Estimating Readmission Rates among Pancreatic Cancer Patients, submitted.
}
\author{
Kyu Ha Lee and Sebastien Haneuse \cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

	### simulating survival data ###
	
	set.seed(123456)
	
	n = 2000
	p = 2

	x	= matrix(0, n, p)	
	x[,1]	= rnorm(n, 0, 2)	
	x[,2]	= sample(c(0, 1), n, replace = TRUE)
			
	# - true value of parameter
	
	beta.true = c(0.5, 0.5)
	
	alpha.true = 1.5		
	kappa.true = 0.02

	cens <- c(30, 40)		

	simData <- simSurv(x, beta.true, alpha.true, kappa.true, cens) 
	
   	summary(simData)			 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ univariate analysis }% __ONLY ONE__ keyword per line
