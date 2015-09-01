\name{simSurv}
\alias{simSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function that simulates independent/cluster-correlated right-censored survival data under Weibull/Weibull-Normal model.
}
\description{
The function to simulate independent/cluster-correlated right-censored survival data under Weibull/Weibull-Normal model.
}
\usage{
simSurv(cluster=NULL, x, beta.true, alpha.true, kappa.true, sigmaV.true=NULL, cens)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cluster}{
	a vector of cluster information for \code{n} subjects. The cluster membership must be set to consecutive positive integers, \eqn{1:J}. 
	Required only when generating clustered data.
}
  \item{x}{
	covariate matrix, \code{n} observations by \code{p} variables.
}
  \item{beta.true}{
	true value for \eqn{\beta}.
}
  \item{alpha.true}{
	true value for \eqn{\alpha}.
}
  \item{kappa.true}{
	true value for \eqn{\kappa}.
}
  \item{sigmaV.true}{
	true value for \eqn{\sigma_V}. 
	Required only when generating clustered data.
}
  \item{cens}{
	a vector with two numeric elements. The right censoring times are generated from Uniform(\eqn{cens[1]}, \eqn{cens[2]}).
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{simSurv} returns a data.frame containing univariate time-to-event outcomes from \code{n} subjects. 
It is of dimension \eqn{n\times 2}: the columns correspond to \eqn{y}, \eqn{\delta}.\cr
\item{y}{a vector of \code{n} times to the event}
\item{delta}{a vector of \code{n} censoring indicators for the event time (1=event occurred, 0=censored)}
}
\author{
Kyu Ha Lee and Sebastien Haneuse \cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}

\examples{
	set.seed(123456)
	
	J = 110
	nj = 50
	n = J * nj

	cluster <- rep(1:J, each = nj)	

	x	= matrix(0, n, 2)	
	x[,1]	= rnorm(n, 0, 2)	
	x[,2]	= sample(c(0, 1), n, replace = TRUE)

	beta.true = c(0.5, 0.5)
	
	alpha.true = 1.5		
	kappa.true = 0.02
	sigmaV.true = 0.1

	cens <- c(30, 40)		

	simData <- simSurv(cluster, x, beta.true, alpha.true, kappa.true, 
				sigmaV.true, cens) 		 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ univariate analysis }% __ONLY ONE__ keyword per line
