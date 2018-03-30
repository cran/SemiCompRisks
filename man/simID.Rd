\name{simID}
\alias{simID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function that simulates independent/cluster-correlated semi-competing risks data under semi-Markov Weibull/Weibull-MVN models.
}
\description{
The function to simulate independent/cluster-correlated semi-competing risks data under semi-Markov Weibull/Weibull-MVN models.
}
\usage{
simID(id=NULL, x1, x2, x3, beta1.true, beta2.true, beta3.true,
			alpha1.true, alpha2.true, alpha3.true, 
			kappa1.true, kappa2.true, kappa3.true, 
			theta.true, SigmaV.true=NULL, cens)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
	a vector of cluster information for \code{n} subjects. The cluster membership must be set to consecutive positive integers, \eqn{1:J}. 
	Required only when generating clustered data.
}
  \item{x1}{
	covariate matrix, \code{n} observations by \code{p1} variables.
}
  \item{x2}{
	covariate matrix, \code{n} observations by \code{p2} variables.
}
  \item{x3}{
	covariate matrix, \code{n} observations by \code{p3} variables.
}
  \item{beta1.true}{
	true value for \eqn{\beta_1}.
}
  \item{beta2.true}{
	true value for \eqn{\beta_2}.
}
  \item{beta3.true}{
	true value for \eqn{\beta_3}.
}
  \item{alpha1.true}{
	true value for \eqn{\alpha_1}.
}
  \item{alpha2.true}{
	true value for \eqn{\alpha_2}.
}
  \item{alpha3.true}{
	true value for \eqn{\alpha_3}.
}
  \item{kappa1.true}{
	true value for \eqn{\kappa_1}.
}
  \item{kappa2.true}{
	true value for \eqn{\kappa_2}.
}
  \item{kappa3.true}{
	true value for \eqn{\kappa_3}.
}
  \item{theta.true}{
	true value for \eqn{\theta}.
}
  \item{SigmaV.true}{
	true value for \eqn{\Sigma_V}. Required only when generating clustered data.
}
  \item{cens}{
	a vector with two numeric elements. The right censoring times are generated from Uniform(\eqn{cens[1]}, \eqn{cens[2]}).
}
}

 
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{simIDcor} returns a data.frame containing semi-competing risks outcomes from \code{n} subjects. 
It is of dimension \eqn{n\times 4}: the columns correspond to \eqn{y_1}, \eqn{\delta_1}, \eqn{y_2}, \eqn{\delta_2}. \cr
\item{y1}{a vector of \code{n} times to the non-terminal event}
\item{y2}{a vector of \code{n} times to the terminal event}
\item{delta1}{a vector of \code{n} censoring indicators for the non-terminal event time (1=event occurred, 0=censored)}
\item{delta2}{a vector of \code{n} censoring indicators for the terminal event time (1=event occurred, 0=censored)}
}
\author{
Kyu Ha Lee and Sebastien Haneuse \cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}

\examples{
	library(MASS)
	set.seed(123456)
	
	J = 110
	nj = 50
	n = J * nj

	id <- rep(1:J, each = nj)

	kappa1.true <- 0.05
	kappa2.true <- 0.01
	kappa3.true <- 0.01
	alpha1.true <- 0.8
	alpha2.true <- 1.1
	alpha3.true <- 0.9
	beta1.true <- c(0.5, 0.8, -0.5)
	beta2.true <- c(0.5, 0.8, -0.5)
	beta3.true <- c(1, 1, -1)
	SigmaV.true <- matrix(0.25,3,3)
	theta.true <- 0.5
	cens <- c(90, 90)

	cov1 <- matrix(rnorm((length(beta1.true)-1)*n, 0, 1), n, length(beta1.true)-1)
	cov2 <- sample(c(0, 1), n, replace = TRUE)
	x1 <- as.data.frame(cbind(cov1, cov2))
	x2 <- as.data.frame(cbind(cov1, cov2))
	x3 <- as.data.frame(cbind(cov1, cov2))
	
	simData <- simID(id, x1, x2, x3, beta1.true, beta2.true, beta3.true, 
				alpha1.true, alpha2.true, alpha3.true, 
				kappa1.true, kappa2.true, kappa3.true, 
				theta.true, SigmaV.true, cens)	
   				    				
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line
