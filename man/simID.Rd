\name{simID}
\alias{simID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to simulate semi-competing risks data under Weibull model
}
\description{
The function to simulate semi-competing risks data under Weibull model
}
\usage{
simID(x1, x2, x3, beta1.true, beta2.true, beta3.true, alpha1.true, alpha2.true, 
		alpha3.true, kappa1.true, kappa2.true, kappa3.true, gamma.true, cens)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
Covariate matrix, \code{n} observations by \code{p1} variables
}
  \item{x2}{
Covariate matrix, \code{n} observations by \code{p2} variables
}
  \item{x3}{
Covariate matrix, \code{n} observations by \code{p3} variables
}
  \item{beta1.true}{
True value for \eqn{\beta_1}
}
  \item{beta2.true}{
True value for \eqn{\beta_2}
}
  \item{beta3.true}{
True value for \eqn{\beta_3}
}
  \item{alpha1.true}{
True value for \eqn{\alpha_1}
}
  \item{alpha2.true}{
True value for \eqn{\alpha_2}
}
  \item{alpha3.true}{
True value for \eqn{\alpha_3}
}
  \item{kappa1.true}{
True value for \eqn{\kappa_1}
}
  \item{kappa2.true}{
True value for \eqn{\kappa_2}
}
  \item{kappa3.true}{
True value for \eqn{\kappa_3}
}
  \item{gamma.true}{
True value for \eqn{\gamma}
}

  \item{cens}{
A vector with two numeric elements. The right censoring times are generated from Uniform(\eqn{cens[1]}, \eqn{cens[2]})
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{simID} returns an object of class \code{simID} \cr
\item{y1}{a vector of \code{n} times to the non-terminal event}
\item{y2}{a vector of \code{n} times to the terminal event}
\item{delta1}{a vector of \code{n} censoring indicators for the non-terminal event time (1=event occurred, 0=censored)}
\item{delta2}{a vector of \code{n} censoring indicators for the terminal event time (1=event occurred, 0=censored)}
}
%%\references{
%%Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2014), 
%%Bayesian semiparametric analysis of semicompeting risks data: 
%%investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, \emph{in press}.\cr \cr
%%Lee, K. H., Dominici, F., Schrag, D., and Haneuse, S.,
%%Hierarchical models for cluster-correlated semicompeting risks data with application to quality of end-of-life care for pancreatic cancer, \emph{submitted}. \cr
%%}
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

	### simulating semi-competing risks data ###
	
	set.seed(123456)
	
	n = 2000
	p1 = 2
	p2 = 2
	p3 = 2		

	x	= matrix(0, n, p1)	
	x[,1]	= rnorm(n, 0, 2)	
	x[,2]	= sample(c(0, 1), n, replace = TRUE)
	
	x1 = as.matrix(x, ncol = p1)
	x2 = as.matrix(x, ncol = p2)
	x3 = as.matrix(x, ncol = p3)
	
	# - true value of parameter
	
	beta1.true = c(0.5, 0.5)
	beta2.true = c(0.5, 0.5)
	beta3.true = c(1, 1)
	
	theta.true = 2
	gamma.true = rgamma(n, shape = 1/theta.true, rate = 1/theta.true)

	alpha1.true = 1.5
	alpha2.true = 0.8
	alpha3.true = 0.5
		
	kappa1.true = 0.02
	kappa2.true = 0.1
	kappa3.true = 0.15
	
	cens <- c(30, 40)		

	simData <- simID(x1, x2, x3, beta1.true, beta2.true, beta3.true, alpha1.true, 
    			alpha2.true, alpha3.true, kappa1.true, kappa2.true, kappa3.true, 
   				 gamma.true, cens) 
   				    				
   	summary(simData)			 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line
