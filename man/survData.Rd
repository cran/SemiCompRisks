\name{survData}
\alias{survData}
\docType{data}
\title{
Simulated survival data using \code{simSurv} function.
}
\description{
Simulated survival data using \code{simSurv} function. See \code{Examples} for the setting of data generation.
}
\usage{data(survData)}
\format{
  A data frame with 2000 observations on the following 4 variables.
  \describe{
    \item{\code{time}}{the time to event}
    \item{\code{event}}{the censoring indicators for the event time; 1=event observed, 0=censored}
    \item{\code{cov1}}{the first column of covariate matrix x}
    \item{\code{cov2}}{the second column of covariate matrix x}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
	
### simulating survival data ###
	
set.seed(123456)
	
n = 2000
p = 2

x		= matrix(0, n, p)	
colnames(x)	= c("cov1", "cov2")	
x[,1]		= rnorm(n, 0, 2)	
x[,2]		= sample(c(0, 1), n, replace = TRUE)
				
# - true value of parameter
	
beta.true = c(0.5, 0.5)	
alpha.true = 1.5		
kappa.true = 0.02

cens <- c(30, 40)		

simData <- simSurv(x, beta.true, alpha.true, kappa.true, cens) 
 
y <- simData$y
delta <- simData$delta
	
survData <- as.data.frame(cbind(y, delta, x))
colnames(survData)[1:2]	<- c("time", "event")

}
\keyword{datasets}
\keyword{ univariate analysis }% __ONLY ONE__ keyword per line