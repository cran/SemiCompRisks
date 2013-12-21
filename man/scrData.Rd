\name{scrData}
\alias{scrData}
\docType{data}
\title{
Simulated semi-competing risks data using \code{simID} function.
}
\description{
Simulated semi-competing risks data using \code{simID} function. See \code{Examples} for the setting of data generation.
}
\usage{data(scrData)}
\format{
  A data frame with 2000 observations on the following 10 variables.
  \describe{
    \item{\code{time1}}{the time to non-terminal event}
    \item{\code{event1}}{the censoring indicators for the non-terminal event time; 1=event observed, 0=censored/truncated}
    \item{\code{time2}}{the time to terminal event}
    \item{\code{event2}}{the censoring indicators for the terminal event time; 1=event observed, 0=censored}
    \item{\code{x1.1}}{the first column of covariate matrix x1}
    \item{\code{x1.2}}{the second column of covariate matrix x1}
    \item{\code{x2.1}}{the first column of covariate matrix x2}
    \item{\code{x2.2}}{the second column of covariate matrix x2}
    \item{\code{x3.1}}{the first column of covariate matrix x3}
    \item{\code{x3.2}}{the second column of covariate matrix x3}
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

### simulating semi-competing risks data ###
	
set.seed(123456)
	
n = 2000
p1 = 2
p2 = 2
p3 = 2		

x		= matrix(0, n, p1)
colnames(x)	= c("cov1", "cov2")	
x[,1]		= rnorm(n, 0, 2)	
x[,2]		= sample(c(0, 1), n, replace = TRUE)
	
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
alpha2.true, alpha3.true, kappa1.true, kappa2.true, kappa3.true, gamma.true, cens)
   				 
y1 <- simData$y1
y2 <- simData$y2
delta1 <- simData$delta1
delta2 <- simData$delta2
   				 
scrData <- as.data.frame(cbind(y1, delta1, y2, delta2, x1, x2, x3))

colnames(scrData)[1:4]	<- c("time1", "event1", "time2", "event2")  
}
\keyword{datasets}
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line