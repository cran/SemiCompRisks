\name{BayesSurv}
\alias{BayesSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to fit Bayesian parametric and semi-parametric regression models to univariate survival data.
}
\description{
The function to fit Bayesian parametric and semi-parametric regression models to univariate survival data.
}
\usage{
BayesSurv(survData, hyperParams, startValues, mcmcParams, numReps, thin, path, 
burninPerc = 0.5, type = "semi-parametric", nChain = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
	The data frame containing univariate time-to-event outcome with covariate matrix from \code{n} subjects. See *Examples*. 
}
  \item{hyperParams}{
	a vector containing hyperparameter values for hierarchical model. See *Details* and *Examples*.
}
  \item{startValues}{
  	a list containing vectors of starting values of the parameters. See *Details* and *Examples*.
}
  \item{mcmcParams}{
	a vector containing a variable required for Metropolis-Hastings (MH) algorithm. See *Details* and *Examples*.
}
  \item{numReps}{
	total number of scans
}
  \item{thin}{
	extent of thinning
}
  \item{path}{
	the name of directory where the results are saved
}
  \item{burninPerc}{
	the proportion of burn-in
}
  \item{type}{
	the type of analysis: either "\code{parametric}" or "\code{semi-parametric}"
}
  \item{nChain}{
	the number of chains
}
}


\details{
  The function \code{BayesSurv} fits Bayesian semi-parametric (a mixture of piecewise constant) and parametric (Weibull) models.
  The covariates are incorporated via Cox proportional hazards model:
  \deqn{h(t_{i} | x_{i}) = h_{0}(t_{i})\exp(\beta x_{i}), t_{i}>0,}
  
  where \eqn{h_0} is an unspecified baseline hazard function and \eqn{\beta} is a vector of \eqn{p} log-hazard ratio regression parameters.
  For `semi-parametric analysis', let \eqn{s_{\max}}$ denote the largest observed event time. 
  Then, consider the finite partition of the relevant time axis into \eqn{J + 1} disjoint intervals: \eqn{0<s_1<s_2<...<s_{J+1} = s_{\max}}. 
  For notational convenience, let \eqn{I_j=(s_{j-1}, s_j]} denote the \eqn{j^{th}} partition. 
  For given a partition, \eqn{s = (s_1, \dots, s_{J + 1})}, we assume the log-baseline hazard functions is piecewise constant:
  \deqn{\lambda_{0}(t)=\log h_{0}(t) = \sum_{j=1}^{J + 1} \lambda_{j} I(t\in I_{j})}

   where \eqn{I(\cdot)} is the indicator function and \eqn{s_0 \equiv 0}. Note, this specification is general in that the partitions of the time axes differ across the three hazard functions.
   In our proposed Bayesian framework, our prior choices are:

	\deqn{\pi(\beta) \propto 1,}
	\deqn{\lambda | J, \mu_{\lambda}, \sigma_{\lambda}^2 \sim \mathcal{N}_{J+1}(\mu_{\lambda}1, \sigma_{\lambda}^2\Sigma_{\lambda}),}
	\deqn{J \sim \mathcal{P}(\alpha),}
	\deqn{\pi(s | J) \propto \frac{(2J+1)! \prod_{j=1}^{J+1}(s_j-s_{j-1})}{(s_{J+1})^{(2J+1)}},}
	\deqn{\pi(\mu_{\lambda}) \propto 1,}
	\deqn{\sigma_{\lambda}^{-2} \sim \mathcal{G}(a, b),}
	
	Note that \eqn{J} and \eqn{s} are treated as random and the priors for \eqn{J} and \eqn{s} jointly form a time-homogeneous Poisson process prior for the partition.
	The number of time splits and their positions are therefore updated within our computational scheme using reversible jump MCMC.
		
	For `parametric analysis', we assume the Weibull model for \eqn{h_{0}(t) = \alpha \kappa t^{\alpha-1}}.
   In our Bayesian framework, our prior choices are:	
 	\deqn{\pi(\beta) \propto 1,}
 	\deqn{\pi(\alpha) \sim \mathcal{G}(a, b),}
 	\deqn{\pi(\kappa) \sim \mathcal{G}(c, d),} 	
}

\value{
\code{BayesSurv} returns an object of class \code{BaysSurv}. \cr 
\code{names(object$chain1)} shows the list of posterior samples of model parameters, the number of acceptance in MH algorithm, etc.
}

\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2013). 
Bayesian Semi-parametric Analysis of Semi-competing Risks Data: 
Estimating Readmission Rates among Pancreatic Cancer Patients, submitted.
}
\author{
Kyu Ha Lee and Sebastien Haneuse \cr
}
%\note{
%If the prespecified value of \code{save} is less than that of \code{num.reps}, the results are saved
%as \code{.Rdata} file under the directory \code{working directory/mcmcOutcome}. 
%Trace plots of posterior samples for the parameters are also provided, with the first half taken as burn-in.
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{print.BayesSurv}}, \code{\link{summary.BayesSurv}}, \code{\link{plot.BayesSurv}}
}


\examples{
	
\dontrun{	
	
# loading simulated data set	
data(survData)

n = dim(survData)[1]
p = dim(survData)[2]-2

###############################
# setting hyperparameter values

a	<- 0.5	# prior parameters for 1/sigma^2
b	<- 0.01
alpha	<- 5	# prior parameter for J
c_lam 	<- 1	# prior parameter for MVN-ICAR specification

hyperParams <- c(a, b, alpha, c_lam)

#########################
# setting starting values

s_max	<- max(survData$time[survData$event == 1])

beta	<- rep(1, p)
s	<- c(seq(8, max(survData$time[survData$event == 1]), 8), s_max);
J	<- length(s) - 1
lambda	<- runif(J+1, -3, 0)
sigSq_lam<- var(lambda)
mu_lam	<- mean(lambda)

# chain 1

startValues <- list()
startValues[[1]] <- as.vector(c(beta, J, mu_lam, sigSq_lam, lambda, s))

# chain 2

beta	<- rep(0.2, p)
lambda	<- runif(J+1, -3, -1)

startValues[[2]] <- as.vector(c(beta, J, mu_lam, sigSq_lam, lambda, s))

##################################################
# setting variable values needed for MHG algorithm

C		<- 0.70
delPert		<- 0.5
s_propBI	<- floor(sort(unique(survData$time[survData$event == 1])))
s_propBI	<- s_propBI[s_propBI < s_max]
num_s_propBI	<- length(s_propBI)
J_max 		<- 30 					
time_lambda	<- seq(1:39)
nTime_lambda	<- length(time_lambda)

mcmcParams <- c(C, delPert, num_s_propBI, J_max, s_max, nTime_lambda, s_propBI, time_lambda)

##################################################
# number of chains

numReps	= 2e6
thin	= 20
burninPerc	= 0.5
path1	= "outcome/"
type = "semi-parametric"
nChain	= 2

# fitting Bayesian semi-parametric regression model to univariate survival data	
# In practice, set 'numReps' to a larger number such as 1000000

fitSurv <- BayesSurv(survData, hyperParams, startValues, mcmcParams, numReps, thin, path1, 
burninPerc, type, nChain)
		
print(fitSurv)
summary(fitSurv)		
		
## plot for estimates of baseline hazard function
plot(fitSurv)		
		
}
				

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ univariate analysis }% __ONLY ONE__ keyword per line

