\name{BayesSurv}
\alias{BayesSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to perform Bayesian parametric and semi-parametric regression analysis for univariate time-to-event data.
}
\description{
Univariate regression analysis of time-to-event data is performed using Cox proportional hazards model in a Bayesian framework. The prior of baseline hazard function can be specified by either parametric Weibull or non-parametric mixture of piecewise exponential models (PEM).
}
\usage{
BayesSurv(survData, hyperParams, startValues, mcmcParams, numReps, thin, path, 
burninPerc = 0.5, hz.type = "Weibull", nChain = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
	The data frame containing univariate survival data with covariate matrix from \code{n} subjects. It is of dimension \eqn{n\times(2+p)}, where \eqn{p} is the number of covariates. The first two columns correspond to \eqn{y}, \eqn{\delta}.
}
  \item{hyperParams}{
	a vector containing hyperparameter values. For PEM model (\code{hz.type}="PEM"), it is of length 4; c(\eqn{a}, \eqn{b}, \eqn{\alpha}, \eqn{c_{\lambda}}). For Weibull model (\code{hz.type}="Weibull"), it is of length 4; c(\eqn{a}, \eqn{b}, \eqn{c}, \eqn{d}).
}
  \item{startValues}{
  	a list containing vectors of starting values of parameters. The length of list must be equal to the number of chains (\code{nChain}). For PEM model (\code{hz.type}="PEM"), each of list components is of length \eqn{p+2(J+1)+3}; c(\eqn{\beta}, \eqn{J}, \eqn{\mu_{\lambda}}, \eqn{\sigma_{\lambda}^2}, \eqn{\lambda}, \eqn{s}). For Weibull model (\code{hz.type}="Weibull"), each of list components is of length \eqn{p+2}; c(\eqn{\beta}, \eqn{\alpha}, \eqn{\kappa}).
}

  \item{mcmcParams}{
	a vector containing variables required for Metropolis-Hastings-Green (MHG) algorithm. \cr
    For PEM model (\code{hz.type}="PEM"), it is of length \eqn{(6+num\_s\_propBI+nTime\_lambda}); c(\eqn{C}, delPert, num_s_propBI, \eqn{J_{max}}, \eqn{s_{max}}, nTime_lambda, s_propBI, time_lambda). The perturbation random variable is drawn from Uniform(0.5-delPert, 0.5+delPert) in the birth/death move for \eqn{h_{0}}. nTime_lambda is \code{length}(time_lambda), where time_lambda is a set of time points where the posterior estimates of \eqn{log(h_{0})} is calculated. s_propBI is a set of candidates for proposal values in birth move of \eqn{h_{0}}. \cr
    For Weibull model (\code{hz.type}="Weibull"), it is the value for the variance of proposal density in Metropolis-Hastings step for \eqn{\alpha}.
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
  \item{hz.type}{
	prior specification for baseline hazard functions : either "\code{Weibull}" or "\code{PEM}"
}
  \item{nChain}{
	the number of chains
}
}


\details{
  The function \code{BayesSurv} implements Bayesian semi-parametric (a mixture of piecewise constant) and parametric (Weibull) models. Let \eqn{t_{i}} denote time to event of interest from subject \eqn{i=1,...,n}. The covariates \eqn{x_{i}} are incorporated via Cox proportional hazards model:
  \deqn{h(t_{i} | x_{i}) = h_{0}(t_{i})\exp(\beta x_{i}), t_{i}>0,}
  where \eqn{h_0} is an unspecified baseline hazard function and \eqn{\beta} is a vector of \eqn{p} log-hazard ratio regression parameters.
  For "semi-parametric analysis", let \eqn{s_{\max}} denote the largest observed event time.
  Then, consider the finite partition of the relevant time axis into \eqn{J + 1} disjoint intervals: \eqn{0<s_1<s_2<...<s_{J+1} = s_{\max}}. 
  For notational convenience, let \eqn{I_j=(s_{j-1}, s_j]} denote the \eqn{j^{th}} partition. 
  For given a partition, \eqn{s = (s_1, \dots, s_{J + 1})}, we assume the log-baseline hazard functions is piecewise constant:
  \deqn{\lambda_{0}(t)=\log h_{0}(t) = \sum_{j=1}^{J + 1} \lambda_{j} I(t\in I_{j}),}

   where \eqn{I(\cdot)} is the indicator function and \eqn{s_0 \equiv 0}. In our proposed Bayesian framework, our prior choices are:

	\deqn{\pi(\beta) \propto 1,}
	\deqn{\lambda | J, \mu_{\lambda}, \sigma_{\lambda}^2 \sim \mathcal{N}_{J+1}(\mu_{\lambda}1, \sigma_{\lambda}^2\Sigma_{\lambda}),}
	\deqn{J \sim \mathcal{P}(\alpha),}
	\deqn{\pi(s | J) \propto \frac{(2J+1)! \prod_{j=1}^{J+1}(s_j-s_{j-1})}{(s_{J+1})^{(2J+1)}},}
	\deqn{\pi(\mu_{\lambda}) \propto 1,}
	\deqn{\sigma_{\lambda}^{-2} \sim \mathcal{G}(a, b).}
	
	Note that \eqn{J} and \eqn{s} are treated as random and the priors for \eqn{J} and \eqn{s} jointly form a time-homogeneous Poisson process prior for the partition.
	The number of time splits and their positions are therefore updated within our computational scheme using reversible jump MCMC.
		
	For "parametric analysis", we assume the Weibull model for \eqn{h_{0}(t) = \alpha \kappa t^{\alpha-1}}.
   In our Bayesian framework, our prior choices are:	
 	\deqn{\pi(\beta) \propto 1,}
 	\deqn{\pi(\alpha) \sim \mathcal{G}(a, b),}
 	\deqn{\pi(\kappa) \sim \mathcal{G}(c, d).} 	
}

\value{
\code{BayesSurv} returns an object of class \code{BayesSurv}. \cr 
\code{names(object$chain1)} shows the list of posterior samples of model parameters, the number of acceptance in MH algorithm, etc.
}

\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2014), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, \emph{in press}.
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
	
# loading a data set	
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
hz.type = "PEM"
nChain	= 2

# fitting Bayesian semi-parametric regression model to univariate survival data	

fitSurv <- BayesSurv(survData, hyperParams, startValues, mcmcParams, numReps, thin, path1, 
burninPerc, hz.type, nChain)
		
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

