\name{BayesSurvcor}
\alias{BayesSurvcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to perform Bayesian parametric and semi-parametric regression analysis for cluster-correlated univariate time-to-event data.
}
\description{
Univariate regression analysis of clustered time-to-event data is performed using Cox proportional hazards frailty model in a Bayesian framework. The prior of baseline hazard function can be specified by either parametric Weibull or non-parametric mixture of piecewise exponential models (PEM). Either parametric normal or non-parametric Dirichlet process mixture of normals (DPM) can be adopted for the prior of cluster-specific random effects distribution.
}
\usage{
BayesSurvcor(survData, hyperParams, startValues, mcmcParams, numReps, thin, path, 
burninPerc = 0.5, hz.type = "Weibull", re.type = "Normal", storeV = TRUE, nChain = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
	The data frame containing univariate survival data with covariate matrix from \code{n} subjects. It is of dimension \eqn{n\times(3+p)}, where \eqn{p} is the number of covariates. The first two columns correspond to \eqn{y}, \eqn{\delta}, cluster membership. Note that the cluster membership must set to consecutive positive integers; 1:J.
}
  \item{hyperParams}{
	a vector containing hyperparameter values for hierarchical model.\cr
    For Weibull-Normal, (\code{hz.type}="Weibull", \code{re.type}="Normal"), it is of length 6; c(\eqn{a}, \eqn{b}, \eqn{c}, \eqn{d}, \eqn{\rho_1}, \eqn{\rho_2}).\cr
    For Weibull-DPM, (\code{hz.type}="Weibull", \code{re.type}="DPM"), it is of length 10; c(\eqn{a}, \eqn{b}, \eqn{c}, \eqn{d}, \eqn{\mu_0}, \eqn{\zeta_0}, \eqn{a_0}, \eqn{b_0}, \eqn{a_{\tau}}, \eqn{b_{\tau}}). \cr
    For PEM-Normal, (\code{hz.type}="PEM", \code{re.type}="Normal"), it is of length 6; c(\eqn{a}, \eqn{b}, \eqn{\alpha}, \eqn{c_{\lambda}}, \eqn{\rho_1}, \eqn{\rho_2}). \cr
    For PEM-DPM, (\code{hz.type}="PEM", \code{re.type}="DPM"), it is of length 10; c(\eqn{a}, \eqn{b}, \eqn{\alpha}, \eqn{c_{\lambda}}, \eqn{\mu_0}, \eqn{\zeta_0}, \eqn{a_0}, \eqn{b_0}, \eqn{a_{\tau}}, \eqn{b_{\tau}}).
}
  \item{startValues}{
  	a list containing vectors of starting values of parameters. The length of list must be equal to the number of chains (\code{nChain}). \cr
    For Weibull-Normal, (\code{hz.type}="Weibull", \code{re.type}="Normal"), each of list components is of length \eqn{p+J+3}; c(\eqn{\beta}, \eqn{\alpha}, \eqn{\kappa}, \eqn{V}, \eqn{\zeta}).\cr
    For Weibull-DPM, (\code{hz.type}="Weibull", \code{re.type}="DPM"), each of list components is of length \eqn{p+2J+3}; c(\eqn{\beta}, \eqn{\alpha}, \eqn{\kappa}, \eqn{V}, \eqn{c_v}, \eqn{\tau}). \cr
    For PEM-Normal, (\code{hz.type}="PEM", \code{re.type}="Normal"), each of list components is of length \eqn{p+2(K+1)+J+4};       c(\eqn{\beta}, \eqn{K}, \eqn{\mu_{\lambda}}, \eqn{\sigma_{\lambda}^2}, \eqn{\lambda}, \eqn{s}, \eqn{V}, \eqn{\zeta}).\cr
    For PEM-DPM, (\code{hz.type}="PEM", \code{re.type}="DPM"), each of list components is of length \eqn{p+2(K+1)+3J+4}; c(\eqn{\beta}, \eqn{K}, \eqn{\mu_{\lambda}}, \eqn{\sigma_{\lambda}^2}, \eqn{\lambda}, \eqn{s}, \eqn{V}, \eqn{c_v}, \eqn{\tau}). 
}
  \item{mcmcParams}{
	a vector containing variables required for Metropolis-Hastings-Green (MHG) algorithm. \cr
      For Weibull models (Weibull-Normal and Weibull-DPM: \code{hz.type}="Weibull"), it is of length 2; c(mhProp_alpha_var, mhProp_V_var). mhProp_alpha_var and mhProp_V_var are the variance of proposal density in Metropolis-Hastings step for \eqn{\alpha} and \eqn{V}, respectively. \cr
      For PEM models (PEM-Normal and PEM-DPM: \code{hz.type}="PEM"), it is of length \eqn{(7+num\_s\_propBI+nTime\_lambda)}; c(\eqn{C}, delPert, num_s_propBI, \eqn{K_{max}}, \eqn{s_{max}}, nTime_lambda, s_propBI, time_lambda, mhProp_V_var). The perturbation random variable is drawn from Uniform(0.5-delPert, 0.5+delPert) in the birth/death move for \eqn{h_{0}}. nTime_lambda is \code{length}(time_lambda), where time_lambda is a set of time points where the posterior estimates of \eqn{log(h_{0})} is calculated. s_propBI is a set of candidates for proposal values in birth move of \eqn{h_{0}}. mhProp_V_var is the variance of proposal density in Metropolis-Hastings step for \eqn{V}.
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
  \item{re.type}{
	prior specification for cluster-specific random effects distribution: either "\code{Normal}" or "\code{DPM}"
}
  \item{storeV}{
	logical value to determine whether all the posterior samples of \eqn{V} are to be stored.
}
  \item{nChain}{
	the number of chains
}
}


\details{
    The function \code{BayesSurvcor} fits Bayesian semi-parametric (piecewise exponential mixture) and parametric (Weibull) models to cluster-correlated univariate time-to-event data. Let \eqn{t_{ji}} denote time to event of interest from subject \eqn{i=1,...,n_j} in cluster \eqn{j=1,...,J}. The covariates \eqn{x_{ji}} are incorporated via Cox proportional hazards model:
    \deqn{h(t_{ji} | x_{ji}) = h_{0}(t_{ji})\exp(\beta x_{ji} + V_{j}), t_{ji}>0,}
    
    where \eqn{h_0} is an unspecified baseline hazard function and \eqn{\beta} is a vector of \eqn{p} log-hazard ratio regression parameters. \eqn{V_j}s are cluster-specific random effects.
    For parametric Normal prior specification for a vector of cluster-specific random effects, we assume \eqn{V} arise as i.i.d. draws from a mean 0 Normal distribution with variance \eqn{\sigma^2}. Specifically, the priors can be written as follows:
    \deqn{V_j \sim Normal(0, \sigma^2),}
 	\deqn{\sigma^2 \sim inverse-Gamma(\rho_1, \rho_2).}
    For DPM prior specification for \eqn{V_j}, we consider non-parametric Dirichlet process mixture of Normal distributions: the \eqn{V_j}'s' are draws from a finite mixture of M Normal distributions, each with their own mean and variance, (\eqn{\mu_m}, \eqn{\sigma_m^2}) for \eqn{m=1,...,M}. Let \eqn{m_j\in\{1,...,M\}} denote the specific component to which the \eqn{j}th cluster belongs. Since the class-specific (\eqn{\mu_m}, \eqn{\sigma_m^2}) are not known they are taken to be draws from some distribution, \eqn{G_0}, often referred to as the centering distribution. Furthermore, since the true class memberships are unknown, we denote the probability that the \eqn{j}th cluster belongs to any given class by the vector \eqn{p=(p_1,..., p_M)} whose components add up to 1.0. In the absence of prior knowledge regarding the distribution of class memberships for the \eqn{J} clusters across the \eqn{M} classes, a natural prior for \eqn{p} is the conjugate symmetric \eqn{Dirichlet(\tau/M,...,\tau/M)} distribution; the hyperparameter, \eqn{\tau}, is often referred to as a the precision parameter. The prior can be represented as follows (\eqn{M} goes to infinity):
	\deqn{V_j | m_j \sim Normal(\mu_{m_j}, \sigma_{m_j}^2),}
	\deqn{(\mu_m, \sigma_m^2) \sim G_{0},~~ for ~m=1,...,M,}
	\deqn{m_j | p \sim Discrete(m_j| p_1,...,p_M),}
	\deqn{p \sim Dirichlet(\tau/M,...,\tau/M),}
    where \eqn{G_0} is taken to be a multivariate Normal/inverse-Gamma (NIG) distribution for which the probability density function is the following product:
	\deqn{f_{NIG}(\mu, \sigma^2 | \mu_0, \zeta_0, a_0, b_0) = f_{Normal}(\mu | 0, 1/\zeta_0^2) \times f_{inv-Gamma}(\sigma^2 | a_0, b_0).}
    In addition, we use \eqn{Gamma(a_{\tau}, b_{\tau})} as the hyperprior for \eqn{\tau}.
    
  For non-parametric prior specification (PEM) for baseline hazard function, let \eqn{s_{\max}} denote the largest observed event time. Then, consider the finite partition of the relevant time axis into \eqn{K + 1} disjoint intervals: \eqn{0<s_1<s_2<...<s_{K+1} = s_{\max}}. For notational convenience, let \eqn{I_k=(s_{k-1}, s_k]} denote the \eqn{k^{th}} partition. For given a partition, \eqn{s = (s_1, \dots, s_{K + 1})}, we assume the log-baseline hazard functions is piecewise constant:
  \deqn{\lambda_{0}(t)=\log h_{0}(t) = \sum_{k=1}^{K + 1} \lambda_{k} I(t\in I_{k}),}

   where \eqn{I(\cdot)} is the indicator function and \eqn{s_0 \equiv 0}. In our proposed Bayesian framework, our prior choices are:

	\deqn{\pi(\beta) \propto 1,}
	\deqn{\lambda | K, \mu_{\lambda}, \sigma_{\lambda}^2 \sim MVN_{K+1}(\mu_{\lambda}1, \sigma_{\lambda}^2\Sigma_{\lambda}),}
	\deqn{K \sim Poisson(\alpha),}
	\deqn{\pi(s | K) \propto \frac{(2K+1)! \prod_{k=1}^{K+1}(s_k-s_{k-1})}{(s_{K+1})^{(2K+1)}},}
	\deqn{\pi(\mu_{\lambda}) \propto 1,}
	\deqn{\sigma_{\lambda}^{-2} \sim Gamma(a, b).}
	
	Note that \eqn{K} and \eqn{s} are treated as random and the priors for \eqn{K} and \eqn{s} jointly form a time-homogeneous Poisson process prior for the partition. The number of time splits and their positions are therefore updated within our computational scheme using reversible jump MCMC.
		
	For parametric Weibull prior specification for baseline hazard function, \eqn{h_{0}(t) = \alpha \kappa t^{\alpha-1}}.
   In our Bayesian framework, our prior choices are:	
 	\deqn{\pi(\beta) \propto 1,}
 	\deqn{\pi(\alpha) \sim Gamma(a, b),}
 	\deqn{\pi(\kappa) \sim Gamma(c, d).}
}

\value{
\code{BayesSurvcor} returns an object of class \code{BayesSurvcor}. \cr 
\code{names(object$chain1)} shows the list of posterior samples of model parameters, the number of acceptance in MH algorithm, etc.
}

\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2014), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, \emph{in press}.\cr \cr
Lee, K. H., Dominici, F., Schrag, D., and Haneuse, S.,
Hierarchical models for cluster-correlated semicompeting risks data with application to quality of end-of-life care for pancreatic cancer, \emph{arXiv:1502.00526}; submitted. \cr
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
data(scrCorData)

survData <- scrCorData[,3:8]
colnames(survData)[c(1,2)] <- c("time", "event")

n = dim(survData)[1]
p = dim(survData)[2]-3
J <- length(unique(survData$cluster))

###############################
# setting hyperparameter values


a		<- 0.7	# prior parameters for 1/sigma_1^2
b		<- 0.7
alpha	<- 10	# prior parameter for K
c_lam 	<- 1	# prior parameter for MVN-ICAR specification

rho1 	<- 0.5 # prior for zeta
rho2 	<- 0.01

hyperParams <- c(a, b, alpha, c_lam, rho1, rho2)



#########################
# setting starting values

s_max		<- max(survData$time[survData$event == 1])

grid <- s_max/50

beta		<- rep(0.1, p)
s			<- unique(sort(c(sample(1:s_max, 10), s_max)))
K			<- length(s) - 1
lambda		<- runif(K+1, -4, -3)
sigSq_lam	<- var(lambda)
mu_lam		<- mean(lambda)
V	<- rep(0, J)
zeta <- 0.5



# chain 1

startValues <- list()
startValues[[1]] <- as.vector(c(beta, K, mu_lam, sigSq_lam, lambda, s, V, zeta))

# chain 2

beta		<- rep(0.2, p)
lambda		<- runif(K+1, -4.1, -3)
sigSq_lam	<- var(lambda)
mu_lam		<- mean(lambda)
V	<- rep(0.1, J)
zeta <- 0.7

startValues[[2]] <- as.vector(c(beta, K, mu_lam, sigSq_lam, lambda, s, V, zeta))

##################################################
# setting variable values needed for MH algorithm

C				<- 0.20
delPert		<- 0.5
s_propBI	<- floor(sort(unique(survData$time[survData$event == 1])))
s_propBI	<- s_propBI[s_propBI < s_max]
num_s_propBI	<- length(s_propBI)
K_max 			<- 50
time_lambda	<- seq(grid, s_max, grid)
nTime_lambda	<- length(time_lambda)
mhProp_V_var <- 0.05


mcmcParams <- c(C, delPert, num_s_propBI, K_max, s_max, nTime_lambda, s_propBI, time_lambda, 
mhProp_V_var)





##################################################
# number of chains

numReps	= 2e6
thin	= 20
burninPerc	= 0.5
path1	= "outcome/"
hz.type = "PEM"
re.type     = "Normal"
storeV      = TRUE
nChain	= 2

# fitting Bayesian semi-parametric regression model to univariate survival data	

fit <- BayesSurvcor(survData, hyperParams, startValues, mcmcParams, numReps, thin, path1, 
burninPerc, hz.type, re.type, storeV, nChain)
		
print(fit)
summary(fit)		
		
## plot for estimates of baseline hazard function
plot(fit)				
		
}
				

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ univariate analysis }% __ONLY ONE__ keyword per line

