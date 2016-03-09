\name{BayesSurv}
\alias{BayesSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to implement Bayesian parametric and semi-parametric regression analyses for univariate time-to-event data.
}
\description{
Independent/cluster-correlated univariate right-censored survival data can be analyzed using hierarchical models. The prior for the baseline hazard function can be specified by either parametric (Weibull) model or non-parametric mixture of piecewise exponential models (PEM).
}
\usage{
BayesSurv(Y, lin.pred, data, cluster=NULL, model="Weibull", 
          hyperParams, startValues, mcmc, path=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
	a data.frame containing univariate time-to-event outcomes from \code{n} subjects. It is of dimension \eqn{n\times 2}: the columns correspond to \eqn{y}, \eqn{\delta}.
}
  \item{lin.pred}{
	a formula object that corresponds to \eqn{h}.
}
  \item{data}{
	a data.frame in which to interpret the variables named in the formula in \code{lin.pred}.
}
  \item{cluster}{
	a vector of cluster information for \code{n} subjects. The cluster membership must be consecutive positive integers, \eqn{1:J}.  
}
  \item{model}{
	a character vector that specifies the type of components in a model. 
	The first element is for the specification of baseline hazard functions: "Weibull" or "PEM". 
	The second element needs to be set only for clustered data and is for the specification of cluster-specific random effects distribution: "Normal" or "DPM".
}
  \item{hyperParams}{
	a list containing lists or vectors for hyperparameter values in hierarchical models. Components include, 
	\code{WB} (a list containing a numeric vector for Weibull hyperparameters: \code{WB.ab}), 
	\code{PEM} (a list containing numeric vectors for PEM hyperparameters: \code{PEM.ab}, \code{PEM.alpha}).
	Models for clustered data require additional components,
	\code{Normal} (a list containing a numeric vector for hyperparameters in a Normal prior: \code{Normal.ab}), 
	\code{DPM} (a list containing numeric vectors for DPM hyperparameters: \code{DPM.ab}, \code{aTau}, \code{bTau}).
	See Details and Examples below.
}
  \item{startValues}{
	a list containing vectors of starting values for model parameters. It can be specified as the object returned by the function \code{\link{initiate.startValues}}.
}
  \item{mcmc}{
	a list containing variables required for MCMC sampling. Components include,
	\code{run} (a list containing numeric values for setting for the overall run: 
	\code{numReps}, total number of scans; \code{thin}, extent of thinning; \code{burninPerc}, the proportion of burn-in).
	\code{storage} (a list containing numeric values for storing posterior samples for cluster-specific random effects: 
	\code{storeV}, a logical value to determine whether all the posterior samples of \eqn{V} are to be stored.)
	\code{tuning} (a list containing numeric values relevant to tuning parameters for specific updates in Metropolis-Hastings-Green (MHG) algorithm: 
	\code{mhProp_V_var}, the variance of proposal density for \eqn{V} in DPM models; 
	\code{mhProp_alpha_var}, the variance of proposal density for \eqn{\alpha} in Weibull models; 
	\code{C}, a numeric value for the proportion that determines the sum of probabilities of choosing the birth and the death moves in PEM models. 
	The value should not exceed 0.8;
	\code{delPert}, the perturbation parameter in the birth update in PEM models. The values must be between 0 and 0.5;
	If \code{rj.scheme}=1, the birth update will draw the proposal time split from \eqn{1:s_{max}}.
	If \code{rj.scheme}=2, the birth update will draw the proposal time split from uniquely ordered failure times in the data. Only required for PEM models;
	\code{K_max}, the maximum number of splits allowed at each iteration in MHG algorithm for PEM models;
	\code{time_lambda} - time points at which the log-hazard function is calculated for \code{\link{plot.Bayes}}, Only required for PEM models).
	See Details and Examples below.
}
  \item{path}{
	the name of directory where the results are saved.
}
}

\details{
  The function \code{BayesSurv} implements Bayesian semi-parametric (piecewise exponential mixture) and parametric (Weibull) models to univariate time-to-event data. Let \eqn{t_{ji}} denote time to event of interest from subject \eqn{i=1,...,n_j} in cluster \eqn{j=1,...,J}. The covariates \eqn{x_{ji}} are incorporated via Cox proportional hazards model:
    \deqn{h(t_{ji} | x_{ji}) = h_{0}(t_{ji})\exp(x_{ji}^{\top}\beta + V_{j}), t_{ji}>0,}
    
    where \eqn{h_0} is an unspecified baseline hazard function and \eqn{\beta} is a vector of \eqn{p} log-hazard ratio regression parameters. \eqn{V_j}'s are cluster-specific random effects.
    For parametric Normal prior specification for a vector of cluster-specific random effects, we assume \eqn{V} arise as i.i.d. draws from a mean 0 Normal distribution with variance \eqn{\sigma^2}. Specifically, the priors can be written as follows:
    \deqn{V_j \sim Normal(0, \sigma^2),}
 	\deqn{\zeta=1/\sigma^2 \sim Gamma(a_{N}, b_{N}).}
    For DPM prior specification for \eqn{V_j}, we consider non-parametric Dirichlet process mixture of Normal distributions: the \eqn{V_j}'s' are draws from a finite mixture of M Normal distributions, each with their own mean and variance, (\eqn{\mu_m}, \eqn{\sigma_m^2}) for \eqn{m=1,...,M}. Let \eqn{m_j\in\{1,...,M\}} denote the specific component to which the \eqn{j}th cluster belongs. Since the class-specific (\eqn{\mu_m}, \eqn{\sigma_m^2}) are not known they are taken to be draws from some distribution, \eqn{G_0}, often referred to as the centering distribution. Furthermore, since the true class memberships are unknown, we denote the probability that the \eqn{j}th cluster belongs to any given class by the vector \eqn{p=(p_1,..., p_M)} whose components add up to 1.0. In the absence of prior knowledge regarding the distribution of class memberships for the \eqn{J} clusters across the \eqn{M} classes, a natural prior for \eqn{p} is the conjugate symmetric \eqn{Dirichlet(\tau/M,...,\tau/M)} distribution; the hyperparameter, \eqn{\tau}, is often referred to as a the precision parameter. The prior can be represented as follows (\eqn{M} goes to infinity):
	\deqn{V_j | m_j \sim Normal(\mu_{m_j}, \sigma_{m_j}^2),}
	\deqn{(\mu_m, \sigma_m^2) \sim G_{0},~~ for ~m=1,...,M,}
	\deqn{m_j | p \sim Discrete(m_j| p_1,...,p_M),}
	\deqn{p \sim Dirichlet(\tau/M,...,\tau/M),}
    where \eqn{G_0} is taken to be a multivariate Normal/inverse-Gamma (NIG) distribution for which the probability density function is the following product:
	\deqn{f_{NIG}(\mu, \sigma^2 | \mu_0, \zeta_0, a_0, b_0) = f_{Normal}(\mu | 0, 1/\zeta_0^2) \times f_{Gamma}(\zeta=1/\sigma^2 | a_0, b_0).}
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
 	 	
	We provide a detailed description of the hierarchical models for cluster-correlated univariate survival data. The models for independent data can be obtained by removing cluster-specific random effects, \eqn{V_j}, and its corresponding prior specification from the description given above. 	
}

\value{
	\code{BayesSurv} returns an object of class \code{Bayes}. \cr 
}

\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2015), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, 64, 2, 253-273.\cr \cr
Lee, K. H., Dominici, F., Schrag, D., and Haneuse, S. (2016),
Hierarchical models for semicompeting risks data with application to quality of end-of-life care for pancreatic cancer, \emph{Journal of the American Statistical Association}, in press. \cr
}
\author{
Kyu Ha Lee and Sebastien Haneuse \cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
\note{
The posterior samples of \eqn{V_g} are saved separately in \code{working directory/path}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{initiate.startValues}}, \code{\link{print.Bayes}}, \code{\link{summary.Bayes}}, \code{\link{plot.Bayes}}
}

\examples{
	
\dontrun{		
# loading a data set	
data(survData)
Y <- survData[,c(1,2)]
cluster <- survData[,3]
lin.pred <- as.formula( ~ cov1 + cov2)

#####################
## Hyperparameters ##
#####################

## Weibull baseline hazard function: alpha1, kappa1
##
WB.ab <- c(0.5, 0.01) # prior parameters for alpha
##
WB.cd <- c(0.5, 0.05) # prior parameters for kappa

## PEM baseline hazard function: 
##
PEM.ab <- c(0.7, 0.7) # prior parameters for 1/sigma^2
##
PEM.alpha <- 10 # prior parameters for K

## Normal cluster-specific random effects
##
Normal.ab 	<- c(0.5, 0.01) 		# prior for zeta

## DPM cluster-specific random effects
##
DPM.ab <- c(0.5, 0.01)
aTau  <- 1.5
bTau  <- 0.0125

##
hyperParams <- list(WB=list(WB.ab=WB.ab, WB.cd=WB.cd),
                    PEM=list(PEM.ab=PEM.ab, PEM.alpha=PEM.alpha),
                    Normal=list(Normal.ab=Normal.ab),
                    DPM=list(DPM.ab=DPM.ab, aTau=aTau, bTau=bTau))
                    
###################
## MCMC SETTINGS ##
###################

## Setting for the overall run
##
numReps    <- 2000
thin       <- 10
burninPerc <- 0.5

## Settings for storage
##
storeV    <- TRUE

## Tuning parameters for specific updates
##
##  - those common to all models
mhProp_V_var     <- 0.05
##
## - those specific to the Weibull specification of the baseline hazard functions
mhProp_alpha_var <- 0.01
##
## - those specific to the PEM specification of the baseline hazard functions
C        <- 0.2
delPert  <- 0.5
rj.scheme <- 1
K_max    <- 50
s_max    <- max(Y$time[Y$event == 1])
time_lambda <- seq(1, s_max, 1)

##
mcmc.WB  <- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
                    storage=list(storeV=storeV),
                    tuning=list(mhProp_alpha_var=mhProp_alpha_var, mhProp_V_var=mhProp_V_var))
##
mcmc.PEM <- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
                    storage=list(storeV=storeV),
                    tuning=list(mhProp_V_var=mhProp_V_var, C=C, delPert=delPert,
                    rj.scheme=rj.scheme, K_max=K_max, time_lambda=time_lambda))

################################################################
## Analysis of Independent Univariate Survival Data ############
################################################################

####################
## WEIBULL ##
####################

##
myModel <- "Weibull"
myPath  <- "Output/01-Results-WB/"

startValues      <- vector("list", 2)
startValues[[1]] <- initiate.startValues(Y, lin.pred, survData, model=myModel)
startValues[[2]] <- initiate.startValues(Y, lin.pred, survData, model=myModel, 
WB.alpha=1.12)
                
##
fit_WB <- BayesSurv(Y, lin.pred, survData, cluster=NULL, model=myModel, 
                  hyperParams, startValues, mcmc.WB, path=myPath)
                  
fit_WB
summ.fit_WB <- summary(fit_WB); names(summ.fit_WB)
summ.fit_WB
plot(fit_WB, tseq=seq(from=0, to=30, by=5))
plot(fit_WB, tseq=seq(from=0, to=30, by=5), plot.est = "BH")
names(fit_WB.plot <- plot(fit_WB, tseq=seq(0, 30, 5), plot=FALSE))                  

#########
## PEM ##
#########
                
##
myModel <- "PEM"
myPath  <- "Output/02-Results-PEM/"

startValues      <- vector("list", 2)
startValues[[1]] <- initiate.startValues(Y, lin.pred, survData, model=myModel)
startValues[[2]] <- initiate.startValues(Y, lin.pred, survData, model=myModel, 
beta=rep(0.1,2))                                
                
##
fit_PEM <- BayesSurv(Y, lin.pred, survData, cluster=NULL, model=myModel,
                   hyperParams, startValues, mcmc.PEM, path=myPath)
                   
fit_PEM
summ.fit_PEM <- summary(fit_PEM); names(summ.fit_PEM)
summ.fit_PEM
plot(fit_PEM)
plot(fit_PEM, plot.est = "BH")
names(fit_PEM.plot <- plot(fit_PEM, plot=FALSE))                   

###############################################################
## Analysis of Correlated Univariate Survival Data ############
###############################################################

####################
## WEIBULL-NORMAL ##
####################

##
myModel <- c("Weibull", "Normal")
myPath  <- "Output/03-Results-WB_Normal/"

startValues      <- vector("list", 2)

startValues[[1]] <- initiate.startValues(Y, lin.pred, survData, model=myModel, cluster)
startValues[[2]] <- initiate.startValues(Y, lin.pred, survData, model=myModel, cluster,
Normal.zeta=0.95)

##
fit_WB_N <- BayesSurv(Y, lin.pred, survData, cluster, model=myModel,
                        hyperParams, startValues, mcmc.WB, path=myPath)
                        
fit_WB_N
summ.fit_WB_N <- summary(fit_WB_N); names(summ.fit_WB_N)
summ.fit_WB_N
plot(fit_WB_N, tseq=seq(from=0, to=30, by=5))
plot(fit_WB_N, tseq=seq(from=0, to=30, by=5), plot.est = "BH")
names(fit_WB_N.plot <- plot(fit_WB_N, tseq=seq(0, 30, 5), plot=FALSE))                        

#################
## WEIBULL-DPM ##
#################

##
myModel <- c("Weibull", "DPM")
myPath  <- "Output/04-Results-WB_DPM/"

startValues      <- vector("list", 2)
startValues[[1]] <- initiate.startValues(Y, lin.pred, survData, model=myModel, cluster)
startValues[[2]] <- initiate.startValues(Y, lin.pred, survData, model=myModel, cluster,
Normal.zeta=0.95)

##
fit_WB_DPM <- BayesSurv(Y, lin.pred, survData, cluster, model=myModel,
                        hyperParams, startValues, mcmc.WB, path=myPath)

fit_WB_DPM
summ.fit_WB_DPM <- summary(fit_WB_DPM); names(summ.fit_WB_DPM)
summ.fit_WB_DPM
plot(fit_WB_DPM, tseq=seq(from=0, to=30, by=5))
plot(fit_WB_DPM, tseq=seq(from=0, to=30, by=5), plot.est = "BH")
names(fit_WB_DPM.plot <- plot(fit_WB_DPM, tseq=seq(0, 30, 5), plot=FALSE))

################
## PEM-NORMAL ##
################

##
myModel <- c("PEM", "Normal")
myPath  <- "Output/05-Results-PEM_Normal/"

startValues      <- vector("list", 2)
startValues[[1]] <- initiate.startValues(Y, lin.pred, survData, model=myModel, cluster)
startValues[[2]] <- initiate.startValues(Y, lin.pred, survData, model=myModel, cluster,
Normal.zeta=0.95)

##
fit_PEM_N <- BayesSurv(Y, lin.pred, survData, cluster, model=myModel,
                            hyperParams, startValues, mcmc.PEM, path=myPath)

fit_PEM_N
summ.fit_PEM_N <- summary(fit_PEM_N); names(summ.fit_PEM_N)
summ.fit_PEM_N
plot(fit_PEM_N)
plot(fit_PEM_N, plot.est = "BH")
names(fit_PEM_N.plot <- plot(fit_PEM_N, plot=FALSE)) 

#############
## PEM-DPM ##
#############

##
myModel <- c("PEM", "DPM")
myPath  <- "Output/06-Results-PEM_DPM/"

startValues      <- vector("list", 2)
startValues[[1]] <- initiate.startValues(Y, lin.pred, survData, model=myModel, cluster)
startValues[[2]] <- initiate.startValues(Y, lin.pred, survData, model=myModel, cluster,
Normal.zeta=0.95)

##
fit_PEM_DPM <- BayesSurv(Y, lin.pred, survData, cluster, model=myModel,
                        hyperParams, startValues, mcmc.PEM, path=myPath)
                        
fit_PEM_DPM
summ.fit_PEM_DPM <- summary(fit_PEM_DPM); names(summ.fit_PEM_DPM)
summ.fit_PEM_DPM
plot(fit_PEM_DPM)
plot(fit_PEM_DPM, plot.est = "BH")
names(fit_PEM_DPM.plot <- plot(fit_PEM_DPM, plot=FALSE))                          
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ univariate analysis }% __ONLY ONE__ keyword per line

