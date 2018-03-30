\name{BayesID_HReg}
\alias{BayesID_HReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to implement Bayesian parametric and semi-parametric analyses for semi-competing risks data in the context of hazard regression (HReg) models.
}
\description{
Independent/cluster-correlated semi-competing risks data can be analyzed using HReg models that have a hierarchical structure. The priors for baseline hazard functions can be specified by either parametric (Weibull) model or non-parametric mixture of piecewise exponential models (PEM). The option to choose between parametric multivariate normal (MVN) and non-parametric Dirichlet process mixture of multivariate normals (DPM) is available for the prior of cluster-specific random effects distribution. The conditional hazard function for time to the terminal event given time to non-terminal event can be modeled based on either Markov (it does not depend on the timing of the non-terminal event) or semi-Markov assumption (it does depend on the timing).
}
\usage{
BayesID_HReg(Formula, data, id=NULL, model=c("semi-Markov", "Weibull"),
        hyperParams, startValues, mcmcParams, path=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Formula}{
    a \code{Formula} object of the form \eqn{y_1}+\eqn{\delta_1} | \eqn{y_2}+\eqn{\delta_2} ~ \eqn{x_1} | \eqn{x_2} | \eqn{x_3}. See Details and Examples below.
}
  \item{data}{
	a data.frame in which to interpret the variables named in \code{Formula}.
}
  \item{id}{
	a vector of cluster information for \code{n} subjects. The cluster membership must be consecutive positive integers, \eqn{1:J}.  
}
  \item{model}{
	a character vector that specifies the type of components in a model. 
	The first element is for the assumption on \eqn{h_3}: "semi-Markov" or "Markov". 
	The second element is for the specification of baseline hazard functions: "Weibull" or "PEM". 
	The third element needs to be set only for clustered semi-competing risks data and is for the specification of cluster-specific random effects distribution: "MVN" or "DPM".
}
  \item{hyperParams}{
	a list containing lists or vectors for hyperparameter values in hierarchical models. Components include, 
	\code{theta} (a numeric vector for hyperparameter in the prior of subject-specific frailty variance component), 
	\code{WB} (a list containing numeric vectors for Weibull hyperparameters: \code{WB.ab1}, \code{WB.ab2}, \code{WB.ab3}, \code{WB.cd1}, \code{WB.cd2}, \code{WB.cd3}), 
	\code{PEM} (a list containing numeric vectors for PEM hyperparameters: \code{PEM.ab1}, \code{PEM.ab2}, \code{PEM.ab3}, \code{PEM.alpha1}, \code{PEM.alpha2}, \code{PEM.alpha3}).
	Models for clustered data require additional components,
	\code{MVN} (a list containing numeric vectors for MVN hyperparameters: \code{Psi_v}, \code{rho_v}), 
	\code{DPM} (a list containing numeric vectors for DPM hyperparameters: \code{Psi0}, \code{rho0}, \code{aTau}, \code{bTau}).
	See Details and Examples below.
}
  \item{startValues}{
  	a list containing vectors of starting values for model parameters. It can be specified as the object returned by the function \code{\link{initiate.startValues_HReg}}.
}
  \item{mcmcParams}{
	a list containing variables required for MCMC sampling. Components include,
	\code{run} (a list containing numeric values for setting for the overall run: 
	\code{numReps}, total number of scans; \code{thin}, extent of thinning; \code{burninPerc}, the proportion of burn-in).
	\code{storage} (a list containing numeric values for storing posterior samples for subject- and cluster-specific random effects: 
	\code{nGam_save}, the number of \eqn{\gamma} to be stored; 
	\code{storeV}, a vector of three logical values to determine whether all the posterior samples of \eqn{V_1}, \eqn{V_2}, \eqn{V_3} are to be stored).
	\code{tuning} (a list containing numeric values relevant to tuning parameters for specific updates in Metropolis-Hastings-Green (MHG) algorithm: 
	\code{mhProp_theta_var}, the variance of proposal density for \eqn{\theta}; \code{mhProp_Vg_var}, the variance of proposal density for \eqn{V_g} in DPM models; 
	\code{mhProp_alphag_var}, the variance of proposal density for \eqn{\alpha_g} in Weibull models; 
	\code{Cg}, a vector of three proportions that determine the sum of probabilities of choosing the birth and the death moves in PEM models. 
	The sum of the three elements should not exceed 0.6;
	\code{delPertg}, the perturbation parameters in the birth update in PEM models. The values must be between 0 and 0.5;
	If \code{rj.scheme}=1, the birth update will draw the proposal time split from \eqn{1:s_{max}}.
	If \code{rj.scheme}=2, the birth update will draw the proposal time split from uniquely ordered failure times in the data. Only required for PEM models;
	\code{Kg_max}, the maximum number of splits allowed at each iteration in MHG algorithm for PEM models;
	\code{time_lambda1}, \code{time_lambda2}, \code{time_lambda3} - time points at which the log-hazard functions are calculated for \code{\link{predict.Bayes_HReg}}, Only required for PEM models).
	See Details and Examples below.
}
  \item{path}{
	the name of directory where the results are saved.
}
}
\details{
   We view the semi-competing risks data as arising from an underlying illness-death model system in which individuals may undergo one or more of three transitions: 1) from some initial condition to non-terminal event, 2) from some initial condition to terminal event, 3) from non-terminal event to terminal event. Let \eqn{t_{ji1}}, \eqn{t_{ji2}} denote time to non-terminal and terminal event from subject \eqn{i=1,...,n_j} in cluster \eqn{j=1,...,J}. The system of transitions is modeled via the specification of three hazard functions:
    \deqn{h_1(t_{ji1} | \gamma_{ji}, x_{ji1}, V_{j1}) = \gamma_{ji} h_{01}(t_{ji1})\exp(x_{ji1}^{\top}\beta_1 +V_{j1}), t_{ji1}>0,}
    \deqn{h_2(t_{ji2} | \gamma_{ji}, x_{ji2}, V_{j2}) = \gamma_{ji} h_{02}(t_{ji2})\exp(x_{ji2}^{\top}\beta_2 +V_{j2}), t_{ji2}>0,}
    \deqn{h_3(t_{ji2} | t_{ji1}, \gamma_{ji}, x_{ji3}, V_{j3}) = \gamma_{ji} h_{03}(t_{ji2})\exp(x_{ji3}^{\top}\beta_3 +V_{j3}), 0<t_{ji1}<t_{ji2},}
    
    where \eqn{\gamma_{ji}} is a subject-specific frailty and \eqn{V_j}=(\eqn{V_{j1}}, \eqn{V_{j2}}, \eqn{V_{j3}}) is a vector of cluster-specific random effects (each specific to one of the three possible transitions), taken to be independent of \eqn{x_{ji1}}, \eqn{x_{ji2}}, and \eqn{x_{ji3}}.
    For \eqn{g \in \{1,2,3\}}, \eqn{h_{0g}} is an unspecified baseline hazard function and \eqn{\beta_g} is a vector of \eqn{p_g} log-hazard ratio regression parameters.
    The \eqn{h_{03}} is assumed to be Markov with respect to \eqn{t_1}. We refer to the model specified by three conditional hazard functions as the Markov model.
    An alternative specification is to model the risk of terminal event following non-terminal event as a function of the sojourn time. Specifically, retaining \eqn{h_1} and \eqn{h_2} as above,
    we consider modeling \eqn{h_3} as follows:
    \deqn{h_3(t_{ji2} | t_{ji1}, \gamma_{ji}, x_{ji3}, V_{j3}) = \gamma_{ji} h_{03}(t_{ji2}-t_{ji1})\exp(x_{ji3}^{\top}\beta_3 +V_{j3}), 0<t_{ji1}<t_{ji2}.}
    We refer to this alternative model as the semi-Markov model. \cr
    For parametric MVN prior specification for a vector of cluster-specific random effects, we assume \eqn{V_j} arise as i.i.d. draws from a mean 0 MVN distribution with variance-covariance matrix \eqn{\Sigma_V}. The diagonal elements of the \eqn{3\times 3} matrix \eqn{\Sigma_V} characterize variation across clusters in risk for non-terminal, terminal and terminal following non-terminal event, respectively, which is not explained by covariates included in the linear predictors. Specifically, the priors can be written as follows:
    \deqn{V_j \sim MVN(0, \Sigma_V),}
 	\deqn{\Sigma_V \sim inverse-Wishart(\Psi_v, \rho_v).}
    For DPM prior specification for \eqn{V_j}, we consider non-parametric Dirichlet process mixture of MVN distributions: the \eqn{V_j}'s are draws from a finite mixture of M multivariate Normal distributions, each with their own mean vector and variance-covariance matrix, (\eqn{\mu_m}, \eqn{\Sigma_m}) for \eqn{m=1,...,M}. Let \eqn{m_j\in\{1,...,M\}} denote the specific component to which the \eqn{j}th cluster belongs. Since the class-specific (\eqn{\mu_m}, \eqn{\Sigma_m}) are unknown they are taken to be draws from some distribution, \eqn{G_0}, often referred to as the centering distribution. Furthermore, since the true class memberships are not known, we denote the probability that the \eqn{j}th cluster belongs to any given class by the vector \eqn{p=(p_1,..., p_M)} whose components add up to 1.0. In the absence of prior knowledge regarding the distribution of class memberships for the \eqn{J} clusters across the \eqn{M} classes, a natural prior for \eqn{p} is the conjugate symmetric \eqn{Dirichlet(\tau/M,...,\tau/M)} distribution; the hyperparameter, \eqn{\tau}, is often referred to as a the precision parameter. The prior can be represented as follows (\eqn{M} goes to infinity):
	\deqn{V_j | m_j \sim MVN(\mu_{m_j}, \Sigma_{m_j}),}
	\deqn{(\mu_m, \Sigma_m) \sim G_{0},~~ for ~m=1,...,M,}
	\deqn{m_j | p \sim Discrete(m_j| p_1,...,p_M),}
	\deqn{p \sim Dirichlet(\tau/M,...,\tau/M),}
    where \eqn{G_0} is taken to be a multivariate Normal/inverse-Wishart (NIW) distribution for which the probability density function is the following product:
	\deqn{f_{NIW}(\mu, \Sigma | \Psi_0, \rho_0) = f_{MVN}(\mu | 0, \Sigma) \times f_{inv-Wishart}(\Sigma | \Psi_0, \rho_0).}
    We consider \eqn{Gamma(a_{\tau}, b_{\tau})} as the prior for concentration parameter \eqn{\tau}.
    
    For non-parametric PEM prior specification for baseline hazard functions, let \eqn{s_{g,\max}} denote the largest observed event time for each transition \eqn{g \in \{1,2,3\}}.
    Then, consider the finite partition of the relevant time axis into \eqn{K_{g} + 1} disjoint intervals: \eqn{0<s_{g,1}<s_{g,2}<...<s_{g, K_g+1} = s_{g, \max}}. For notational convenience, let \eqn{I_{g,k}=(s_{g, k-1}, s_{g, k}]} denote the \eqn{k^{th}} partition. For given a partition, \eqn{s_g = (s_{g,1}, \dots, s_{g, K_g + 1})}, we assume the log-baseline hazard functions is piecewise constant:
    \deqn{\lambda_{0g}(t)=\log h_{0g}(t) = \sum_{k=1}^{K_g + 1} \lambda_{g,k} I(t\in I_{g,k}),}
    where \eqn{I(\cdot)} is the indicator function and \eqn{s_{g,0} \equiv 0}. Note, this specification is general in that the partitions of the time axes differ across the three hazard functions. our prior choices are, for \eqn{g\in\{1,2,3\}}:
	\deqn{\lambda_g | K_g, \mu_{\lambda_g}, \sigma_{\lambda_g}^2 \sim MVN_{K_g+1}(\mu_{\lambda_g}1, \sigma_{\lambda_g}^2\Sigma_{\lambda_g}),}
	\deqn{K_g \sim Poisson(\alpha_g),}
	\deqn{\pi(s_g | K_g) \propto \frac{(2K_g+1)! \prod_{k=1}^{K_g+1}(s_{g,k}-s_{g,k-1})}{(s_{g,K_g+1})^{(2K_g+1)}},}
	\deqn{\pi(\mu_{\lambda_g}) \propto 1,}
	\deqn{\sigma_{\lambda_g}^{-2} \sim Gamma(a_g, b_g).}
	Note that \eqn{K_g} and \eqn{s_g} are treated as random and the priors for \eqn{K_g} and \eqn{s_g} jointly form a time-homogeneous Poisson process prior for the partition. The number of time splits and their positions are therefore updated within our computational scheme using reversible jump MCMC.
	
	For parametric Weibull prior specification for baseline hazard functions, \eqn{h_{0g}(t) = \alpha_g \kappa_g t^{\alpha_g-1}}. In our Bayesian framework, our prior choices are, for \eqn{g\in\{1,2,3\}}:
    \deqn{\pi(\alpha_g) \sim Gamma(a_g, b_g),}
 	\deqn{\pi(\kappa_g) \sim Gamma(c_g, d_g).}
    Our prior choice for remaining model parameters in all of four models (Weibull-MVN, Weibull-DPM, PEM-MVN, PEM-DPM) is given as follows:
    \deqn{\pi(\beta_g) \propto 1,}
	\deqn{\gamma_{ji}|\theta \sim Gamma(\theta^{-1}, \theta^{-1}),}
	\deqn{\theta^{-1} \sim Gamma(\psi, \omega).}  
	
	We provide a detailed description of the hierarchical models for cluster-correlated semi-competing risks data. The models for independent semi-competing risks data can be obtained by removing cluster-specific random effects, \eqn{V_j}, and its corresponding prior	specification from the description given above.
}


\value{
\code{BayesID_HReg} returns an object of class \code{Bayes_HReg}. \cr
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2015), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, 64, 2, 253-273.\cr \cr
Lee, K. H., Dominici, F., Schrag, D., and Haneuse, S. (2016),
Hierarchical models for semicompeting risks data with application to quality of end-of-life care for pancreatic cancer, \emph{Journal of the American Statistical Association}, 111, 515, 1075-1095. \cr \cr
Alvares, D., Haneuse, S., Lee, C., Lee, K. H. (2018+),
SemiCompRisks: an R package for independent and cluster-correlated analyses of semi-competing risks data, \emph{submitted}, arXiv:1801.03567. \cr
}
\author{
Kyu Ha Lee and Sebastien Haneuse\cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
\note{
The posterior samples of \eqn{\gamma} and \eqn{V_g} are saved separately in \code{working directory/path}. 
For a dataset with large \eqn{n}, \code{nGam_save} should be carefully specified considering the system memory and the storage capacity.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{initiate.startValues_HReg}}, \code{\link{print.Bayes_HReg}}, \code{\link{summary.Bayes_HReg}}, \code{\link{predict.Bayes_HReg}}
}
\examples{

\dontrun{   
# loading a data set
data(scrData)
id=scrData$cluster

form <- Formula(time1 + event1 | time2 + event2 ~ x1 + x2 + x3 | x1 + x2 | x1 + x2)

#####################
## Hyperparameters ##
#####################

## Subject-specific frailty variance component
##  - prior parameters for 1/theta
##
theta.ab <- c(0.7, 0.7)

## Weibull baseline hazard function: alphas, kappas
##
WB.ab1 <- c(0.5, 0.01) # prior parameters for alpha1
WB.ab2 <- c(0.5, 0.01) # prior parameters for alpha2
WB.ab3 <- c(0.5, 0.01) # prior parameters for alpha3
##
WB.cd1 <- c(0.5, 0.05) # prior parameters for kappa1
WB.cd2 <- c(0.5, 0.05) # prior parameters for kappa2
WB.cd3 <- c(0.5, 0.05) # prior parameters for kappa3

## PEM baseline hazard function
##
PEM.ab1 <- c(0.7, 0.7) # prior parameters for 1/sigma_1^2
PEM.ab2 <- c(0.7, 0.7) # prior parameters for 1/sigma_2^2
PEM.ab3 <- c(0.7, 0.7) # prior parameters for 1/sigma_3^2
##
PEM.alpha1 <- 10 # prior parameters for K1
PEM.alpha2 <- 10 # prior parameters for K2
PEM.alpha3 <- 10 # prior parameters for K3

## MVN cluster-specific random effects
##
Psi_v <- diag(1, 3)
rho_v <- 100

## DPM cluster-specific random effects
##
Psi0  <- diag(1, 3)
rho0  <- 10
aTau  <- 1.5
bTau  <- 0.0125

##
hyperParams <- list(theta=theta.ab,
                WB=list(WB.ab1=WB.ab1, WB.ab2=WB.ab2, WB.ab3=WB.ab3,
                       WB.cd1=WB.cd1, WB.cd2=WB.cd2, WB.cd3=WB.cd3),
                   PEM=list(PEM.ab1=PEM.ab1, PEM.ab2=PEM.ab2, PEM.ab3=PEM.ab3,
                       PEM.alpha1=PEM.alpha1, PEM.alpha2=PEM.alpha2, PEM.alpha3=PEM.alpha3),
                   MVN=list(Psi_v=Psi_v, rho_v=rho_v),
                   DPM=list(Psi0=Psi0, rho0=rho0, aTau=aTau, bTau=bTau))
                    
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
nGam_save <- 0
storeV    <- rep(TRUE, 3)

## Tuning parameters for specific updates
##
##  - those common to all models
mhProp_theta_var  <- 0.05
mhProp_Vg_var     <- c(0.05, 0.05, 0.05)
##
## - those specific to the Weibull specification of the baseline hazard functions
mhProp_alphag_var <- c(0.01, 0.01, 0.01)
##
## - those specific to the PEM specification of the baseline hazard functions
Cg        <- c(0.2, 0.2, 0.2)
delPertg  <- c(0.5, 0.5, 0.5)
rj.scheme <- 1
Kg_max    <- c(50, 50, 50)
sg_max    <- c(max(scrData$time1[scrData$event1 == 1]),
               max(scrData$time2[scrData$event1 == 0 & scrData$event2 == 1]),
               max(scrData$time2[scrData$event1 == 1 & scrData$event2 == 1]))

time_lambda1 <- seq(1, sg_max[1], 1)
time_lambda2 <- seq(1, sg_max[2], 1)
time_lambda3 <- seq(1, sg_max[3], 1)               

##
mcmc.WB  <- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
                storage=list(nGam_save=nGam_save, storeV=storeV),
                tuning=list(mhProp_theta_var=mhProp_theta_var,
                mhProp_Vg_var=mhProp_Vg_var, mhProp_alphag_var=mhProp_alphag_var))

##
mcmc.PEM <- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
                storage=list(nGam_save=nGam_save, storeV=storeV),
                tuning=list(mhProp_theta_var=mhProp_theta_var,
                mhProp_Vg_var=mhProp_Vg_var, Cg=Cg, delPertg=delPertg,
                rj.scheme=rj.scheme, Kg_max=Kg_max,
                time_lambda1=time_lambda1, time_lambda2=time_lambda2,
                time_lambda3=time_lambda3))
    
#####################
## Starting Values ##
#####################

##
Sigma_V <- diag(0.1, 3)
Sigma_V[1,2] <- Sigma_V[2,1] <- -0.05
Sigma_V[1,3] <- Sigma_V[3,1] <- -0.06
Sigma_V[2,3] <- Sigma_V[3,2] <- 0.07

#################################################################
## Analysis of Independent Semi-Competing Risks Data ############
#################################################################

#############
## WEIBULL ##
#############

##
myModel <- c("semi-Markov", "Weibull")
myPath  <- "Output/01-Results-WB/"

startValues      <- initiate.startValues_HReg(form, scrData, model=myModel, nChain=2)

##
fit_WB <- BayesID_HReg(form, scrData, id=NULL, model=myModel,
                hyperParams, startValues, mcmc.WB, path=myPath)
				
fit_WB
summ.fit_WB <- summary(fit_WB); names(summ.fit_WB)
summ.fit_WB
pred_WB <- predict(fit_WB, tseq=seq(from=0, to=30, by=5))
plot(pred_WB, plot.est="Haz")
plot(pred_WB, plot.est="Surv")

#########
## PEM ##
#########

##						
myModel <- c("semi-Markov", "PEM")
myPath  <- "Output/02-Results-PEM/"

startValues      <- initiate.startValues_HReg(form, scrData, model=myModel, nChain=2)

##
fit_PEM <- BayesID_HReg(form, scrData, id=NULL, model=myModel,
                 hyperParams, startValues, mcmc.PEM, path=myPath)
				
fit_PEM
summ.fit_PEM <- summary(fit_PEM); names(summ.fit_PEM)
summ.fit_PEM
pred_PEM <- predict(fit_PEM)
plot(pred_PEM, plot.est="Haz")
plot(pred_PEM, plot.est="Surv")
					
#################################################################
## Analysis of Correlated Semi-Competing Risks Data #############
#################################################################

#################
## WEIBULL-MVN ##
#################

##
myModel <- c("semi-Markov", "Weibull", "MVN")
myPath  <- "Output/03-Results-WB_MVN/"

startValues      <- initiate.startValues_HReg(form, scrData, model=myModel, id, nChain=2)

##
fit_WB_MVN <- BayesID_HReg(form, scrData, id, model=myModel,
                    hyperParams, startValues, mcmc.WB, path=myPath)
                    
fit_WB_MVN
summ.fit_WB_MVN <- summary(fit_WB_MVN); names(summ.fit_WB_MVN)
summ.fit_WB_MVN
pred_WB_MVN <- predict(fit_WB_MVN, tseq=seq(from=0, to=30, by=5))
plot(pred_WB_MVN, plot.est="Haz")
plot(pred_WB_MVN, plot.est="Surv")


#################
## WEIBULL-DPM ##
#################

##
myModel <- c("semi-Markov", "Weibull", "DPM")
myPath  <- "Output/04-Results-WB_DPM/"

startValues      <- initiate.startValues_HReg(form, scrData, model=myModel, id, nChain=2)

##
fit_WB_DPM <- BayesID_HReg(form, scrData, id, model=myModel,
                    hyperParams, startValues, mcmc.WB, path=myPath)

fit_WB_DPM
summ.fit_WB_DPM <- summary(fit_WB_DPM); names(summ.fit_WB_DPM)
summ.fit_WB_DPM
pred_WB_DPM <- predict(fit_WB_MVN, tseq=seq(from=0, to=30, by=5))
plot(pred_WB_DPM, plot.est="Haz")
plot(pred_WB_DPM, plot.est="Surv")

#############
## PEM-MVN ##
#############

##
myModel <- c("semi-Markov", "PEM", "MVN")
myPath  <- "Output/05-Results-PEM_MVN/"

startValues      <- initiate.startValues_HReg(form, scrData, model=myModel, id, nChain=2)

##
fit_PEM_MVN <- BayesID_HReg(form, scrData, id, model=myModel,
                    hyperParams, startValues, mcmc.PEM, path=myPath)
                    
fit_PEM_MVN
summ.fit_PEM_MVN <- summary(fit_PEM_MVN); names(summ.fit_PEM_MVN)
summ.fit_PEM_MVN
pred_PEM_MVN <- predict(fit_PEM_MVN)
plot(pred_PEM_MVN, plot.est="Haz")
plot(pred_PEM_MVN, plot.est="Surv")

#############
## PEM-DPM ##
#############

##
myModel <- c("semi-Markov", "PEM", "DPM")
myPath  <- "Output/06-Results-PEM_DPM/"

startValues      <- initiate.startValues_HReg(form, scrData, model=myModel, id, nChain=2)
                    
##
fit_PEM_DPM <- BayesID_HReg(form, scrData, id, model=myModel,
                    hyperParams, startValues, mcmc.PEM, path=myPath)
                    
fit_PEM_DPM
summ.fit_PEM_DPM <- summary(fit_PEM_DPM); names(summ.fit_PEM_DPM)
summ.fit_PEM_DPM
pred_PEM_DPM <- predict(fit_PEM_DPM)
plot(pred_PEM_DPM, plot.est="Haz")
plot(pred_PEM_DPM, plot.est="Surv")
                    
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ semi-competing risks analysis }
\keyword{ hazard regression models }% __ONLY ONE__ keyword per line

