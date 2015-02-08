\name{BayesIDcor}
\alias{BayesIDcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to fit parametric and semi-parametric hierarchical models to cluster-correlated semi-competing risks data.
}
\description{
Cluster-correlated semi-competing risks data can be analyzed using hierarchical models. The prior of baseline hazard functions can be specified by either parametric (Weibull) model or non-parametric mixture of piecewise exponential models (PEM). The option to choose between parametric multivariate normal (MVN) and non-parametric Dirichlet process mixture of multivariate normals (DPM) is available for the prior of cluster-specific random effects distribution. The conditional hazard function for time to the terminal event given time to non-terminal can be modeled as Markov (it does not depend on the timing of the non-terminal event) or semi-Markov (it does depend on the timing)
}
\usage{
BayesIDcor(survData, nCov, hyperParams, startValues, mcmcParams, nGam_save, numReps, 
thin, path = "results/", burninPerc = 0.5, hz.type = "Weibull", re.type = "MVN",
model = "Markov", storeV= rep(TRUE, 3), nChain = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
	The data frame containing semi-competing risks data with cluster information and covariate matrices from \code{n} subjects. It is of dimension \eqn{n\times(5+p_1+p_2+p_3)}, where \eqn{p_1},\eqn{p_2},\eqn{p_3} are the number of covariates included in the conditional hazard functions \eqn{h_1}, \eqn{h_2}, \eqn{h_3}, repectively. The first five columns correspond to \eqn{y_1}, \eqn{\delta_1}, \eqn{y_2}, \eqn{\delta_2}, cluster membership. Note that the cluster membership must set to consecutive positive integers; 1:J.  
}
  \item{nCov}{
	a numeric vector that contains dimensions of covariate space: c(\eqn{p_1},\eqn{p_2},\eqn{p_3})
}
  \item{hyperParams}{
	a vector containing hyperparameter values for hierarchical model.\cr
      For Weibull-MVN, (\code{hz.type}="Weibull", \code{re.type}="MVN"), it is of length 24; c(\eqn{a_1}, \eqn{b_1}, \eqn{a_2}, \eqn{b_2}, \eqn{a_3}, \eqn{b_3}, \eqn{c_1}, \eqn{d_1}, \eqn{c_2}, \eqn{d_2}, \eqn{c_3}, \eqn{d_3},\eqn{\psi}, \eqn{\omega}, \eqn{\rho_v}, \eqn{\Psi_v}). Note the matrix \eqn{\Psi_v} must be read by columns.\cr
    For Weibull-DPM, (\code{hz.type}="Weibull", \code{re.type}="DPM"), it is of length 30; c(\eqn{a_1}, \eqn{b_1}, \eqn{a_2}, \eqn{b_2}, \eqn{a_3}, \eqn{b_3}, \eqn{c_1}, \eqn{d_1}, \eqn{c_2}, \eqn{d_2}, \eqn{c_3}, \eqn{d_3},\eqn{\psi}, \eqn{\omega}, \eqn{\mu_0}, \eqn{\zeta_0}, \eqn{\Psi_0}, \eqn{\rho_0}, \eqn{a_{\tau}}, \eqn{b_{\tau}}). Note the matrix \eqn{\Psi_0} must be read by columns.\cr
    For PEM-MVN, (\code{hz.type}="PEM", \code{re.type}="MVN"), it is of length 24; c(\eqn{a_1}, \eqn{b_1}, \eqn{a_2}, \eqn{b_2}, \eqn{a_3}, \eqn{b_3}, \eqn{\alpha_1}, \eqn{\alpha_2}, \eqn{\alpha_3}, \eqn{c_{\lambda1}}, \eqn{c_{\lambda2}}, \eqn{c_{\lambda3}}, \eqn{\psi}, \eqn{\omega}, \eqn{\rho_v}, \eqn{\Psi_v}). Note the matrix \eqn{\Psi_v} must be read by columns. \cr
    For PEM-DPM, (\code{hz.type}="PEM", \code{re.type}="DPM"), it is of length 30; c(\eqn{a_1}, \eqn{b_1}, \eqn{a_2}, \eqn{b_2}, \eqn{a_3}, \eqn{b_3}, \eqn{\alpha_1}, \eqn{\alpha_2}, \eqn{\alpha_3}, \eqn{c_{\lambda1}}, \eqn{c_{\lambda2}}, \eqn{c_{\lambda3}}, \eqn{\psi}, \eqn{\omega}, \eqn{\mu_0}, \eqn{\zeta_0}, \eqn{\Psi_0}, \eqn{\rho_0}, \eqn{a_{\tau}}, \eqn{b_{\tau}}). Note the matrix \eqn{\Psi_0} must be read by columns.
}
  \item{startValues}{
  	a list containing vectors of starting values of parameters. The length of list must be equal to the number of chains (\code{nChain}). \cr
    For Weibull-MVN, (\code{hz.type}="Weibull", \code{re.type}="MVN"), each of list components is of length \eqn{p_1+p_2+p_3+n+3J+16}; c(\eqn{\beta_1}, \eqn{\beta_2}, \eqn{\beta_3}, \eqn{\alpha_1}, \eqn{\alpha_2}, \eqn{\alpha_3}, \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3}, \eqn{\theta}, \eqn{\gamma}, \eqn{V_1}, \eqn{V_2}, \eqn{V_3}, \eqn{\Sigma_V}). Note the matrix \eqn{\Sigma_V} must be read by columns.\cr
    For Weibull-DPM, (\code{hz.type}="Weibull", \code{re.type}="DPM"), each of list components is of length \eqn{p_1+p_2+p_3+n+4J+8}; c(\eqn{\beta_1}, \eqn{\beta_2}, \eqn{\beta_3}, \eqn{\alpha_1}, \eqn{\alpha_2}, \eqn{\alpha_3}, \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3}, \eqn{\theta}, \eqn{\gamma}, \eqn{V_1}, \eqn{V_2}, \eqn{V_3}, \eqn{c}, \eqn{\tau}). \cr
    For PEM-MVN, (\code{hz.type}="PEM", \code{re.type}="MVN"), each of list components is of length \eqn{p_1+p_2+p_3+10+n+2(K_1+1)+2(K_2+1)+2(K_3+1)+2+3J+9};       c(\eqn{\beta_1}, \eqn{\beta_2}, \eqn{\beta_3}, \eqn{K_1}, \eqn{K_2}, \eqn{K_3}, \eqn{\mu_{\lambda_1}}, \eqn{\mu_{\lambda_2}}, \eqn{\mu_{\lambda_3}}, \eqn{\sigma_{\lambda_1}^2}, \eqn{\sigma_{\lambda_2}^2}, \eqn{\sigma_{\lambda_3}^2}, \eqn{\theta}, \eqn{\gamma}, \eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\lambda_3}, \eqn{s_1}, \eqn{s_2}, \eqn{s_3}, \eqn{V_1}, \eqn{V_2}, \eqn{V_3}, \eqn{\Sigma_V}). Note the matrix \eqn{\Sigma_V} must be read by columns.
    For PEM-DPM, (\code{hz.type}="PEM", \code{re.type}="DPM"), each of list components is of length \eqn{p_1+p_2+p_3+10+n+2(K_1+1)+2(K_2+1)+2(K_3+1)+2+4J+1}; c(\eqn{\beta_1}, \eqn{\beta_2}, \eqn{\beta_3}, \eqn{K_1}, \eqn{K_2}, \eqn{K_3}, \eqn{\mu_{\lambda_1}}, \eqn{\mu_{\lambda_2}}, \eqn{\mu_{\lambda_3}}, \eqn{\sigma_{\lambda_1}^2}, \eqn{\sigma_{\lambda_2}^2}, \eqn{\sigma_{\lambda_3}^2}, \eqn{\theta}, \eqn{\gamma}, \eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\lambda_3}, \eqn{s_1}, \eqn{s_2}, \eqn{s_3}, \eqn{V_1}, \eqn{V_2}, \eqn{V_3}, \eqn{c}, \eqn{\tau}). 
      
    
}
  \item{mcmcParams}{
	a vector containing variables required for Metropolis-Hastings-Green (MHG) algorithm. \cr
      For Weibull models (Weibull-MVN and Weibull-DPM: \code{hz.type}="Weibull"), it is of length 7; c(mhProp_alpha1_var, mhProp_alpha2_var, mhProp_alpha3_var, mhProp_theta_var, mhProp_V1_var, mhProp_V2_var, mhProp_V3_var). mhProp_alpha1_var and mhProp_V1_var are the variance of proposal density in Metropolis-Hastings step for \eqn{\alpha_1} and \eqn{V_1}, respectively. \cr
      For PEM models (PEM-MVN and PEM-DPM: \code{hz.type}="PEM"), it is of length \eqn{(18+num\_s\_propBI1+num\_s\_propBI2+num\_s\_propBI3+nTime\_lambda1+nTime\_lambda2+nTime\_lambda3+4}); c(\eqn{C_1}, \eqn{C_2}, \eqn{C_3}, delPert1, delPert2, delPert3, num_s_propBI1, num_s_propBI2, num_s_propBI3, \eqn{K_{1max}}, \eqn{K_{2max}}, \eqn{K_{3max}}, \eqn{s_{1max}}, \eqn{s_{2max}}, \eqn{s_{3max}}, nTime_lambda1, nTime_lambda2, nTime_lambda3, s_propBI1, s_propBI2, s_propBI3, time_lambda1,time_lambda2, time_lambda3, mhProp_theta_var, mhProp_V1_var, mhProp_V2_var, mhProp_V3_var). The perturbation random variable is drawn from Uniform(0.5-delPert1, 0.5+delPert1) in the birth/death move for \eqn{h_{01}}. nTime_lambda1 is \code{length}(time_lambda1), where time_lambda1 is a set of time points where the posterior estimates of \eqn{log(h_{01})} is calculated. s_propBI1 is a set of candidates for proposal values in birth move of \eqn{h_{01}}. mhProp_theta_var and mhProp_V1_var are the variance of proposal density in Metropolis-Hastings step for \eqn{\theta} and \eqn{V_1}, respectively.
}
  \item{nGam_save}{
	the number of \eqn{\gamma} to be stored
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
	prior specification for cluster-specific random effects distribution: either "\code{MVN}" or "\code{DPM}"
}
  \item{model}{
	assumptions on \eqn{h_{03}}: either "\code{Markov}" or "\code{semi-Markov}"
}
  \item{storeV}{
	logical values to determine whether all the posterior samples of \eqn{V_1}, \eqn{V_2}, \eqn{V_3} are to be stored.
}
  \item{nChain}{
	the number of chains
}
}

\details{
    We view the semi-competing risks data as arising from an underlying illness-death model system in which individuals may undergo one or more of three transitions: 1) no event to non-terminal event, 2) no event to terminal event, 3) non-terminal event to terminal event. Let \eqn{t_{ji1}}, \eqn{t_{ji2}} denote time to non-terminal and terminal event from subject \eqn{i=1,...,n_j} in cluster \eqn{1,...,J}. The system of transitions is modeled via the specification of three hazard functions:
    \deqn{h_1(t_{ji1} | \gamma_{ji}, V_{j1}, x_{ji1}) = \gamma_{ji} h_{01}(t_{ji1})\exp(\beta_1 x_{ji1}+V_{j1}), t_{ji1}>0,}
    \deqn{h_2(t_{ji2} | \gamma_{ji}, V_{j2}, x_{ji2}) = \gamma_{ji} h_{02}(t_{ji2})\exp(\beta_2 x_{ji2}+V_{j2}), t_{ji2}>0,}
    \deqn{h_3(t_{ji2} | t_{ji1}, \gamma_{ji}, V_{j3}, x_{ji3}) = \gamma_{ji} h_{03}(t_{ji2})\exp(\beta_3 x_{ji3}+V_{j3}), 0<t_{ji1}<t_{ji2},}
    
    where \eqn{\gamma_{ji}} is a subject-specific frailty and \eqn{V_j}=(\eqn{V_{j1}}, \eqn{V_{j2}}, \eqn{V_{j3}}) is a vector of cluster-specific random effects (each specific to one of the three possible transitions), taken to be distributed independently of \eqn{x_{ji1}}, \eqn{x_{ji2}}, and \eqn{x_{ji3}}.
    For \eqn{g \in \{1,2,3\}}, \eqn{h_{0g}} is an unspecified baseline hazard function and \eqn{\beta_g} is a vector of \eqn{p_g} log-hazard ratio regression parameters.
    The \eqn{h_{03}} is assumed to be Markov with respect to \eqn{t_1}. We refer to the model specified by three conditional hazard functions as the Markov model.
    An alternative specification is to model the risk of terminal event following non-terminal event as a function of the sojourn time. Specifically, retaining \eqn{h_1} and \eqn{h_2} as above,
    we consider modeling \eqn{h_3} as follows:
    \deqn{h_3(t_{ji2} | t_{ji1}, \gamma_{ji}, V_{j3}, x_{ji3}) = \gamma_{ji} h_{03}(t_{ji2}-t_{ji1})\exp(\beta_3 x_{ji3}+V_{j3}), 0<t_{ji1}<t_{ji2}.}
    We refer to this alternative model as the semi-Markov model. \cr
    For parametric MVN prior specification for a vector of cluster-specific random effects, we assume \eqn{V_j} arise as i.i.d. draws from a mean 0 multivariate Normal distribution with variance-covariance matrix \eqn{\Sigma_V}. The diagonal elements of the \eqn{3\times 3} matrix \eqn{\Sigma_V} characterize variation across clusters in risk for non-terminal, terminal and terminal following non-terminal event, respectively, which is not explained by covariates included in the linear predictors. Specifically, the priors can be written as follows:
    \deqn{V_j \sim MVN(0, \Sigma_V),}
 	\deqn{\Sigma_V \sim inverse-Wishart(\Psi_v, \rho_v).}
    For DPM prior specification for \eqn{V_j}, we consider non-parametric Dirichlet process mixture of MVN distributions: the \eqn{V_j}'s' are draws from a finite mixture of M multivariate Normal distributions, each with their own mean vector and variance-covariance matrix, (\eqn{\mu_m}, \eqn{\Sigma_m}) for \eqn{m=1,...,M}. Let \eqn{m_j\in\{1,...,M\}} denote the specific component to which the \eqn{j}th cluster belongs. Since the class-specific (\eqn{\mu_m}, \eqn{\Sigma_m}) are unknown they are taken to be draws from some distribution, \eqn{G_0}, often referred to as the centering distribution. Furthermore, since the true class memberships are not known, we denote the probability that the \eqn{j}th cluster belongs to any given class by the vector \eqn{p=(p_1,..., p_M)} whose components add up to 1.0. In the absence of prior knowledge regarding the distribution of class memberships for the \eqn{J} clusters across the \eqn{M} classes, a natural prior for \eqn{p} is the conjugate symmetric \eqn{Dirichlet(\tau/M,...,\tau/M)} distribution; the hyperparameter, \eqn{\tau}, is often referred to as a the precision parameter. The prior can be represented as follows (\eqn{M} goes to infinity):
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
    Our prior choice for remaining model parameters in all of four models (Weibull-MVN, Weibull-DPM, PEM-MVN, PEM-DPM) is as follows:
    \deqn{\pi(\beta_g) \propto 1,}
	\deqn{\gamma_{ji}|\theta \sim Gamma(\theta^{-1}, \theta^{-1}),}
	\deqn{\theta^{-1} \sim Gamma(\psi, \omega).}
}



\value{
\code{BayesIDcor} returns an object of class \code{BayesIDcor}. \cr
\code{names(object$chain1)} shows the list of posterior samples of model parameters, the number of acceptance in MHG algorithm, etc.
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2014), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, \emph{in press}.\cr \cr
Lee, K. H., Dominici, F., Schrag, D., and Haneuse, S.,
Hierarchical models for cluster-correlated semicompeting risks data with application to quality of end-of-life care for pancreatic cancer, \emph{arXiv:1502.00526}; submitted. \cr
}
\author{
Kyu Ha Lee and Sebastien Haneuse\cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
\note{
The posterior samples of \eqn{\gamma} and \eqn{V_j} are saved separately in \code{working directory/path}.
For a dataset with large \eqn{n}, \code{nGam_save} should be carefully specified considering the system memory and the storage capacity.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{BayesID}}, \code{\link{print.BayesID}}, \code{\link{summary.BayesID}}, \code{\link{plot.BayesID}}
}
\examples{

\dontrun{

# loading a data set
data(scrCorData)
    
scrData <- scrCorData
	
n = dim(scrData)[1]
p1 = 3
p2 = 3
p3 = 3
nCov <- c(p1, p2, p3)
J = length(unique(scrData[,5]))

###################################################################################
# analysis of clustered semi-competing risks data using semi-Markov PEM-MVN model #
###################################################################################
	
a1		<- 0.7	# prior parameters for 1/sigma_1^2
b1		<- 0.7
a2		<- 0.7	# prior parameters for 1/sigma_2^2
b2		<- 0.7
a3		<- 0.7	# prior parameters for 1/sigma_3^2
b3		<- 0.7
alpha1	<- 10	# prior parameter for K1
alpha2	<- 10	# prior parameter for K2
alpha3	<- 10	# prior parameter for K3
c_lam1 	<- 1	# prior parameter for MVN-ICAR specification
c_lam2 	<- 1
c_lam3 	<- 1
psi		<- 0.7	# prior parameters for 1/theta
omega 	<- 0.7

rho_v 	<- 100
Psi_v	<- diag(1, 3)
    
nGam_save <- 10 # the number of subjects whose gamma parameters are being saved
    
hyperParams <- c(a1, b1, a2, b2, a3, b3, alpha1, alpha2, alpha3, c_lam1, c_lam2, c_lam3,
psi, omega, rho_v, Psi_v)


#########################
# setting starting values
			
s1_max		<- max(scrData$time1[scrData$event1 == 1])
s2_max		<- max(scrData$time2[scrData$event1 == 0 & scrData$event2 == 1])
s3_max		<- max(scrData$time2[scrData$event1 == 1 & scrData$event2 == 1])
    
grid1 <- s1_max/50
grid2 <- s2_max/50
grid3 <- s3_max/50

# chain 1    
    
beta1		<- rep(0.1, p1)
beta2		<- rep(0.1, p2)
beta3		<- rep(0.1, p3)
    
s1			<- unique(sort(c(sample(1:s1_max, 10), s1_max)))
s2			<- unique(sort(c(sample(1:s2_max, 10), s2_max)))
s3			<- unique(sort(c(sample(1:s3_max, 10), s3_max)))
K1			<- length(s1) - 1
K2			<- length(s2) - 1
K3			<- length(s3) - 1
lambda1		<- runif(K1+1, -3, -2)
lambda2		<- runif(K2+1, -3, -2)
lambda3		<- runif(K3+1, -3, -2)
    
sigSq_lam1	<- var(lambda1)
sigSq_lam2	<- var(lambda2)
sigSq_lam3	<- var(lambda3)
mu_lam1		<- mean(lambda1)
mu_lam2		<- mean(lambda2)
mu_lam3		<- mean(lambda3)
theta		<- 1
gamma		<- rep(1, n)
    
V1	<- rep(0, J)
V2	<- rep(0, J)
V3	<- rep(0, J)
    
Sigma_V <- diag(1, 3)

startValues <- list()
startValues[[1]] <- as.vector(c(beta1, beta2, beta3, K1, K2, K3, mu_lam1, mu_lam2, mu_lam3,
sigSq_lam1, sigSq_lam2, sigSq_lam3, theta, gamma, lambda1, lambda2, lambda3, s1, s2, s3,
V1, V2, V3, Sigma_V))
    
# chain 2    

beta1		<- rep(-0.1, p1)
beta2		<- rep(-0.1, p2)
beta3		<- rep(-0.1, p3)
    
K1			<- length(s1) - 1
K2			<- length(s2) - 1
K3			<- length(s3) - 1
lambda1		<- runif(K1+1, -3, -2)
lambda2		<- runif(K2+1, -3, -2)
lambda3		<- runif(K3+1, -3, -2)
    
sigSq_lam1	<- var(lambda1)
sigSq_lam2	<- var(lambda2)
sigSq_lam3	<- var(lambda3)
mu_lam1		<- mean(lambda1)
mu_lam2		<- mean(lambda2)
mu_lam3		<- mean(lambda3)
theta		<- 0.5
gamma		<- rep(1.1, n)
    
V1	<- rep(0.1, J)
V2	<- rep(0.1, J)
V3	<- rep(0.1, J)
    
Sigma_V <- diag(1.1, 3)
    
startValues[[2]] <- as.vector(c(beta1, beta2, beta3, K1, K2, K3, mu_lam1, mu_lam2, mu_lam3,
sigSq_lam1, sigSq_lam2, sigSq_lam3, theta, gamma, lambda1, lambda2, lambda3, s1, s2, s3,
V1, V2, V3, Sigma_V))
    
    
    
    
    

##################################################
# setting variable values needed for MHG algorithm

C1				<- 0.20
C2				<- 0.20
C3				<- 0.20
delPert1		<- 0.5
delPert2		<- 0.5
delPert3		<- 0.5
s_propBI1		<- seq(1, s1_max, 1)
s_propBI1		<- s_propBI1[s_propBI1 < s1_max]
s_propBI2		<- seq(1, s2_max, 1)
s_propBI2		<- s_propBI2[s_propBI2 < s2_max]
s_propBI3		<- seq(1, s3_max, 1)
s_propBI3		<- s_propBI3[s_propBI3 < s3_max]
    
num_s_propBI1	<- length(s_propBI1)
num_s_propBI2	<- length(s_propBI2)
num_s_propBI3	<- length(s_propBI3)
K1_max 			<- 50
K2_max 			<- 50
K3_max 			<- 50
    
time_lambda1	<- seq(grid1, s1_max, grid1)
time_lambda2	<- seq(grid2, s2_max, grid2)
time_lambda3	<- seq(grid3, s3_max, grid3)
nTime_lambda1	<- length(time_lambda1)
nTime_lambda2	<- length(time_lambda2)
nTime_lambda3	<- length(time_lambda3)
    
mhProp_theta_var <- 0.05
mhProp_V1_var <- 0.05
mhProp_V2_var <- 0.05
mhProp_V3_var <- 0.05
    
mcmcParams <- c(C1, C2, C3, delPert1, delPert2, delPert3, num_s_propBI1, num_s_propBI2,
num_s_propBI3, K1_max, K2_max, K3_max, s1_max, s2_max, s3_max, nTime_lambda1, nTime_lambda2,
nTime_lambda3, s_propBI1, s_propBI2, s_propBI3, time_lambda1, time_lambda2, time_lambda3,
mhProp_theta_var,  mhProp_V1_var, mhProp_V2_var, mhProp_V3_var)
    
    


##################################################
# number of chains

numReps		= 2e6
thin		= 20
burninPerc	= 0.5
path        = "results/"
hz.type 	= "PEM"
re.type     = "MVN"
model 		= "semi-Markov"
storeV      = rep(TRUE, 3)
nChain		= 2

# fitting semi-Markov PEM-MVN model to clustered semi-competing risks data	

fit <- BayesIDcor(scrData, nCov, hyperParams, startValues, mcmcParams, nGam_save, numReps,
    thin, path, burninPerc, hz.type, re.type, model, storeV, nChain)
		
print(fit)
summary(fit)

## plot for estimates of baseline hazard functions
plot(fit)	
	
}

		
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line

