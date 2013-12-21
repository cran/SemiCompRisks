\name{BayesID}
\alias{BayesID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to fit Bayesian parametric and semi-parametric regression models to semi-competing risks data.
}
\description{
The function to fit Bayesian illness-death models to semi-competing risks data. 
In the model, the conditional hazard for the terminal event given that the non-terminal event has occurred (\eqn{h_3}) is assumed 
to be either Markov or semi-Markov with respect to the timing of the non-terminal event.
}
\usage{
BayesID(survData, nCov, hyperParams, startValues, mcmcParams, nGam_save, numReps, 
thin, path, burninPerc = 0.5, type = "semi-parametric", model = "Markov", nChain = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
	The data frame containing semi-competing risks data with covariate matrices from \code{n} subjects. See *Examples*. 
}
  \item{nCov}{
	a numeric vector that contains dimension of covariate space: c(\code{p1}, \code{p2}, \code{p3})
}
  \item{hyperParams}{
	a vector containing hyperparameter values for hierarchical model. See *Details* and *Examples*.
}
  \item{startValues}{
  	a list containing vectors of starting values of the parameters. See *Details* and *Examples*.
}
  \item{mcmcParams}{
	a vector containing variables required for Metropolis-Hastings-Green (MGH) algorithm. See *Details* and *Examples*.
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
  \item{type}{
	the type of analysis: either "\code{parametric}" or "\code{semi-parametric}"
}
  \item{model}{
	assumptions on \eqn{h_{03}}: either "\code{Markov}" or "\code{semi-Markov}"
}
  \item{nChain}{
	the number of chains
}
}
\details{
  In order to analyze semi-competing risks data, we view the data as arising from an underlying illness-death model system in which 
  individuals may undergo one or more of three transition: 1) no event to non-terminal event, 2) no event to terminal event, 3) non-terminal event to terminal event.
  The system of transitions is modeled via the specification of three hazard functions: 
  \deqn{h_1(t_{1i} | \gamma_i, x_{1i}) = \gamma_i h_{01}(t_{1i})\exp(\beta_1 x_{1i}), t_{1i}>0,}
  \deqn{h_2(t_{2i} | \gamma_i, x_{2i}) = \gamma_i h_{02}(t_{2i})\exp(\beta_2 x_{2i}), t_{2i}>0,}
  \deqn{h_3(t_{2i} | t_{1i}, \gamma_i, x_{3i}) = \gamma_i h_{03}(t_{2i})\exp(\beta_3 x_{3i}), 0<t_{1i}<t_{2i},}
  
  where \eqn{\gamma_i} is a subject-specific frailty, taken to be distributed independently of \eqn{x_{1i}}, \eqn{x_{2i}}, and \eqn{x_{3i}}.
  For \eqn{g \in \{1,2,3\}}, \eqn{h_{0g}} is an unspecified baseline hazard function and \eqn{\beta_g} is a vector of \eqn{p_g} log-hazard ratio regression parameters.
  The \eqn{h_{03}} is assumed to Markov with respect to \eqn{t_1}. we refer to the model specified by three conditional hazard functions as the `Markov model'.
  An alternative specification is to model the risk of terminal event following non-terminal event as a function of the sojourn time. Specifically, retaining \eqn{h_1} and \eqn{h_2} as above,
  we consider modeling \eqn{h_3} as follows:
  \deqn{h_3(t_{2i} | t_{1i}, \gamma_i, x_{3i}) = \gamma_i h_{03}(t_{2i}-t_{1i})\exp(\beta_3 x_{3i}), 0<t_{1i}<t_{2i}.}
  We refer to this alternative model as the `semi-Markov' model.
  For `semi-parametric analysis', let \eqn{s_{g,\max}}$ denote the largest observed event time for each transition \eqn{g \in \{1,2,3\}}. 
  Then, consider the finite partition of the relevant time axis into \eqn{J_{g} + 1} disjoint intervals: \eqn{0<s_{g,1}<s_{g,2}<...<s_{g, J_g+1} = s_{g, \max}}. 
  For notational convenience, let \eqn{I_{g,j}=(s_{g, j-1}, s_{g, j}]} denote the \eqn{j^{th}} partition. 
  For given a partition, \eqn{s_g = (s_{g,1}, \dots, s_{g, J_g + 1})}, we assume the log-baseline hazard functions is piecewise constant:
  \deqn{\lambda_{0g}(t)=\log h_{0g}(t) = \sum_{j=1}^{J_g + 1} \lambda_{g,j} I(t\in I_{g,j})}

   where \eqn{I(\cdot)} is the indicator function and \eqn{s_{g,0} \equiv 0}. Note, this specification is general in that the partitions of the time axes differ across the three hazard functions.
   In our proposed Bayesian framework, our prior choices are, for \eqn{g\in\{1,2,3\}}:

	\deqn{\pi(\beta_g) \propto 1,}
	\deqn{\lambda_g | J_g, \mu_{\lambda_g}, \sigma_{\lambda_g}^2 \sim \mathcal{N}_{J_g+1}(\mu_{\lambda_g}1, \sigma_{\lambda_g}^2\Sigma_{\lambda_g}),}
	\deqn{J_g \sim \mathcal{P}(\alpha_g),}
	\deqn{\pi(s_g | J_g) \propto \frac{(2J_g+1)! \prod_{j=1}^{J_g+1}(s_j-s_{j-1})}{(s_{g,J_g+1})^{(2J_g+1)}},}
	\deqn{\pi(\mu_{\lambda_g}) \propto1,}
	\deqn{\sigma_{\lambda_g}^{-2} \sim \mathcal{G}(a_g, b_g),}
	\deqn{\gamma_i|\theta \sim \mathcal{G}(\theta^{-1}, \theta^{-1}),}
	\deqn{\theta^{-1} \sim \mathcal{G}(\psi, \omega).}
	
	Note that \eqn{J_g} and \eqn{s_g} are treated as random and the priors for \eqn{J_g} and \eqn{s_g} jointly form a time-homogeneous Poisson process prior for the partition.
	The number of time splits and their positions are therefore updated within our computational scheme using reversible jump MCMC.
	
	For `parametric analysis', we assume the Weibull model for \eqn{h_{0g}(t) = \alpha_g \kappa_g t^{\alpha_g-1}}.
   In our Bayesian framework, our prior choices are, for \eqn{g\in\{1,2,3\}}:	
 	\deqn{\pi(\beta_g) \propto 1,}
 	\deqn{\pi(\alpha_g) \sim \mathcal{G}(a_g, b_g),}
 	\deqn{\pi(\kappa_g) \sim \mathcal{G}(c_g, d_g),} 	
	\deqn{\gamma_i|\theta \sim \mathcal{G}(\theta^{-1}, \theta^{-1}),}
	\deqn{\theta^{-1} \sim \mathcal{G}(\psi, \omega).}  

}


\value{
\code{BayesID} returns an object of class \code{BaysID}. \cr
\code{names(object$chain1)} shows the list of posterior samples of model parameters, the number of acceptance in MGH algorithm, etc.
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2013). 
Bayesian Semi-parametric Analysis of Semi-competing Risks Data: 
Estimating Readmission Rates among Pancreatic Cancer Patients, submitted.
}
\author{
Kyu Ha Lee and Sebastien Haneuse\cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
\note{
The posterior samples of \eqn{\gamma} are saved separately in \code{working directory/path}. 
For a dataset with large \eqn{n}, \code{nGam_save} should be carefully specified considering the system memory and the storage capacity.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{print.BayesID}}, \code{\link{summary.BayesID}}, \code{\link{plot.BayesID}}, \code{\link{ehr}}
}
\examples{

\dontrun{

# loading simulated data set
data(scrData)
	
n = dim(scrData)[1]	
p1 = 2	
p2 = 2
p3 = 2
nCov <- c(p1, p2, p3) 	

###############################
# setting hyperparameter values for semi-parametric analysis
	
a1		<- 0.7	# prior parameters for 1/sigma_1^2
b1		<- 0.7
a2		<- 0.7	# prior parameters for 1/sigma_2^2
b2		<- 0.7
a3		<- 0.7	# prior parameters for 1/sigma_3^2
b3		<- 0.7
alpha1	<- 10	# prior parameter for J1
alpha2	<- 10	# prior parameter for J2
alpha3	<- 10	# prior parameter for J3
c_lam1 	<- 1	# prior parameter for MVN-ICAR specification
c_lam2 	<- 1	
c_lam3 	<- 1
psi	<- 0.7	# prior parameters for 1/theta
omega	<- 0.7	

nGam_save <- 5  # the number of subjects whose gamma parameters are being saved

hyperParams <- c(a1, b1, a2, b2, a3, b3, alpha1, alpha2, alpha3, c_lam1, c_lam2, c_lam3, 
psi, omega)

#########################
# setting starting values
			
s1_max		<- max(scrData$time1[scrData$event1 == 1])
s2_max		<- max(scrData$time2[scrData$event1 == 0 & scrData$event2 == 1])
s3_max		<- max(scrData$time2[scrData$event1 == 1 & scrData$event2 == 1])

beta1		<- c(0.1, 0.1)
beta2		<- c(0.1, 0.1)
beta3		<- c(0.1, 0.1)
s1		<- c(seq(4, 36, 4), s1_max)
s2		<- c(seq(4, 36, 4), s2_max)
s3		<- c(seq(4, 36, 4), s3_max)
J1		<- length(s1) - 1
J2		<- length(s2) - 1
J3		<- length(s3) - 1
lambda1		<- runif(J1+1, -2, -1)
lambda2		<- runif(J2+1, -2, -1)
lambda3		<- runif(J3+1, -2, -1)

sigSq_lam1	<- var(lambda1)
sigSq_lam2	<- var(lambda2)
sigSq_lam3	<- var(lambda3)
mu_lam1		<- mean(lambda1)
mu_lam2		<- mean(lambda2)
mu_lam3		<- mean(lambda3)
theta		<- 3.3
gamma		<- rgamma(n, 1/theta, 1/theta)

# chain 1

startValues <- list()
startValues[[1]] <- as.vector(c(beta1, beta2, beta3, J1, J2, J3, mu_lam1, mu_lam2, mu_lam3,
 sigSq_lam1, sigSq_lam2, sigSq_lam3, theta, gamma, lambda1, lambda2, lambda3, s1, s2, s3))

beta1		<- c(0.2, 0.3)
beta2		<- c(0.2, 0.3)
beta3		<- c(0.2, 0.3)
lambda1		<- runif(J1+1, -2, -1)
lambda2		<- runif(J2+1, -2, -1)
lambda3		<- runif(J3+1, -2, -1)

# chain 2

startValues[[2]] <- as.vector(c(beta1, beta2, beta3, J1, J2, J3, mu_lam1, mu_lam2, mu_lam3, 
sigSq_lam1, sigSq_lam2, sigSq_lam3, theta, gamma, lambda1, lambda2, lambda3, s1, s2, s3))

##################################################
# setting variable values needed for MHG algorithm

C1		<- 0.20
C2		<- 0.20
C3		<- 0.20
delPert1	<- 0.5
delPert2	<- 0.5
delPert3	<- 0.5
s_propBI1	<- floor(sort(unique(scrData$time1[scrData$event1 == 1])))
s_propBI1	<- unique(s_propBI1[s_propBI1 < s1_max])
s_propBI2	<- floor(sort(unique(scrData$time2[scrData$event1 == 0 & scrData$event2 == 1])))
s_propBI2	<- unique(s_propBI2[s_propBI2 < s2_max])
s_propBI3	<- floor(sort(unique(scrData$time2[scrData$event1 == 1 & scrData$event2 == 1])))
s_propBI3	<- unique(s_propBI3[s_propBI3 < s3_max])

num_s_propBI1	<- length(s_propBI1)
num_s_propBI2	<- length(s_propBI2)
num_s_propBI3	<- length(s_propBI3)
J1_max 		<- 50
J2_max 		<- 50
J3_max 		<- 50
time_lambda1	<- 1:36
time_lambda2	<- 1:36
time_lambda3	<- 1:36
nTime_lambda1	<- length(time_lambda1)
nTime_lambda2	<- length(time_lambda2)
nTime_lambda3	<- length(time_lambda3)

mhProp_theta_var <- 0.05

mcmcParams <- c(C1, C2, C3, delPert1, delPert2, delPert3, num_s_propBI1, num_s_propBI2, 
num_s_propBI3, J1_max, J2_max, J3_max, s1_max, s2_max, s3_max, nTime_lambda1, 
nTime_lambda2, nTime_lambda3, s_propBI1, s_propBI2, s_propBI3, time_lambda1, 
time_lambda2, time_lambda3, mhProp_theta_var)

##################################################
# number of chains

numReps		= 2e6
thin		= 20
burninPerc	= 0.5
path1		= "outcome/"
type 		= "semi-parametric"
model 		= "Markov"
nChain		= 2

# fitting Bayesian semi-parametric regression model to semi-competing risks data	
# In practice, set 'numReps' to a larger number such as 1000000

fitScr <- BayesID(scrData, nCov, hyperParams, startValues, mcmcParams, nGam_save, 
numReps, thin, path1, burninPerc, type, model, nChain)
		
print(fitScr)
summary(fitScr)

## plot for estimates of baseline hazard functions
plot(fitScr)	
	
## calculating conditional explanatory hazard ratio		

vals <- ehr(c(1.5, 0), c(2, 1), fitScr)
plot(vals)

}

		
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line

