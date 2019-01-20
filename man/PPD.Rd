\name{PPD}
\alias{PPD}
	%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Function to predict the joint probability involving two event times in Bayesian illness-death models
}
\description{
	\code{PPD} is a function to predict the joint probability involving two event times in Bayesian illness-death models.
}
\usage{
PPD(fit, x1, x2, x3, t1, t2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
    an object of class \code{Bayes_HReg}. Currently, the function is available for PEM illness-death models.
}
  \item{x1}{
	a vector of covariates for \eqn{h_1} with which to predict.
}
  \item{x2}{
	a vector of covariates for \eqn{h_2} with which to predict.
}
  \item{x3}{
	a vector of covariates for \eqn{h_3} with which to predict.
}
  \item{t1}{
	time to non-terminal event for which the joint probability is calculated.
}
  \item{t2}{
	time to terminal event for which the joint probability is calculated.
}
}
\details{
Using the posterior predictive density, given (\eqn{x_1}, \eqn{x_2}, \eqn{x_3}), one can predict any joint probability involving the two event times such as \eqn{P(T_1<t_1, T_2<t_2| x_1, x_2, x_3)} for \eqn{0<t_1\le t_2} and \eqn{P(T_1=\infty, T_2<t_2| x_1, x_2, x_3)} for \eqn{t_2>0}.
}
\value{
\item{F_u}{Predicted \eqn{P(T_1\le t_1, T_2\le t_2| x_1, x_2, x_3)} in the upper wedge of the support of \eqn{(T_1, T_2)}.}
\item{F_l}{Predicted \eqn{P(T_1=\infty, T_2\le t_2| x_1, x_2, x_3)} in the lower wedge of the support of \eqn{(t1, t2)}.}
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2015), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, 64, 2, 253-273.\cr
}
\author{
Kyu Ha Lee and Sebastien Haneuse\cr
Maintainer: Kyu Ha Lee <klee15239@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BayesID_HReg}}
}
\examples{
\dontrun{   
# loading a data set
data(scrData)
id=scrData$cluster

form <- Formula(time1 + event1 | time2 + event2 ~ x1 + x2 | x1 + x2 | x1 + x2)

#####################
## Hyperparameters ##
#####################

## Subject-specific frailty variance component
##  - prior parameters for 1/theta
##
theta.ab <- c(0.7, 0.7)

## PEM baseline hazard function
##
PEM.ab1 <- c(0.7, 0.7) # prior parameters for 1/sigma_1^2
PEM.ab2 <- c(0.7, 0.7) # prior parameters for 1/sigma_2^2
PEM.ab3 <- c(0.7, 0.7) # prior parameters for 1/sigma_3^2
##
PEM.alpha1 <- 10 # prior parameters for K1
PEM.alpha2 <- 10 # prior parameters for K2
PEM.alpha3 <- 10 # prior parameters for K3

##
hyperParams <- list(theta=theta.ab,
                   PEM=list(PEM.ab1=PEM.ab1, PEM.ab2=PEM.ab2, PEM.ab3=PEM.ab3,
                    PEM.alpha1=PEM.alpha1, PEM.alpha2=PEM.alpha2,
                    PEM.alpha3=PEM.alpha3))
                    
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

## Tuning parameters for specific updates
##
##  - those common to all models
mhProp_theta_var  <- 0.05
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
mcmc.PEM <- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
                storage=list(nGam_save=nGam_save),
                tuning=list(mhProp_theta_var=mhProp_theta_var,
                Cg=Cg, delPertg=delPertg,
                rj.scheme=rj.scheme, Kg_max=Kg_max,
                time_lambda1=time_lambda1, time_lambda2=time_lambda2,
                time_lambda3=time_lambda3))
    
##						
myModel <- c("semi-Markov", "PEM")
myPath  <- "Output/02-Results-PEM/"

startValues      <- initiate.startValues_HReg(form, scrData, model=myModel, nChain=2)

##
fit_PEM <- BayesID_HReg(form, scrData, id=NULL, model=myModel,
                 hyperParams, startValues, mcmc.PEM, path=myPath)
				
PPD(fit_PEM, x1=c(1,1), x2=c(1,1), x3=c(1,1), t1=3, t2=6)
		
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ semi-competing risks analysis }
\keyword{ hazard regression models }
