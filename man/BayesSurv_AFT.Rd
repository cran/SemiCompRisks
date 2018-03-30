\name{BayesSurv_AFT}
\alias{BayesSurv_AFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to implement Bayesian parametric and semi-parametric analyses for univariate survival data in the context of accelerated failure time (AFT) models.
}
\description{
Independent univariate survival data can be analyzed using AFT models that have a hierarchical structure. The proposed models can accomodate left-truncated and/or interval-censored data. An efficient computational algorithm that gives users the flexibility to adopt either a fully parametric (log-Normal) or a semi-parametric (Dirichlet process mixture) model specification is developed.
}
\usage{
BayesSurv_AFT(Formula, data, model = "LN", hyperParams, startValues,
                mcmcParams, path = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Formula}{
    a \code{Formula} object of the form \eqn{L} | \eqn{y_{L}}+\eqn{y_{U}} ~ \eqn{x}. See Details and Examples below.
}
  \item{data}{
    a data.frame in which to interpret the variables named in \code{Formula}.
}
  \item{model}{
    The specification of baseline survival distribution: "LN" or "DPM".
}
  \item{hyperParams}{
a list containing lists or vectors for hyperparameter values in hierarchical models. Components include,
\code{LN} (a list containing numeric vectors for log-Normal hyperparameters: \code{LN.ab}),
\code{DPM} (a list containing numeric vectors for DPM hyperparameters: \code{DPM.mu}, \code{DPM.sigSq}, \code{DPM.ab}, \code{Tau.ab}).
See Details and Examples below.
}
  \item{startValues}{
  	a list containing vectors of starting values for model parameters. It can be specified as the object returned by the function \code{\link{initiate.startValues_AFT}}.
}
  \item{mcmcParams}{
a list containing variables required for MCMC sampling. Components include,
\code{run} (a list containing numeric values for setting for the overall run:
\code{numReps}, total number of scans; \code{thin}, extent of thinning; \code{burninPerc}, the proportion of burn-in).
\code{tuning} (a list containing numeric values relevant to tuning parameters for specific updates in Metropolis-Hastings (MH) algorithm:
\code{beta.prop.var}, the variance of proposal density for \eqn{\beta}; \code{mu.prop.var}, the variance of proposal density for \eqn{\mu}; \code{zeta.prop.var}, the variance of proposal density for \eqn{1/\sigma^2}).

}
  \item{path}{
    the name of directory where the results are saved.
}
}
\details{
The function \code{BayesSurv_AFT} implements Bayesian semi-parametric (DPM) and parametric (log-Normal) models to univariate time-to-event data in the presence of left-truncation and/or interval-censoring. Consider a univariate AFT model that relates the covariate \eqn{x_i} to survival time \eqn{T_i} for the \eqn{i^{\textrm{th}}} subject:
\deqn{\log(T_i) =  x_i^{\top}\beta + \epsilon_i,}
where \eqn{\epsilon_i} is a random variable whose distribution determines that of \eqn{T_i} and \eqn{\beta} is a vector of regression parameters. Considering the interval censoring, the time to the event for the \eqn{i^{\textrm{th}}} subject satisfies \eqn{c_{ij}\leq T_i <c_{ij+1}}. Let \eqn{L_i} denote the left-truncation time.
For the Bayesian parametric analysis, we take  \eqn{\epsilon_i} to follow the Normal(\eqn{\mu}, \eqn{\sigma^2}) distribution for \eqn{\epsilon_i}. The following prior distributions are adopted for the model parameters:
\deqn{\pi(\beta, \mu) \propto 1,}
\deqn{\sigma^2 \sim \textrm{Inverse-Gamma}(a_{\sigma}, b_{\sigma}).}

For the Bayesian semi-parametric analysis, we assume that \eqn{\epsilon_i} is taken as draws from the DPM of normal distributions:
\deqn{\epsilon\sim DPM(G_0, \tau).}
We refer readers to \code{\link{print.Bayes_AFT}} for a detailed illustration of DPM specification. We adopt a non-informative flat prior on the real line for the regression parameters \eqn{\beta} and a Gamma(\eqn{a_{\tau}}, \eqn{b_{\tau}}) hyperprior for the precision parameter \eqn{\tau}.
}
\value{
\code{BayesSurv_AFT} returns an object of class \code{Bayes_AFT}. \cr
}
\references{
Lee, K. H., Rondeau, V., and Haneuse, S. (2017),
Accelerated failure time models for semicompeting risks data in the presence of complex censoring, \emph{Biometrics}, 73, 4, 1401-1412. \cr \cr
Alvares, D., Haneuse, S., Lee, C., Lee, K. H. (2018+),
SemiCompRisks: an R package for independent and cluster-correlated analyses of semi-competing risks data, \emph{submitted}, arXiv:1801.03567. \cr
}
\author{
Kyu Ha Lee and Sebastien Haneuse\cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
%\note{
%TBU
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{initiate.startValues_AFT}}, \code{\link{print.Bayes_AFT}}, \code{\link{summary.Bayes_AFT}}, \code{\link{predict.Bayes_AFT}}
}
\examples{

\dontrun{
# loading a data set
data(survData)
survData$yL <- survData$yU <- survData[,1]
survData$yU[which(survData[,2] == 0)] <- Inf
survData$LT <- rep(0, dim(survData)[1])

form <- Formula(LT | yL + yU ~ cov1 + cov2)

#####################
## Hyperparameters ##
#####################

## log-Normal model
##
LN.ab <- c(0.3, 0.3)

## DPM model
##
DPM.mu <- log(12)
DPM.sigSq <- 100
DPM.ab <-  c(2, 1)
Tau.ab <- c(1.5, 0.0125)

##
hyperParams <- list(LN=list(LN.ab=LN.ab),
DPM=list(DPM.mu=DPM.mu, DPM.sigSq=DPM.sigSq, DPM.ab=DPM.ab, Tau.ab=Tau.ab))

###################
## MCMC SETTINGS ##
###################

## Setting for the overall run
##
numReps    <- 100
thin       <- 1
burninPerc <- 0.5

## Tuning parameters for specific updates
##
##  - those common to all models
beta.prop.var	<- 0.01
mu.prop.var	<- 0.1
zeta.prop.var	<- 0.1

##
mcmcParams	<- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
tuning=list(beta.prop.var=beta.prop.var, mu.prop.var=mu.prop.var,
zeta.prop.var=zeta.prop.var))

################################################################
## Analysis of Independent univariate survival data ############
################################################################

###############
## logNormal ##
###############

##
myModel <- "LN"
myPath  <- "Output/01-Results-LN/"

startValues      <- initiate.startValues_AFT(form, survData, model=myModel, nChain=2)

##
fit_LN <- BayesSurv_AFT(form, survData, model=myModel, hyperParams,
startValues, mcmcParams, path=myPath)

fit_LN
summ.fit_LN <- summary(fit_LN); names(summ.fit_LN)
summ.fit_LN
pred_LN <- predict(fit_LN, time = seq(0, 35, 1), tseq=seq(from=0, to=30, by=5))
plot(pred_LN, plot.est="Haz")
plot(pred_LN, plot.est="Surv")

#########
## DPM ##
#########

##
myModel <- "DPM"
myPath  <- "Output/02-Results-DPM/"

startValues      <- initiate.startValues_AFT(form, survData, model=myModel, nChain=2)

##
fit_DPM <- BayesSurv_AFT(form, survData, model=myModel, hyperParams,
startValues, mcmcParams, path=myPath)

fit_DPM
summ.fit_DPM <- summary(fit_DPM); names(summ.fit_DPM)
summ.fit_DPM
pred_DPM <- predict(fit_DPM, time = seq(0, 35, 1), tseq=seq(from=0, to=30, by=5))
plot(pred_DPM, plot.est="Haz")
plot(pred_DPM, plot.est="Surv")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ univariate analysis }
\keyword{ accelerated failure time models }% __ONLY ONE__ keyword per line

