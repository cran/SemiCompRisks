\name{BayesID_AFT}
\alias{BayesID_AFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to implement Bayesian parametric and semi-parametric analyses for semi-competing risks data in the context of accelerated failure time (AFT) models.
}
\description{
Independent semi-competing risks data can be analyzed using AFT models that have a hierarchical structure. The proposed models can accomodate left-truncated and/or interval-censored data. An efficient computational algorithm that gives users the flexibility to adopt either a fully parametric (log-Normal) or a semi-parametric (Dirichlet process mixture) model specification is developed.
}
\usage{
BayesID_AFT(Y, lin.pred, data, model = "LN", hyperParams, startValues,
mcmcParams, path=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Y}{
a data.frame containing (interval-censored and/or left-truncated) semi-competing risks outcomes from \code{n} subjects. It is of dimension \eqn{n\times 5}: the columns correspond to \eqn{c_{j}}, \eqn{c_{j+1}}, \eqn{c_{k}}, \eqn{c_{k+1}}, \eqn{L}. See Details and Examples below.
}
\item{lin.pred}{
a list containing three formula objects whose right-hand side specifies the covariate terms for the transition \eqn{g}=1,2,3.
}
\item{data}{
a data.frame in which to interpret the variables named in the formulas in \code{lin.pred}.
}
\item{model}{
The specification of baseline survival distribution: "LN" or "DPM".
}
\item{hyperParams}{
a list containing lists or vectors for hyperparameter values in hierarchical models. Components include,
\code{theta} (a numeric vector for hyperparameter in the prior of subject-specific frailty variance component),
\code{LN} (a list containing numeric vectors for log-Normal hyperparameters: \code{LN.ab1}, \code{LN.ab2}, \code{LN.ab3}),
\code{DPM} (a list containing numeric vectors for DPM hyperparameters: \code{DPM.mu1}, \code{DPM.mu2}, \code{DPM.mu3}, \code{DPM.sigSq1}, \code{DPM.sigSq2}, \code{DPM.sigSq3}, \code{DPM.ab1}, \code{DPM.ab2}, \code{DPM.ab3}, \code{Tau.ab1}, \code{Tau.ab2}, \code{Tau.ab3}).
See Details and Examples below.
}
\item{startValues}{
a list containing vectors of starting values for model parameters. It can be specified as the object returned by the function \code{\link{initiate.startValues_AFT}}.
}
\item{mcmcParams}{
a list containing variables required for MCMC sampling. Components include,
\code{run} (a list containing numeric values for setting for the overall run:
\code{numReps}, total number of scans; \code{thin}, extent of thinning; \code{burninPerc}, the proportion of burn-in).
\code{storage} (a list containing numeric values for storing posterior samples for subject- and cluster-specific random effects:
\code{nGam_save}, the number of \eqn{\gamma} to be stored; \code{nY1_save}, the number of \eqn{y1} to be stored; \code{nY2_save}, the number of \eqn{y2} to be stored; \code{nY1.NA_save}, the number of \eqn{y1.NA} to be stored).
\code{tuning} (a list containing numeric values relevant to tuning parameters for specific updates in Metropolis-Hastings (MH) algorithm:
\code{mhProp_theta_var}, the variance of proposal density for \eqn{\theta}; \code{betag.prop.var}, the variance of proposal density for \eqn{\beta_g}; \code{mug.prop.var}, the variance of proposal density for \eqn{\mu_{g}}; \code{zetag.prop.var}, the variance of proposal density for \eqn{1/\sigma_g^2};
\code{gamma.prop.var}, the variance of proposal density for \eqn{\gamma}).
See Details and Examples below.
}
\item{path}{
the name of directory where the results are saved.
}
}
\details{
We view the semi-competing risks data as arising from an underlying illness-death model system in which individuals may undergo one or more of three transitions: 1) from some initial condition to non-terminal event, 2) from some initial condition to terminal event, 3) from non-terminal event to terminal event. Let \eqn{T_{i1}}, \eqn{T_{i2}} denote time to non-terminal and terminal event from subject \eqn{i=1,...,n}.  We propose to directly model the times of the events via the following AFT model specification:

\deqn{\log(T_{i1}) = x_{i1}^\top\beta_1 + \gamma_i + \epsilon_{i1},  T_{i1} > 0,}
\deqn{\log(T_{i2}) = x_{i2}^\top\beta_2 + \gamma_i + \epsilon_{i2},  T_{i2} > 0,}
\deqn{\log(T_{i2} - T_{i1}) = x_{i3}^\top\beta_3 + \gamma_i + \epsilon_{i3},  T_{i2} >  T_{i1},}
where \eqn{x_{ig}} is a vector of transition-specific covariates, \eqn{\beta_g} is a corresponding vector of transition-specific regression parameters and \eqn{\epsilon_{ig}} is a transition-specific random variable whose distribution determines that of the corresponding transition time, \eqn{g \in \{1,2,3\}}. \eqn{\gamma_i} is a study participant-specific random effect that induces positive dependence between the two event times, thereby performing a role analogous to that performed by frailties in models for the hazard function.
Let \eqn{L_{i}} denote the time at study entry (i.e. the left-truncation time). Furthermore, suppose that study participant \eqn{i} was observed at follow-up times \eqn{\{c_{i1},\ldots, c_{im_i}\}} and let \eqn{c_i^*} denote the time to the end of study or to administrative right-censoring. Considering interval-censoring for both events, the times to non-terminal and terminal event for the \eqn{i^{th}} study participant satisfy \eqn{c_{ij}\leq T_{i1}< c_{ij+1}} for some \eqn{j} and \eqn{c_{ik}\leq T_{i2}< c_{ik+1}} for some \eqn{k}, respectively. Then the observed outcomes for the \eqn{i^{th}} study participant can be succinctly denoted by \eqn{\{c_{ij}, c_{ij+1}, c_{ik}, c_{ik+1}, L_{i}\}}.

For the Bayesian semi-parametric analysis, we proceed by adopting independent DPM of normal distributions for each \eqn{\epsilon_{ig}}. More precisely, \eqn{\epsilon_{ig}} is taken to be an independent draw from a mixture of \eqn{M_g} normal distributions with means and variances (\eqn{\mu_{gr}}, \eqn{\sigma_{gr}^2}), for \eqn{r \in \{1,\ldots,M_g\}}. Since the class-specific \eqn{(\mu_{gr}, \sigma_{gr}^2)} are not known, they are taken to be draws from some common distribution, \eqn{G_{g0}}, often referred to as the centering distribution. Furthermore, since the `true' class membership for any given study participant is not known, we let \eqn{p_{gr}} denote the probability of belonging to the \eqn{r^{th}} class for transition \eqn{g} and \eqn{p_g} = \eqn{(p_{g1}, \ldots, p_{gM_g})} the collection of such probabilities. Note, \eqn{p_g} is defined at the level of the population (i.e. is not study participant-specific) and its components add up to 1.0. In the absence of prior knowledge regarding the distribution of class memberships for the \eqn{n} individuals across the \eqn{M_g} classes, \eqn{p_g} is assumed to follow a conjugate symmetric Dirichlet\eqn{(\tau_g/M_g,\ldots,\tau_g/M_g)} distribution, where \eqn{\tau_g} is referred to as the precision parameter. The finite mixture distribution can then be succinctly represented as:
\deqn{\epsilon_{ig} | r_{i} \sim Normal(\mu_{r_{i}}, \sigma_{r_{i}}^2),}
\deqn{(\mu_{gr}, \sigma_{gr}^2) \sim G_{g0}, ~~for~ r=1,\ldots,M_g,}
\deqn{r_{i}| p_g \sim Discrete(r_{i} | p_{g1},\ldots,p_{gM_g}),}
\deqn{p_g \sim Dirichlet(\tau_g/M_g, \ldots, \tau_g/M_g).}
Letting \eqn{M_g} approach infinity, this specification is referred to as a DPM of normal distributions. In our proposed framework, we specify a Gamma(\eqn{a_{\tau_g}}, \eqn{b_{\tau_g}}) hyperprior for \eqn{\tau_g}. For regression parameters, we adopt non-informative flat priors on the real line. For \eqn{\gamma}=\eqn{\{\gamma_1, \ldots, \gamma_n\}}, we assume that each \eqn{\gamma_i} is an independent random draw from a Normal(0, \eqn{\theta}) distribution. In the absence of prior knowledge on the variance component \eqn{\theta}, we adopt a conjugate inverse-Gamma hyperprior, IG(\eqn{a_\theta}, \eqn{b_\theta}). Finally, We take the \eqn{G_{g0}} as a normal distribution centered at \eqn{\mu_{g0}} with a variance \eqn{\sigma_{g0}^2} for \eqn{\mu_{gr}} and an IG(\eqn{a_{\sigma_g}}, \eqn{b_{\sigma_g}}) for \eqn{\sigma_{gr}^2}.

For the Bayesian parametric analysis, we build on the log-Normal formulation and take the \eqn{\epsilon_{ig}} to follow independent Normal(\eqn{\mu_g}, \eqn{\sigma_g^2}) distributions, \eqn{g}=1,2,3. For location parameters \eqn{\{\mu_1, \mu_2, \mu_3\}}, we adopt non-informative flat priors on the real line. For \eqn{\{\sigma_1^2, \sigma_2^2, \sigma_3^2\}}, we adopt independent inverse Gamma distributions, denoted IG(\eqn{a_{\sigma g}}, \eqn{b_{\sigma g}}). For \eqn{\beta_g}, \eqn{\gamma}, and \eqn{\theta}, we adopt the same priors as those adopted for the DPM model.
}
\value{
\code{BayesID_AFT} returns an object of class \code{Bayes_AFT}. \cr
}
\references{
Lee, K. H., Rondeau, V., and Haneuse, S. (2016),
Accelerated failure time models for semicompeting risks data in the presence of complex censoring, \emph{under review}. \cr}
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
\code{\link{initiate.startValues_AFT}}, \code{\link{print.Bayes_AFT}}, \code{\link{summary.Bayes_AFT}}, \code{\link{plot.Bayes_AFT}}
}
\examples{

\dontrun{
# loading a data set
data(scrData)
Y <- matrix(NA, dim(scrData)[1], 5)
Y[,1] <- Y[,2] <- scrData[,1]
Y[,3] <- Y[,4] <- scrData[,3]
Y[which(scrData[,2] == 0),2] <- Inf
Y[which(scrData[,4] == 0),4] <- Inf
Y[,5] <- rep(0, dim(scrData)[1])

form1 <- as.formula( ~ x1 + x2 + x3)
form2 <- as.formula( ~ x1 + x2)
form3 <- as.formula( ~ x1 + x2)
lin.pred <- list(form1, form2, form3)

#####################
## Hyperparameters ##
#####################

## Subject-specific random effects variance component
##
theta.ab <- c(0.5, 0.05)

## log-Normal model
##
LN.ab1 <- c(0.3, 0.3)
LN.ab2 <- c(0.3, 0.3)
LN.ab3 <- c(0.3, 0.3)

## DPM model
##
DPM.mu1 <- log(12)
DPM.mu2 <- log(12)
DPM.mu3 <- log(12)

DPM.sigSq1 <- 100
DPM.sigSq2 <- 100
DPM.sigSq3 <- 100

DPM.ab1 <-  c(2, 1)
DPM.ab2 <-  c(2, 1)
DPM.ab3 <-  c(2, 1)

Tau.ab1 <- c(1.5, 0.0125)
Tau.ab2 <- c(1.5, 0.0125)
Tau.ab3 <- c(1.5, 0.0125)

##
hyperParams <- list(theta=theta.ab,
LN=list(LN.ab1=LN.ab1, LN.ab2=LN.ab2, LN.ab3=LN.ab3),
DPM=list(DPM.mu1=DPM.mu1, DPM.mu2=DPM.mu2, DPM.mu3=DPM.mu3, DPM.sigSq1=DPM.sigSq1,
DPM.sigSq2=DPM.sigSq2, DPM.sigSq3=DPM.sigSq3, DPM.ab1=DPM.ab1, DPM.ab2=DPM.ab2,
DPM.ab3=DPM.ab3, Tau.ab1=Tau.ab1, Tau.ab2=Tau.ab2, Tau.ab3=Tau.ab3))

###################
## MCMC SETTINGS ##
###################

## Setting for the overall run
##
numReps    <- 300
thin       <- 3
burninPerc <- 0.5

## Setting for storage
##
nGam_save <- 10
nY1_save <- 10
nY2_save <- 10
nY1.NA_save <- 10

## Tuning parameters for specific updates
##
##  - those common to all models
betag.prop.var	<- c(0.01,0.01,0.01)
mug.prop.var	<- c(0.1,0.1,0.1)
zetag.prop.var	<- c(0.1,0.1,0.1)
gamma.prop.var	<- 0.01

##
mcmcParams	<- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
storage=list(nGam_save=nGam_save, nY1_save=nY1_save, nY2_save=nY2_save, nY1.NA_save=nY1.NA_save),
tuning=list(betag.prop.var=betag.prop.var, mug.prop.var=mug.prop.var,
zetag.prop.var=zetag.prop.var, gamma.prop.var=gamma.prop.var))

#################################################################
## Analysis of Independent Semi-competing risks data ############
#################################################################

###############
## logNormal ##
###############

##
myModel <- "LN"
myPath  <- "Output/01-Results-LN/"

startValues      <- vector("list", 2)
startValues[[1]] <- initiate.startValues_AFT(Y, lin.pred, scrData, model=myModel)
startValues[[2]] <- initiate.startValues_AFT(Y, lin.pred, scrData, model=myModel, theta = 0.20)

##
fit_LN <- BayesID_AFT(Y, lin.pred, scrData, model=myModel, hyperParams,
startValues, mcmcParams, path=myPath)

fit_LN
summ.fit_LN <- summary(fit_LN); names(summ.fit_LN)
summ.fit_LN
plot(fit_LN, time = seq(0, 35, 1), tseq=seq(from=0, to=30, by=5))
plot(fit_LN, time = seq(0, 35, 1), tseq=seq(from=0, to=30, by=5), plot.est = "BH")
names(fit_LN.plot <- plot(fit_LN, time = seq(0, 35, 1),  tseq=seq(0, 30, 5), plot=FALSE))


#########
## DPM ##
#########

##
myModel <- "DPM"
myPath  <- "Output/02-Results-DPM/"

startValues      <- vector("list", 2)
startValues[[1]] <- initiate.startValues_AFT(Y, lin.pred, scrData, model=myModel)
startValues[[2]] <- initiate.startValues_AFT(Y, lin.pred, scrData, model=myModel, theta = 0.23)

##
fit_DPM <- BayesID_AFT(Y, lin.pred, scrData, model=myModel, hyperParams,
startValues, mcmcParams, path=myPath)

fit_DPM
summ.fit_DPM <- summary(fit_DPM); names(summ.fit_DPM)
summ.fit_DPM
plot(fit_DPM, time = seq(0, 35, 1),  tseq=seq(from=0, to=30, by=5))
plot(fit_DPM, time = seq(0, 35, 1),  tseq=seq(from=0, to=30, by=5), plot.est = "BH")
names(fit_DPM.plot <- plot(fit_DPM, time = seq(0, 35, 1), tseq=seq(0, 30, 5), plot=FALSE))
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ semi-competing risks analysis }
\keyword{ accelerated failure time models }% __ONLY ONE__ keyword per line
