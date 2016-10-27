\name{initiate.startValues_HReg}
\alias{initiate.startValues_HReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function that initiates starting values for a single chain.
}
\description{
The function initiates starting values for a single chain for hazard regression (HReg) models. Users are allowed to set some non-null values to starting values for a set of parameters. The function will automatically generate starting values for any parameters whose values are not specified.
}
\usage{
initiate.startValues_HReg(Y, lin.pred, data, model, cluster = NULL, 
                   beta1 = NULL, beta2 = NULL, beta3 = NULL, beta = NULL,
                   gamma.ji = NULL, theta = NULL,
                   V.j1 = NULL, V.j2 = NULL, V.j3 = NULL, V.j = NULL,
                   WB.alpha = NULL, WB.kappa = NULL, 
                   PEM.lambda1=NULL, PEM.lambda2=NULL, PEM.lambda3=NULL, PEM.lambda=NULL,
                   PEM.s1=NULL, PEM.s2=NULL, PEM.s3=NULL, PEM.s=NULL,
                   PEM.mu_lam=NULL, PEM.sigSq_lam=NULL,
                   MVN.SigmaV = NULL, Normal.zeta = NULL, 
                   DPM.class = NULL, DPM.tau = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
	For \code{BayesID_HReg}, it is a data.frame containing semi-competing risks outcomes from \code{n} subjects.
	For \code{BayesSurv_HReg}, it is a data.frame containing univariate time-to-event outcomes from \code{n} subjects.
}
  \item{lin.pred}{
	For \code{BayesID_HReg}, it is a list containing three formula objects that correspond to \eqn{h_g()}, \eqn{g}=1,2,3.
	For \code{BayesSurv_HReg}, it is a formula object that corresponds to \eqn{h()}.
}
  \item{data}{
	a data.frame in which to interpret the variables named in the formula(s) in \code{lin.pred}.
}
  \item{model}{
	a character vector that specifies the type of components in a model. Check \code{\link{BayesID_HReg}} and \code{\link{BayesSurv_HReg}}.
}
  \item{cluster}{
	a vector of cluster information for \code{n} subjects. The cluster membership must be set to consecutive positive integers, \eqn{1:J}.  
}
  \item{beta1}{
	starting values of \eqn{\beta_1} for \code{BayesID_HReg}.
}
  \item{beta2}{
	starting values of \eqn{\beta_2} for \code{BayesID_HReg}.
}
  \item{beta3}{
	starting values of \eqn{\beta_3} for \code{BayesID_HReg}.
}
  \item{beta}{
	starting values of \eqn{\beta} for \code{BayesSurv_HReg}.
}
  \item{gamma.ji}{
	starting values of \eqn{\gamma} for \code{BayesID_HReg}.
}
  \item{theta}{
	starting values of \eqn{\theta} for \code{BayesID_HReg}.
}
  \item{V.j1}{
	starting values of \eqn{V_{j1}} for \code{BayesID_HReg}.
}
  \item{V.j2}{
	starting values of \eqn{V_{j2}} for \code{BayesID_HReg}.
}
  \item{V.j3}{
	starting values of \eqn{V_{j3}} for \code{BayesID_HReg}.
}
  \item{V.j}{
	starting values of \eqn{V_{j}} for \code{BayesSurv_HReg}.
}
  \item{WB.alpha}{
	starting values of the Weibull parameters, \eqn{\alpha_g} for \code{BayesID_HReg}.
	starting values of the Weibull parameter, \eqn{\alpha} for \code{BayesSurv_HReg}.
}
  \item{WB.kappa}{
	starting values of the Weibull parameters, \eqn{\kappa_g} for \code{BayesID_HReg}.
	starting values of the Weibull parameter, \eqn{\kappa} for \code{BayesSurv_HReg}.
}
  \item{PEM.lambda1}{
	starting values of the PEM parameters, \eqn{\lambda_1} for \code{BayesID_HReg}.
}
  \item{PEM.lambda2}{
	starting values of the PEM parameters, \eqn{\lambda_2} for \code{BayesID_HReg}.
}
  \item{PEM.lambda3}{
	starting values of the PEM parameters, \eqn{\lambda_3} for \code{BayesID_HReg}.
}
  \item{PEM.lambda}{
	starting values of \eqn{\lambda} for \code{BayesSurv_HReg}.
}
  \item{PEM.s1}{
	starting values of the PEM parameters, \eqn{s_1} for \code{BayesID_HReg}.
}
  \item{PEM.s2}{
	starting values of the PEM parameters, \eqn{s_2} for \code{BayesID_HReg}.
}
  \item{PEM.s3}{
	starting values of the PEM parameters, \eqn{s_3} for \code{BayesID_HReg}.
}
  \item{PEM.s}{
	starting values of \eqn{s} for \code{BayesSurv_HReg}.
}
  \item{PEM.mu_lam}{
	starting values of the PEM parameters, \eqn{\mu_{\lambda,g}} for \code{BayesID_HReg}.
	starting values of the PEM parameter, \eqn{\mu_{\lambda}} for \code{BayesSurv_HReg}.
}
  \item{PEM.sigSq_lam}{
	starting values of the PEM parameters, \eqn{\sigma_{\lambda,g}^2} for \code{BayesID_HReg}.
	starting values of the PEM parameter, \eqn{\sigma_{\lambda}^2} for \code{BayesSurv_HReg}.
}
  \item{MVN.SigmaV}{
	starting values of \eqn{\Sigma_V} in DPM models for \code{BayesID_HReg}.
}
  \item{Normal.zeta}{
	starting values of \eqn{\zeta} in DPM models for \code{BayesSurv_HReg}.
}
  \item{DPM.class}{
	starting values of the class membership in DPM models for \code{BayesID_HReg} and \code{BayesSurv_HReg}.
}
  \item{DPM.tau}{
	starting values of \eqn{\tau} in DPM models for \code{BayesID_HReg} and \code{BayesSurv_HReg}.
}
}
%\details{
%	
%}
\value{
	\code{initiate.startValues_HReg} returns a list containing starting values for a sigle chain that can be used for \code{BayesID_HReg} and \code{BayesSurv_HReg}.
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2015), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, 64, 2, 253-273.\cr \cr
Lee, K. H., Dominici, F., Schrag, D., and Haneuse, S. (2016),
Hierarchical models for semicompeting risks data with application to quality of end-of-life care for pancreatic cancer, \emph{Journal of the American Statistical Association}, 111, 515, 1075-1095. \cr
}
\author{
Sebastien Haneuse and Kyu Ha Lee\cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
%\note{
%Function that initiates starting values for a single chain
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{BayesID_HReg}}, \code{\link{BayesSurv_HReg}}
}
\examples{
## See Examples in \code{\link{BayesID_HReg}} and \code{\link{BayesSurv_HReg}}.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line
\keyword{ univariate analysis }
\keyword{ hazard regression models }% __ONLY ONE__ keyword per line
