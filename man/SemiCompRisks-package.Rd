\name{SemiCompRisks-package}
\alias{SemiCompRisks-package}
\alias{SemiCompRisks}
\docType{package}
\title{
Algorithms for fitting parametric and semi-parametric models to semi-competing risks data / univariate survival data.
}
\description{
    The package provides a Bayesian framework to perform the analysis of semi-competing risks or univariate survival data. The framework is flexible in the sense that:\cr    
    1) it can handle cluster-correlated or independent data; \cr
    2) the option to choose between parametric (Weibull) and semi-parametric (mixture of piecewise exponential) prior specification for baseline hazard function(s) is available;  \cr
    3) for clustered data, the option to choose between parametric (multivariate Normal for semicompeting risks data, Normal for univariate survival data) and semi-parametric (Dirichlet process mixture) prior specification for random effects distribution is available;  \cr
    4) for semi-competing risks data, the option to choose between Makov and semi-Makov model is available.  \cr
}
\details{
The package includes following functions:
\tabular{ll}{
\code{BayesID} \tab Bayesian analysis of independent semi-competing risks data  \cr
\code{BayesIDcor} \tab Bayesian analysis of cluster-correlated semi-competing risks data  \cr
\code{BayesSurv} \tab Bayesian analysis of independent univariate survival data \cr
\code{BayesSurvcor} \tab Bayesian analysis of cluster-correlated univariate survival data \cr
\code{ehr} \tab  The function to calculate the conditional explanatory hazard ratio (EHR) \cr
\code{simID} \tab The function to simulate semi-competing risks data under Weibull model \cr
\code{simSurv} \tab The function to simulate right censored survival data under Weibull model \cr
}
\tabular{ll}{
Package: \tab SemiCompRisks\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-1-21\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Kyu Ha Lee and Sebastien Haneuse \cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2014), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, \emph{in press}.\cr \cr
Lee, K. H., Dominici, F., Schrag, D., and Haneuse, S.,
Hierarchical models for cluster-correlated semicompeting risks data with application to quality of end-of-life care for pancreatic cancer, \emph{arXiv:1502.00526}; submitted. \cr
}
%Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%documentation directory ~~
\keyword{ package }
%\seealso{
%Optional links to other man pages, e.g. ~~
%\code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%
%}
