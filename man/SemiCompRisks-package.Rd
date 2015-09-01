\name{SemiCompRisks-package}
\alias{SemiCompRisks-package}
\alias{SemiCompRisks}
\docType{package}
\title{
Algorithms for fitting parametric and semi-parametric models to semi-competing risks data / univariate survival data.
}
\description{
    The package provides functions to perform the analysis of semi-competing risks or univariate survival data. The framework is flexible in the sense that:\cr    
    1) it can handle cluster-correlated or independent data; \cr
    2) the option to choose between parametric (Weibull) and semi-parametric (mixture of piecewise exponential) specification for baseline hazard function(s) is available;  \cr
    3) for clustered data, the option to choose between parametric (multivariate Normal for semicompeting risks data, Normal for univariate survival data) and semi-parametric (Dirichlet process mixture) specification for random effects distribution is available;  \cr
    4) for semi-competing risks data, the option to choose between Makov and semi-Makov model is available.  \cr
}
\details{
The package includes following functions:
\tabular{ll}{
\code{BayesID} \tab Bayesian analysis of semi-competing risks data  \cr
\code{BayesSurv} \tab Bayesian analysis of univariate survival data \cr
\code{FreqID} \tab Frequentist analysis of semi-competing risks data  \cr
\code{FreqSurv} \tab Frequentist analysis of univariate survival data \cr
\code{initiate.startValues} \tab Initiating starting values for Bayesian estimations \cr
\code{simID} \tab Simulating semi-competing risks data under Weibull/Weibull-MVN model \cr
\code{simSurv} \tab Simulating survival data under Weibull/Weibull-Normal model \cr
}
\tabular{ll}{
Package: \tab SemiCompRisks\cr
Type: \tab Package\cr
Version: \tab 2.2\cr
Date: \tab 2015-8-31\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Kyu Ha Lee and Sebastien Haneuse \cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2015), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, 64, 2, 253-273.\cr \cr
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
