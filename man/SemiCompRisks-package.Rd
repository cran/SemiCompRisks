\name{SemiCompRisks-package}
\alias{SemiCompRisks-package}
\alias{SemiCompRisks}
\docType{package}
\title{
Algorithms for fitting parametric and semi-parametric models to semi-competing risks data / univariate survival data.
}
\description{
    The package provides functions to perform the analysis of semi-competing risks or univariate survival data with either hazard regression (HReg) models or accelerated failure time (AFT) models. The framework is flexible in the sense that:\cr
    1) it can handle cluster-correlated or independent data; \cr
    2) the option to choose between parametric (Weibull) and semi-parametric (mixture of piecewise exponential) specification for baseline hazard function(s) is available;  \cr
    3) for clustered data, the option to choose between parametric (multivariate Normal for semicompeting risks data, Normal for univariate survival data) and semi-parametric (Dirichlet process mixture) specification for random effects distribution is available;  \cr
    4) for semi-competing risks data, the option to choose between Makov and semi-Makov model is available.  \cr
}
\details{
The package includes following functions:
\tabular{ll}{
\code{BayesID_HReg} \tab Bayesian analysis of semi-competing risks data using HReg models \cr
\code{BayesID_AFT} \tab Bayesian analysis of semi-competing risks data using AFT models \cr
\code{BayesSurv_HReg} \tab Bayesian analysis of univariate survival data using HReg models \cr
\code{BayesSurv_AFT} \tab Bayesian analysis of univariate survival data using AFT models \cr
\code{FreqID_HReg} \tab Frequentist analysis of semi-competing risks data using HReg models \cr
\code{FreqSurv_HReg} \tab Frequentist analysis of univariate survival data using HReg models\cr
\code{initiate.startValues_HReg} \tab Initiating starting values for Bayesian estimations with HReg models \cr
\code{initiate.startValues_AFT} \tab Initiating starting values for Bayesian estimations with AFT models \cr
\code{simID} \tab Simulating semi-competing risks data under Weibull/Weibull-MVN model \cr
\code{simSurv} \tab Simulating survival data under Weibull/Weibull-Normal model \cr
}
\tabular{ll}{
Package: \tab SemiCompRisks\cr
Type: \tab Package\cr
Version: \tab 3.3\cr
Date: \tab 2019-1-27\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Kyu Ha Lee, Catherine Lee, Danilo Alvares, and Sebastien Haneuse \cr
Maintainer: Kyu Ha Lee <klee15239@gmail.com>
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2015), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, 64, 2, 253-273.\cr \cr
Lee, K. H., Dominici, F., Schrag, D., and Haneuse, S. (2016),
Hierarchical models for semicompeting risks data with application to quality of end-of-life care for pancreatic cancer, \emph{Journal of the American Statistical Association}, 111, 515, 1075-1095. \cr \cr
Lee, K. H., Rondeau, V., and Haneuse, S. (2017),
Accelerated failure time models for semicompeting risks data in the presence of complex censoring, \emph{Biometrics}, 73, 4, 1401-1412. \cr \cr
Alvares, D., Haneuse, S., Lee, C., Lee, K. H. (2018+),
SemiCompRisks: an R package for independent and cluster-correlated analyses of semi-competing risks data, \emph{submitted}, arXiv:1801.03567. \cr
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
