\name{SemiCompRisks-package}
\alias{SemiCompRisks-package}
\alias{SemiCompRisks}
\docType{package}
\title{
Algorithms for fitting parametric and semi-parametric models to semi-competing risks data / univariate survival data.
}
\description{
The package provides following Bayesian frameworks to analyze semi-competing risks/univariate survival data: \cr \cr
1) Bayesian semi-parametric model for analysis of semi-competing risks data, \cr
2) Bayesian parametric model for analysis of semi-competing risks data, \cr
3) Bayesian semi-parametric model for analysis of univariate right censored survival data,\cr
4) Bayesian parametric model for analysis of univariate right censored survival data \cr
}
\details{
The package includes following functions:
\tabular{ll}{
\code{BayesID} \tab The function to fit Bayesian illness-death models to semi-competing risks data  \cr
\code{BayesSurv} \tab The function to fit Bayesian parametric and semi-parametric regression models \cr
				\tab to univariate survival data \cr
\code{ehr} \tab  The function to calculate the conditional explanatory hazard ratio (EHR) \cr
\code{simID} \tab The function to simulate semi-competing risks data under Weibull model \cr
\code{simSurv} \tab The function to simulate right censored survival data under Weibull model \cr
}
\tabular{ll}{
Package: \tab SemiCompRisks\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-21\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Kyu Ha Lee and Sebastien Haneuse \cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2013). 
Bayesian Semi-parametric Analysis of Semi-competing Risks Data: 
Estimating Readmission Rates among Pancreatic Cancer Patients, submitted.
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
