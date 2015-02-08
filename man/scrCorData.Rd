\name{scrCorData}
\alias{scrCorData}
\docType{data}
\title{
A simulated clustered semi-competing risks data set
}
\description{
Simulated semi-competing risks data
}
\usage{data(scrCorData)}
\format{
  A data frame with 2000 observations on the following 14 variables.
  \describe{
    \item{\code{time1}}{the time to non-terminal event}
    \item{\code{event1}}{the censoring indicators for the non-terminal event time; 1=event observed, 0=censored/truncated}
    \item{\code{time2}}{the time to terminal event}
    \item{\code{event2}}{the censoring indicators for the terminal event time; 1=event observed, 0=censored}
    \item{\code{cluster}}{cluster numbers}
    \item{\code{x1.1}}{the first column of covariate matrix x1}
    \item{\code{x1.2}}{the second column of covariate matrix x1}
    \item{\code{x1.3}}{the third column of covariate matrix x1}
    \item{\code{x2.1}}{the first column of covariate matrix x2}
    \item{\code{x2.2}}{the second column of covariate matrix x2}
    \item{\code{x2.3}}{the third column of covariate matrix x2}
    \item{\code{x3.1}}{the first column of covariate matrix x3}
    \item{\code{x3.2}}{the second column of covariate matrix x3}
    \item{\code{x3.3}}{the third column of covariate matrix x3}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(scrCorData)
}
\keyword{datasets}
\keyword{semi-competing risks analysis }
