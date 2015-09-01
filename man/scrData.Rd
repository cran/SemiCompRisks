\name{scrData}
\alias{scrData}
\docType{data}
\title{
A simulated clustered semi-competing risks data set
}
\description{
Simulated semi-competing risks data
}
\usage{data(scrData)}
\format{
  a data frame with 2000 observations on the following 14 variables.
  \describe{
    \item{\code{time1}}{the time to non-terminal event}
    \item{\code{event1}}{the censoring indicators for the non-terminal event time; 1=event observed, 0=censored/truncated}
    \item{\code{time2}}{the time to terminal event}
    \item{\code{event2}}{the censoring indicators for the terminal event time; 1=event observed, 0=censored}
    \item{\code{cluster}}{cluster numbers}
    \item{\code{x1}}{a vector of continuous covarate}
    \item{\code{x2}}{a vector of continuous covarate}
    \item{\code{x3}}{a vector of continuous covarate}
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
data(scrData)
}
\keyword{datasets}
\keyword{semi-competing risks analysis }
