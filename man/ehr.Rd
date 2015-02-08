\name{ehr}
\alias{ehr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to calculate the conditional explanatory hazard ratio (EHR)
}
\description{
The function to calculate the conditional explanatory hazard ratio (EHR)
}
\usage{
ehr(Cov2 = NULL, Cov3 = NULL, x, xlim = c(0, 10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Cov2}{
	a covariate vector of length \code{p2} that corresponds to transition "from no event to terminal event"
}
  \item{Cov3}{
	a covariate vector of length \code{p3} that corresponds to transition "from non-terminal event to terminal event"
}
  \item{x}{
	an object of class \code{BayesID}
}
  \item{xlim}{
	the time (x-axis) limits (x1, x2) of the plot.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{ehr} returns an object of class \code{ehr}. \cr
\item{ehr}{estimated conditional EHR (posterior median) at \code{time} points}
\item{time}{time points at which EHR is estimated}
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2014), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, \emph{in press}.
}
\author{
Kyu Ha Lee and Sebastien Haneuse \cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{BayesID}}, \code{\link{plot.ehr}}
}




% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line
