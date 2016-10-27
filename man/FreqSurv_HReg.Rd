\name{FreqSurv_HReg}
\alias{FreqSurv_HReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to fit parametric Weibull models for the frequentist analysis of univariate survival data.
}
\description{
Independent univariate right-censored survival data can be analyzed using hierarchical models.
}
\usage{
FreqSurv_HReg(Y, lin.pred, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
	a data.frame containing univariate time-to-event outcomes from \code{n} subjects. It is of dimension \eqn{n\times 2}: the columns correspond to \eqn{y}, \eqn{\delta}.
}
  \item{lin.pred}{
	a formula object that corresponds to \eqn{h}.
}
  \item{data}{
	a data.frame in which to interpret the variables named in the formula in \code{lin.pred}.
}
}
\details{
	See \code{\link{BayesSurv_HReg}} for a detailed description of the models.
}
\value{
	\code{FreqSurv_HReg} returns an object of class \code{Freq_HReg}. \cr 
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2015), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, 64, 2, 253-273.\cr \cr
}
\author{
Sebastien Haneuse and Kyu Ha Lee\cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
%\note{
%further notes
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{print.Freq_HReg}}, \code{\link{summary.Freq_HReg}}, \code{\link{plot.Freq_HReg}}, \code{\link{BayesSurv_HReg}}.
}
\examples{
\dontrun{	
# loading a data set	
data(survData)
Y <- survData[,c(1,2)]
lin.pred <- as.formula( ~ cov1 + cov2)

fit_WB <- FreqSurv_HReg(Y, lin.pred, data=survData)
fit_WB
summ.fit_WB <- summary(fit_WB); names(summ.fit_WB)
summ.fit_WB
plot(fit_WB, tseq=seq(from=0, to=30, by=5))
plot(fit_WB, tseq=seq(from=0, to=30, by=5), plot.est = "BH")
names(fit_WB.plot <- plot(fit_WB, tseq=seq(0, 30, 5), plot=FALSE))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ frequentist framework }
\keyword{ univariate analysis }
\keyword{ hazard regression models }% __ONLY ONE__ keyword per line
