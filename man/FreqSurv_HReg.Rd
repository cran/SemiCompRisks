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
FreqSurv_HReg(Formula, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Formula}{
    a \code{Formula} object of the form \eqn{y}+\eqn{\delta} ~ \eqn{x}. See Details and Examples below.
}
\item{data}{
    a data.frame in which to interpret the variables named in \code{Formula}.
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
Alvares, D., Haneuse, S., Lee, C., Lee, K. H. (2018+),
SemiCompRisks: an R package for independent and cluster-correlated analyses of semi-competing risks data, \emph{submitted}, arXiv:1801.03567. \cr
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
	\code{\link{print.Freq_HReg}}, \code{\link{summary.Freq_HReg}}, \code{\link{predict.Freq_HReg}}, \code{\link{BayesSurv_HReg}}.
}
\examples{
\dontrun{	
# loading a data set	
data(survData)

form <- Formula(time + event ~ cov1 + cov2)

fit_WB <- FreqSurv_HReg(form, data=survData)
fit_WB
summ.fit_WB <- summary(fit_WB); names(summ.fit_WB)
summ.fit_WB
pred_WB <- predict(fit_WB, tseq=seq(from=0, to=30, by=5))
plot(pred_WB, plot.est="Haz")
plot(pred_WB, plot.est="Surv")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ frequentist framework }
\keyword{ univariate analysis }
\keyword{ hazard regression models }% __ONLY ONE__ keyword per line
