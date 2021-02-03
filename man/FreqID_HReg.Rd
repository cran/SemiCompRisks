\name{FreqID_HReg}
\alias{FreqID_HReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to fit parametric Weibull models for the frequentist anlaysis of semi-competing risks data.
}
\description{
Independent semi-competing risks data can be analyzed using hierarchical models. Markov or semi-Markov assumption can be adopted for the conditional hazard function for time to the terminal event given time to non-terminal event.
}
\usage{
FreqID_HReg(Formula, data, model="semi-Markov", frailty = TRUE, na.action = "na.fail",
subset=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Formula}{
    a \code{Formula} object, with the outcome on the left of a \eqn{\sim}, and covariates on the right. It is of the form, time to non-terminal event + corresponding censoring indicator | time to terminal event + corresponding censoring indicator \eqn{\sim} covariates for \eqn{h_1} | covariates for \eqn{h_2} | covariates for \eqn{h_3}: i.e., \eqn{y_1}+\eqn{\delta_1} | \eqn{y_2}+\eqn{\delta_2} ~ \eqn{x_1} | \eqn{x_2} | \eqn{x_3}. 
}
  \item{data}{
	a data.frame in which to interpret the variables named in \code{Formula}.
}
\item{model}{
	a character value that specifies the type of a model based on the assumption on \eqn{h_3}: "semi-Markov" or "Markov". 
}
  \item{frailty}{
	a logical value to determine whether to include the subject-specific shared frailty term, \eqn{\gamma}, into the model.
}
\item{na.action}{
how NAs are treated. See \code{model.frame}.
}
\item{subset}{
a specification of the rows to be used: defaults to all rows. See \code{model.frame}.
}
}
\details{
	See \code{\link{BayesID_HReg}} for a detailed description of the models.
}
\value{
	\code{FreqID_HReg} returns an object of class \code{Freq_HReg}. \cr
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2015), 
Bayesian semiparametric analysis of semicompeting risks data: 
investigating hospital readmission after a pancreatic cancer diagnosis, \emph{Journal of the Royal Statistical Society: Series C}, 64, 2, 253-273.\cr \cr
Alvares, D., Haneuse, S., Lee, C., Lee, K. H. (2019),
SemiCompRisks: An R package for the analysis of independent and cluster-correlated semi-competing risks data, \emph{The R Journal}, 11, 1, 376-400. \cr
}
\author{
Sebastien Haneuse and Kyu Ha Lee\cr
Maintainer: Kyu Ha Lee <klee15239@gmail.com>
}
%%\note{
%%further notes
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{print.Freq_HReg}}, \code{\link{summary.Freq_HReg}}, \code{\link{predict.Freq_HReg}}, \code{\link{BayesID_HReg}}.
}
\examples{
\dontrun{
# loading a data set
data(scrData)

form <- Formula(time1 + event1 | time2 + event2 ~ x1 + x2 + x3 | x1 + x2 | x1 + x2)

fit_WB	<- FreqID_HReg(form, data=scrData, model="semi-Markov")

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
\keyword{ semi-competing risks analysis }
\keyword{ hazard regression models }% __ONLY ONE__ keyword per line
