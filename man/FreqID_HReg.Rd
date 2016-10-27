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
FreqID_HReg(Y, lin.pred, data, model="semi-Markov", frailty = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
	a data.frame containing semi-competing risks outcomes from \code{n} subjects. It is of dimension \eqn{n\times 4}: the columns correspond to \eqn{y_1}, \eqn{\delta_1}, \eqn{y_2}, \eqn{\delta_2}.
}
  \item{lin.pred}{
	a list containing three formula objects that correspond to \eqn{h_g}, \eqn{g}=1,2,3.
}
  \item{data}{
	a data.frame in which to interpret the variables named in the formulas in \code{lin.pred}.
}
\item{model}{
	a character value that specifies the type of a model based on the assumption on \eqn{h_3}: "semi-Markov" or "Markov". 
}
  \item{frailty}{
	a logical value to determine whether to include the subject-specific shared frailty term, \eqn{\gamma}, into the model.
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

}
\author{
Sebastien Haneuse and Kyu Ha Lee\cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
%%\note{
%%further notes
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{print.Freq_HReg}}, \code{\link{summary.Freq_HReg}}, \code{\link{plot.Freq_HReg}}, \code{\link{BayesID_HReg}}.
}
\examples{
\dontrun{
# loading a data set
data(scrData)
Y <- scrData[,c(1,2,3,4)]
form1 <- as.formula( ~ x1 + x2 + x3)
form2 <- as.formula( ~ x1 + x2)
form3 <- as.formula( ~ x1 + x2)
lin.pred <- list(form1, form2, form3)

fit_WB	<- FreqID_HReg(Y, lin.pred, data=scrData, model="semi-Markov")

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
\keyword{ semi-competing risks analysis }
\keyword{ hazard regression models }% __ONLY ONE__ keyword per line
