\name{initiate.startValues_AFT}
\alias{initiate.startValues_AFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function that initiates starting values for a single chain.
}
\description{
The function initiates starting values for a single chain for accelrated failture time (AFT) models. Users are allowed to set some non-null values to starting values for a set of parameters. The function will automatically generate starting values for any parameters whose values are not specified.
}
\usage{
    initiate.startValues_AFT(Formula, data, model, nChain=1,
                            beta1=NULL, beta2=NULL, beta3=NULL, beta=NULL,
                            gamma=NULL, theta=NULL,
                            y1=NULL, y2=NULL, y=NULL,
                            LN.mu=NULL, LN.sigSq=NULL,
                            DPM.class1=NULL, DPM.class2=NULL, DPM.class3=NULL,
                            DPM.class=NULL, DPM.mu1=NULL, DPM.mu2=NULL,
                            DPM.mu3=NULL, DPM.mu=NULL, DPM.zeta1=NULL,
                            DPM.zeta2=NULL, DPM.zeta3=NULL, DPM.zeta=NULL,
                            DPM.tau=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Formula}{
	For \code{BayesID_AFT}, it is a data.frame containing semi-competing risks outcomes from \code{n} subjects. See \code{BayesID_AFT}.
	For \code{BayesSurv_AFT}, it is a data.frame containing univariate time-to-event outcomes from \code{n} subjects. See \code{BayesSurv_AFT}.
	For \code{BayesID_AFT}, it is a list containing three formula objects that correspond to the transition \eqn{g}=1,2,3.
	For \code{BayesSurv_AFT}, it is a formula object that corresponds to \eqn{log(t)}.
}
  \item{data}{
	a data.frame in which to interpret the variables named in the formula(s) in \code{lin.pred}.
}
  \item{model}{
	a character vector that specifies the type of components in a model. Check \code{\link{BayesID_AFT}} and \code{\link{BayesSurv_AFT}}.
}
\item{nChain}{
The number of chains.
}
  \item{beta1}{
	starting values of \eqn{\beta_1} for \code{BayesID_AFT}.
  }
  \item{beta2}{
	starting values of \eqn{\beta_2} for \code{BayesID_AFT}.
  }
  \item{beta3}{
	starting values of \eqn{\beta_3} for \code{BayesID_AFT}.
  }
  \item{beta}{
	starting values of \eqn{\beta} for \code{BayesSurv_AFT}.
  }
  \item{gamma}{
	starting values of \eqn{\gamma} for \code{BayesID_AFT}.
  }
  \item{theta}{
	starting values of \eqn{\theta} for \code{BayesID_AFT}.
  }
  \item{y1}{
	starting values of \eqn{log(t_1)} for \code{BayesID_AFT}.
  }
  \item{y2}{
	starting values of \eqn{log(t_2)} for \code{BayesID_AFT}.
  }
  \item{y}{
      starting values of \eqn{log(t)} for \code{BayesSurv_AFT}.
  }
  \item{LN.mu}{
	starting values of \eqn{\beta_0} in logNormal models for \code{BayesID_AFT} and \code{BayesSurv_AFT}.
  }
  \item{LN.sigSq}{
	starting values of \eqn{\sigma^2} in logNormal models for \code{BayesID_AFT} and \code{BayesSurv_AFT}.
  }
  \item{DPM.class1}{
      starting values of the class membership for transition 1 in DPM models for \code{BayesID_AFT}.
  }
  \item{DPM.class2}{
      starting values of the class membership for transition 2 in DPM models for \code{BayesID_AFT}.
  }
  \item{DPM.class3}{
      starting values of the class membership for transition 3 in DPM models for \code{BayesID_AFT}.
  }
  \item{DPM.class}{
      starting values of the class membership in DPM models for \code{BayesSurv_AFT}.
  }
  \item{DPM.mu1}{
      starting values of \eqn{\mu_1} in DPM models for \code{BayesID_AFT}.
  }
  \item{DPM.mu2}{
      starting values of \eqn{\mu_2} in DPM models for \code{BayesID_AFT}.
  }
  \item{DPM.mu3}{
      starting values of \eqn{\mu_3} in DPM models for \code{BayesID_AFT}.
  }
  \item{DPM.mu}{
      starting values of \eqn{\mu} in DPM models for \code{BayesSurv_AFT}.
  }
  \item{DPM.zeta1}{
      starting values of \eqn{\zeta_{1}} in DPM models for \code{BayesID_AFT}.
  }
  \item{DPM.zeta2}{
      starting values of \eqn{\zeta_{2}} in DPM models for \code{BayesID_AFT}.
  }
  \item{DPM.zeta3}{
      starting values of \eqn{\zeta_{3}} in DPM models for \code{BayesID_AFT}.
  }
  \item{DPM.zeta}{
      starting values of \eqn{\zeta} in DPM models for \code{BayesSurv_AFT}.
  }
  \item{DPM.tau}{
	starting values of \eqn{\tau} in DPM models for \code{BayesID_AFT} and \code{BayesSurv_AFT}.
  }
}
%\details{
%	
%}
\value{
	\code{initiate.startValues_AFT} returns a list containing starting values for a sigle chain that can be used for \code{BayesID_AFT} and \code{BayesSurv_AFT}.
}
\references{
Lee, K. H., Rondeau, V., and Haneuse, S. (2017),
Accelerated failure time models for semicompeting risks data in the presence of complex censoring, \emph{Biometrics}, 73, 4, 1401-1412. \cr \cr
Alvares, D., Haneuse, S., Lee, C., Lee, K. H. (2018+),
SemiCompRisks: an R package for independent and cluster-correlated analyses of semi-competing risks data, \emph{submitted}, arXiv:1801.03567. \cr
}
\author{
Sebastien Haneuse and Kyu Ha Lee\cr
Maintainer: Kyu Ha Lee <klee15239@gmail.com>
}
%\note{
%Function that initiates starting values for a single chain
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{BayesID_AFT}}, \code{\link{BayesSurv_AFT}}
}
\examples{
## See Examples in \code{\link{BayesID_AFT}} and \code{\link{BayesSurv_AFT}}.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian framework }
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line
\keyword{ univariate analysis }
\keyword{ accelerated failure time models }% __ONLY ONE__ keyword per line
