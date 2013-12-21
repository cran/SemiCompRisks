\name{BMT}
\alias{BMT}
\docType{data}
\title{
Data on 137 Bone Marrow Transplant Patients
}
\description{
Data on 137 Bone Marrow Transplant Patients
}
\usage{data(BMT)}
\format{
  A data frame with 137 observations on the following 22 variables.
  \describe{
    \item{\code{g}}{Disease group; 1-ALL, 2-AML low-risk, 3-high-risk}
    \item{\code{T1}}{Time (in days) to death or on study time}
    \item{\code{T2}}{Disease-Free survival time (time to relapse, death or end of study)}
    \item{\code{delta1}}{Death indicator; 1-Dead, 0-Alive}
    \item{\code{delta2}}{Relapse indicator; 1-Relapsed, 0-Disease-Free}
    \item{\code{delta3}}{Disease-Free survival indicator; 1-Dead or relapsed, 0-Alive disease-free}
    \item{\code{TA}}{Time (in days) to acute graft-versus-host disease}
    \item{\code{deltaA}}{Acute graft-versus-host disease indicator; 1-Developed acute graft-versus-host disease, 0-Never developed acute graft-versus-host disease}
    \item{\code{TC}}{Time (in days) to chronic graft-versus-host disease}
    \item{\code{deltaC}}{Chronic graft-versus-host disease indicator; 1-Developed chronic graft-versus-host disease, 0-Never developed chronic graft-versus-host disease}
    \item{\code{TP}}{Time (in days) to return of platelets to normal levels}
    \item{\code{deltaP}}{Platelet recovery indicator; 1-Platelets returned to normal levels, 0-Platelets never returned to normal levels}
    \item{\code{Z1}}{Patient age in years}
    \item{\code{Z2}}{Donor age in years}
    \item{\code{Z3}}{Patient sex; 1-Male, 2-Female}
    \item{\code{Z4}}{Donor sex; 1-Male, 2-Female}
    \item{\code{Z5}}{Patient CMV status; 1-CMV positive, 0-CMV negative}
    \item{\code{Z6}}{Donor CMV status; 1-CMV positive, 0-CMV negative}
    \item{\code{Z7}}{Waiting time to transplant in days}
    \item{\code{Z8}}{FAB; 1-FAB Grade 4 or 5 and AML, 0-Otherwise}
    \item{\code{Z9}}{Hospital; 1-The Ohio State University, 2-Alfred, 3-St. Vincent, 4-Hahnemann}
    \item{\code{Z10}}{MTX used as a graft-versus-host-prophylactic; 1-Yes, 0-No}
  }
}
\source{
1. R package "\code{KMsurv}". \cr
2. Klein, J. P. and Moeschberger M. L. (2005). 
Survival Analysis: Techniques for Censored and Truncated Data.
}
\references{
Klein, J. P. and Moeschberger M. L. (2005). 
Survival Analysis: Techniques for Censored and Truncated Data.
}
\examples{
data(BMT)
}
\keyword{datasets}
