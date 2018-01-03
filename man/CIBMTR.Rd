\name{CIBMTR}
\alias{CIBMTR}
\docType{data}
\title{
Center for International Blood and Bone Marrow Transplant Research (CIBMTR) data
}
\description{
We provide a dataset with five covariates from a study of acute graft-versus-host (GVHD) disease with 9651 patients who underwent first allogeneic hematopoietic cell transplant. We also provide an algorithm to simulate semi-competing risks outcome data. 
}
\usage{data("CIBMTR")}
\format{
  A data frame with 9651 observations on the following 5 variables.
  \describe{
    \item{\code{sexP}}{patient sex: \code{M}-Male, \code{F}-Female}
    \item{\code{ageP}}{patient age: \code{LessThan10}, \code{10to19}, \code{20to29}, \code{30to39}, \code{40to49}, \code{50to59}, \code{60plus}}
    \item{\code{dType}}{disease type: \code{AML}-Acute Myeloid Leukemia, \code{ALL}-Acute Lymphoblastic Leukemia, \code{CML}-Chronic Myeloid Leukemia, \code{MDS}-Myelodysplastic Syndrome}
    \item{\code{dStatus}}{disease stage: \code{Early}-early, \code{Int}-intermediate, \code{Adv}-advanced}
    \item{\code{donorGrp}}{human leukocyte antigen compatibility: \code{HLA_Id_Sib}-identical sibling, \code{8_8}-8/8, \code{7_8}-7/8}
  }
}
\details{
See Examples below for an algorithm to simulate semi-competing risks outcome data. 
}
\source{
Center for International Blood and Bone Marrow Transplant Research
}
\references{
Lee, C., Lee, S.J., Haneuse, S. (2017+). Time-to-event analysis when the event is defined on a finite time interval. \emph{under review}.
}
\seealso{
	\code{\link{CIBMTR_Params}}
}
\examples{
data(CIBMTR_Params)
data(CIBMTR)

## CREATING DUMMY VARIABLES ##

# Sex (M: reference)
CIBMTR$sexP <- as.numeric(CIBMTR$sexP)-1

# Age (LessThan10: reference)
CIBMTR$ageP20to29 <- as.numeric(CIBMTR$ageP=="20to29")
CIBMTR$ageP30to39 <- as.numeric(CIBMTR$ageP=="30to39")
CIBMTR$ageP40to49 <- as.numeric(CIBMTR$ageP=="40to49")
CIBMTR$ageP50to59 <- as.numeric(CIBMTR$ageP=="50to59")
CIBMTR$ageP60plus <- as.numeric(CIBMTR$ageP=="60plus")

# Disease type (AML: reference)
CIBMTR$dTypeALL <- as.numeric(CIBMTR$dType=="ALL")
CIBMTR$dTypeCML <- as.numeric(CIBMTR$dType=="CML")
CIBMTR$dTypeMDS <- as.numeric(CIBMTR$dType=="MDS")

# Disease status (Early: reference)
CIBMTR$dStatusInt <- as.numeric(CIBMTR$dStatus=="Int")
CIBMTR$dStatusAdv <- as.numeric(CIBMTR$dStatus=="Adv")

# HLA compatibility (HLA_Id_Sib: reference)
CIBMTR$donorGrp8_8 <- as.numeric(CIBMTR$donorGrp=="8_8")
CIBMTR$donorGrp7_8 <- as.numeric(CIBMTR$donorGrp=="7_8")

# Covariate matrix
x <- CIBMTR[,c("sexP","ageP20to29","ageP30to39","ageP40to49","ageP50to59","ageP60plus",
"dTypeALL","dTypeCML","dTypeMDS","dStatusInt","dStatusAdv","donorGrp8_8","donorGrp7_8")]

# Set the parameter values
beta1 <- CIBMTR_Params$beta1.true
beta2 <- CIBMTR_Params$beta2.true
beta3 <- CIBMTR_Params$beta3.true
alpha1 <- CIBMTR_Params$alpha1.true
alpha2 <- CIBMTR_Params$alpha2.true
alpha3 <- CIBMTR_Params$alpha3.true
kappa1 <- CIBMTR_Params$kappa1.true
kappa2 <- CIBMTR_Params$kappa2.true
kappa3 <- CIBMTR_Params$kappa3.true
theta <- CIBMTR_Params$theta.true

set.seed(1405)
simCIBMTR <- simID(cluster=NULL, x, x, x, beta1, beta2, beta3, alpha1, alpha2, alpha3,
          kappa1, kappa2, kappa3, theta, SigmaV.true=NULL, cens=c(365,365))
          
names(simCIBMTR) <- c("time1", "event1", "time2", "event2")
CIBMTR <- cbind(simCIBMTR, CIBMTR) 
head(CIBMTR)
}
\keyword{datasets}
