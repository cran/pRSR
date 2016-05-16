\name{pRSR-package}
\alias{pRSR-package}
\alias{pRSR}
\docType{package}
\title{
It tests periodicity for any series using response surface regression (RSR). 
}
\description{
It tests periodicity for any series using response surface regression (RSR).
Whole response surface is integrated in the package. Therefore, one can easily get the pvalue
for testing periodicity for any series length. However, the main focus of this algorithm is for
short series. Plot for response surface regression for different quantile values can also be obtained.  
}
\details{                 
\tabular{ll}{
Package: pRSR
Type: Package
Title:  Test of Periodicity using Response Surface Regression
Version: 3.1.1
Date: 2016-05-12
Author: M. S. Islam
Maintainer: M. S. Islam <shahed-sta@sust.edu>
Depends: R (>= 2.15.1)
Description: Tests periodicity in short time series using response surface regression. Moreover,
some useful graphs can be produced.
License: GPL (>= 2)
LazyLoad: yes
}

}
\author{
M. S. Islam
Maintainer: M. S. Islam <shahed-sta@sust.edu>
}
\references{
Islam, M.S. (2008). Peridocity, Change Detection and Prediction in Microarrays. 
Ph.D. Thesis, The University of Western Ontario. 

MacKinnon, J. G. (2002). Computing numerical distribution functions in econometrics. 
In proceedings of High Performance Computing Systems and Applications, 
edited by Pollard, A., Mewhort, D. J. and Weaver, D. F. Springer US. Vol. 451, 455-471.
}
\keyword{ package }

\examples{
# #Testing periodicity
z<-SimulateHReg(20, f=2.5/20, 1, 2)
pvalrsr(z)  # finding p-value using RSR

# For comparing with Fisher's g test
# library(GeneCycle)
# fisher.g.test(z) # Fisher's g test

# Plot for 75%, 90% and 95% quantiles.
plotrsr(n=10:50, q=c(75,90,95))
}
