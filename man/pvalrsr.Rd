\name{pvalrsr}
\alias{pvalrsr}
\title{
Finds pvalue for testing periodicity using RSR.
}
\description{
It finds pvalue for tesing periodicity for any series using RSR. 
}
\usage{
pvalrsr(x,  t=1:length(x), nf=150, Numpq = 11)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{series to be tested for periodicity}
  \item{t}{vector of corresponding time points}
  \item{nf}{number of frequencies to enumerate}
  \item{Numpq}{Numebr of indices for interpolation in RSR}
}
\details{
A full RSR is integral part of the package. This was done using likelihood ratio statistic
for simulated series from white noise process. For more information about the procedure, please see
the first reference. 
}
\value{
pvalue
}
\references{
Islam, M.S. (2008). Peridocity, Change Detection and Prediction in Microarrays. 
Ph.D. Thesis, The University of Western Ontario. 

MacKinnon, J. G. (2001). Computing numerical distribution functions in econometrics. 
In proceedings of High Performance Computing Systems and Applications, 
edited by Pollard, A., Mewhort, D. J. and Weaver, D. F. Springer US. Vol. 451, 455-471.
}
\author{
M. S. Islam
}
\examples{
# Non-Fourier frequency
z<-SimulateHReg(20, f=2.5/20, 1, 2)
pvalrsr(z)  # finding p-value using RSR

# For comparing with Fisher's g test
# library(GeneCycle)
# fisher.g.test(z) # Fisher's g test


# Fourier frequency
y<-SimulateHReg(20, f=2/20, 1, 2)
pvalrsr(y)  # finding p-value using RSR
# For comparing with Fisher's g test
# library(GeneCycle)
# fisher.g.test(z) # Fisher's g test
}
\keyword{ts}