\name{FitHReg}
\alias{FitHReg}
\title{ Fits Three Parameter Harmonic Regression }
\description{
Estimates A, B and f in the harmonic regression,
y(t)=mu+A*cos(2*pi*f*t)+B*sin(2*pi*f*t)+e(t)
using LS.
}

\usage{FitHReg(y, t = 1:length(y), nf=150)}

\arguments{
  \item{y}{ series }
  \item{t}{ time points }
  \item{nf}{ nf, number of frequencies to enumerate }
}

\details{
Program is interfaced to C for efficient computation.
}

\value{
Object of class "HReg" produced.
This is a list with components:
'coefficients', 'residuals', 'Rsq', 'fstatistic', 'sigma', 'freq', 'LRStat'
corresponding to the 3 regression coefficients, residuals, R-squared,
F-statistic, residual sd, optimal frequency and LR-test statistic for
null hypothesis white noise.
}

\author{M.S. Islam}

\references{
Islam, M.S. (2008).  
Peridocity, Change Detection and Prediction in Microarrays. 
Ph.D. Thesis, The University of Western Ontario. 
}


\seealso{ \code{\link{GetFitHReg}} }

\examples{
z<-SimulateHReg(10, f=2.5/10, 1, 2)
FitHReg(z)
}

\keyword{ ts }
