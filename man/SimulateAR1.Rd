\name{SimulateAR1}
\alias{SimulateAR1}
\title{ Simulate AR(1) series }
\description{
An AR(1) series with mean zero and variance 1 and
with autocorrelation paramater phi is simulated. 
}
\usage{
SimulateAR1(n, phi)
}

\arguments{
  \item{n}{ length of series }
  \item{phi}{autocorrelation parameter }
}
\details{
The model equation is:
z[t] = phi*z[t-1]+a[t],
where z[1] is N(0,1) and a[t] are NID(0, siga),
\eqn{siga=\sqrt(1/(1-phi^2))}.
}
\value{
autocorrelated time series of length n
}
\author{M.S. Islam}

\seealso{  \code{\link{FitHReg}}, \code{\link{SimulateHReg}}}

\examples{
e<-SimulateAR1(10^4, phi=0.8)
mean(e)
sd(e)
acf(e, lag.max=5, plot=FALSE)
}

\keyword{ ts }
