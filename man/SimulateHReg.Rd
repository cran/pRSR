\name{SimulateHReg}
\alias{SimulateHReg}
\title{ Simulate Harmonic Regression}
\description{
Simulates a harmonic regression.  Possible error distributions are
normal, t(5), t(5,6), AR1.
}
\usage{
SimulateHReg(n, f, A, B, simPlot = FALSE, Dist = "n", phi = 0.3)
}
\arguments{
  \item{n}{length of series }
  \item{f}{frequency }
  \item{A}{cosine amplitude }
  \item{B}{sine amplitude }
  \item{simPlot}{ plot simulated series }
  \item{Dist}{ one of "n" for normal, "t" for t(5), "s" skewed t(5,6), "a"
  autocorrelated AR1 with parameter phi }
  \item{phi}{ only used if AR1 error distribution is selected }
}
\details{
This will simulate harmonic hegression
}
\value{
vector of length n, simulated harmonic series
}

\author{M.S. Islam}


\seealso{  \code{\link{FitHReg}}}
\examples{
z<-SimulateHReg(10, f=2.5/10, 1, 2)
FitHReg(z)
}
\keyword{ ts }
