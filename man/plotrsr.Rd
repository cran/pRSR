\name{plotrsr}
\alias{plotrsr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of response surface regression for different quantile
}
\description{
It produces a single plot of response surface regression for different quantile values. 
}
\usage{
plotrsr(n = 20:50, q = c(90, 95),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sequence of series size for the plot}
  \item{q}{Vector or single quantile value  }
  \item{\dots}{Further arguments for function lines() }
}

\value{
Plot of RSR
}
\references{
MacKinnon, J. G. (2002). Computing numerical distribution functions in econometrics. 
In proceedings of High Performance Computing Systems and Applications, 
edited by Pollard, A., Mewhort, D. J. and Weaver, D. F. Springer US. Vol. 451, 455-471.
}
\author{
M. S. Islam
}

\seealso{
 \code{\link{pvalrsr}}
}
\examples{
# Plot for 75%, 90% and 95% quantiles.
plotrsr(n=10:50, q=c(75,90,95))
# Plot for 80%, 90%, 95% and 99% quantiles.
# We use color red and dashed line
plotrsr(n=10:50, q=c(80,90,95, 99), col=2, lty=3)
}
\keyword{ts}