\name{GetFitHReg}
\alias{GetFitHReg}

\title{Compute loglikelihood ratio test statistic }

\description{
The loglikelihood ratio test statistic is computed for testing for periodicity
}

\usage{
GetFitHReg(y, t, nf=150)
}
\arguments{
  \item{y}{ vector containing the series }
  \item{t}{ vector of corresponding time points }
  \item{nf}{ nf, number of frequencies to enumerate }
}
\details{
This function interfaces with C code for fast evaluation.
}
\value{
LR statistic and 
estimated frequency
}

\references{
Islam, M.S. (2008).  
Peridocity, Change Detection and Prediction in Microarrays. 
Ph.D. Thesis, The University of Western Ontario. 
}


\examples{
#Simple Examples
z<-SimulateHReg(10, f=2.5/10, 1, 2)
GetFitHReg(z)
t<-seq(2,20,2)
GetFitHReg(y=z, t=t)
GetFitHReg(z, nf=25)
}

\keyword{ ts }
