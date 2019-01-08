\name{bb.mix2}
\alias{bb.mix2}
\title{
  mixture distribution of two beta-binomial distributions
}
\description{
  Maximum likelihood estimation of model parameters for a mixture of two 
  beta-binomial distributions
}
\usage{
bb.mix2(nA, nTotal, pi1=0.01, rho1=0.01, pi2=0.5, rho2=0.1, maxIt = 100, 
min.pi = 1e-04, min.rho = 1e-04, converged = 1e-05, reEst = 1, traceIt = 0)
}
\arguments{
  \item{nA}{
  a vector of the number of sucesses across n experiments, or the number of reads with alternative alleles across n loci. 
}
  \item{nTotal}{
  a vector of the number of trials across n experiments, or the total number of reads (i.e., read-depth) across n loci. 
}
  \item{pi1}{
  expected proportion of success for the first mixture component.
}
  \item{rho1}{
  over-dispersion parameter for the first mixture component.
}
  \item{pi2}{
  expected proportion of success for the second mixture component.
}
  \item{rho2}{
  over-dispersion parameter for the second mixture component.
}
  \item{maxIt}{
  maximum iteration of EM algorithm.
}
  \item{min.pi}{
  minimum value of pi1 or pi2.
}
  \item{min.rho}{
  minimum value of rho1 or rho2.
}
  \item{converged}{
  threshold for coveregnce checking. 
}
  \item{reEst}{
  if reEst=1, re-estimate the the paramerers of the mixture componenets. 
  If reEst=0, only output he likelihood
}
  \item{traceIt}{
  whether to print out some tracing information. 
}
}
\details{
  Consider a generative model for a beta-binomial distribution. Y given Z follows a binomial distribution with sucess probability rho, and rho follows a beta distribution with parameters a and b. Then in our notation pi = a/(a + b), and rho = 1/(1 + a + b). 
}
\value{
A list of five components
\item{posteriorP }{posterior probablity that one observation belong to one of the two mixture components}
\item{nA}{the input values of nA. 
}
\item{nTotal}{the input values of nTotal. 
}
\item{theta }{a vector of length five, giving the values of the parameter estimates for pi1, rho1, pi2, rho2, and tau, mixture proportion for the 1st component}
\item{logLik }{log likelihood across iterations}

}
\author{ Wei Sun weisun2013@gmail.com }
\keyword{ methods }

