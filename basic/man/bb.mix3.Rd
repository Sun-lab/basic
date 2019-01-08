\name{bb.mix3}
\alias{bb.mix3}
\title{
  mixture distribution of three beta-binomial distributions
}
\description{
  Maximum likelihood estimation of model parameters for a mixture of three 
  beta-binomial distributions
}
\usage{
bb.mix3(nA, nTotal, theta = c(0.01, 0.01, 0.5, 0.1,0.99, 0.01, 0.3, 0.4), 
maxIt = 100, range.pi1 = c(1e-04, 0.15), range.pi2 = c(0.25, 0.75), 
range.pi3 = c(0.85, 0.9999), range.rho1 = c(0.01, 0.5), 
range.rho2 = c(0.01, 0.5), range.rho3 = c(0.01, 0.5), 
min.psi = 1e-04, converged = 1e-05, reEst = 2, ws = NULL, 
traceIt = 0)
}
\arguments{
  \item{nA}{
  a vector of the number of sucesses across n experiments, or the number of reads with alternative alleles across n loci. 
}
  \item{nTotal}{
  a vector of the number of trials across n experiments, or the total number of reads (i.e., read-depth) across n loci. 
}
  \item{theta}{
  initial values for 8 parameters in the order of pi1, rho1, pi2, rho2, pi3, rho3, psi1, and psi2. pi's and rho's are the expected proportion of success and over-dispersion for the three mixture components. psi1 and psi2 are the mixture proportions for the 1st and 2nd components. 
}
  \item{maxIt}{
  maximum iteration of EM algorithm.
}
  \item{range.pi1, range.pi2, range.pi3, range.rho1, range.rho2, range.rho3}{
  allowed ranges for the parameters pi's and rho's.
}
  \item{min.psi}{
  minimum value of any mixture proportion.
}
  \item{converged}{
  threshold for coveregnce checking. 
}
  \item{reEst}{
  if reEst=0, only update the posterior probability and output the likelihood. If reEst=1, re-estimate mixture proportion. If reEst=2, re-estimate all the parameters. 
}
  \item{ws}{
  weights of each the input data, a vector of the same lenght as nA. 
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
\item{posteriorP }{posterior probablity that one observation belong to one of the three mixture components}
\item{nA}{the input values of nA. }
\item{nTotal}{the input values of nTotal. }
\item{theta }{a vector of length eight, giving the values of the parameter estimates for pi1, rho1, pi2, rho2, pi3, rho3, psi1, and pi2. }
\item{logLik }{log likelihood across iterations}

}
\author{ Wei Sun weisun2013@gmail.com }
\keyword{ methods }

