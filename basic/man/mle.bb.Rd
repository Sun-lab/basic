\name{mle.bb}
\alias{mle.bb}
\title{
weighted MLE of beta-binomial distribution
}
\description{
weighted Maximum Likelihood Estimates of the mean and over-dispersion parameter for a beta-binomial distribution. 
}
\usage{
mle.bb(nA, nTotal, pi = 0.5, rho = 0.1, min.pi=1e-4, max.pi=0.9999, 
         min.rho=1e-4, max.rho=0.9999, ws = NULL)
}
\arguments{
  \item{nA}{
   number of success, or the number of alternative reads in sequencing data 
   applications. Given nTotal, nA is the random variable that follows a 
   beta-binomial distribution. nA should be a vector of non-negative inteters. 
}
  \item{nTotal}{
  number of trials, or the total number of reads (i.e., the number of 
  reference or alternative reads) in  data applications. nTotal should be 
  a vector of non-negative inteters of the same length as nA. In addition, 
  nTotal[i] >= nA[i]. 
}
  \item{pi}{
  initial value of \code{pi}, see details. 
}
  \item{rho}{
  initial value of \code{rho}, see details. 
}
  \item{min.pi, max.pi}{
  minimum and maximum value of \code{pi}
}
  \item{min.rho, max.rho}{
  minimum and maximum value of \code{rho}
}
  \item{ws}{
  the weights for each observatoin. A non-negative vector of the same length 
  as nA or nTotal.
}
}
\details{
The beta-binomial distribution extends the binomial distribution by allowing 
the probability of success at each trial following a beta distribution, with 
parameters \eqn{\alpha} and \eqn{\beta}. There are a few different ways to 
parameterize a beta-binomial distribution. One is to use \eqn{\alpha} and \eqn{\beta} as parameters. We use an alternative parameterization 
\eqn{\pi = \alpha/(\alpha + \beta) }, which is the mean parameter, and 
\eqn{\rho = 1/(1 + \alpha + \beta) }, which is known as the ``intra class'' 
or ``intra cluster'' correlation.  
}
\value{
A list of three entries
  \item{pi }{MLE of \eqn{\pi}. }
  \item{theta }{MLE of \eqn{\theta}. }
  \item{logLik }{log likelihood. }
}
\examples{

set.seed(1999)
nTotal = rnbinom(n=8000, size=1.1, mu=40)
nA     = rep(0, length(nTotal))
prob1  = 0.001; rho1 = 0.005

for(i in 1:length(nTotal)){
    nA[i] = rbetabinom(1, size=nTotal[i], prob=prob1, rho=rho1)
}

date()
fit0 = mle.bb(nA, nTotal)
date()
fit1 = vglm(cbind(nA, nTotal-nA) ~ 1, betabinomial)
date()

fit0
fit1

}
\author{ Wei Sun weisun2013@gmail.com }
\keyword{ methods }
