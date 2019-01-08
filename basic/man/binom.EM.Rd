\name{binom.EM}
\alias{binom.EM}
\title{
  clustering mutations based on bulk sequencing data and copy number calls
}
\description{
  Using an EM algorithm, clustering mutations based on VAF esitmated from bulk sequencing data, and account for the copy number calls for each mutation. 
}
\usage{
binom.EM(y, z, cn, k, tau0 = NULL, converge = 1e-05, maxIt = 100, tau.max=0.5-1e-6)
}
\arguments{
  \item{y}{
  a vector of the number of alternative reads for each mutation
}
  \item{z}{
  a vector of the total number of reads for each mutation
}
  \item{cn}{
  a vector of copy number calls for each mutation
}
  \item{k}{
  the number of mutations
}
  \item{tau0}{
  initial values of cluster means. Since we expect somatic mutations occur in one of the two alleles, cluster means should be within the range of (0, 0.5). If tau0 is NULL, cluster means will be choosen uniformly from the interval of (0, 0.5). 
}
  \item{converge}{
  the EM algorithm is declared coverged if the change of cluster means and mixure proportions are smaller than coverge. 
}
  \item{maxIt}{
  the maximum number of iterations. 
}
  \item{tau.max}{
  the maximum value of tau. 
}

}
\details{
  Given copy number for the i-th observation, denoted by cn[i], the cluster mean of binomial dsitribution is adjusted by tau*2/cn[i]. 
}
\value{
A list of three entries
  \item{pi1 }{MLE of mixture proportions \eqn{\pi}. }
  \item{tau1 }{MLE of cluster means \eqn{\tau}. }
  \item{postP }{a matrix of posterior probabilities. }
}
\author{ Wei Sun weisun2013@gmail.com }
\keyword{ methods }
