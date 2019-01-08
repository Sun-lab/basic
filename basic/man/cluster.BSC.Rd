\name{cluster.BSC}
\alias{cluster.BSC}
\title{
  clustering somatic mutations using both bulk tumor and single cell sequencing data
}
\description{
  clustering somatic mutations using mutation calls, read-depth, and mutant read count from both bulk tumor and multiple single tumor cells. 
}
\usage{
cluster.BSC(b.cn, b.asrec, b.trec, mutations, asrec, trec, nclusters,  
         straify.by.maf=TRUE, stratify.n.intervals=20, stratify.window=2,  
         epsilon=0.001, alpha1=0.02, beta1=0.8, alpha2=0.05, beta2=0.6, 
         theta = c(0.01, 0.01, 0.5, 0.1,0.99, 0.01, 0.3, 0.4), 
         maxIt=1000, converged=1e-5)
}
\arguments{
  \item{b.cn}{
 a vector of length n, giving the copy number for the n loci of interest. 
}
  \item{b.asrec}{
 a vector of allele-specific read count from bulk tumor sample 
}
  \item{b.trec}{
 a vector of total read count from bulk tumor sample 
}
  \item{mutations}{
 an n x m matrix of mutation calls for n loci in m single cells. Values 0 or 1 indicate the absence/presence of somatic mutation and missing value indate fail to pass the filteres for mutation call (e.g., low read-depth).  
}
  \item{asrec}{
 an n x m matrix of the number of reads with mutant allele, for n loci in m  cells. The rows and columns of \code{asrec} match the rows and columns of \code{mutations}. 
}
  \item{trec}{
an n x m matrix of read-depth for n loci in m  cells. The rows and columns of \code{trec} match the rows and columns of \code{mutations}. 
}
  \item{nclusters}{
the number of clusers. It could be a vector, e.g., 2:5, and then the results of 2 to 5 clusers will be reported. 
}
  \item{straify.by.maf}{
 we fit a mixture of beta-binomial distribution to decribe read counts from single cells. If straify.by.maf=TRUE, we fit stratif somatic mutations based MAF from bulk tumor and fit such mixture distributions for each strata separately.  
}
  \item{stratify.n.intervals, stratify.window}{
 If straify.by.maf=TRUE, we cover somatic mutations (ordered by their MAF) by sliding windows of size stratify.window*1/stratify.n.intervals, and slide the window by 1/stratify.n.intervals each time. Because these statum have overlaps, one mutations may be used more than once to fit the mixture distribution. After fitting the model aross all stratum, a somatic mutation will use the mixture distirbution fitted by the closest strata. The distance between a somatic mutation and a strata is defined by the distnace between the sotmatic mtuations's MAF and the center of the strata. 
}
  \item{epsilon}{
we assume the underlying true mutation status of each locus in a cell follows a Bernoulli distribution. The success probablity of the Bernoulli distriubtion is bounded by [epsilon, 1 - epsilon]. 
}
  \item{alpha1, alpha2}{
initial values of false positive rates of mutation calls for the loci that pass the mutation calling filters (\code{alpha1}) or not (\code{alpha2}). 
}
  \item{beta1, beta2}{
initial values of 1 - false negative rate of mutation calls for the loci that pass the mutation calling filters (\code{beta1}) or not (\code{beta2}). 
}
  \item{theta}{
a vector of the initial values for pi1, rho1, pi2, rho2, pi3, rho3, psi1, and psi2. We model the number of mutant reads for those loci that do not pass the mutation calling filer (e.g., with low read-depth) by a mixture of three beta-binomial distributions with paramters \code{pik} and \code{rhok} for the k-th mixture component, and mixture proportion psi1, psi2, and 1-psi1-psi2.
}
  \item{maxIt}{
  the maximum number of iterations.
}
  \item{converged}{
  threshold for coveregnce checking. The algorithm is converged if the difference of any parameter between the current and previous iteration is smaller than \code{converged}.
}
}
\details{
  The mutations matrix have all three values 0, 1, and missing. When we fit the mixtrue distribution of three beta-binomials. We first use all the data, and then re-estimate mixture proportions using the data with mising mutation calls. 
}
\value{
A list of the results for each number of clusters, as well as two additional elements: nclust.by.BIC and summary.cluster.
\item{nclust.by.BIC }{number of clusters selected by BIC. For example, "n2" means 2 clusters. }
\item{summary.cluster }{a data frame that summarize the BIC, number of iterations, and gap: the maximum difference of parameter values between the last two iterations. If gap < convergence, the algorithm is converged. Otherwise, the update was terminated because it reaches the maximum number of iterations. }

For each number of cluster, the results is saved as a list with the folllowing componenets. 
\item{Ehats }{ estimate of E. }
\item{theta }{ estimate of theta. }
\item{gammas }{ a maxtir of nclust rows and m columns, with the (i,j)-th elment giving the success probalbity that the underlying true mutation occur for the i-th cluster in the j-th cell. }
\item{taus }{ estimate of taus: mean values for the binomial distributions that model the MAF of bulk tumor sample. }
\item{nus }{ estimate of nus: cluster mixture propotions. }
\item{postPs }{ a maxtir of dimension n x nclust, the posterior probablity each mutation belong to one the clusters. }
\item{logLik }{ log likelihood. }
\item{BIC }{ BIC. }
\item{nIt }{ the number of iterations. }
\item{gap }{ the maximum difference of parameter values between the last two iterations.}
}
\author{ Wei Sun weisun2013@gmail.com }
\keyword{ methods }
