\name{dist.binom}
\alias{dist.binom}
\title{
calculate distance for two binomial variables
}
\description{
calculate distance for two binomial variables, using truncated likelihood ratio test statistic. Here the two variables correspond to two loci
}
\usage{
dist.binom(y1, z1, y2, z2, cn1=NULL, cn2=NULL)
}
\arguments{
  \item{y1}{
  number of success for the 1st variable/locus
}
  \item{z1}{
  number of trials for the 1st variable/locus
}
  \item{y2}{
  number of success for the 2nd variable/locus
}
  \item{z2}{
  number of trials for the 2nd variable/locus
}
  \item{cn1}{
  copy number for the 1st locus
}
  \item{cn2}{
  copy number for the 2nd locus
}

}
\value{
likelihood ratio statistic.
}

\author{ Wei Sun weisun2013@gmail.com }
\keyword{ methods }

