\name{logsumexp}
\alias{logsumexp}
\title{
  calculate log(sum(exp()))
}
\description{
  calculate log(sum(exp())) for values saved in log scale to avoid overflow or underflow in non-log scale. 
}
\usage{
logsumexp(v)
}
\arguments{
  \item{v}{
  a vector of input values. 
}
}
\value{
a scalar of the log(sum(exp(v)))
}
\author{ Wei Sun weisun2013@gmail.com }
\keyword{ utilities }
