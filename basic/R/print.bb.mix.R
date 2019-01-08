print.bb.mix <-
function(x, ...)
{

  nA     = x$nA
  nTotal = x$nTotal
  
  cat(sprintf("nA[1:%d]:\n", length(nA)))
  print(nA[1:min(5, length(nA))])
  if(length(nA) > 5){ cat(" ...") }
  cat("\n")
  
  cat(sprintf("nTotal[1:%d]:\n", length(nTotal)))
  print(nTotal[1:min(5, length(nTotal))])
  if(length(nTotal) > 5){ cat(" ...") }
  cat("\n")
  
  cat("theta:\n", sep="")
  print(x$theta)
  cat("\n")
  
}

