mle.bb <-
function(nA, nTotal, pi=0.5, rho=0.1, min.pi=1e-4, max.pi=0.9999, 
         min.rho=1e-4, max.rho=0.9999, ws=NULL)
{
  
  if(!is.numeric(nA)){
    stop("nA must be a numerical vector\n")
  }
  
  if(!is.numeric(nTotal)){
    stop("nTotal must be a numerical vector\n")
  }
  
  if(length(nA) != length(nTotal)){
    stop("nA and nTotal have different lengths\n")  
  }
  
  if(any(nA > nTotal)){
    stop("nA must be smaller than nTotal\n")  
  }
  
  if(is.null(ws)){ ws = rep(1, length(nA)) }
  
  # -------------------------------------------------------------
  # find MLE
  # -------------------------------------------------------------
  
  logLik = 0.0
  
  W = .C("mle_bb", as.integer(length(nA)), as.double(nA), as.double(nTotal), 
         as.double(ws), pi=as.double(pi), rho=as.double(rho), 
         as.double(min.pi), as.double(max.pi), as.double(min.rho), 
         as.double(max.rho), logLik=as.double(logLik), PACKAGE="basic")
  
  list(pi=W$pi, rho=W$rho, logLik=W$logLik)
}
