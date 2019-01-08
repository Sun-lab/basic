binom.EM <-
function(y, z, cn, k, tau0=NULL, converge=1e-5, maxIt=100, tau.max=0.5-1e-6){
  
  if(length(y) != length(z)){
    stop("length of y and z should be the same\n")
  }
  
  if(any(y != round(y)) || any(y < 0)){
    stop("y must be non-negative integers\n")
  }
  
  if(any(z != round(z)) || any(z <= 0)){
    stop("z must be positive integers\n")
  }
  
  if(any(y > z)){
    stop("y must be equal to or smaller than z\n")
  }
  
  if(is.null(tau0)){
    tau0 = seq(0, 0.5, length.out=k+2)
    tau0 = tau0[2:(k+1)]
  }
  
  if(!is.null(tau0) && length(tau0) != k){
    stop("length of tau0 should be the same as k\n")
  }
  
  nn   = length(y)
  pi0  = rep(1/k, k) ## cluster priors
  
  for(it in 1:maxIt){
    # cat(it, date(), "\n")
    # cat(tau1, "\n")
    logP = matrix(NA, nrow=nn, ncol=k)
    
    for(j in 1:k){
      
      logP[,j]  = dbinom(y, z, tau0[j]*2/cn, log=TRUE)
    }
    
    log.postP.num = t(t(logP) + log(pi0))
    dim(log.postP.num)
    log.postP.num[1:2,]
    
    log.postP = log.postP.num - apply(log.postP.num, 1, logsumexp)
    postP     = exp(log.postP)
    
    tau1 = colSums(y*postP)/colSums(2*z/cn*postP)
    tau1[which(tau1 > tau.max)] = tau.max
    pi1  = colSums(postP)
    pi1  = pi1/sum(pi1)
    
    if(all(abs(c(pi1, tau1) - c(pi0, tau0)) < converge)){
      break
    }
    
    pi0 = pi1
    tau0 = tau1
  }
  
  list(pi=pi1, tau=tau1, postP=postP)
}
