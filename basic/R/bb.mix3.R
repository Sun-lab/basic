bb.mix3 <- 
function(nA, nTotal, theta = c(0.01, 0.01, 0.5, 0.1,0.99, 0.01, 0.3, 0.4), 
         maxIt=100, range.pi1=c(1e-4, 0.15), range.pi2=c(0.25, 0.75), 
         range.pi3=c(0.85, 0.9999), range.rho1=c(0.01, 0.5), 
         range.rho2=c(0.01, 0.5), range.rho3=c(0.01, 0.5), 
         min.psi=1e-4, converged=1e-5, reEst=2, ws=NULL, traceIt=0)
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
  
  logLik = rep(NA, maxIt)
  theta0 = theta
  
  for(it in 1:maxIt){
    if(traceIt){
      print(sprintf("%d: %s", it, date()))
    }
    
    # ----------------------------------------------------------------
    # E-step: estimate class probability
    # ----------------------------------------------------------------
    logL.Mx = matrix(NA, nrow=length(nTotal), ncol=3)
    
    for(i in 1:length(nTotal)){
      logL.Mx[i,1] = dbetabinom(nA[i], size=nTotal[i], prob=theta[1], 
                                rho=theta[2], log=TRUE)
      
      logL.Mx[i,2] = dbetabinom(nA[i], size=nTotal[i], prob=theta[3], 
                                rho=theta[4], log=TRUE)
      
      logL.Mx[i,3] = dbetabinom(nA[i], size=nTotal[i], prob=theta[5], 
                                rho=theta[6], log=TRUE)
    }
    
    logL.Mx[,1] = logL.Mx[,1] + log(theta[7])
    logL.Mx[,2] = logL.Mx[,2] + log(theta[8])
    logL.Mx[,3] = logL.Mx[,3] + log(1 - theta[7] - theta[8])
    
    logLik[it] = sum(apply(logL.Mx, 1, logsumexp))
    
    # calculate posterior prob while keeping numbers in log scale
    postPs1  = 1/(1+exp(logL.Mx[,2]-logL.Mx[,1])+exp(logL.Mx[,3]-logL.Mx[,1]))
    postPs2  = 1/(1+exp(logL.Mx[,1]-logL.Mx[,2])+exp(logL.Mx[,3]-logL.Mx[,2]))
    pp = cbind(postPs1, postPs2, 1 - postPs1 - postPs2)
    
    if(reEst == 0){ 
      if(traceIt){
        print("reEst=0")
      }
      
      break 
    }
    
    # re-estimate psi1
    theta[7]  = sum(postPs1)/length(postPs1)
    theta[8]  = sum(postPs2)/length(postPs2)
    
    if(reEst == 2){ 
      # re-estimate other parameters
      fit1  = mle.bb(nA, nTotal, pi=theta[1], rho=theta[2], 
                     min.pi=range.pi1[1], max.pi=range.pi1[2], 
                     min.rho=range.rho1[1], max.rho=range.rho1[2], 
                     ws=pp[,1])
      
      fit2  = mle.bb(nA, nTotal, pi=theta[3], rho=theta[4], 
                     min.pi=range.pi2[1], max.pi=range.pi2[2], 
                     min.rho=range.rho2[1], max.rho=range.rho2[2], 
                     ws=pp[,2])
      
      fit3  = mle.bb(nA, nTotal, pi=theta[5], rho=theta[6], 
                     min.pi=range.pi3[1], max.pi=range.pi3[2], 
                     min.rho=range.rho3[1], max.rho=range.rho3[2], 
                     ws=pp[,3])
      
      theta[1] = fit1$pi
      theta[2] = fit1$rho
      
      theta[3] = fit2$pi
      theta[4] = fit2$rho
      
      theta[5] = fit3$pi
      theta[6] = fit3$rho
    }
    
    if(max(abs(theta - theta0)) < converged){
      break
    }else{
      theta0 = theta
    }
    
    if(traceIt){
      print(theta)
    }
  }
  
  colnames(pp) = c("c1", "c2", "c3")
  names(theta) = c("pi1", "rho1", "pi2", "rho2", "pi3", "rho3", "psi1", "psi2")
  
  mix1 = list(posteriorP=pp, nA=nA, nTotal=nTotal, theta=theta, logLik=logLik)
  class(mix1) = "bb.mix"
  mix1
}
