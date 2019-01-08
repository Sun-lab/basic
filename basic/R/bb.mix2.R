bb.mix2 <-
function(nA, nTotal, pi1=0.01, rho1=0.01, pi2=0.5, rho2=0.1, maxIt=100, 
         min.pi=1e-4, min.rho=1e-4, converged=1e-5, reEst=1, traceIt=0)
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
  
  tau1   = 0.5
  theta0 = c(pi1, rho1, pi2, rho2, tau1)
  logLik = rep(NA, maxIt)
  
  for(it in 1:maxIt){
    if(traceIt){
      print(sprintf("%d: %s", it, date()))
    }
    
    # estimate class probability
    logL.Mx = matrix(NA, nrow=length(nTotal), ncol=2)
    
    for(i in 1:length(nTotal)){
      logL.Mx[i,1] = dbetabinom(nA[i], size=nTotal[i], prob=pi1, 
                                rho=rho1, log=TRUE)
      logL.Mx[i,2] = dbetabinom(nA[i], size=nTotal[i], prob=pi2, 
                                rho=rho2, log=TRUE)
    }
    
    logL.Mx[,1] = logL.Mx[,1] + log(tau1)
    logL.Mx[,2] = logL.Mx[,2] + log(1- tau1)
    
    logLik[it] = sum(apply(logL.Mx, 1, logsumexp))
    
    # calculate posterior prob while keeping numbers in log scale
    postPs  = 1 / (1 + exp(logL.Mx[,2] - logL.Mx[,1]))
    pp = cbind(postPs, 1 - postPs)
    
    if(reEst == 0){ 
      if(traceIt){
        print("reEst=0")
      }
      
      break 
    }
    
    # re-estimate tau1
    tau1  = sum(postPs)/length(postPs)
    
    # re-estimate other parameters
    fit1  = mle.bb(nA, nTotal, pi1, rho1, min.pi, min.rho, ws=pp[,1])
    fit2  = mle.bb(nA, nTotal, pi2, rho2, min.pi, min.rho, ws=pp[,2])
    
    pi1  = fit1$pi
    rho1 = fit1$rho
    
    pi2  = fit2$pi
    rho2 = fit2$rho
    
    theta = c(pi1, rho1, pi2, rho2, tau1)
    
    if(max(abs(theta - theta0)) < converged){
      break
    }else{
      theta0 = theta
    }
    
    if(traceIt){
      print(theta)
    }
  }
  
  colnames(pp) = c("c1", "c2")
  
  names(theta) = c("pi1", "rho1", "pi2", "rho2", "tau1")
  mix1 = list(posteriorP=pp, nA=nA, nTotal=nTotal, theta=theta, logLik=logLik)
  class(mix1) = "bb.mix"
  
  mix1
}
