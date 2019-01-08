bb.mix3.plot <- function(fit0){
  
  nA = fit0$nA
  nTotal = fit0$nTotal
  
  y1 = y2 = y3 = rep(NA, length(nTotal))
  
  for(i in 1:length(y1)){
    y1[i] = rbetabinom(1, nTotal[i], prob=fit0$theta["pi1"], 
                       rho=fit0$theta["rho1"])
    
    y2[i] = rbetabinom(1, nTotal[i], prob=fit0$theta["pi2"], 
                       rho=fit0$theta["rho2"])
    
    y3[i] = rbetabinom(1, nTotal[i], prob=fit0$theta["pi3"], 
                       rho=fit0$theta["rho3"])
  }
  
  w2kp = which(nTotal > 0)
  
  d1 = density(y1[w2kp]/nTotal[w2kp], bw=0.05, from=-0.1, to=1.1)
  d2 = density(y2[w2kp]/nTotal[w2kp], bw=0.05, from=-0.1, to=1.1)
  d3 = density(y3[w2kp]/nTotal[w2kp], bw=0.05, from=-0.1, to=1.1)
  
  if(!all(d1$x == d2$x)){
    stop("mismatch of x values\n.")
  }
  
  psi1 = fit0$theta["psi1"]
  psi2 = fit0$theta["psi2"]
  psi3 = 1 - psi1 - psi2
  
  plot(density(nA[w2kp]/nTotal[w2kp], bw=0.05), main="", type="l", lwd=2, 
       col="grey", xlab="MAF of loci w/ missing call")
  lines(d1$x, d1$y*psi1, col="red", lwd=2)
  lines(d2$x, d2$y*psi2, col="darkblue", lwd=2)
  lines(d3$x, d3$y*psi3, col="green", lwd=2)
  lines(d1$x, d1$y*psi1 + d2$y*psi2 + d3$y*psi3, col="orange", lwd=1)
  legend("topright", legend=c("observed MAF", "fitted MAF"), lty=c(1,1), 
         lwd=c(2,1), col=c("grey", "orange"), bty="n")
  
}

