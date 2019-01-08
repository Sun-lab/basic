dist.binom <-
function(y1, z1, y2, z2, cn1=NULL, cn2=NULL){
  
  if(is.null(cn1)){ cn1 = rep(2, length(y1)) }
  
  if(is.null(cn2)){ cn2 = rep(2, length(y2)) }
  
  r1 = (y1/z1)*cn1/2
  r2 = (y2/z2)*cn2/2
  
  b   = (z1*cn2 + y1*cn1 + z2*cn1 + y2*cn2)/2
  D   = b*b - cn1*cn2*(y1 + y2)*(z1 + z2)
  r12 = (b - sqrt(D))/(2*(z1 + z2))
  
  if(any(r1 < 0 | r1 > 1)){ stop("invalid values for y1 or z1") }
  if(any(r2 < 0 | r2 > 1)){ stop("invalid values for y2 or z2") }

  r1[which(r1 <= 1e-5)] = 1e-5
  r1[which(r1 >= 1-1e-5)] = 1-1e-5
  
  r2[which(r2 <= 1e-5)] = 1e-5
  r2[which(r2 >= 1-1e-5)] = 1-1e-5
  
  r12[which(r12 <= 1e-5)] = 1e-5
  r12[which(r12 >= 1-1e-5)] = 1-1e-5
  
  logL = y1*log(r1/r12) + (z1-y1)*log((1-2*r1/cn1)/(1-2*r12/cn1))
  logL = logL + y2*log(r2/r12) + (z2-y2)*log((1-2*r2/cn2)/(1-2*r12/cn2))
  
  if(any(logL < -1e-12)){ stop("we expect logL >= 0") }
  logL[which(logL < 1e-12)] = 1e-12
  
  logL
}
