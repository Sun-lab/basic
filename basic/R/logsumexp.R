logsumexp <-
function(v){
  ## log(sum(exp(v)))
  v  = v[is.finite(v)]
  lv = length(v)
  if(lv==0){ return(-Inf) }
  if(lv==1){ return(v[1]) }
  
  w1  = which.max(v)
  res = sum(exp(v[-w1] - v[w1]))
  lse = v[w1] + log1p(res)
  lse
}
