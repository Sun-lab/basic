cluster.SCS <-
function(mutations, asrec, trec, nclusters, epsilon=0.001, 
         alpha1=0.02, beta1=0.8, alpha2=0.05, beta2=0.6, 
         theta = c(0.01, 0.01, 0.5, 0.1, 0.99, 0.01, 0.3, 0.4),
         maxIt=1000, converged=1e-5){
  
  # ---------------------------------------------------------
  # check the input data for single cells
  # ---------------------------------------------------------
  
  umuts = unique(c(mutations))
  
  if(! setequal(umuts, c(NA, 0, 1))){
    stop("mutations must take three values 0, 1, or NA")
  }
  
  if(! all(dim(mutations) == dim(asrec)) ){
    stop("dimensions of mutations and asrec do not match\n")
  }
  
  if(! all(dim(mutations) == dim(trec)) ){
    stop("dimensions of mutations and trec do not match\n")
  }
  
  if(any(c(asrec) > c(trec)) ){
    stop("asrec must be no larger than trec\n")
  }
  
  if(any(c(asrec) < 0) || max(abs(c(asrec) - round(c(asrec)))) > 1e-10 ){
    stop("asrec must non-negative integers\n")
  }
  
  if(any(c(trec) < 0) || max(abs(c(trec) - round(c(trec)))) > 1e-10 ){
    stop("trec must non-negative integers\n")
  }
  
  nn = nrow(mutations) # number of mutations
  mm = ncol(mutations) # number of cells
  
  # -----------------------------------------------------------------
  # Initialize Ehats. Fit a mixture distribution of three 
  # beta-binomials, with cluster 1 being no mutations, and clusters
  # 2 and 3 being cells with mutations in one or both alleles. 
  #
  # we fit the mixture distributions using all the data, 
  # and then re-estimate mixture proportion using those with 
  # missing genotype calls
  # -----------------------------------------------------------------
  
  nA     = c(asrec)
  nTotal = c(trec)
  is.NAs = which(c(is.na(mutations)))
  
  fit0 = bb.mix3(nA, nTotal, theta)
  fit1 = bb.mix3(nA[is.NAs], nTotal[is.NAs], fit0$theta, reEst=1)
  
  pi1  = fit1$theta["pi1"]
  rho1 = fit1$theta["rho1"]
  pi2  = fit1$theta["pi2"]
  rho2 = fit1$theta["rho2"]
  pi3  = fit1$theta["pi3"]
  rho3 = fit1$theta["rho3"]
  psi1 = fit1$theta["psi1"]
  psi2 = fit1$theta["psi2"]
  psi3 = 1 - psi1 - psi2
  
  psi2p = psi2/(psi2 + psi3)
  psi3p = psi3/(psi2 + psi3)
  
  ## impute Ds
  Ds = mutations
  pR.D0 = pR.D1 = matrix(NA, nrow=nrow(Ds), ncol=ncol(Ds))
  
  for(i in 1:nn){
    for(j in 1:mm){
      if(is.na(mutations[i,j])){
        if(trec[i,j] >= 3){
          f1  = dbetabinom(asrec[i,j], trec[i,j], prob=pi1, rho=rho1)
          f2  = dbetabinom(asrec[i,j], trec[i,j], prob=pi2, rho=rho2)
          f3  = dbetabinom(asrec[i,j], trec[i,j], prob=pi3, rho=rho3)
          Ds[i,j] = (psi2*f2 + psi3*f3)/(psi1*f1 + psi2*f2 + psi3*f3)
        }
        pR.D0[i,j] = f1
        pR.D1[i,j] = psi2p*f2 + psi3p*f3
      }
    }
  }

  for(i in 1:nn){
    w2fill = which(is.na(mutations[i,]) & trec[i,] < 3)
    if(length(w2fill) > 0){
      if(length(w2fill) == ncol(mutations)){
        Ds[i,w2fill] = 0.5
      }else{
        Ds[i,w2fill] = mean(Ds[i,-w2fill])
      }
    }
  }
  
  # ---------------------------------------------------------
  # perform clustering
  # ---------------------------------------------------------
  
  d1 = dist(Ds, method="manhattan")
  h1 = hclust(d1)
  
  # ========================================================================
  # iterate through the number of clusters
  # ========================================================================
  
  results = list()
  alpha.betas = c(alpha1, beta1, alpha2, beta2)
  
  for(nclust in nclusters){
    
    message(sprintf("nclust = %d, time: %s\n", nclust, date()))
    
    Ehats = Ds
    c1    = cutree(h1, k=nclust)
    
    alpha1 = alpha.betas[1]
    beta1  = alpha.betas[2]
    alpha2 = alpha.betas[3]
    beta2  = alpha.betas[4]
    
    # ---------------------------------------------------------
    # estimate the cluster means for single cell data: Es
    # ---------------------------------------------------------
    
    gammas = matrix(NA, nrow=mm, ncol=nclust)
    rownames(gammas) = colnames(mutations)
    
    for(j in 1:nclust){
      gammas[,j] = colMeans(Ehats[which(c1==j),])
    }
    gammas[which(gammas > 1-epsilon)] = 1 - epsilon
    gammas[which(gammas < epsilon)]   = epsilon
    
    nu = rep(1/nclust, nclust)
    
    # --------------------------------------------------------------
    # iteratively update parameters
    # --------------------------------------------------------------
    
    for(it in 1:maxIt){
      #     theta = c(alpha1, beta1, alpha2, beta2)
      #     names(theta) = c("alpha1", "beta1", "alpha2", "beta2")
      #     
      #     image(t(Ehats), main=it)
      #     cat(it, date(), "\n")
      #     print(cor(gammas))
      #     cat("\n")
      
      if(it %% 100 == 0){
        message(sprintf("iteration %s, time: %s\n", it, date()))
      }
      
      # cluster mutations
      log.postPs = matrix(NA, nrow=nn, ncol=nclust)
      
      for(i in 1:nn){
        for(k in 1:nclust){
          # likelihood for single cell data
          Lk.SC = Ehats[i,]*gammas[,k] + (1 - Ehats[i,])*(1-gammas[,k])
          log.postPs[i,k] = log(nu[k]) + sum(log(Lk.SC))
        }
      }
      
      log.postPs = log.postPs - apply(log.postPs, 1, logsumexp)
      postPs     = exp(log.postPs)
      
      nu = colSums(postPs)
      nu = nu/sum(nu)
      
      # estimate Es
      for(i in 1:nn){
        for(j in 1:mm){
          
          ## pD.E0 = p(observed data | E=0), pD.E1 = p(observed data | E=1)
          if(is.na(mutations[i,j])){
            if(trec[i,j] == 0){
              pD.E0 = pD.E1 = 1
            }else{
              pD.E0 = (1 - alpha2)*pR.D0[i,j] + alpha2*pR.D1[i,j]
              pD.E1 = (1 - beta2)*pR.D0[i,j]  + beta2*pR.D1[i,j]
            }
          }else{
            # here Ds[i,j] are observed value of 0 or 1
            pD.E0 = (1 - alpha1)*(1 - Ds[i,j]) + alpha1*Ds[i,j]
            pD.E1 = (1 - beta1)*(1 - Ds[i,j])  + beta1*Ds[i,j]
          }
          
          ## pE0 = p(E=0), pE1 = p(E=1)
          pE1 = sum(postPs[i,]*gammas[j,])
          pE0 = sum(postPs[i,]*(1 - gammas[j,]))
          
          Ehats[i,j] = pE1*pD.E1/(pE1*pD.E1 + pE0*pD.E0)
        }
      }
      
      theta.ab0 = c(alpha1, beta1, alpha2, beta2)
      names(theta.ab0) = c("alpha1", "beta1", "alpha2", "beta2")
      
      D.obs  = c(Ds[which(!is.na(mutations))])
      E.obs  = c(Ehats[which(!is.na(mutations))])
      
      alpha1  = sum(D.obs*(1 - E.obs))/sum(1-E.obs)
      beta1   = sum(D.obs*E.obs)/sum(E.obs)
      
      D.miss = c(Ds[which(is.na(mutations))])
      E.miss = c(Ehats[which(is.na(mutations))])
      
      alpha2  = sum(D.miss*(1 - E.miss))/sum(1-E.miss)
      beta2   = sum(D.miss*E.miss)/sum(E.miss)
      
      # restimate gammas
      for(j in 1:nclust){
        gammas[,j] = colSums(Ehats*postPs[,j])/sum(postPs[,j])
      }
      gammas[which(gammas > 1-epsilon)] = 1 - epsilon
      gammas[which(gammas < epsilon)]   = epsilon
      
      theta.ab = c(alpha1, beta1, alpha2, beta2)
      names(theta.ab) = c("alpha1", "beta1", "alpha2", "beta2")
      
      gap   = max(abs(theta.ab - theta.ab0))
      
      if(gap < converged){ break }
      
    }
    
    # --------------------------------------------------------------
    # calculate likelihood
    # --------------------------------------------------------------
    
    logLikE = 0
    
    # log-likelihood for E
    for(i in 1:nn){
      for(k in 1:nclust){
        # likelihood for single cell data
        Lk.SC   = Ehats[i,]*gammas[,k] + (1 - Ehats[i,])*(1-gammas[,k])
        logLikE = logLikE + log(nu[k]) + sum(log(Lk.SC))
      }
    }
    
    # log-likelihood for D or R
    logLikDR = 0
    
    for(i in 1:nn){
      for(j in 1:mm){
        
        ## pD.E0 = p(observed data | E=0), pD.E1 = p(observed data | E=1)
        if(is.na(mutations[i,j])){
          if(trec[i,j] == 0){
            pD.E0 = pD.E1 = 1
          }else{
            pD.E0 = (1 - alpha2)*pR.D0[i,j] + alpha2*pR.D1[i,j]
            pD.E1 = (1 - beta2)*pR.D0[i,j]  + beta2*pR.D1[i,j]
          }
        }else{
          # here Ds[i,j] are observed value of 0 or 1
          pD.E0 = (1 - alpha1)*(1 - Ds[i,j]) + alpha1*Ds[i,j]
          pD.E1 = (1 - beta1)*(1 - Ds[i,j])  + beta1*Ds[i,j]
        }
        
        logLikDR = logLikDR + log(pD.E0*(1 - Ehats[i,j]) + pD.E1*Ehats[i,j])
      }
    }
    
    logLikE
    logLikDR
    
    logL = logLikE + logLikDR
    # df = (nclust*mm) [for gamma] + nclust-1 [for nu]
    BIC  = -2*logL + log(nn*mm)*(nclust*mm + nclust - 1)
    
    # --------------------------------------------------------------
    # prepare output
    # --------------------------------------------------------------
    
    theta = c(alpha1, beta1, alpha2, beta2, pi1, rho1, pi2, rho2, pi3, rho3, 
              psi1, psi2)
    
    names(theta) = c("alpha1", "beta1", "alpha2", "beta2", "pi1", "rho1", 
                     "pi2", "rho2", "pi3", "rho3", "psi1", "psi2")
    
    lab1 = paste0("n", nclust)
    results[[lab1]] = list(E=Ehats, theta=theta, gammas=gammas, nu=nu, 
                           postPs=postPs, logLik=logL, BIC=BIC, nIt=it, 
                           gap=gap)
  }
  
  getItem <- function(v, label){
    if(is.null(v)){
      bic = NA
    }else{
      bic=v[[label]]
    }
    bic
  }
  
  sBIC = sapply(results,  getItem, label="BIC")
  sIt  = sapply(results,  getItem, label="nIt")
  sGap = sapply(results,  getItem, label="gap")
  
  sumClust = data.frame(nclust=nclusters, BIC=sBIC, nIt=sIt, gap=sGap)

  results$nclust.by.BIC = names(results)[which.min(sBIC)]
  results$summary.cluster = sumClust
  
  results
}

