
cluster.BSC <-
function(b.cn, b.asrec, b.trec, mutations, asrec, trec, nclusters,  
         straify.by.maf=TRUE, stratify.n.intervals=20, stratify.window=2,  
         epsilon=0.001, alpha1=0.02, beta1=0.8, alpha2=0.05, beta2=0.6, 
         theta = c(0.01, 0.01, 0.5, 0.1,0.99, 0.01, 0.3, 0.4), 
         maxIt=1000, converged=1e-5){
  
  # ---------------------------------------------------------
  # check whether b.cn is a valid vector
  # ---------------------------------------------------------
  if(length(b.cn) != length(b.asrec)){
    stop("b.cn and b.asrec do not have the same length.\n")
  }
  
  if(any(b.cn <= 0) || max(abs(b.cn - round(b.cn))) > 1e-10 ){
    stop("b.cn must positive integers\n")
  }
  
  # ---------------------------------------------------------
  # check whether b.asrec and b.trec are valid vectors
  # ---------------------------------------------------------
  if(length(b.asrec) != length(b.trec)){
    stop("b.asrec and b.trec do not have the same length.\n")
  }
  
  if(length(b.asrec) != nrow(mutations)){
    stop("dimensions of b.asrec and mutations do not match.\n")
  }
  
  if(any(b.asrec > b.trec) ){
    stop("b.asrec must be no larger than b.trec\n")
  }
  
  if(any(b.asrec <= 0) || max(abs(b.asrec - round(b.asrec))) > 1e-10 ){
    stop("asrec must positive integers.\n")
  }
  
  if(any(b.trec <= 0) || max(abs(b.trec - round(b.trec))) > 1e-10 ){
    stop("trec must non-negative integers.\n")
  }
  
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
  
  # ---------------------------------------------------------
  # check parameter values
  # ---------------------------------------------------------
  
  rd1 = round(stratify.n.intervals)
  
  if(stratify.n.intervals < 0 || stratify.n.intervals != rd1){
    stop("stratify.n.intervals must be positive integers\n")
  }
  
  if(stratify.window < 0 || stratify.window >= stratify.n.intervals){
    stop("please choose stratify.step between 0 and stratify.n.intervals\n")
  }
  
  if(stratify.window != round(stratify.window)){
    stop("stratify.step must be integers\n")
  }
  
  # -----------------------------------------------------------------
  # Initialize Ehats. Fit a mixture distribution of three 
  # beta-binomials, with cluster 1 being no mutations, and clusters
  # 2 and 3 being cells with mutations in one or both alleles. 
  # There are two stategies:
  # (1) fit the mixture distribution using all the data
  # (2) fit the mixture distribution for each group of data
  #
  # we fit the mixture distributions using all the data,
  # and then re-estimate mixture proportion using those with
  # missing genotype calls
  # -----------------------------------------------------------------
  
  message(sprintf("initialize mutation calls in single cells\n"))
  
  mafs = (b.asrec/b.trec)*(b.cn/2)
  
  if(straify.by.maf){
    s1 = 1/stratify.n.intervals
    qs.start = c(0, seq(s1, 1-stratify.window*s1, by=s1))
    qs.end   = seq(stratify.window*s1, 1, by=s1)
    qs.center = (qs.start + qs.end)/2
    
    # calculate the center of each sliding window
    starts  = c(0, seq(s1))
    centers = quantile(mafs, probs=qs.center)
    dist1   = outer(mafs, centers, FUN="-")
    grps    = apply(abs(dist1), 1, which.min)
    table(grps)
    
    qs = c(0, quantile(mafs, probs=seq(s1, 1, by=s1)))
    
    thetas = matrix(NA, nrow=stratify.n.intervals-stratify.window+1, ncol=8)
    fits = list()
    
    message(sprintf("fit mixture of beta-binomials, stratified by maf\n"))
    
    for(j in 1:nrow(thetas)){
      message(".", appendLF=FALSE)
      
      wwj  = which(mafs > qs[j] & mafs <= qs[j+stratify.window])
      
      asrec.j = c(asrec[wwj,])
      trec.j  = c(trec[wwj,])
      is.na.j = which(c(is.na(mutations[wwj,])))
      
      # use all the data to estimate mixture components, and then
      # use data with missing mutation to esimate mixture proportions
      fit1 = bb.mix3(asrec.j, trec.j, theta)
      fit2 = bb.mix3(asrec.j[is.na.j], trec.j[is.na.j], fit1$theta, reEst=1)
      
      thetas[j,] = fit2$theta
      fits[[j]]  = fit2
    }
    message("")
    colnames(thetas) = names(fit1$theta)
    thetas
    
  }else{
    nA     = c(asrec)
    nTotal = c(trec)
    is.NAs = which(c(is.na(mutations)))
    
    fit0a = bb.mix3(nA, nTotal, theta)
    fit0b = bb.mix3(nA[is.NAs], nTotal[is.NAs], fit0a$theta, reEst=1)
    
  }

#   fit0b$theta
#   
#   par(mfrow=c(3,3), mar=c(2,4,1,1))
#   bb.mix3.plot(fit0b)
#   for(k in 1:8){
#     bb.mix3.plot(fits[[k]])
#   }
#   
#   par(mfrow=c(4,2), mar=c(2,4,1,1))
#   for(k in 1:8){
#     plot(thetas[,k], type="b", ylab=colnames(thetas)[k], ylim=c(0,1))
#     abline(h=fit0b$theta[k])
#   }

  # ---------------------------------------------------------
  # impute somatic mutation calls
  # ---------------------------------------------------------
  
  Ds = mutations
  pR.D0 = pR.D1 = matrix(NA, nrow=nrow(Ds), ncol=ncol(Ds))
  
  for(i in 1:nn){
    
    if(straify.by.maf){
      fit1 = fits[[grps[i]]]
    }else{
      fit1 = fit0b
    }

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
  # calculate distance of bulk tumor samples
  # ---------------------------------------------------------
  
  dist.bulk = matrix(NA, nrow=nn, ncol=nn)
  
  # y1=b.asrec[i]; z1=b.trec[i]; y2=b.asrec; z2=b.trec;
  # cn1=b.cn[i]; cn2=b.cn
  
  for(i in 1:(nn-1)){
    dist.bulk[i,] = dist.binom(b.asrec[i], b.trec[i], b.asrec, b.trec, 
                               b.cn[i], b.cn)
    dist.bulk[,i] = dist.bulk[i,]
  }
  
  diag(dist.bulk) = 0
  
  dim(dist.bulk)
  dist.bulk[1:5,1:5]
  
  d2 = dist.bulk
  d2[which(d2 < 1e-6)] = 1e-6
  d2 = log10(d2) + 6
  
  d1 = dist(Ds, method="manhattan")
  d1 = as.matrix(d1)
  summary(c(d1))
  summary(c(d2))
  
  # ---------------------------------------------------------
  # perform clustering
  # ---------------------------------------------------------
  
  d12 = d1 + d2
  h12 = hclust(as.dist(d12))
  
  # ========================================================================
  # iterate through the number of clusters
  # ========================================================================

  results = list()
  alpha.betas = c(alpha1, beta1, alpha2, beta2)

  for(nclust in nclusters){
    
    message(sprintf("nclust = %d, time: %s\n", nclust, date()))

    Ehats = Ds
    c12   = cutree(h12, k=nclust)
    
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
      gammas[,j] = colMeans(Ehats[which(c12==j),])
    }
    gammas[which(gammas > 1-epsilon)] = 1 - epsilon
    gammas[which(gammas < epsilon)]   = epsilon
    
    # ---------------------------------------------------------
    # estimate cluster means bulk tumors
    # ---------------------------------------------------------
    
    taus = rep(NA, nclust)
    
    # the derivative of the likelihood function
    ftau <- function(tau, ys, zs, cns, postPs){
      logL = sum(postPs*(ys*log(2*tau/cns) + (zs - ys)*log(1 - 2*tau/cns)))
      logL
    }
    
    for(j in 1:nclust){
      wwj = which(c12==j)
      ys  = b.asrec[wwj]
      zs  = b.trec[wwj]
      cns = b.cn[wwj]
      ps  = rep(1, length(wwj))
      
      uj = optimize(ftau, c(1e-6,1/2-1e-6), ys=ys, zs=zs, cns=cns, postPs=ps, 
                    maximum=TRUE)
      taus[j] = uj$maximum
    }
    
    taus0 = taus
    
    # ---------------------------------------------------------
    # iteratively update parameters
    # ---------------------------------------------------------
    
    nu = rep(1/nclust, nclust)
    
    for(it in 1:maxIt){
      #     theta = c(alpha1, beta1, alpha2, beta2)
      #     names(theta) = c("alpha1", "beta1", "alpha2", "beta2")
      #     
      #     image(t(Ehats), main=it)
      #     cat(it, date(), "\n")
      #     
      #     print(theta)
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
          # likelihood for bulk tumor data
          lk.B  = dbinom(b.asrec[i], b.trec[i], prob=2*taus[k]/b.cn[i], log=TRUE)
          
          log.postPs[i,k] = log(nu[k]) + sum(log(Lk.SC)) + lk.B
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
      
      # restimate taus
      for(j in 1:nclust){
        uj = optimize(ftau, c(1e-6,0.5-1e-6), ys=b.asrec, zs=b.trec, cns=b.cn, 
                      postPs=postPs[,j], maximum=TRUE)
        taus[j] = uj$maximum
      }
      
      theta.ab = c(alpha1, beta1, alpha2, beta2)
      names(theta.ab) = c("alpha1", "beta1", "alpha2", "beta2")
      
      gap   = max(abs(theta.ab - theta.ab0), abs(taus - taus0))
      if(gap < converged){ break }
      taus0 = taus
    }
    
    # --------------------------------------------------------------
    # calculate likelihood
    # --------------------------------------------------------------
    
    logLikE = logLikB = 0
    Lk.SC = Lk.BK = rep(NA, nclust)
    
    # log-likelihood for E of single cells and Bulk tumor
    for(i in 1:nn){
      for(k in 1:nclust){
        # likelihood for single cell data
        Lk.SC[k] = prod(Ehats[i,]*gammas[,k] + (1 - Ehats[i,])*(1-gammas[,k]))
        Lk.BK[k] = dbinom(b.asrec[i], b.trec[i], prob=2*taus[k]/b.cn[i])
      }
      logLikE = logLikE + log(sum(nu * Lk.SC))
      logLikB = logLikB + log(sum(nu * Lk.BK))
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
    
    logLikB
    logLikE
    logLikDR
    
    logL = logLikB + logLikE + logLikDR
    # df = (nclust*mm) [for gamma] + nclust-1 [for nu] + nclust-1 [for tau]
    BIC  = -2*logL + log(nn*(mm+1))*(nclust*mm + 2*nclust - 2)
    
    theta = c(alpha1, beta1, alpha2, beta2, pi1, rho1, pi2, rho2, pi3, rho3, 
              psi1, psi2)
    
    names(theta) = c("alpha1", "beta1", "alpha2", "beta2", "pi1", "rho1", 
                     "pi2", "rho2", "pi3", "rho3", "psi1", "psi2")
    
    lab1 = paste0("n", nclust)
    results[[lab1]] = list(E=Ehats, theta=theta, gammas=gammas, taus=taus, 
                           nu=nu, postPs=postPs, logLik=logL, BIC=BIC,
                           nIt=it, gap=gap)
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


