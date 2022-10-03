SWAG_GS= function(S,burnin = round(S*.1),thin = 10,
                  save_all = 1,
                  mh.delta.star = .1){
  
  ###########################
  ## set hyper params
  S0 = eye(p); S0.inv = solve(S0) 
  V0 = S0; V0.inv = S0.inv
  
  eta0 = p + 2
  S1 = eye(p1); S1.inv = solve(S1)
  S2 = eye(p2); S2.inv= solve(S2)
  eta1 = p1 + 2; eta2 = p2 + 2
  
  nu.domain = c((p+2):(p+50)); ND = length(nu.domain)
  
  # priors
  N1 = eye(p1)
  N2 = eye(p2)
  xi1 = p1 + 2
  xi2 =  p2 + 2
  
  ## intialize values
  # level 1
  U = lapply(1:g,function(j)matrix(rnorm(ns[j]*p),ncol=p))
  Sig = Sig.inv = Psi = Psi.inv = lapply(1:g,function(j)eye(p))
  # level 2 and 3
  nu0 = p+2
  nu = p+2
  V0 = V0.inv = eye(p)
  V1 = lapply(1:g,function(j)rwish(S1,eta1))
  V1.inv = lapply(V1,solve)
  V2 = lapply(1:g,function(j)rwish(S2,eta2))
  V2.inv = lapply(V2,solve)
  # level 3
  K1 = K1.inv = eye(p1)
  K2 = K2.inv = eye(p2)
  S0 = S0.inv = eye(p)
  eta0 = p+2
  
  pis = .5

  len.out = length(seq(from = burnin+1, to = S, by=thin))

  ## create storage
  if (save_all == 0){
    # layer 1
    cov.out = cov.inv = array(NA,dim = c(len.out,p*p*g))
    nu.out = array(NA,dim = c(len.out,2))
    eta0.out = array(NA,dim = c(len.out,1))
    pis.out = array(NA,dim=c(len.out,1))
    acc = 0
    index = 1
  } else {
    # layer 1
    Psi.out = array(NA,dim = c(S,p*p*g))
    Sig.out = array(NA,dim = c(S,p*p*g))
    cov.inv = array(NA,dim = c(S,p*p*g))
    U.out   = array(NA,dim = c(S,sum(ns*p)))
    pis.out = array(NA,dim=c(S,1))
    acc = 0
    # layer 2
    nu.out = array(NA,dim = c(S,2))
    V0.out = array(NA,dim = c(S,p*p))
    V1.out = array(NA,dim = c(S,p1*p1*g))
    V2.out = array(NA,dim = c(S,p2*p2*g))
    # layer 3
    K1.out = array(NA,dim = c(S,p1*p1))
    K2.out = array(NA,dim = c(S,p2*p2))
    eta0.out = array(NA,dim = c(S,1))
  }
  
  ## GS
  for( s in 1:S ){
    
    
    # sample gamma/ nu corresponding to kronecker structure
    # sample nu kron - uniform on nu.domain - same as shrink to pooled procedure
    d.sig = matrix(NA,nrow=g,ncol=ND)
    for ( j in 1:g ){
      Vj.inv = kronecker(V2.inv[[j]],V1.inv[[j]])
      Vj = kronecker(V2[[j]],V1[[j]])
      d.sig[j,] = sapply(nu.domain,function(k)
        dmatT.propto(Y.list[[j]],k-p+1,pis^(1/2)*U[[j]],(1-pis) * Vj * (k-p-1),k,
                     Omega.inv = Vj.inv / (k-p-1) / (1-pis) ) )
    }
    # same nu for each sigj
    d.sig = apply(d.sig,2,sum)
    d.sig = exp(d.sig-max(d.sig))
    probs = d.sig/sum(d.sig)
    nu = sample(nu.domain,size = 1,prob = probs)
    
    ## sample xSigj
    for( j in 1:g ){
      V = kronecker(V2[[j]],V1[[j]])
      Y.tilde = Y.list[[j]] - pis^(1/2)*U[[j]]
      M = solve(V * (nu - p - 1) + t(Y.tilde) %*% Y.tilde / (1-pis))
      Sig.inv[[j]] = rwish(M,nu+ns[j]-1)
      Sig[[j]] = solve(Sig.inv[[j]])
    }
    
    ## propose pis.star from symmetric proposal
    if (s<1000){
      mh.delta = .5
      acc = 0
    } else {
      mh.delta = mh.delta.star
    }
    pi.hat = pis[1]
    pis.star = MH_sym_proposal_01(pi.hat, mh.delta)
    ## compute acceptance ratio, joint distn' of all variables and pi
    d.star = sapply(1:g,function(j)dmatnorm(Y.list[[j]],0,
                                            eye(ns[j]),
                                            pis.star * Psi[[j]] + (1-pis.star) * Sig[[j]],
                                            if_log = TRUE))
    
    d.s = sapply(1:g,function(j)dmatnorm(Y.list[[j]],0,
                                         eye(ns[j]),
                                         pi.hat * Psi[[j]] + (1-pi.hat) * Sig[[j]],
                                         if_log = TRUE))
    R = sum(d.star) - sum(d.s) +
      dbeta(pis.star,1/2,1/2,log=T) - dbeta(pi.hat,1/2,1/2,log=T) 
    ## if u<r, set pis to be pis.star
    if (log(runif(1))<R){
      pis = pis.star
      acc = acc + 1
    }
    
    
    ## sample Uj
    for ( j in 1:g ){
      W = solve(pis/(1-pis)*Sig.inv[[j]] + Psi.inv[[j]])
      M = pis^(1/2)/(1-pis) * Y.list[[j]] %*% Sig.inv[[j]] %*% W
      U[[j]] = rmatnorm(M,eye(ns[j]),W)
    }
    
    
    # sample nu0 - uniform on nu.domain - same as shrink to pooled procedure
    # integrate out psij
    d.sig.nu.domain = array(NA,dim=ND)
    for ( k in 1:length(nu.domain)){
      nu.star = nu.domain[k]
      d.sig = sapply(1:g,function(j) dmatT.propto(U[[j]],nu.star-p+1,0,V0*(nu.star - p - 1),nu.star,
                                                  Omega.inv = V0.inv/(nu.star - p - 1)) )
      d.sig.nu.domain[k] = sum(unlist(d.sig))
    }
    d.temp=exp(d.sig.nu.domain-max(d.sig.nu.domain))
    probs = d.temp/sum(d.temp)
    nu0 = sample(nu.domain,size = 1,prob = probs)
    
    
    
    ## sample Psijs
    for ( j in 1:g){
      M = solve(V0*(nu0 - p - 1) + t(U[[j]]) %*% U[[j]])
      Psi.inv[[j]] = rwish(M,nu0 + ns[j] - 1)
      Psi[[j]] = solve(Psi.inv[[j]])
    }
    
    
    ### wishart means
    
    ## sample V0
    M0 = solve(Reduce('+',Psi.inv)*(nu0-p-1) + S0.inv*eta0)
    V0 = rwish(M0,eta0 + nu0*g)
    V0.inv = solve(V0)
    
    ## sample g V1s
    Sig.inv.chol = lapply(Sig.inv,function(l)t(chol(l))) ### flipped so chol(S) = LL^T
    for ( j in 1:g ){
      L = array(Sig.inv.chol[[j]],dim=c(p1,p2,p)) # checked- should be correct
      helper = lapply(1:p,function(k) L[,,k] %*% t(V2[[j]]) %*% t(L[,,k]))
      M = solve(Reduce('+',helper)*(nu-p-1) + S1.inv * eta1)
      V1[[j]] = rwish(M,eta1 + nu * p2)
      V1.inv[[j]] = solve(V1[[j]])
    }
    
    ## sample g V2s 
    for ( j in 1:g ){
      L = array(Sig.inv.chol[[j]],dim=c(p1,p2,p))
      helper = lapply(1:p,function(k) t(L[,,k]) %*% t(V1[[j]]) %*% L[,,k])
      M = solve(Reduce('+',helper)*(nu-p-1) + S2.inv * eta2)
      V2[[j]] = rwish(M,eta2 + nu * p1)
      V2.inv[[j]] = solve(V2[[j]])
    }
    
    ### add layer to shrink homo to kron
    
    # sample K1
    V0.chol = t(chol(V0)) ### flipped so chol(S) = LL^T
    L = array(V0.chol,dim=c(p1,p2,p)) 
    helper = lapply(1:p,function(k) (L[,,k]) %*% t(K2.inv) %*% t(L[,,k]))
    M = solve(Reduce('+',helper)*eta0 + N1 * (xi1 - p1 - 1))
    K1.inv = rwish(M, p2*eta0 + xi1)
    K1 = solve(K1.inv)
    
    # sample K2
    helper = lapply(1:p,function(k) t(L[,,k]) %*% t(K1.inv) %*% (L[,,k]))
    M = solve(Reduce('+',helper)*eta0 + N2 * (xi2 - p2 - 1))
    K2.inv = rwish(M, p1*eta0 + xi2)
    K2 = solve(K2.inv)
    
    # format to be S0 for code
    S0 = kronecker(K2,K1)
    S0.inv = kronecker(K2.inv,K1.inv)
    
    # sample eta0
    d.v0 = sapply(nu.domain,function(k)
      dwish(V0,S0/k,k,if_log = T) )
    d.v0 = exp(d.v0-max(d.v0))
    probs = d.v0/sum(d.v0)
    eta0 = sample(nu.domain,size = 1,prob = probs)
    
    
    cov.inv.temp = cov.temp = list()
    for ( j in 1:g ){
      cov.temp[[j]] = (1-pis)*Sig[[j]]+pis*Psi[[j]]
      cov.inv.temp[[j]] = solve(cov.temp[[j]])
    }
    
    ## store output
    if (save_all == 0){
       if((s>burnin)&((s %% thin)==0)){
	  # layer 1
      	  cov.out[index,] = unlist(cov.temp)
      	  cov.inv[index,] = unlist(cov.inv.temp)
      	  nu.out[index,] = c(nu0,nu)
      	  eta0.out[index,] = eta0
      	  pis.out[index,] = pis
	  index = index + 1
      }
    } else {
      # layer 1
      Sig.out[s,]  = unlist(Sig)
      Psi.out[s,] = unlist(Psi)
      cov.inv[s,] = unlist(cov.inv.temp)
      U.out[s,] = unlist(U)
      pis.out[s,] = pis
      # layer 2
      V0.out[s,] = c(V0)
      V1.out[s,] = unlist(V1)
      V2.out[s,] = unlist(V2)
      nu.out[s,] = c(nu0,nu)
      # layer 3- shrink homo to kron
      K1.out[s,] = c(K1)
      K2.out[s,] = c(K2)
      eta0.out[s,] = eta0
    }
    
    
  }
  
  tosave.ind = seq(from =burnin, to = S, by = thin)
  if (save_all == 0){
    return(list("cov.out" = cov.out,
                "cov.inv" = cov.inv,"eta0" = eta0.out,
                "nu" = nu.out,"pis" = pis.out))
  } else {
    return(list("Sig" = Sig.out[tosave.ind,], "Psi" = Psi.out[tosave.ind,],
                "U" = U.out[tosave.ind,],
                "V0" = V0.out[tosave.ind,], "V1" = V1.out[tosave.ind,],
                "V2" = V2.out[tosave.ind,],"nu" = nu.out[tosave.ind,],
                "K1" = K1.out[tosave.ind,], "K2" = K2.out[tosave.ind,],
                "eta0" = eta0.out[tosave.ind,],"cov.inv" = cov.inv[tosave.ind,],
                "pis" = pis.out[tosave.ind],"acc" = acc/(S-500),
                "nu.domain" = nu.domain))
  }
  
}

####################################
########################################
cov.func = function(X){
  ## sample covariance matrix
  ## same as base r cov
  ## X: n x p matrix that we want a pxp column sample covariance matrix from
  
  X = as.matrix(X)
  N = nrow(X)
  
  X.temp = t(apply(X,1,function(KK) KK - colMeans(X)))
  
  t(X.temp) %*% X.temp / (N - 1)
}
cov.mle = function(X){
  ## MLE covariance matrix
  ## X: n x p matrix that we want a pxp column sample covariance matrix from
  
  X = as.matrix(X)
  N = nrow(X)
  
  X.temp = t(apply(X,1,function(KK) KK - colMeans(X)))
  
  t(X.temp) %*% X.temp / (N)
}
cov.pool = function(X,group){
  # sample pooled covariance matrix
  
  X = as.matrix(X)
  group = as.factor(group)
  N = length(group)
  
  sum.cov = tapply(seq_len(N),group,function(KK)
    length(KK) * cov.mle(X[KK,]))
  
  Reduce("+",sum.cov)/(N-length(unique(group)))
  
}
cov.pool.mle = function(X,group){
  # MLE pooled covariance matrix
  
  X = as.matrix(X)
  group = as.factor(group)
  N = length(group)
  
  sum.cov = tapply(seq_len(N),group,function(KK)
    length(KK) * cov.mle(X[KK,])) ## sum up Y^TY
  
  Reduce("+",sum.cov)/(N) ## MLE divides by total N- see MKB book
  
}
cov.shrink.pm = function(S,n,S0,nu,
                         de.meaned = F){
  ## Evaluate posterior mean from following hierarchical model:
  # Y ~ N_(nxp)(0,Sig)
  # Therefore, S = Y'Y ~ W(Sig,n)
  # Sig ~ IW(S0^(-1)/(nu-p-1),nu)
  # if de.meaned = true, use n-1 instead of n in degrees of freedom on S
  
  nu.star = nu + n
  if ( de.meaned == T ){
    nu.star = nu.star - 1
  }
  
  out = (S0 * (nu - p - 1) + S)/(nu.star - p - 1)
  
  return(out)
  
}
cov.shrink.steinBE = function(S,n,S0,nu,
                         de.meaned = F){
  ## Evaluate inverse of posterior mean of precision matrix from following hierarchical model
  ## (i.e. Bayes estimator under Stein loss):
  # Y ~ N_(nxp)(0,Sig)
  # Therefore, S = Y'Y ~ W(Sig,n)
  # Sig ~ IW(S0^(-1)/(nu-p-1),nu)
  # if de.meaned = true, use n-1 instead of n in degrees of freedom on S
  
  nu.star = nu + n
  if ( de.meaned == T ){
    nu.star = nu.star - 1
  }
  
  out = (S0 * (nu - p - 1) + S) / (nu.star)
  
  return(out)
  
}
cov.kron.mle = function(X,itmax = 100,eps = 1e-5,
                        de.meaned = F){
  ## block coordinate descent algorithm
  ## X: p1 x p2 x n array, de-meaned
  ## output: Psi p1 x p1 row covariance
  ## output: Sig p2 x p2 column covariance
  # if de.meaned = true, use n-1 instead of n in degrees of freedom on observed data sampling distn
  
  
  # 
  p1 = dim(X)[1]
  p2 = dim(X)[2]
  N = dim(X)[3] # since removing mean
  
  if (de.meaned == T){
    N = N-1
  }
  
  # initialize sig tilde
  Sig.tilde = matrix(rowMeans(apply(X,3,cov.mle)),ncol = p2)
  Sig.tilde.inv = qr.solve(Sig.tilde)
  
  # initialize stopping checks
  Psi.tilde = Psi.tilde.old = matrix(0,ncol = p1,nrow = p1)
  Sig.tilde.old = matrix(0,ncol = p2, nrow = p2)
  check = F
  it = 0
  while((check == FALSE) & (it < itmax)){
    
    Psi.tilde = matrix(rowSums(apply(X,3,function(KK)KK %*% Sig.tilde.inv %*% t(KK))),
                       ncol = p1)/(N*p2)
    Psi.tilde.inv = qr.solve(Psi.tilde)
    
    Sig.tilde = matrix(rowSums(apply(X,3,function(KK)t(KK) %*% Psi.tilde.inv %*% KK)),
                       ncol = p2)/(N*p1)
    Sig.tilde.inv = qr.solve(Sig.tilde)
    
    if (all(abs(Sig.tilde.old-Sig.tilde)<eps) &
        all(abs(Psi.tilde.old-Psi.tilde)<eps)){
      check = TRUE
    }
    
    # update for next iteration in while loop
    it = it+1 
    Psi.tilde.old = Psi.tilde
    Sig.tilde.old = Sig.tilde
    
  }
  
  return(list("Psi" = Psi.tilde,"Sigma" = Sig.tilde,
              "it" = it))
  
}

cov.kron.pool.mle = function(X,group,itmax = 100,eps = 1e-5,
                             de.meaned = F){
  ## block coordinate descent algorithm
  ## X: p1 x p2 x n array
  ## output: Psi p1 x p1 row covariance
  ## output: Sig p2 x p2 column covariance
  # if de.meaned = true, use n-1 instead of n in degrees of freedom on observed data sampling distn
  
  # params
  group = as.factor(group)
  
  p1 = dim(X)[1]
  p2 = dim(X)[2]
  N = dim(X)[3]; n = N
  
  if (de.meaned == TRUE){
    N = N - length(unique(group))
  }
  
  # initialize sig tilde
  Sig.tilde = matrix(rowMeans(apply(X,3,cov.mle)),ncol = p2)
  Sig.tilde.inv = qr.solve(Sig.tilde)
  
  # initialize stopping checks
  Psi.tilde = matrix(0,ncol = p1,nrow = p1)
  check = F
  it = 0
  Psi.tilde.old = matrix(0,ncol = p1,nrow = p1)
  Sig.tilde.old = matrix(0,ncol = p2, nrow = p2)
  while((check == FALSE) & (it < itmax)){
    
    Psi.tilde = tapply(seq_len(n),group,function(KK)
      matrix(rowSums(apply(X[,,KK],3,function(MM)MM %*% Sig.tilde.inv %*% t(MM))),
             ncol = p1)
      )
    Psi.tilde = Reduce("+",Psi.tilde)/(N*p2)
    Psi.tilde.inv = qr.solve(Psi.tilde)
    
    Sig.tilde = tapply(seq_len(n),group,function(KK)
      matrix(rowSums(apply(X[,,KK],3,function(MM)t(MM) %*% Psi.tilde.inv %*% MM)),
             ncol = p2)
    )
    Sig.tilde = Reduce("+",Sig.tilde)/(N*p1)
    Sig.tilde.inv = qr.solve(Sig.tilde)
    
    if (all(abs(Sig.tilde.old-Sig.tilde)<eps) &
        all(abs(Psi.tilde.old-Psi.tilde)<eps)){
      check = TRUE
    }
    
    # update for next iteration in while loop
    it = it+1 
    Psi.tilde.old = Psi.tilde
    Sig.tilde.old = Sig.tilde
    
  }
  
  return(list("Psi" = Psi.tilde,"Sigma" = Sig.tilde))
  
}
