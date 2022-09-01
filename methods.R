multiple_shrinkage_GS = function(S,burnin = round(S*.1),thin = 10,
                                 save_all = 1,
                                 mh.delta = .1,
                                 single.weight = 1){
  
  ###########################
  ## set fixed hyper params
  ##### MAKE THIS EASIER FOR USER TO EDIT
  
  # priors on layer 3 hetero kronecker Wisharts
  S1 = eye(p1); S1.inv = solve(S1)
  S2 = eye(p2); S2.inv= solve(S2)
  eta1 = p1 + 2; eta2 = p2 + 2
  
  # priors on layer 4 homo kronecker Wisharts
  N1 = eye(p1) 
  N2 = eye(p2)
  xi1 = p1 + 2; xi2 =  p2 + 2
  
  # domain for degrees of freedom
  df.domain = c((p+2):(p+50)); DFD = length(df.domain)
  
  
  ## intialize values for GS
  # level 1
  U = lapply(1:g,function(j)matrix(rnorm(ns[j]*p),ncol=p))
  Sig = Sig.inv = Psi = Psi.inv = lapply(1:g,function(j)eye(p))
  # level 2 and 3
  nu0 = p+2
  gam = p+2
  V0 = V0.inv = eye(p)
  V1 = lapply(1:g,function(j)rwish(S1,eta1))
  V1.inv = lapply(V1,solve)
  V2 = lapply(1:g,function(j)rwish(S2,eta2))
  V2.inv = lapply(V2,solve)
  # level 4
  K1 = K1.inv = eye(p1)
  K2 = K2.inv = eye(p2)
  S0 = S0.inv = eye(p)
  eta0 = p+2
  
  w = rep(.5,g)
  
  
  ## create storage
  if (save_all == 0){
    # layer 1
    Psi.out = array(NA,dim = c(S,p*p*g))
    Sig.out = array(NA,dim = c(S,p*p*g))
    cov.inv = array(NA,dim = c(S,p*p*g))
    df.out = array(NA,dim = c(S,2))
    eta0.out = array(NA,dim = c(S,1))
    w.out = array(NA,dim=c(S,g))
  } else {
    # layer 1
    Psi.out = array(NA,dim = c(S,p*p*g))
    Sig.out = array(NA,dim = c(S,p*p*g))
    cov.inv = array(NA,dim = c(S,p*p*g))
    U.out   = array(NA,dim = c(S,sum(ns*p)))
    w.out = array(NA,dim=c(S,g))
    acc = rep(0,g)
    # layer 2
    df.out = array(NA,dim = c(S,2))
    V0.out = array(NA,dim = c(S,p*p))
    V1.out = array(NA,dim = c(S,p1*p1*g))
    V2.out = array(NA,dim = c(S,p2*p2*g))
    # layer 3
    K1.out = array(NA,dim = c(S,p1*p1))
    K2.out = array(NA,dim = c(S,p2*p2))
    eta0.out = array(NA,dim = c(S,1))
  }
  
  # helpers
  acc = rep(0,g)
  
  ## GS
  for( s in 1:S ){
    
    ## sample (gamma, {Sigj}) jointly
    # sample gamma; Sigjs are marginalized out
    d.sig = matrix(NA,nrow=g,ncol=DFD)
    for ( j in 1:g ){
      Vj.inv = kronecker(V2.inv[[j]],V1.inv[[j]])
      Vj = kronecker(V2[[j]],V1[[j]])
      d.sig[j,] = sapply(df.domain,function(k)
        dmatT.propto(Y.list[[j]],k-p+1,w[j]^(1/2)*U[[j]],(1-w[j]) * Vj * (k-p-1),k,
                     Omega.inv = Vj.inv / (k-p-1) / (1-w[j]) ) )
    }
    d.sig = apply(d.sig,2,sum)
    d.sig = exp(d.sig-max(d.sig))
    probs = d.sig/sum(d.sig)
    gam = sample(df.domain,size = 1,prob = probs)
    
    # sample Sigj
    for( j in 1:g ){
      V = kronecker(V2[[j]],V1[[j]])
      Y.tilde = Y.list[[j]] - w[j]^(1/2)*U[[j]]
      M = solve(V * (gam - p - 1) + t(Y.tilde) %*% Y.tilde / (1-w[j]))
      Sig.inv[[j]] = rwish(M,gam+ns[j]-1)
      Sig[[j]] = solve(Sig.inv[[j]])
    }
    
    ## sample (w,{Uj}) jointly
    # sample w; U is marginalized out
    # propose w.star from symmetric proposal
    if ( single.weight == 1 ){
      w.hat = w[1]
      w.star = MH_sym_proposal_01(w.hat, mh.delta)
      # compute acceptance ratio, joint distn' of all variables and w; U is marginalized out
      d.star = sapply(1:g,function(j)dmatnorm(Y.list[[j]],0,
                                              eye(ns[j]),
                                              w.star * Psi[[j]] + (1-w.star) * Sig[[j]],
                                              if_log = TRUE))
      
      d.s = sapply(1:g,function(j)dmatnorm(Y.list[[j]],0,
                                           eye(ns[j]),
                                           w.hat * Psi[[j]] + (1-w.hat) * Sig[[j]],
                                           if_log = TRUE))
      R = sum(d.star) - sum(d.s) +
        dbeta(w.star,1/2,1/2,log=T) - dbeta(w.hat,1/2,1/2,log=T) 
      # if u<r, set w to be w.star
      if (log(runif(1))<R){
        w = rep(w.star,g)
        acc[1] = acc[1] + 1
      }
    } else {
      for ( j in 1:g ){
        wj.hat = w[j]
        wj.star = MH_sym_proposal_01(wj.hat, mh.delta)
        # compute acceptance ratio, joint distn' of all variables and w; U is marginalized out
        d.star = dmatnorm(Y.list[[j]],0,
                          eye(ns[j]),
                          wj.star * Psi[[j]] + (1-wj.star) * Sig[[j]],
                          TRUE)
        
        d.s = dmatnorm(Y.list[[j]],0,
                       eye(ns[j]),
                       wj.hat * Psi[[j]] + (1-wj.hat) * Sig[[j]],
                       TRUE)
        R = (d.star) - (d.s) +
          dbeta(wj.star,1/2,1/2,log=T) - dbeta(wj.hat,1/2,1/2,log=T) 
        # if u<r, set wj to be wj.star
        if (log(runif(1))<R){
          w[j] = wj.star
          acc[j] = acc[j] + 1
        }
      }
      
    }
    
    ## sample Uj
    for ( j in 1:g ){
      W = solve(w[j]/(1-w[j])*Sig.inv[[j]] + Psi.inv[[j]])
      M = w[j]^(1/2)/(1-w[j]) * Y.list[[j]] %*% Sig.inv[[j]] %*% W
      U[[j]] = rmatnorm(M,eye(ns[j]),W)
    }
    
    
    ## sample (nu0,{Psij}) jointly
    # sample nu0; Psij is marginalized out
    d.sig.df.domain = array(NA,dim=DFD)
    for ( k in 1:length(df.domain)){
      nu.star = df.domain[k]
      d.sig = sapply(1:g,function(j) dmatT.propto(U[[j]],nu.star-p+1,0,V0*(nu.star - p - 1),nu.star,
                                                  Omega.inv = V0.inv/(nu.star - p - 1)) )
      d.sig.df.domain[k] = sum(unlist(d.sig))
    }
    d.temp=exp(d.sig.df.domain-max(d.sig.df.domain))
    probs = d.temp/sum(d.temp)
    nu0 = sample(df.domain,size = 1,prob = probs)
    
    # sample Psijs
    for ( j in 1:g){
      M = solve(V0*(nu0 - p - 1) + t(U[[j]]) %*% U[[j]])
      Psi.inv[[j]] = rwish(M,nu0 + ns[j] - 1)
      Psi[[j]] = solve(Psi.inv[[j]])
    }
    
    
    ### layer 3
    
    ## sample V0
    M0 = solve(Reduce('+',Psi.inv)*(nu0-p-1) + S0.inv*eta0)
    V0 = rwish(M0,eta0 + nu0*g)
    V0.inv = solve(V0)
    
    ## sample g V1s
    Sig.inv.chol = lapply(Sig.inv,function(l)t(chol(l))) ### flipped so chol(S) = LL^T
    for ( j in 1:g ){
      L = array(Sig.inv.chol[[j]],dim=c(p1,p2,p)) # checked- should be correct
      helper = lapply(1:p,function(k) L[,,k] %*% t(V2[[j]]) %*% t(L[,,k]))
      M = solve(Reduce('+',helper)*(gam-p-1) + S1.inv * eta1)
      V1[[j]] = rwish(M,eta1 + gam * p2)
      V1.inv[[j]] = solve(V1[[j]])
    }
    
    ## sample g V2s 
    for ( j in 1:g ){
      L = array(Sig.inv.chol[[j]],dim=c(p1,p2,p))
      helper = lapply(1:p,function(k) t(L[,,k]) %*% t(V1[[j]]) %*% L[,,k])
      M = solve(Reduce('+',helper)*(gam-p-1) + S2.inv * eta2)
      V2[[j]] = rwish(M,eta2 + gam * p1)
      V2.inv[[j]] = solve(V2[[j]])
    }
    
    ### layer 4
    
    ## sample K1
    V0.chol = t(chol(V0)) ### flipped so chol(S) = LL^T
    L = array(V0.chol,dim=c(p1,p2,p)) 
    helper = lapply(1:p,function(k) (L[,,k]) %*% t(K2.inv) %*% t(L[,,k]))
    M = solve(Reduce('+',helper)*eta0 + N1 * (xi1 - p1 - 1))
    K1.inv = rwish(M, p2*eta0 + xi1)
    K1 = solve(K1.inv)
    
    ## sample K2
    helper = lapply(1:p,function(k) t(L[,,k]) %*% t(K1.inv) %*% (L[,,k]))
    M = solve(Reduce('+',helper)*eta0 + N2 * (xi2 - p2 - 1))
    K2.inv = rwish(M, p1*eta0 + xi2)
    K2 = solve(K2.inv)
    
    ## format to be S0 for code
    S0 = kronecker(K2,K1)
    S0.inv = kronecker(K2.inv,K1.inv)
    
    ## sample eta0
    d.v0 = sapply(df.domain,function(k)
      dwish(V0,S0/k,k,if_log = T) )
    d.v0 = exp(d.v0-max(d.v0))
    probs = d.v0/sum(d.v0)
    eta0 = sample(df.domain,size = 1,prob = probs)
    
    # helper for output
    cov.inv.temp = list()
    for ( k in 1:g ){
      cov.inv.temp[[k]] = solve(w[j]*Psi[[k]] + (1-w[j])*Sig[[k]])
    }
    
    ## store output
    if (save_all == 0){
      # layer 1
      Sig.out[s,]  = unlist(Sig)
      Psi.out[s,] = unlist(Psi)
      cov.inv[s,] = unlist(cov.inv.temp)
      df.out[s,] = c(nu0,gam)
      eta0.out[s,] = eta0
      w.out[s,] = w
    } else {
      # layer 1
      Sig.out[s,]  = unlist(Sig)
      Psi.out[s,] = unlist(Psi)
      cov.inv[s,] = unlist(cov.inv.temp)
      U.out[s,] = unlist(U)
      w.out[s,] = w
      # layer 2
      V0.out[s,] = c(V0)
      V1.out[s,] = unlist(V1)
      V2.out[s,] = unlist(V2)
      df.out[s,] = c(nu0,gam)
      # layer 3- shrink homo to kron
      K1.out[s,] = c(K1)
      K2.out[s,] = c(K2)
      eta0.out[s,] = eta0
    }
    
    
  }
  
  tosave.ind = seq(from = burnin,to=S,by=thin)
  if (save_all == 0){
    return(list("Sig" = Sig.out[tosave.ind,], "Psi" = Psi.out[tosave.ind,],
                "cov.inv" = cov.inv[tosave.ind,],"eta0" = eta0.out[tosave.ind,],
                "df" = df.out[tosave.ind,],"w" = w.out[tosave.ind,]))
  } else {
    return(list("Sig" = Sig.out[tosave.ind,], "Psi" = Psi.out[tosave.ind,],
                "U" = U.out[tosave.ind,],
                "V0" = V0.out[tosave.ind,], "V1" = V1.out[tosave.ind,],
                "V2" = V2.out[tosave.ind,],"df" = df.out[tosave.ind,],
                "K1" = K1.out[tosave.ind,], "K2" = K2.out[tosave.ind,],
                "eta0" = eta0.out[tosave.ind,],"cov.inv" = cov.inv[tosave.ind,],
                "w" = w.out[tosave.ind,],"acc" = acc/(S-500),
                "df.domain" = df.domain))
  }
  
}


dmatT.propto = function(X,nu,M,Omega,nu0,
                        Omega.inv = qr.solve(Omega)){
  # output is log output
  
  p = ncol(X)
  n = nrow(X)
  
  # if X is a vector not a matrix
  if (is.null(p)){
    p = length(X)
    n=1
  }
  
  # help with huge numbers
  main.det = determinant(eye(n) + (X-M) %*% Omega.inv %*% t(X-M))$modulus
  # equal to main.det, probably faster but introduces complex numbers in case of numerical error: 
  # main.det = sum(log(eigen((X-M) %*% Omega.inv %*% t(X-M))$val + 1))
  gamma.div = gamma.mv(p,(nu+n+p-1)/2,T) - gamma.mv(p,(nu+p-1)/2,T)
  
  out = gamma.div  + (-p*n/2) * log(nu0-p-1) + 
    (-(nu+n+p-1)/2) * (main.det)
  
  return(out)
  
}
##########################################################
## Matrix Normal Density
# X: nxp matrix you want to evaluate density at
# M: nxp mean of normal density
# U: nxn row covariance matrix
# V: pxp column covariance matrix
##########################################################
dmatnorm = function(X,
                    M = matrix(0,ncol=ncol(X),nrow=nrow(X)),
                    U = eye(nrow(X)),
                    V = eye(ncol(X)),
                    U.inv = qr.solve(U),
                    V.inv = qr.solve(V),
                    if_log = FALSE){
  # function for density of X_{ n\times p} \sim N(M, V kron U)
  N = nrow(X)
  P = ncol(X)
  out = (2*pi)^(-N*P/2) * det(U.inv)^(P/2) * det(V.inv)^(N/2)*
    exp(-mat_trace( V.inv %*% t(X-M) %*% U.inv %*% (X-M) )/2)
  if (if_log == TRUE){
    out = -N*P*log(2*pi)/2 - P*log(det(U))/2 - N*log(det(V))/2 - 
      mat_trace( V.inv %*% t(X-M) %*% U.inv %*% (X-M) )/2
  }
  return(out)
}
##########################################################
## Multivariate Gamma Function
# p: dimension
# a: object
# output: Gamma_p(a)
##########################################################
gamma.mv = function(p,a,log.out = FALSE){
  out = pi^(p*(p-1)/4) * prod(gamma(a+(1-c(1:p))/2))
  if (log.out == TRUE){
    out = (p*(p-1)/4) * log(pi) + sum(lgamma(a+(1-c(1:p))/2))
  }
  return(out)
}
##########################################################
## Symmetric proposal on domain (0,1) for Metropolis algorithm
# given in Hoff page 190
# center: the previous value, where the proposal uniform distribution is centered at
# delta: neighborhood around "center" to draw uniformly from
# out: draw from symmetric proposal on domain (0,1)
##########################################################
MH_sym_proposal_01 = function(center, delta){
  out = runif(1,center-delta,center+delta)
  if (out<0){
    out = abs(out)
  } else if (out>1){
    out = 2-out
  }
  return(out)
}
##########################################################
## Wishart Density
# X: pxp matrix you want to evaluate density at
# M: pxp scale matrix
# nu: scaler degrees of freedom
##########################################################
dwish = function(X,M,nu, if_log = FALSE) {
  # function for density of X \sim W(M, nu)
  P = ncol(X)
  gamma.p = pi^(P*(P-1)/4) * prod(gamma( (nu+1-c(1:P))/2 ))
  out = 1/(2^(nu*P/2) * gamma.p * det(M)^(nu/2) ) * det(X)^((nu-P-1)/2) * 
    exp(mat_trace(-qr.solve(M)%*%X/2))
  if ( if_log == TRUE ){
    out = (-1) * ((nu*P/2)*log(2) + (P*(P-1)/4)*log(pi) + sum(lgamma( (nu+1-c(1:P))/2 )) +
                    (nu/2)*determinant(M)$mod) + ((nu-P-1)/2)*determinant(X)$mod +
      mat_trace(-qr.solve(M)%*%X/2)
  }
  return(out)
}