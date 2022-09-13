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
