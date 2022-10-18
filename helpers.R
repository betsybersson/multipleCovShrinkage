#############################################
## general
#############################################

##########################################################
## Identity
# N: length of diagonal
##########################################################
eye = function(N){
  diag(rep(1,N))
}
##########################################################
## Trace of matrix
# X: square matrix
##########################################################
mat_trace = function(X){
  sum(diag(X))
}
##########################################################
## Covariance Stein Loss
# A: Estimator
# B: Truth
##########################################################
loss_stein = function(A,B){
  helper = A %*% qr.solve(B)
  ld_help = determinant(helper,logarithm = TRUE)$modulus[1]
  mat_trace(helper) - ld_help - nrow(A)
}
##########################################################
## Squared Error Loss
# A: Estimator
# B: Truth
##########################################################
loss_sq = function(A,B){
  mat_trace(crossprod(A-B))
}
##########################################################
## Covariance Quadratic Loss
# A: Estimator
# B: Truth
##########################################################
loss_quad = function(A,B){
  helper = A %*% qr.solve(B) - eye(nrow(A))
  mat_trace(helper %*% helper) 
}


#############################################
## distribution-related
#############################################
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
## Random Multivariate or Matrix-variate Normal Sample
# M: nxp mean of normal density
# U: nxn row covariance matrix
# V: pxp column covariance matrix
##########################################################
rmatnorm = function(M,U = eye(nrow(M)),V = eye(ncol(M))) {
  # function for density of Y_{ n\times p} \sim N(M, V kron U)
  N = nrow(M)
  P = ncol(M)
  
  E = matrix(rnorm(N*P),ncol=P)
  
  out = M + t(chol(U)) %*% E %*% chol(V)
  
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
    out = -N*P*log(2*pi)/2 - P*determinant(U)$mod[1]/2 - N*determinant(V)$mod[1]/2 - 
      mat_trace( V.inv %*% t(X-M) %*% U.inv %*% (X-M) )/2
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


#############################################
## SWAG specific
#############################################
vec.array = function(X){
  ## X: p1 x p2 x n array
  ## output: x: n x p1p2 matrix
  
  t(apply(X,3,c))
  
}

vec.inv.array = function(x,p1,p2){
  ## input: x: n x p1p2 matrix
  ## output: X: p1 x p2 x n array 
  
  N = nrow(x)
  XX = array(NA,dim=c(p1,p2,N))
  
  for ( nn in 1:N ){
    XX[,,nn] = matrix(x[nn,],ncol=p2,nrow=p1)
  }
  
  return(XX)
  
}
## not finished below- probably want some kind of function like this so I don't have to keep typing it out
list.to.3d.array = function(L){
  # input: list of length J containing matrices of dimension p1 x p2
  # output: array of dimension p1 x p2 x J 
  
  J = length(L)
  p1 = nrow(L[[1]])
  p2 = ncol(L[[1]])
  
  array.out = array(NA,dim=c(p1,p2,J))
  for ( jj in 1:J ){
    array.out[,,jj] = L[[jj]]
  }
  
  return(array.out)
  
}
rep.array = function(A,p){
  ## repeat matrix A p times
  ## output: X : array of dimension (dim(A),p)
  p1 = nrow(A)
  p2 = ncol(A)
  out.array = array(NA,dim=c(p1,p2,p))
  for(j in seq_len(p)){
    out.array[,,j] = A
  }
  return(out.array)
}
rep.vec = function(A,p){
  ## repeat vector A p times, cbind
  ## output: X : length(A),p
  n=length(A)
  out.mat = array(NA,dim=c(n,p))
  for(j in seq_len(p)){
    out.mat[,j] = A
  }
  return(out.mat)
}
lbind = function(L){
  # stack (rbind) matrices of dim njxp in list L of length g
  
  p = ncol(L[[1]])
  ns = unlist(lapply(L,nrow)); N = sum(ns)
  g = length(ns)
  
  mat.out = matrix(NA,ncol = p,nrow = N)
  group = rep(1:g,times =  ns)
  
  for ( j in 1:g){
    mat.out[group == j,] = as.matrix(L[[j]])
  }
  
  return(list("mat" = mat.out,"group" = group))
}

## copied directly from hoff amen
rwish <-function(S0,nu=dim(S0)[1]+2){
    # sample from a Wishart distribution 
    # with expected value nu*S0 
    sS0<-chol(S0)
    Z<-matrix(rnorm(nu*dim(S0)[1]),nu,dim(S0)[1])%*%sS0
    t(Z)%*%Z
}
