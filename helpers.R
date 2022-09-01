### formatting scripts
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

##########################################################
#### generic
##########################################################

##########################################################
## Identity
# N: length of diagonal
##########################################################
eye = function(N){
  diag(rep(1,N))
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
## Covariance Quadratic Loss
# A: Estimator
# B: Truth
##########################################################
loss_quad = function(A,B){
  helper = A %*% qr.solve(B) - eye(nrow(A))
  mat_trace(helper %*% helper) 
}
##########################################################
## Squared Error Loss
# A: Estimator
# B: Truth
##########################################################
loss_sq = function(A,B){
  mat_trace(crossprod(A-B))
}
