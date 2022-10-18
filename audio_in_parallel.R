source("./cov_functions.R")
source("./helpers.R")
library(foreach)
library(doParallel)

## read data
df = readRDS("./temp/speech-MFCCs.rds")

## dimensions
ns = unlist(lapply(df,function(j)dim(j)[1]))
p1 = dim(df[[1]])[2]
p2 = dim(df[[1]])[3]
p = p1*p2
J = length(df)


# ### limit J
# J = 3
# ns = ns[1:J]
# df = df[1:3]


n.test = 100
set.seed(123)
tosave = sapply(ns,function(j)sample(j,n.test))
set.seed(Sys.time())

vec.3d = function(XX){
  ## input: XX of dim (n,p1,p2)
  ## output: matrix of dim (n,p1*p2); 
  ## (columns p2 stacked on top of each other)
  t(apply(XX,1,c))
}

Y.test = list()
Y.list = list()
for( j in 1:J){
  Y.test[[j]] = vec.3d(df[[j]][tosave[,j],,])
  Y.list[[j]] = vec.3d(df[[j]][-tosave[,j],,])
}

# update
ns = ns - n.test

Y.list = lapply(Y.list,scale)

g = J

## options
S = 3000
burnin = 500
thin = 4
mh.delta.star = .1


###########################
## set hyperx params
S0 = eye(p); S0.inv = solve(S0) 
V0 = S0; V0.inv = S0.inv

eta0 = p + 2
S1 = eye(p1); S1.inv = solve(S1)
S2 = eye(p2); S2.inv= solve(S2)
eta1 = p1 + 2; eta2 = p2 + 2

# nu.domain = c((p+2):(p+50))
nu.domain = round(seq(from = (p+2), to = (2*p), length.out = 10) )
ND = length(nu.domain)

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
# layer 1
cov.out = cov.inv = array(NA,dim = c(len.out,p*p*g))
nu.out = array(NA,dim = c(len.out,2))
eta0.out = array(NA,dim = c(len.out,1))
pis.out = array(NA,dim=c(len.out,1))
index = 1


# helper
S.list = list()
for ( j in 1:g ){
  S.list[[j]] = t(Y.list[[j]]) %*% Y.list[[j]]
}

mh.delta = .5


# parallel shit
#setup parallel backend to use many processors
cores=detectCores()


## GS
tic1 = Sys.time()
for( s in 1:S ){
  
  tic = Sys.time()
  
  
  # sample gamma/ nu corresponding to kronecker structure
  # sample nu kron - uniform on nu.domain - same as shrink to pooled procedure
  ## parallel step
  cl <- makeCluster(min(cores[1] - 1,g))  # dontoverload your computer
  registerDoParallel(cl)
  
  out.temp = array(NA,dim=ND)
  
  d.sig <- foreach(j = 1:g, .combine='rbind') %dopar% {
    
    Vj.inv = kronecker(V2.inv[[j]],V1.inv[[j]])
    Vj = kronecker(V2[[j]],V1[[j]])
    
    helper = (Y.list[[j]] - pis^(1/2)*U[[j]]) %*% Vj.inv %*% 
      t(Y.list[[j]] - pis^(1/2)*U[[j]]) / (1-pis)
    det.helper = det(helper)
    trace.helper = mat_trace(helper)
    trace.helper.sq = mat_trace(helper %*% t(helper))

    for ( k in 1:ND){
      
      NU = nu.domain[k]-p+1
      NU0 = nu.domain[k]
      
      # main.det = determinant(eye(ns[j]) + helper/ (NU0-p-1))$modulus
      
      if ((NU0-p-1)<8){
        main.det = determinant(eye(ns[j]) + helper/ (NU0-p-1))$modulus
      } else { # else, use approximation that only requires det of helper
        main.det = log(1 + det.helper + trace.helper/(NU0-p-1) +
                         trace.helper^2/2/(NU0-p-1)^2 -
                         trace.helper.sq/2/(NU0-p-1)^2)
      }
      
      gamma.div = gamma.mv(p,(NU+ns[j]+p-1)/2,T) - gamma.mv(p,(NU+p-1)/2,T)
      
      out.temp[k] = gamma.div  + (-p*ns[j]/2) * log(NU0-p-1) + 
        (-(NU+ns[j]+p-1)/2) * (main.det)

    }
    
    out.temp
    
  }
  
  #stop cluster
  stopCluster(cl)

  # same nu for each sigj
  d.sig = apply(d.sig,2,sum)
  d.sig = exp(d.sig-max(d.sig))
  probs = d.sig/sum(d.sig)
  nu = sample(nu.domain,size = 1,prob = probs)
  
  
  ## sample Sigj
  cl <- makeCluster(min(cores[1]-1,g))  # dontoverload your computer
  registerDoParallel(cl)
  
  finalMatrix <- foreach(j = 1:g, .combine='cbind') %dopar% {
    V = kronecker(V2[[j]],V1[[j]])
    Y.tilde = Y.list[[j]] - pis^(1/2)*U[[j]]
    M = qr.solve(V * (nu - p - 1) + t(Y.tilde) %*% Y.tilde / (1-pis))
    temp1 = rwish(M,nu+ns[j]-1)
    temp2 = qr.solve(temp1)
    
    list(temp1,temp2)
  }
  
  #stop cluster
  stopCluster(cl)
  
  Sig.inv = finalMatrix[1,]
  Sig = finalMatrix[2,]
  

  
  
  ## propose pis.star from symmetric proposal
  if (s>1000){
    mh.delta = mh.delta.star
  }
  pi.hat = pis[1]
  pis.star = MH_sym_proposal_01(pi.hat, mh.delta)
  ## compute acceptance ratio, joint distn' of all variables and pi
  
  # start cluster
  cl <- makeCluster(min(cores[1]-1,g))  # dontoverload your computer
  registerDoParallel(cl)
  
  finalMatrix <- foreach(j = 1:g, .combine='cbind') %dopar% {
    
    d.star = dmatnorm_edited(S.list[[j]],ns[j],
                             pis.star * Psi[[j]] + (1-pis.star) * Sig[[j]])
    d.s = dmatnorm_edited(S.list[[j]],ns[j],
                          pi.hat * Psi[[j]] + (1-pi.hat) * Sig[[j]])
    c(d.star,d.s)
  }
  
  #stop cluster
  stopCluster(cl)
  
  d.star = finalMatrix[1,]
  d.s = finalMatrix[2,]
  
  R = sum(d.star) - sum(d.s) +
    dbeta(pis.star,1/2,1/2,log=T) - dbeta(pi.hat,1/2,1/2,log=T) 
  ## if u<r, set pis to be pis.star
  if (log(runif(1))<R){
    pis = pis.star
  }
  

  
  ## sample Uj

  # start cluster
  cl <- makeCluster(min(cores[1]-1,g))  # dontoverload your computer
  registerDoParallel(cl)
  
  U <- foreach(j = 1:g, .combine='c',.packages = c('amen')) %dopar% {
    
    W = solve(pis/(1-pis)*Sig.inv[[j]] + Psi.inv[[j]])
    M = pis^(1/2)/(1-pis) * Y.list[[j]] %*% Sig.inv[[j]] %*% W
    
    out = rmatnorm(M,eye(ns[j]),W)
    
    list(out)
  }
  
  #stop cluster
  stopCluster(cl)


  
  
  # sample nu0 - uniform on nu.domain - same as shrink to pooled procedure
  # integrate out psij
  # start cluster
  cl <- makeCluster(min(cores[1]-1,g))  # dontoverload your computer
  registerDoParallel(cl)
  
  d.sig.nu.domain <- foreach(j = 1:g, .combine='rbind') %dopar% {

    helper = U[[j]] %*% V0.inv %*% t(U[[j]])
    det.helper = det(helper)
    trace.helper = mat_trace(helper)
    trace.helper.sq = mat_trace(helper %*% t(helper))

    for ( k in 1:ND){
      
      NU0 = nu.domain[k]
      NU = NU0-p+1
      
      # main.det = determinant(eye(ns[j]) + helper/(NU0 - p - 1))$modulus
      if ((NU0-p-1)<8){
        main.det = determinant(eye(ns[j]) + helper/ (NU0-p-1))$modulus
      } else { # else, use approximation that only requires det of helper
        main.det = log(1 + det.helper + trace.helper/(NU0-p-1) +
                         trace.helper^2/2/(NU0-p-1)^2 -
                         trace.helper.sq/2/(NU0-p-1)^2)
      }
      
      gamma.div = gamma.mv(p,(NU+ns[j]+p-1)/2,T) - gamma.mv(p,(NU+p-1)/2,T)
      
      out.temp[k] = gamma.div  + (-p*ns[j]/2) * log(NU0-p-1) + 
        (-(NU+ns[j]+p-1)/2) * (main.det)
      
    }
    
    out.temp
    
  }
  
  #stop cluster
  stopCluster(cl)
  
  d.sig.nu.domain = colSums(d.sig.nu.domain)
  
  d.temp=exp(d.sig.nu.domain-max(d.sig.nu.domain))
  probs = d.temp/sum(d.temp)
  nu0 = sample(nu.domain,size = 1,prob = probs)
  

  
  
  cl <- makeCluster(min(cores[1]-1,g))  # dontoverload your computer
  registerDoParallel(cl)
  
  finalMatrix <- foreach(j = 1:g, .combine='cbind') %dopar% {
    M = qr.solve(V0*(nu0 - p - 1) + t(U[[j]]) %*% U[[j]])

    temp1 = rwish(M,nu0 + ns[j] - 1)
    temp2 = qr.solve(temp1)
    
    list(temp1,temp2)
  }
  
  #stop cluster
  stopCluster(cl)
  
  Psi.inv = finalMatrix[1,]
  Psi = finalMatrix[2,]

  
  
  ### wishart means
  
  ## sample V0
  M0 = solve(Reduce('+',Psi.inv)*(nu0-p-1) + S0.inv*eta0)
  V0 = rwish(M0,eta0 + nu0*g)
  V0.inv = qr.solve(V0)
  
  ## sample g V1s
  for ( j in 1:g ){
    
    Sig.inv.chol = t(chol(Sig.inv[[j]]))
    
    L = array(Sig.inv.chol,dim=c(p1,p2,p)) # checked- should be correct
    
    helper = lapply(1:p,function(k) L[,,k] %*% t(V2[[j]]) %*% t(L[,,k]))
    M = qr.solve(Reduce('+',helper)*(nu-p-1) + S1.inv * eta1)
    V1[[j]] = rwish(M,eta1 + nu * p2)
    V1.inv[[j]] = qr.solve(V1[[j]])
  
    helper = lapply(1:p,function(k) t(L[,,k]) %*% t(V1[[j]]) %*% L[,,k])
    M = qr.solve(Reduce('+',helper)*(nu-p-1) + S2.inv * eta2)
    V2[[j]] = rwish(M,eta2 + nu * p1)
    V2.inv[[j]] = qr.solve(V2[[j]])
  }
  
  ### add layer to shrink homo to kron
  
  # sample K1
  V0.chol = t(chol(V0)) ### flipped so chol(S) = LL^T
  L = array(V0.chol,dim=c(p1,p2,p)) 
  
  helper = lapply(1:p,function(k) (L[,,k]) %*% t(K2.inv) %*% t(L[,,k]))
  M = qr.solve(Reduce('+',helper)*eta0 + N1 * (xi1 - p1 - 1))
  K1.inv = rwish(M, p2*eta0 + xi1)
  K1 = qr.solve(K1.inv)
  
  # sample K2
  helper = lapply(1:p,function(k) t(L[,,k]) %*% t(K1.inv) %*% (L[,,k]))
  M = qr.solve(Reduce('+',helper)*eta0 + N2 * (xi2 - p2 - 1))
  K2.inv = rwish(M, p1*eta0 + xi2)
  K2 = qr.solve(K2.inv)
  
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
    cov.inv.temp[[j]] = qr.solve(cov.temp[[j]])
  }
  
  ## store output
  if((s>burnin)&((s %% thin)==0)){
    # layer 1
    cov.out[index,] = unlist(cov.temp)
    cov.inv[index,] = unlist(cov.inv.temp)
    nu.out[index,] = c(nu0,nu)
    eta0.out[index,] = eta0
    pis.out[index,] = pis
    index = index + 1
  }
  
  toc = Sys.time() - tic
  
  if ((s%%thin)==0){
    print(paste0("Iteration ,",s," complete."))
    print(paste0("Time for last iteration:",toc))
    print("---------------")
    print(paste0(cov.temp[[2]][1:3,1:3]))
    print("---------------")
  }
  
  
}

toc1 = Sys.time() - tic1


# put output in list
model = list("cov.out" = cov.out,
             "cov.inv" = cov.inv,"eta0" = eta0.out,
             "nu" = nu.out,"pis" = pis.out)


### save output
saveRDS(model,file = "./output/audio_S3000.RDS")

