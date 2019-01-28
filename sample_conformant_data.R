rm(list=ls())

library(readr)
library(MCMCpack)
library(mvtnorm)
library(BDgraph)

########################### Problem Dimensions
n = 100
p = 150
s = 10

############################ Select a set of n x p (SNPs) covariates

## The synthetic data in the paper use a subset of the real SNPs as covariates, 
# but as the NFBC66 dataset is confidential we'll use scrime to sample similar data
library(scrime)

x = simulateSNPs(n, p, c(3, 2),prop.explain = c(0.9, 0.95))$data
x = cbind(rep(1,n),x)

####################################################################

graph_pattern = 2 # in 2,3,4

snr = 25  # in 5,15,25
  
corr_param = 0.9 # in 0.3 , 0.6 , 0.9
    

### Create the underlying graph
if(graph_pattern==1){
  ### 1) Random but full
  G = matrix(1,s,s)
  Prime = list(c(1:s))
  Res = Prime
  Sep = list()
  
}else if(graph_pattern==2){
  
  ### 2) Block Diagonal structure
  Prime = list(c(1:floor(s*2/3)),
               c((floor(s*2/3)+1):(ceiling(s*4/5)-1)),
               c(ceiling(s*4/5):s))
  
  Res = Prime
  Sep = lapply(Res,function(x) which(x==-99))
  
  G = matrix(0,s,s)
  for(i in Prime){
    G[i,i] = 1
  }
  
}else if(graph_pattern==3){
  
  ### 3) Decomposable model
  Prime = list(c(1:floor(s*5/12),ceiling(s*9/10):s),
               c(floor(s*2/9):(ceiling(s*2/3)-1)),
               c(ceiling(s*2/3):(ceiling(s*4/5)-1)),
               c(ceiling(s*4/5):s))
  
  Sep = list(); H=list()
  for( i in 2:length(Prime)){
    H = union(H,Prime[[i-1]])
    Sep[[i-1]] = intersect( H,Prime[[i]])
  }
  
  Res = list()
  Res[[1]] = Prime[[1]]
  for( i in 2:length(Prime)){
    Res[[i]] = setdiff( Prime[[i]],Sep[[i-1]])
  }
  
  G = matrix(0,s,s)
  for(i in Prime)
    G[i,i] = 1
  
  ## decomp check
  dimnames(G) = list(1:s,1:s)
  length( gRbase::mcsMAT(G - diag(s)) ) > 0
  
  
}else if(graph_pattern==4){
  
  ### 4) Non-decomposable model
  nblocks = 5
  nElemPerBlock =c(floor(s/4),floor(s/2)-1-floor(s/4),
                   ceiling(s*2/3)-1-floor(s/2),7)
  nElemPerBlock = c(nElemPerBlock , s-sum(nElemPerBlock))
  res = 1:s; blockIdx = list()
  for(i in 1:nblocks){
    # blockIdx[[i]] = sample(res,nElemPerBlock[i])
    blockIdx[[i]] = res[1:nElemPerBlock[i]]
    res = setdiff(res,blockIdx[[i]])
  }
  
  
  G = matrix(0,s,s)
  ## add diagonal
  for(i in 1:nblocks)
    G[blockIdx[[i]],blockIdx[[i]]] = 1
  ## add cycle
  G[blockIdx[[1]],blockIdx[[2]]] = 1 ; G[blockIdx[[2]],blockIdx[[1]]] = 1
  G[blockIdx[[1]],blockIdx[[5]]] = 1 ; G[blockIdx[[5]],blockIdx[[1]]] = 1
  G[blockIdx[[2]],blockIdx[[3]]] = 1 ; G[blockIdx[[3]],blockIdx[[2]]] = 1
  G[blockIdx[[3]],blockIdx[[5]]] = 1 ; G[blockIdx[[5]],blockIdx[[3]]] = 1
  
  ## decomp check
  dimnames(G) = list(1:s,1:s)
  length( gRbase::mcsMAT(G -diag(s) ) ) > 0
  
  # Prime = blockIdx
  Res = blockIdx ## this is not correct but not used in the non-decomp case
  
}

### Gamma Pattern
gamma = matrix(0,p+1,s)
gamma[1,] = 1


### 2) Extra Patterns

## outcomes (correlated in the decomp model) have some predictors in common
gamma[6:10,6:9] = 1

## outcomes (correlated in the decomp model) have some predictors in common
gamma[16:20,14:15] = 1

## outcomes (sort-of correlated [pair-wise] in the decomp model)
# have predictors in common
gamma[26:30,6:15] = 1

## outcomes (NOT correlated in the decomp model) have predictors in common
gamma[36:40,c(3:5,16:17)] =1

## these predictors are associated with ALL the outcomes
gamma[46:50,] = 1

### ---

## get for every correlated bunch in the decomposable model,

if(graph_pattern<4){
  # a different set of predictors
  for(i in 1:length(Prime))
    gamma[6:10 + (i+6) * 10, Prime[[i]]] = 1   ## for each Prime component
  
  ## for every Residual instead
  for(i in 1:length(Res))
    gamma[6:10 + (i+11) * 10, Res[[i]]] = 1
  
}else{
  
  for(i in 1:length(Prime))
    gamma[6:10 + (i+4) * 10, Prime[[i]]] = 1   ## for each Prime component
  
  ## for every Residual instead
  for(i in 1:length(Res))
    gamma[6:10 + (i+9) * 10, Res[[i]]] = 1
}

#### Sample the betas
sd_b = 1
b = matrix(rnorm((p+1)*s,5,sd_b),p+1,s)


xb = matrix(NA,n,s)

for(i in 1:s){
  if(sum(gamma[,i])>1){
    xb[,i] = x[,gamma[,i]==1] %*% b[gamma[,i]==1,i]
  }else{
    xb[,i] = rep(1,n) * b[1,i]
  }
}

##Sample the variance
v_r = mean(diag(var(xb))) / snr

nu = s+1

M = matrix(corr_param,s,s)
diag(M) = rep(1,s)

P = BDgraph::rgwish(n=1,adj.g=G,b=3,D=v_r*M)

var = solve(P)

factor = 10 ; factor_min = 0.01; factor_max = 1000
count = 0 ; maxit = 10000

factor_prev = 1

repeat{
  
  var = var / factor * factor_prev
  
  ### Sample the errors and the Ys
  cVar = chol(as.matrix(var))
  err = matrix(rnorm(n*s),n,s) %*% cVar
  y = xb+err
  

  ## Reparametrisation ( assuming PEO is 1:s )
  cVar = t(cVar) # make it lower-tri
  S = diag(diag(cVar))
  sigma = S*S
  L = cVar %*% solve(S)
  rho =  diag(s) - solve(L)
  
  ### S/N Ratio
  emp_snr = mean( diag( var(xb) %*% solve(sigma) ))
  emp_g_snr = mean( diag( var( (err)%*%t(rho) ) %*% solve(sigma) ))
  
  ##############
  
  if( abs(emp_snr - snr) < (snr/10) | count > maxit ){
    break
  }else{
    if( emp_snr < snr ){ # increase factor
      factor_min = factor
    }else{ # decrease factor
      factor_max = factor
    }
    factor_prev = factor
    factor = (factor_min + factor_max)/2
  }
  count = count+1
}

#################
x = x[,-1]   # leave out the intercept because is coded inside already (variable standardised)
data = cbind(y,x)

####################################################################

## Write to file
name = paste("simul_",snr,"_",corr_param,"_",graph_pattern,".txt",sep="")
write.table(x=data,file=name,na = "NAN",col.names=FALSE,row.names=FALSE)

blockList = list(1:ncol(y),ncol(y)+(1:ncol(x)))
blockLabels = c(rep(0,ncol(y)),rep(1,ncol(x)))
write.table(x=blockLabels,file="blocks.txt",na = "NAN",col.names=FALSE,row.names=FALSE)

structureGraph = matrix(c(0,0,
                          1,0),2,2,byrow=TRUE)
write.table(x=structureGraph,file="structureGraph.txt",na = "NAN",col.names=FALSE,row.names=FALSE)

name = paste("simul_",snr,"_",corr_param,"_",graph_pattern,".RData",sep="")
save.image(name)

justName = paste("simul_",snr,"_",corr_param,"_",graph_pattern,sep="")

dir.create("results")
write("make\n",file="call.sh")

command = paste(sep="","./BVS_Reg --method SSUR --dataFile ",justName,".txt --blockFile blocks.txt --structureGraphFile structureGraph.txt --outFilePath results/ --nChains 2 --nIter 50000")
write(command,file="call.sh",append = TRUE)

command = paste(sep="","./BVS_Reg --method dSUR --dataFile ",justName,".txt --blockFile blocks.txt --structureGraphFile structureGraph.txt  --outFilePath results/ --nChains 2 --nIter 50000 # this implies dense residual covariance")
write(command,file="call.sh",append = TRUE)

command = paste(sep="","./BVS_Reg --method HESS --dataFile ",justName,".txt --blockFile blocks.txt --structureGraphFile structureGraph.txt  --outFilePath results/ --nChains 2 --nIter 50000\n")
write(paste(command,"\n"),file="call.sh",append = TRUE)

cat("Now run \nchmod +x call.sh && ./call.sh \n")
