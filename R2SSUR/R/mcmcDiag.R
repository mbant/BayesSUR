#' R2SSUR -- Bayesian Sparse Seemingly Unrelated Regression
#' Diagnose the convergence of the MCMC iterations
#' @export
mcmcDiag <- function(x){
  
  logP <- getEst(x, "logP")
  nChain <- ncol(logP)
  model_size <- getEst(x, "model_size")
  par(mfrow=c(1,2))
  
  matplot(logP, xlab="Iteration", ylab="Log posterior", typ="l", lty=1, col=2:(nChain+1))
  legend("topleft", legend=paste("Chain",1:nChain), lty=1, col=2:(nChain+1), cex=0.5)
  
  matplot(model_size, xlab="Iteration", ylab="Model size", typ="l", lty=1, col=2:(nChain+1))
  legend("topleft", legend=paste("Chain",1:nChain), lty=1, col=2:(nChain+1), cex=0.5)
   
  par(mfrow=c(1,1))
  
  # we might also need other diagnostic plots here together, like temperature, etc.
  
}