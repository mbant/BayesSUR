#' R2SSUR -- Bayesian Sparse Seemingly Unrelated Regression
#' Extract results from the object of fitted Bayesian Sparse Seemingly Unrelated Regression
#' @export
getEst <- function(object, est){
  
  if( est == "G0" ) Est <- as.matrix( read.table(object$output$G) )
  if( est == "gamma" ) Est <- as.matrix( read.table(object$output$gamma) )
  if( est == "beta" ) Est <- as.matrix( read.table(object$output$beta) )
  if( est == "sigmaRho" ) Est <- as.matrix( read.table(object$output$sigmaRho) )
  if( est == "logP" ) Est <- t( as.matrix( read.table(object$output$logP) ) )
  if( est == "model_size" ) Est <- t( as.matrix( read.table(object$output$model_size) ) )
  
  return(Est)
  
}