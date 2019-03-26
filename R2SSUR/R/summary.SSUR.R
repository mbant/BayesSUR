#' R2SSUR -- Bayesian Sparse Seemingly Unrelated Regression
#' I may be better to use this object instead of current "fit" from runRSSUR, and include the input filenames as well
#' @export
summary <- function(object){
  
  ans <- list(status=object$status)
  
  library(Matrix)
  gamma <- as.matrix( read.table(object$output$gamma)>0.5 )
  ans$df <- sum(gamma)
  ans$gamma <- Matrix( gamma, sparse=TRUE )
  ans$beta <- Matrix( as.matrix(read.table(object$output$beta))[gamma], sparse=TRUE )
  ans$sigmaRho <- as.matrix( read.table(object$output$sigmaRho) )
  ans$logP <- t( as.matrix( read.table(object$output$logP) ) )
  
  ans$chainParameters <- object$input[1:3]
  ans$modelParameters <- object$input[3:9]
  ans$hyperParameters <- object$hyperparameters
  ans$outputFiles <- object$output
 
  return(ans) 
}