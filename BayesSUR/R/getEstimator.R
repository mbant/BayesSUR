#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title getEstimator
#' @description
#' Extract results from the object of fitted Bayesian Seemingly Unrelated Regression
#' @name getEstimator
#' @param object fitted "runSUR" model
#' @param estimator the name of one estimator. Default is the latent indicator estimator "gamma"
#' @export
getEstimator <- function(object, estimator="gamma"){
  
  if( sum(estimator %in% c("gamma","beta","G0","logP","model_size"))<1 )
    stop("Please specify a correct estimator!")
  if( estimator == "gamma" ) Est <- as.matrix( read.table(object$output$gamma) )
  if( estimator == "beta" ) Est <- as.matrix( read.table(object$output$beta) )
  if( estimator == "G0" ) Est <- as.matrix( read.table(object$output$G) )
  #if( est == "logP" ) Est <- t( as.matrix( read.table(object$output$logP) ) )
  #if( est == "model_size" ) Est <- rowSums(as.matrix( read.table(object$output$model_size) )) 
  
  return(Est)
  
}