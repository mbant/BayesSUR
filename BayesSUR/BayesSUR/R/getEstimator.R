#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title getEstimator
#' @description
#' Extract results from the object of fitted Bayesian Seemingly Unrelated Regression
#' @name getEstimator
#' @param object fitted "runSUR" model
#' @param estimator the name of one estimator. Default is the latent indicator estimator "gamma". Other options "beta" and "G0" correspond the posterior means of coefficient matrix and response graph, respectively 
#' @export
getEstimator <- function(object, estimator="gamma"){
  
  devAskNewPage(FALSE)
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  if( sum(estimator %in% c("gamma","beta","G0"))<1 )
    stop("Please specify a correct estimator!")
  if( estimator == "gamma" ) Est <- as.matrix( read.table(object$output$gamma) )
  if( estimator == "beta" ) Est <- as.matrix( read.table(object$output$beta) )
  if( estimator == "G0" ) Est <- as.matrix( read.table(object$output$G) )
  #if( est == "logP" ) Est <- t( as.matrix( read.table(object$output$logP) ) )
  #if( est == "model_size" ) Est <- rowSums(as.matrix( read.table(object$output$model_size) )) 
  
  return(Est)
  
}