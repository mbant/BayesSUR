#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title getEstimator
#' @description
#' Summarize the results from fitted model
#' @name summary.BayesSUR
#' @param object fitted "runSUR" model
#' @param PmaxCovariate cutpoint for thresholding the estimated latent indicator variable. Default is 0.5
#' @export
summary.BayesSUR <- function(object, PmaxCovariate=0.5){
  
  ans <- list(status=object$status)
  
  library(Matrix)
  gamma <- as.matrix( read.table(object$output$gamma)>PmaxCovariate )
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