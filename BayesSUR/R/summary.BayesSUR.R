#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title summary
#' @description
#' Summarize the results from fitted model
#' @importFrom Matrix Matrix
#' @name summary.BayesSUR
#' @param object fitted "runSUR" model
#' @param PmaxPredictor cutpoint for thresholding the estimated latent indicator variable. Default is 0.5
#' @param ... Other parameters in the function \code{summary.default.R} file
#' @export
summary.BayesSUR <- function(object, PmaxPredictor=0.5, ...){
  
  ans <- list(status=object$status)
  
  gamma <- as.matrix( read.table(paste(object$output$outFilePath,object$output$gamma,sep=""))>PmaxPredictor )
  ans$df <- sum(gamma)
  ans$gamma <- Matrix( gamma, sparse=TRUE )
  ans$beta <- Matrix( as.matrix(read.table(paste(object$output$outFilePath,object$output$beta,sep="")))[gamma], sparse=TRUE )
  ans$logP <- t( as.matrix( read.table(paste(object$output$outFilePath,object$output$logP,sep="")) ) )
  
  ans$chainParameters <- object$input[1:3]
  ans$modelParameters <- object$input[3:9]
  ans$hyperParameters <- object$hyperparameters
  ans$outputFiles <- object$output
  ans$outputFiles["outFilePath"] <- NULL
 
  return(ans) 
}