#' @title extract the posterior mean of the parameters
#' @description
#' Extract the posterior mean of the parameters of a "BayesSUR" class object.
#' @name getEstimator
#' @param object an object of class "BayesSUR"
#' @param estimator the name of one estimator. Default is the latent indicator estimator "\code{gamma}". Other options "\code{beta}", "\code{Gy}" and "\code{CPO}" 
#' correspond the posterior means of coefficient matrix, response graph and conditional predictive ordinate (CPO) respectively 
#' @param Pmax threshold that truncate the estimator. Default is 0. If the estimator is beta, then beta is truncated based on the latent indicator matrix shresholding at \code{Pmax} 
#' 
#' @return Return the one estimator from an object of class "BayesSUR". It is the posterior mean of the latent indicator variable if \code{estimator="gamma"}, posterior mean of the regression coefficients
#' if \code{estimator="beta"}, posterior mean of the response graph if \code{estimator="Gy"} and the CPO if \code{estimator="CPO"},
#' 
#' @examples
#' data("example_eQTL", package = "BayesSUR")
#' hyperpar <- list( a_w = 2 , b_w = 5 )
#' 
#' fit <- BayesSUR(Y = example_eQTL[["blockList"]][[1]], 
#'                 X = example_eQTL[["blockList"]][[2]],
#'                 data = example_eQTL[["data"]], outFilePath = tempdir(),
#'                 nIter = 100, burnin = 50, nChains = 2, gammaPrior = "hotspot",
#'                 hyperpar = hyperpar, tmpFolder = "tmp/" )
#' 
#' ## check output
#' # extract the posterior mean of the coefficients matrix
#' beta_hat <- getEstimator(fit, estimator="beta")
#' 
#' @export
getEstimator <- function(object, estimator="gamma", Pmax=0){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  if( sum(estimator %in% c("gamma","beta","Gy","CPO"))<1 )
    stop("Please specify a correct estimator!")
  
  if( estimator == "gamma" ){
    Est <- as.matrix( read.table(object$output$gamma) )
    Est[Est<=Pmax] <- 0
  } 
  
  if( estimator == "beta" ){
    Est <- as.matrix( read.table(object$output$beta) )
    Est[as.matrix( read.table(object$output$gamma) )<=Pmax] <- 0
  } 
  
  if( estimator == "Gy" ){
    Est <- as.matrix( read.table(object$output$G) )
    Est[Est<=Pmax] <- 0
  } 
  
  if( estimator == "CPO" ) Est <- as.matrix( read.table(object$output$CPO) )
  
  return(Est)
  
}