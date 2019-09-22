#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title getEstimator
#' @description
#' Extract results from the object of fitted Bayesian Seemingly Unrelated Regression
#' @name getEstimator
#' @param object fitted \code{runSUR} model
#' @param estimator the name of one estimator. Default is the latent indicator estimator "\code{gamma}". Other options "\code{beta}" and "\code{Gy}" correspond the posterior means of coefficient matrix, response graph and CPO, respectively 
#' @param Pmax truncate the estimator. Default is 0. If the estimator is beta, then beta is truncated based on the latent indicator matrix shresholding at \code{Pmax} 
#' 
#' @examples
#' \donttest{
#' data(example_eQTL, package = "BayesSUR")
#' hyperpar <- list( a_w = 2 , b_w = 5 )
#' 
#' fit <- runSUR(example_eQTL[["data"]], outFilePath = "results/",
#'                      Y = example_eQTL[["blockList"]][[1]],
#'                      X = example_eQTL[["blockList"]][[2]],
#'                      nIter = 1000, nChains = 2, gammaPrior = "hotspot",
#'                      hyperpar = hyperpar, tmpFolder="tmp/" )
#' 
#' ## check output
#' # extract the posterior mean of the coefficients matrix
#' beta_hat <- getEstimator(fit, estimator="beta")
#' }
#' 
#' @export
getEstimator <- function(object, estimator="gamma", Pmax=0){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  if( sum(estimator %in% c("gamma","beta","Gy","CPO"))<1 )
    stop("Please specify a correct estimator!")
  
  if( estimator == "gamma" ){
    Est <- as.matrix( read.table(object$output$gamma) )
    Est[Est<=0] <- 0
  } 
  
  if( estimator == "beta" ){
    Est <- as.matrix( read.table(object$output$beta) )
    Est[as.matrix( read.table(object$output$gamma) )<=0] <- 0
  } 
  
  if( estimator == "Gy" ){
    Est <- as.matrix( read.table(object$output$G) )
    Est[Est<=0] <- 0
  } 
  
  if( estimator == "CPO" ) Est <- as.matrix( read.table(object$output$CPO) )
  
  return(Est)
  
}