#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title Extract Model Fitted Values
#' @description
#' Extract model fitted values from the object of fitted Bayesian Seemingly Unrelated Regression
#' @name fitted
#' @param object fitted \code{runSUR} model
#' @param PmaxPredictor truncate the estimated coefficients based on thresholding the estimated latent indicator variable at 0 by default
#' @param ... other arguments
#' 
#' @return Fitted values extracted from the object \code{object}. If the \code{runSUR} specified data standardization, the fitted values are base based on standardized data.
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
#' ## check fitted values
#' fitted.val <- fitted(fit)
#' }
#' 
#' @export
fitted.BayesSUR <- function(object, PmaxPredictor=0, ...){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  
  gamma_hat <- as.matrix( read.table(object$output$gamma) )
  beta_hat <- as.matrix( read.table(object$output$beta) )
  X <- as.matrix( read.table(object$output$X,header=T) )
  
  if( "X0" %in% names(object$output) ){
    X0 <- as.matrix( read.table(object$output$X0) )
    beta_hat[-c(1:ncol(X0)),][gamma_hat<=PmaxPredictor] <- 0
    y.pred <- cbind(X0, X) %*% beta_hat
  }else{
    beta_hat[gamma_hat<=PmaxPredictor] <- 0
    y.pred <- X %*% beta_hat
  }
    
  return(y.pred)
  
}