#' @title fitted response values corresponds to the posterior mean estimates
#' @description
#' Return the fitted response values that correspond to the posterior mean estimates from a "BayesSUR" class object.
#' @name fitted.BayesSUR
#' @param object an object of class "BayesSUR"
#' @param Pmax threshold that truncates the estimated coefficients based on thresholding the estimated latent indicator variable. Default is 0.
#' @param ... other arguments
#' 
#' @return Fitted values extracted from an object of class "BayesSUR". If the \code{BayesSUR} specified data standardization, the fitted values are base based on standardized data.
#' 
#' @examples
#' data("example_eQTL", package = "BayesSUR")
#' hyperpar <- list( a_w = 2 , b_w = 5 )
#' 
#' set.seed(9173)
#' fit <- BayesSUR(Y = example_eQTL[["blockList"]][[1]], 
#'                 X = example_eQTL[["blockList"]][[2]],
#'                 data = example_eQTL[["data"]], outFilePath = tempdir(),
#'                 nIter = 100, burnin = 50, nChains = 2, gammaPrior = "hotspot",
#'                 hyperpar = hyperpar, tmpFolder = "tmp/" )
#' 
#' ## check fitted values
#' fitted.val <- fitted(fit)
#' 
#' @export
fitted.BayesSUR <- function(object, Pmax=0, ...){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  
  gamma_hat <- as.matrix( read.table(object$output$gamma) )
  beta_hat <- as.matrix( read.table(object$output$beta) )
  X <- as.matrix( read.table(object$output$X,header=T) )
  
  if( "X0" %in% names(object$output) ){
    X0 <- as.matrix( read.table(object$output$X0) )
    beta_hat[-c(1:ncol(X0)),][gamma_hat<=Pmax] <- 0
    y.pred <- cbind(X0, X) %*% beta_hat
  }else{
    beta_hat[gamma_hat<=Pmax] <- 0
    y.pred <- X %*% beta_hat
  }
    
  return(y.pred)
  
}