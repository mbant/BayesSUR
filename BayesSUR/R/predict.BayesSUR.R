#' @title predict responses corresponding to the posterior mean of the coefficients, return posterior mean of coefficients or indices of nonzero coefficients
#' @description
#' Predict responses corresponding to the posterior mean of the coefficients, return posterior mean of coefficients or indices of nonzero coefficients of a "BayesSUR" class object.
#' @name predict.BayesSUR
#' 
#' @param object an object of class "BayesSUR"
#' @param newx Matrix of new values for x at which predictions are to be made. Must be a matrix
#' @param type Type of prediction required. Type "response" gives the fitted responses. Type "coefficients" computes the coefficients 
#' truncated the estimated coefficients based on thresholding the estimated latent indicator variable at \code{Pmax}. 
#' Type "nonzero" returns a list of the indices of the nonzero coefficients corresponding to the estimated latent indicator variable thresholding at \code{Pmax}
#' @param Pmax threshold that truncates the estimated coefficients based on thresholding the estimated latent indicator variable. Default is 0.
#' @param ... other arguments
#' 
#' @return Predicted values extracted from an object of class "BayesSUR". If the \code{BayesSUR} specified data standardization, the fitted values are base based on standardized data.
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
#' ## check prediction
#' predict.val <- predict(fit, newx=example_eQTL[["blockList"]][[2]])
#' 
#' @export
predict.BayesSUR <- function(object, newx, type=c("response", "coefficients", "nonzero"), Pmax=0, ...){
  
  type = match.arg(type)
  if (missing(newx) & type=="response"){
      stop("You need to supply a value for 'newx'!")
  }
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  
  gamma_hat <- as.matrix( read.table(object$output$gamma) )
  beta_hat <- as.matrix( read.table(object$output$beta) )
  X <- as.matrix( read.table(object$output$X,header=T) )
  
  if( "X0" %in% names(object$output) ){
    
    X0 <- as.matrix( read.table(object$output$X0) )
    beta_hat[-c(1:ncol(X0)),][gamma_hat<=Pmax] <- 0
    
    if( missing(newx) ){
      y.pred <- cbind(X0, X) %*% beta_hat
    }else{
      y.pred <- newx %*% beta_hat
    }
   
  }else{
    beta_hat[gamma_hat<=Pmax] <- 0
    if( missing(newx) ){
      y.pred <- X %*% beta_hat
    }else{
      y.pred <- newx %*% beta_hat
    }
  }
    
  if(type == "response")
    return(y.pred)
  if(type == "coefficients")
    return(beta_hat)
  if(type == "nonzero")
    return( which(gamma_hat<=Pmax, arr.ind=TRUE) )
  
}
