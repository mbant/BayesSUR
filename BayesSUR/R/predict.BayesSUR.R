#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title make predictions from a "runSUR" object.
#' @description
#' Similar to other predict methods. This function predicts fitted values, coefficients and indexes of nonzero coefficients.
#' @name predict
#' 
#' @param object fitted \code{runSUR} model
#' @param newx Matrix of new values for x at which predictions are to be made. Must be a matrix
#' @param type Type of prediction required. Type "response" gives the fitted responses. Type "coefficients" computes the coefficients 
#' truncated the estimated coefficients based on thresholding the estimated latent indicator variable at \code{Pmax}. 
#' Type "nonzero" returns a list of the indices of the nonzero coefficients corresponding to the estimated latent indicator variable thresholding at \code{Pmax}
#' @param Pmax truncate the estimated coefficients based on thresholding the estimated latent indicator variable at 0 by default
#' @param ... other arguments
#' 
#' @return Predicted values extracted from the object \code{object}. If the \code{runSUR} specified data standardization, the fitted values are base based on standardized data.
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
#'                      hyperpar = hyperpar, tmpFolder = "tmp/" )
#' 
#' ## check prediction
#' predict.val <- predict(fit, newx=example_eQTL[["blockList"]][[2]])
#' }
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
