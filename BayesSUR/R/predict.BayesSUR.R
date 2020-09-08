#' @title predict responses corresponding to the posterior mean of the coefficients, return posterior mean of coefficients or indices of nonzero coefficients
#' @description
#' Predict responses corresponding to the posterior mean of the coefficients, return posterior mean of coefficients or indices of nonzero coefficients of a \code{BayesSUR} class object.
#' @name predict.BayesSUR
#' 
#' @param object an object of class \code{BayesSUR}
#' @param newx Matrix of new values for x at which predictions are to be made. Must be a matrix
#' @param type Type of prediction required. \code{type="response"} gives the fitted responses; \code{type="coefficients"} computes the coefficients 
#' truncated the estimated coefficients based on thresholding the estimated latent indicator variable at \code{Pmax};
#' \code{type="nonzero"} returns a list of the indices of the nonzero coefficients corresponding to the estimated latent indicator variable thresholding at \code{Pmax}
#' @param Pmax threshold that truncates the estimated coefficients based on thresholding the estimated latent indicator variable. Default is 0.
#' @param beta.type the type of estimated coefficients beta for prediction. Default is \code{marginal}, giving marginal beta estimation. If \code{beta.type="conditional"}, it gives conditional beta estimation
#' @param ... other arguments
#' 
#' @return Predicted values extracted from an object of class \code{BayesSUR}. If the \code{BayesSUR} specified data standardization, the fitted values are base based on standardized data.
#' 
#' @examples
#' data("exampleEQTL", package = "BayesSUR")
#' hyperpar <- list( a_w = 2 , b_w = 5 )
#' 
#' set.seed(9173)
#' fit <- BayesSUR(Y = exampleEQTL[["blockList"]][[1]], 
#'                 X = exampleEQTL[["blockList"]][[2]],
#'                 data = exampleEQTL[["data"]], outFilePath = tempdir(),
#'                 nIter = 100, burnin = 50, nChains = 2, gammaPrior = "hotspot",
#'                 hyperpar = hyperpar, tmpFolder = "tmp/" )
#' 
#' ## check prediction
#' predict.val <- predict(fit, newx=exampleEQTL[["blockList"]][[2]])
#' 
#' @export
predict.BayesSUR <- function(object, newx, type="response", Pmax=0, beta.type="marginal", ...){
  
  # type <- match.arg(type)
  # if (missing(newx) & type=="response"){
  #     stop("You need to supply a value for 'newx'!")
  # }
  if( length(type) > 1 ){
    warning("'type' has length > 1 and only the first element will be used")
    type <- type[1]
  }
  if( !(type %in% c("response", "coefficients", "nonzero")) )
    stop("Please specify correct 'type'!")
  if( Pmax<0 | Pmax>1 )
    stop("Please specify correct argument 'Pmax' in [0,1]!")
  if( !(beta.type %in% c("marginal", "conditional")) )
    stop("Please specify acorrect 'beta.type'!")
  
  if( (type %in% c("response", "coefficients")) & (Pmax>0) & (beta.type=="marginal"))
    stop("Pmax > 0 is valid only if the arguments type='coefficients' and beta.type='conditional'!")
  
  gamma_hat <- getEstimator( object, estimator="gamma", Pmax = Pmax, ...)
  beta_hat <- getEstimator( object, estimator="beta", Pmax = Pmax, beta.type=beta.type, ...)
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  X <- as.matrix( read.table(object$output$X,header=T) )
  if( "X0" %in% names(object$output) ){
    X0 <- as.matrix( read.table(object$output$X0) )
  }else{
    X0 <- NULL
  }
    
  if( missing(newx) ){
    y.pred <- cbind(X0, X) %*% beta_hat
  }else{
    y.pred <- newx %*% beta_hat
  }
  
  gamma_out <-  which(gamma_hat==1, arr.ind=TRUE)
  colnames(gamma_out) <- c("predictors", "response")
    
  if(type == "response")
    return( y.pred )
  if(type == "coefficients")
    return( beta_hat )
  if(type == "nonzero")
    return( gamma_out )
  
}
