#' @title extract the posterior mean of the parameters
#' @description
#' Extract the posterior mean of the parameters of a \code{BayesSUR} class object.
#' @name get.estimator
#' @param object an object of class \code{BayesSUR}
#' @param estimator the name of one estimator. Default is the latent indicator estimator "\code{gamma}". Other options "\code{beta}", "\code{Gy}", "\code{CPO}" and "\code{logP}" 
#' correspond the marginal (conditional) coefficient matrix if \code{beta.type="marginal"}(\code{"conditional"}), response graph and conditional predictive ordinate (CPO) respectively 
#' @param Pmax threshold that truncate the estimator "\code{gamma}" or "\code{Gy}". Default is \code{0}
#' @param beta.type the type of output beta. Default is \code{marginal}, giving marginal beta estimation. If \code{beta.type="conditional"}, it gives conditional beta estimation
#' 
#' @return Return the estimator from an object of class \code{BayesSUR}. It is a matrix if the length of argument \code{marginal} is greater than 1. Otherwise, it is a list
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
#' ## check output
#' # extract the posterior mean of the coefficients matrix
#' beta_hat <- get.estimator(fit, estimator="beta")
#' 
#' @export
get.estimator <- function(object, estimator = "gamma", Pmax = 0, beta.type = "marginal"){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  if( sum(!estimator %in% c("gamma","beta","Gy","CPO","logP"))>0 ){
    stop("Please specify correct estimator!")
  }else{
    ret <- rep(list(NULL), length(estimator))
    names(ret) <- estimator
  }
  
  if( "gamma" %in% estimator ){
    ret$gamma <- as.matrix( read.table(object$output$gamma) )
    if(Pmax > 0)
      ret$gamma[ret$gamma<=Pmax] <- 0
    rownames(ret$gamma) <- colnames(read.table(object$output$X,header=T))
    colnames(ret$gamma) <- colnames(read.table(object$output$Y,header=T))
  } 
    
    if( "beta" %in% estimator ){
      ret$beta <- as.matrix( read.table(object$output$beta) )
      
      if( sum(beta.type %in% c("marginal", "conditional"))>0 ){
        if( beta.type == "conditional" ){
          gammas <- as.matrix( read.table(object$output$gamma) )
          ret$beta <- (gammas>=Pmax)*ret$beta/gammas
          ret$beta[is.na(ret$beta)] <- 0
        }
      }else{
        stop("Please specify correct beta.type!")
      }
      
      colnames(ret$beta) <- colnames(read.table(object$output$Y,header=T))
      if("X_0" %in% names(object$output)){
        rownames(ret$beta) <- c(colnames(read.table(object$output$X_0,header=T)), colnames(read.table(object$output$X,header=T)))
      }else{
        rownames(ret$beta) <- colnames(read.table(object$output$X,header=T))
      }
    } 
    
    if( "Gy" %in% estimator ){
      
      covariancePrior <- object$input$covariancePrior
      if(covariancePrior == "HIW"){
        ret$Gy <- as.matrix( read.table(object$output$G) )
      }else{
        stop("Gy is only estimated with hyper-inverse Wishart prior for the covariance matrix of responses!")
      }
      if( Pmax > 0)
        ret$Gy[ret$Gy<=Pmax] <- 0
      rownames(ret$Gy) <- colnames(ret$Gy) <- names(read.table(object$output$Y,header=T))
      
    } 
    
    if( "CPO" %in% estimator ){
      if(is.null(object$output$CPO))
        stop("Please specify argument output_CPO in BayesSUR()!")
      
      ret$CPO <- as.matrix( read.table(object$output$CPO) )
      
      rownames(ret$CPO) <- rownames(as.matrix( read.table(object$output$Y,header=T) ))
      colnames(ret$CPO) <- colnames(as.matrix( read.table(object$output$Y,header=T) ))
      
    } 
    
    if( "logP" %in% estimator ){
      ret$logP <- t( as.matrix( read.table(object$output$logP) ) )
    } 
  
  if(length(estimator)>1){
    return(ret)
  }else{
    return(ret[[1]])
  }
}