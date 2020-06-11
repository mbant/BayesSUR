#' @title extract the posterior mean of the parameters
#' @description
#' Extract the posterior mean of the parameters of a "BayesSUR" class object.
#' @name get.estimator
#' @param object an object of class "BayesSUR"
#' @param estimator the name of one estimator. Default is the latent indicator estimator "\code{gamma}". Other options "\code{beta}", "\code{Gy}", "\code{CPO}" and "\code{logP}" 
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
get.estimator <- function(object, estimator="gamma", Pmax=0){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  if( sum(estimator %in% c("gamma","beta","Gy","CPO","logP"))<1 )
    stop("Please specify a correct estimator!")
  
  if( length(estimator) == 1 ){if( estimator == "gamma" ){
    Est <- as.matrix( read.table(object$output$gamma) )
    Est[Est<=Pmax] <- 0
    rownames(Est) <- colnames(read.table(object$output$X,header=T))
    
    ## Create the return object
    ret <-  Est
    class(ret) <- "Manhattan"
    class(plot) <- c(class(plot), "Manhattan")
    return(ret)
  } 
    
    if( estimator == "beta" ){
      Est <- as.matrix( read.table(object$output$beta) )
      Est[as.matrix( read.table(object$output$gamma) )<=Pmax] <- 0
      
      return(Est)
    } 
    
    if( estimator == "Gy" ){
      
      covariancePrior <- object$input$covariancePrior
      if(covariancePrior == "HIW"){
        Est <- as.matrix( read.table(object$output$G) )
      }else{
        stop("Gy is only estimated with hyper-inverse Wishart prior for the covariance matrix of responses!")
      }
      Est[Est<=Pmax] <- 0
      rownames(Est) <- colnames(Est) <- names(read.table(object$output$Y,header=T))
      
      ## Create the return object
      ret <-  Est
      class(ret) <- "response.graph"
      class(plot) <- c(class(plot), "response.graph")
      return(ret)
    } 
    
    if( estimator == "CPO" ){
      Est <- as.matrix( read.table(object$output$CPO) )
      
      if(is.null(object$output$CPO))
        stop("Please specify argument output_CPO in BayesSUR()!")
      rownames(Est) <- rownames(as.matrix( read.table(object$output$Y,header=T) ))
      colnames(Est) <- colnames(as.matrix( read.table(object$output$Y,header=T) ))
      
      ## Create the return object
      ret <-  Est
      class(ret) <- "CPO"
      class(plot) <- c(class(plot), "CPO")
      return(ret)
    } 
    
    if( estimator == "logP" ){
      logP <- t( as.matrix( read.table(object$output$logP) ) )
      model_size <- as.matrix( read.table(object$output$model_size) )
      ncol_Y <- ncol(read.table(object$output$gamma))
      nIter <- object$input$nIter
      
      covariancePrior <- object$input$covariancePrior
      if(covariancePrior=="HIW" & is.null(object$output$Gvisit)){
        Gvisit <- as.matrix( read.table(object$output$Gvisit) )
        ret <- list(logP=logP, model_size=model_size, Gvisit=Gvisit, 
                    ncol_Y=ncol_Y, nIter=nIter, covariancePrior=covariancePrior)
      }else{
        ret <- list(logP=logP, model_size=model_size,
                    ncol_Y=ncol_Y, nIter=nIter, covariancePrior=covariancePrior)
      }
      
      ## Create the return class
      class(ret) <- "MCMCdiag"
      class(plot) <- c(class(plot), "MCMCdiag")
      return(ret)
    } 
  }
  
  if( length(estimator)==2 & sum(estimator %in% c("gamma", "Gy"))==2 ){
    covariancePrior <- object$input$covariancePrior
    if(covariancePrior == "HIW"){
      Gy <- as.matrix( read.table(object$output$G) )
    }else{
      stop("Gy is only estimated with hyper-inverse Wishart prior for the covariance matrix of responses!")
    }
    gamma <- as.matrix( read.table(object$output$gamma) )
    colnames(gamma) <- names(read.table(object$output$Y,header=T))
    rownames(gamma) <- names(read.table(object$output$X,header=T))
    
    ## Create the return object
    ret <- list(gamma=gamma, Gy=Gy, covariancePrior=covariancePrior)
    class(ret) <- "network"
    class(plot) <- c(class(plot), "network")
    return(ret)
  } 
  
  if( (length(estimator)==3 & sum(estimator %in% c("beta", "gamma", "Gy"))==3)
      | (length(estimator)==2 & sum(estimator %in% c("beta", "gamma"))==2)){
    beta <- as.matrix( read.table(object$output$beta) )
    gamma <- as.matrix( read.table(object$output$gamma) )
    #colnames(beta) <- colnames(gamma) <- names(read.table(object$output$Y,header=T))
    #rownames(gamma) <- names(read.table(object$output$X,header=T))
    nonpen <- nrow(beta) - nrow(gamma)
    if(nonpen > 0){
      rownames(beta) <- c(names(read.table(object$output$X0,header=T)), names(read.table(object$output$X,header=T)))
    }else{
      #rownames(beta) <- names(read.table(object$output$X,header=T))
    }
    covariancePrior <- object$input$covariancePrior
    if( (covariancePrior == "HIW") & ("Gy" %in% estimator) ){
      Gy <- as.matrix( read.table(object$output$G) )
      colnames(Gy) <- rownames(Gy) <- names(read.table(object$output$Y,header=T))
      ret <- list(beta=beta, gamma=gamma, Gy=Gy, covariancePrior=covariancePrior)
    }else{
      ret <- list(beta=beta, gamma=gamma, covariancePrior=covariancePrior)
    }
    
    ## Create the return object
    class(ret) <- "estimator"
    class(plot) <- c(class(plot), "estimator")
    return(ret)
  } 
  
}