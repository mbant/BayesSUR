#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title summary
#' @description
#' Summarize the results from \code{runSUR} fitted model
#' @importFrom Matrix Matrix
#' @name summary.BayesSUR
#' @param object fitted \code{runSUR} model
#' @param PmaxPredictor cutpoint for thresholding the marginal posterior probability of inclusion (mPPI). Default is 0.5
#' @param ... Other parameters in the function \code{summary.default.R} file
#' 
#' @examples
#' \donttest{
#' data(example_eQTL, package = "BayesSUR")
#' hyperpar = list( a_w = 2 , b_w = 5 )
#' 
#' fit = runSUR(example_eQTL[["data"]], outFilePath = "results/",
#'                      Y = example_eQTL[["blockList"]][[1]],
#'                      X = example_eQTL[["blockList"]][[2]],
#'                      nIter = 1000, nChains = 2, gammaPrior = "hotspot",
#'                      hyperpar = hyperpar, tmpFolder="tmp/" )
#' 
#' ## check output
#' # show the summary information
#' summary(fit)
#' }
#' 
#' @export
summary.BayesSUR <- function(object, PmaxPredictor=0.5, ...){
  
  ans <- list(status=object$status)
  ans$elpd <- c(elpd(object, method="LOO"), elpd(object, method="WAIC"))
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  
  ans$CPO <- summary.default(as.vector(as.matrix(read.table(object$output$CPO))))[-4]
  
  gamma <- as.matrix( read.table(object$output$gamma) )
  ans$df <- sum(gamma > PmaxPredictor)
  
  # extract top 10 genomic variables based on average mPPI across all responses
  mean.predictors <- rowMeans(gamma)
  top10.predictors <- mean.predictors[sort.list(mean.predictors, decreasing=TRUE)[1:min(10,nrow(gamma))]]
  names(top10.predictors) <- names(read.table(object$output$X,header=T))[sort.list(mean.predictors, decreasing=TRUE)[1:min(10,nrow(gamma))]]
  
  # extract top 10 response variables based on average mPPI across all responses
  mean.responses <- colMeans(gamma)
  top10.responses <- mean.responses[sort.list(mean.responses, decreasing=TRUE)[1:min(10,ncol(gamma))]]
  names(top10.responses) <- names(read.table(object$output$Y,header=T))[sort.list(mean.responses, decreasing=TRUE)[1:min(10,ncol(gamma))]]
  
  ans$chainParameters <- object$input[1:3]
  ans$modelParameters <- object$input[4:9]
  ans$hyperParameters <- object$input$hyperparameters
  ans$outputFiles <- object$output
  ans$outputFiles["outFilePath"] <- NULL
 
  #cat("\nCall:\n", deparse(object$call)[1], "\n", sep="")
  cat("\nCall:\n  ", paste(unlist(strsplit(deparse(object$call), ","))[1:3],c(",",",",", ...)"),sep="",collapse=""), "\n", sep="")
  cat("\nCPOs:\n"); print(ans$CPO)
  cat("\nNumber of identified predictors (mPPI > ", PmaxPredictor, "): ", sum(gamma > PmaxPredictor), " of ", ncol(gamma), "x", nrow(gamma), "\n", sep="")
  cat("\nTop", min(10,nrow(gamma)), "predictors on average mPPI across all responses:\n")
  print(top10.predictors)
  
  cat("\nTop", min(10,ncol(gamma)), "responses on average mPPI across all predictors:\n")
  print(top10.responses)
  
  cat("\nExpected log pointwise predictive density (elpd):\n", "  elpd.LOO = ", ans$elpd[1], ",  elpd.WAIC = ", ans$elpd[2], "\n", sep="")
  cat("\nMCMC specification:\n", "  iterations = ",  ans$chainParameters$nIter, ",  burn-in = ",  ans$chainParameters$burnin, ",  chains = ",  ans$chainParameters$nChains, "\n", sep="")
  cat("\nModel specification:\n", "  covariance prior: ",  ans$modelParameters$covariancePrior, "\n  gamma prior: ",  ans$modelParameters$gammaPrior, "\n  gamma sampler: ",  
      ans$modelParameters$gammaSampler, "\n  gamma initialisation: ",  ans$modelParameters$gammaInit, "\n", sep="")
  
  if(is.null(ans$hyperParameters)){
    cat("\nHyper-parameters:\n")
    print(unlist(object$input$hyperParameters))
  }
  cat("\n")
  #cat("\noutputFiles:\n", paste(ans$outputFiles))
  
  invisible(ans)
}