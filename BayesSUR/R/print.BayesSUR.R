#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title print a BayesSUR object
#' @description
#' Print a summary of the \code{runSUR} object.
#' @name print
#' @param x object of the fitted \code{runSUR} object.
#' @param PmaxPredictor cutpoint for thresholding the marginal posterior probability of inclusion (mPPI). Default is 0.5
#' @param ... further arguments passed to or from other methods.
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
#' # show the print information
#' print(fit)
#' }
#' 
#' @export
print.BayesSUR <- function(x, PmaxPredictor=0.5, ...){
  
  gamma <- as.matrix( read.table(paste(x$output$outFilePath,x$output$gamma,sep="")) )
  
  #cat("\nCall:\n", deparse(x$call)[1], "\n", sep="")
  cat("\nCall:\n ", paste(unlist(strsplit(deparse(x$call), ","))[1:3],c(",",",",", ...)"),sep="",collapse=""), "\n", sep="")
  cat("\nNumber of identified predictors (mPPI > ", PmaxPredictor, "): ", sum(gamma > PmaxPredictor), " of ", ncol(gamma), "x", nrow(gamma), "\n", sep="")
  cat("\nExpected log pointwise predictive density (elpd):\n", " elpd.LOO = ", elpd(x, method="LOO"), ",  elpd.WAIC = ", elpd(x, method="WAIC"), "\n\n", sep="")

}