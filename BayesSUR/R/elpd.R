#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title prediction accuracy
#' @description
#' Measure the prediction accuracy by the elpd (expected log pointwise predictive density). One method to estimate elpd is by the Bayesian LOO estimate of out-of-sample predictive fit. The other method the WAIC (widely applicable information criterion).
#' @name elpd
#' @param object the object from the runSUR
#' @param method the name of the prediction accuracy index. Default is the "\code{LOO}". The other index is the "\code{WAIC}".
#' @references Vehtari, A., Gelman, A., Gabry, J. (2017). \emph{Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC.} Statistics and Computing, 27(5): 1413–1432.
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
#' # print the prediction accuracy elpd (expected log pointwise predictive density) 
#' # by the Bayesian LOO estimate of out-of-sample predictive fit
#' elpd(fit, method="LOO")
#' }
#' 
#' @export
elpd <- function(object, method="LOO"){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  
  if(toupper(method) == "LOO"){
    elpd <- sum(log(read.table(object$output$CPO)))
    names(elpd) <- "elpd.loo"
  }else if(toupper(method) == "WAIC"){
    elpd <- sum(read.table(object$output$WAIC))
    names(elpd) <- "elpd.waic"
  }else{
    stop("Please give the correct method name!")
  }
  
  return(elpd)
}