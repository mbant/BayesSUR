#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title extract estimated coefficients from a "runSUR" object.
#' @description
#' Extract model estimated coefficients from the object of fitted Bayesian Seemingly Unrelated Regression
#' @name coef
#' 
#' @param object fitted \code{runSUR} model
#' @param Pmax truncate the estimated coefficients based on thresholding the estimated latent indicator variable at 0 by default
#' @param ... other arguments
#' 
#' @return Estimated coefficients are from the object \code{object}. If the \code{runSUR} specified data standardization, the fitted values are base based on standardized data.
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
#' beta.hat <- coef(fit)
#' }
#' 
#' @export
coef.BayesSUR <- function(object, Pmax=0, ...){
  
  predict.BayesSUR(object, type="coefficients", Pmax=Pmax, ...)
    
}
