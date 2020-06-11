#' @title extract the posterior mean of the coefficients of a "BayesSUR" class object
#' @description
#' Extract the posterior mean of the coefficients of a "BayesSUR" class object
#' @name coef.BayesSUR
#' 
#' @param object an object of class "BayesSUR"
#' @param Pmax threshold that truncates the estimated coefficients based on thresholding the estimated latent indicator variable. Default is 0.
#' @param ... other arguments
#' 
#' @return Estimated coefficients are from an object of class "BayesSUR". If the \code{BayesSUR} specified data standardization, the fitted values are base based on standardized data.
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
#' ## check prediction
#' beta.hat <- coef(fit)
#' 
#' @export
coef.BayesSUR <- function(object, Pmax=0, ...){
  
  predict.BayesSUR(object, type="coefficients", Pmax=Pmax, ...)
    
}
