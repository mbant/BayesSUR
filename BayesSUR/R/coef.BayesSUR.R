#' @title coef method for class \code{BayesSUR}
#' @description
#' Extract the posterior mean of the coefficients of a \code{BayesSUR} class object
#' @name coef.BayesSUR
#'
#' @param object an object of class \code{BayesSUR}
#' @param beta.type type of output beta. Default is \code{marginal}, giving 
#' marginal beta estimation. If \code{beta.type="conditional"}, it gives beta 
#' estimation conditional on gamma=1.
#' @param Pmax If \code{Pmax=0.5} and \code{beta.type="conditional"}, it gives 
#' median probability model betas. Default is 0.
#' @param ... other arguments
#'
#' @return Estimated coefficients are from an object of class \code{BayesSUR}. 
#' If the \code{BayesSUR} specified data standardization, the fitted values 
#' are base based on standardized data.
#'
#' @examples
#' data("exampleQTL", package = "BayesSUR")
#' hyperpar <- list(a_w = 2, b_w = 5)
#'
#' set.seed(9173)
#' fit <- BayesSUR(
#'   Y = exampleEQTL[["blockList"]][[1]],
#'   X = exampleEQTL[["blockList"]][[2]],
#'   data = exampleEQTL[["data"]], outFilePath = tempdir(),
#'   nIter = 10, burnin = 0, nChains = 1, gammaPrior = "hotspot",
#'   hyperpar = hyperpar, tmpFolder = "tmp/"
#' )
#'
#' ## check prediction
#' beta.hat <- coef(fit)
#'
#' @export
coef.BayesSUR <- function(object, beta.type = "marginal", Pmax = 0, ...) {
  if (!(beta.type %in% c("marginal", "conditional"))) {
    stop("Please specify correct 'beta.type'!")
  }
  if (Pmax < 0 || Pmax > 1) {
    stop("Please specify a correct argument 'Pmax' in [0,1]!")
  }
  if ((Pmax > 0) && (beta.type == "marginal")) {
    stop("Pmax > 0 is valid only if the argument beta.type = 'conditional'!")
  }
  getEstimator(object, estimator = "beta", Pmax = Pmax, 
               beta.type = beta.type, ...)
}
