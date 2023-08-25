#' @title get fitted responses
#' @description
#' Return the fitted response values that correspond to the posterior mean 
#' estimates from a \code{BayesSUR} class object.
#' @name fitted.BayesSUR
#' 
#' @param object an object of class \code{BayesSUR}
#' @param beta.type type of estimated beta for the fitted model. Default is 
#' \code{marginal}, giving marginal beta estimation. If 
#' \code{beta.type="conditional"}, it gives beta estimation conditional 
#' on gamma=1
#' @param Pmax valid if \code{beta.type="conditional"}. If 
#' \code{beta.type="conditional"} and \code{Pmax=0.5}, it gives median 
#' probability model betas. Default is 0
#' @param ... other arguments
#'
#' @return Fitted values extracted from an object of class \code{BayesSUR}. If 
#' the \code{BayesSUR} specified data standardization, the fitted values are 
#' base based on standardized data.
#'
#' @examples
#' data("exampleEQTL", package = "BayesSUR")
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
#' ## check fitted values
#' fitted.val <- fitted(fit)
#'
#' @export
fitted.BayesSUR <- function(object, Pmax = 0, beta.type = "marginal", ...) {
  if (Pmax < 0 || Pmax > 1) {
    stop("Please specify correct argument 'Pmax' in [0,1]!")
  }
  if ((Pmax > 0) && (beta.type == "marginal")) {
    stop("Pmax > 0 is valid only if the argument beta.type = 'conditional'!")
  }

  beta_hat <- getEstimator(object, estimator = "beta", Pmax = Pmax, 
                           beta.type = beta.type, ...)

  object$output[-1] <- 
    paste(object$output$outFilePath, object$output[-1], sep = "")
  X <- as.matrix(read.table(object$output$X, header = TRUE))

  if ("X0" %in% names(object$output)) {
    X0 <- as.matrix(read.table(object$output$X0))
  } else {
    X0 <- NULL
  }

  y.pred <- cbind(X0, X) %*% beta_hat

  return(y.pred)
}
