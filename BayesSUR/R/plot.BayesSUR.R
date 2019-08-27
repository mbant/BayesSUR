#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title plotEstimator
#' @description
#' Six available plots from the object of fitted Bayesian Seemingly Unrelated Regression
#' @importFrom grDevices dev.hold dev.flush
#' @name plot
#' @param x object of the fitted "runSUR" model.
#' @param which if a subset of the plots is required, specify a subset of the numbers 1:5 which are plots of estimators, response graph, network, manhattan and MCMC diagnosis, respectively.
#' @param ... Other parameters in the function \code{plot.default.R} file
#' @export
plot.BayesSUR <- function(x, which = c(1L:5L), ...){
  
  if (!inherits(x, "BayesSUR")) 
    stop("use only with \"BayesSUR\" objects")
  if (!is.numeric(which) || any(which < 1) || any(which > 5)) 
    stop("'which' must be in 1:5")
  
  show <- rep(FALSE, 5)
  show[which] <- TRUE
  
  devAskNewPage(TRUE)
  
  if (show[1L]) {
    dev.hold()
    plotEstimator(x, ...)
    dev.flush()
  }
  if (show[2L]) {
    dev.hold()
    plotResponseGraph(x, ...)
    dev.flush()
  }
  if (show[3L]) {
    dev.hold()
    plotNetwork(x, ...)
    dev.flush()
  }
  if (show[4L]) {
    dev.hold()
    plotManhattan(x, ...)
    dev.flush()
  }
  if (show[5L]) {
    dev.hold()
    plotMCMCdiag(x, ...)
    dev.flush()
  }
  
  devAskNewPage(options("device.ask.default")[[1]])

}