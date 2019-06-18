#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title plotEstimator
#' @description
#' Six available plots from the object of fitted Bayesian Seemingly Unrelated Regression
#' @name plot
#' @param object fitted "runSUR" model.
#' @param which if a subset of the plots is required, specify a subset of the numbers 1:5 which are plots of estimators, response graph, network, manhattan and MCMC diagnosis, respectively.
#' @export
plot.BayesSUR <- function(object, which = c(1L:5L), ask = prod(par("mfcol")) < length(which) && dev.interactive()){
  
  if (!inherits(object, "BayesSUR")) 
    stop("use only with \"BayesSUR\" objects")
  if (!is.numeric(which) || any(which < 1) || any(which > 5)) 
    stop("'which' must be in 1:5")
  
  show <- rep(FALSE, 5)
  show[which] <- TRUE
  
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  
  if (show[1L]) {
    dev.hold()
    plotEstimator(object)
    dev.flush()
  }
  if (show[2L]) {
    dev.hold()
    plotResponseGraph(object)
    dev.flush()
  }
  if (show[3L]) {
    dev.hold()
    plotNetwork(object)
    dev.flush()
  }
  if (show[4L]) {
    dev.hold()
    plotManhattan(object)
    dev.flush()
  }
  if (show[5L]) {
    dev.hold()
    plotMCMCdiag(object)
    dev.flush()
  }
  
  invisible()

}