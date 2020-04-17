#' @title create a selection of plots for a "BayesSUR" class object
#' @description
#' Convenience function to create a selection of plots for a "BayesSUR" class object. They are plots of estimators, response graph, network, manhattan and MCMC diagnosis indexed by numbers 1:5.
#' @importFrom grDevices dev.hold dev.flush devAskNewPage
#' @name plot.BayesSUR
#' @param x an object of class "BayesSUR".
#' @param which if a subset of the plots is required, specify a subset of the numbers 1:5 which are plots of estimators, response graph, network, manhattan and MCMC diagnosis, respectively. Default is \code{c(1L:4L)} Only \code{c(1,4,5)} is valid for the HRR models.
#' @param ... other arguments
#' 
#' @examples
#' data("example_eQTL", package = "BayesSUR")
#' hyperpar = list( a_w = 2 , b_w = 5 )
#' 
#' set.seed(9173)
#' fit <- BayesSUR(Y = example_eQTL[["blockList"]][[1]], 
#'                 X = example_eQTL[["blockList"]][[2]],
#'                 data = example_eQTL[["data"]], outFilePath = tempdir(),
#'                 nIter = 100, burnin = 0, nChains = 2, gammaPrior = "hotspot",
#'                 hyperpar = hyperpar, tmpFolder = "tmp/" )
#' 
#' ## check output
#' # Show the interactive plots. Note that it needs at least 2000*(nbloc+1) iterations 
#' # for the diagnosis plots where nbloc=3 by default 
#' \donttest{
#' plot.BayesSUR(fit)
#' }
#' 
#' @export
plot.BayesSUR <- function(x, which = c(1L:4L), ...){
  
  if (!inherits(x, "BayesSUR")) 
    stop("use only with \"BayesSUR\" objects")
  if (!is.numeric(which) || any(which < 1) || any(which > 5)) 
    stop("'which' must be in 1:5")
  
  show <- rep(FALSE, 5)
  show[which] <- TRUE
  
  if(x$input$covariancePrior=="IG")
    show[ 2:3 ] <- FALSE
  
  devAskNewPage(TRUE)
  
  if (show[1L]) {
    dev.hold()
    estimators <- getEstimator(x, estimator = c("beta","gamma","Gy"))
    plot.Estimator(estimators, header="\nEstimators", ...)
    dev.flush()
  }
  if (show[2L]) {
    dev.hold()
    Gy <- getEstimator(x, estimator = "Gy")
    plot.ResponseGraph(Gy, ...)
    dev.flush()
  }
  if (show[3L]) {
    dev.hold()
    network <- getEstimator(x, estimator = c("gamma","Gy"))
    plot.Network(network, header="\n\nNetwork respresentation", ...)
    dev.flush()
  }
  if (show[4L]) {
    dev.hold()
    gamma <- getEstimator(x, estimator = "gamma")
    plot.Manhattan(gamma, header="\n\nManhattan-like plots", ...)
    dev.flush()
  }
  if (show[5L]) {
    dev.hold()
    MCMCdiag <- getEstimator(x, estimator = "logP")
    plot.MCMCdiag(MCMCdiag, header="\nMCMC diagnostic plots", ...)
    dev.flush()
  }
  
  devAskNewPage(options("device.ask.default")[[1]])
}