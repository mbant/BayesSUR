#' @title create a selection of plots for a \code{BayesSUR} class object
#' @description
#' Convenience function to create a selection of plots for a \code{BayesSUR} class object. They are plots of estimators, estimated graph Gy, network, Manhattan, MCMC diagnostics and conditional predictive ordinates.
#' @importFrom grDevices dev.hold dev.flush devAskNewPage
#' @name plot.BayesSUR
#' @param x an object of class \code{BayesSUR}
#' @param estimator It is in \code{c(NA, 'beta', 'gamma', 'Gy', 'logP', 'CPO')} and works by combining with argument \code{type}.
#' \itemize{
#'   \item If \code{estimator} is in \code{c("beta", "gamma", "Gy")} and argument \code{type="heatmap"}, it prints heatmaps of the specified estimator in \code{estimator}, and refers to function \code{plotEstimator()} for more other arguments.
#'   \item If \code{estimator="Gy"} and argument \code{type="graph"}, it prints a structure graph of \code{"Gy"}, and refers to function \code{plotGraph()} for more other arguments.
#'   \item If \code{estimator=c("gamma", "Gy")} and argument \code{type="network"}, it prints the estimated network between the response variables and predictors with nonzero coefficients, and refers to function \code{plotMCMCdiag()} for more other arguments.
#'   \item If \code{estimator=NA} (default) and \code{type="NA"} (default), it interactively prints the plots of estimators (i.e., beta, gamma and (or) Gy), response graph Gy, network, Manhattan and MCMC diagnostics.
#' }
#' @param type It is one of \code{"NA"}, \code{"heatmap"}, \code{"graph"}, \code{"network"}, \code{"Manhattan"} and \code{"diagnostics"}, and works by combining with argument \code{estimator}.
#' \itemize{
#'   \item If \code{type="Manhattan"} and argument \code{estimator="gamma"}, it prints Manhattan-like plots for marginal posterior inclusion probabilities (mPIP) and numbers of associated response variables for individual predictors, and refers to function \code{plotManhattan()} for more other arguments.
#'   \item If \code{type="diagnostics"} and argument \code{estimator="logP"} it shows trace plots and diagnostic density plots of a fitted model, and refers to function \code{plotMCMCdiag()} for more other arguments.
#'   \item If \code{type="diagnostics"} and argument \code{estimator="CPO"}, it shows the conditional predictive ordinate (CPO) for each individual of a fitted model, and refers to function \code{plotCPO()} for more other arguments.
#' }
#' @param ... other arguments, see functions \code{plotEstimator()}, \code{plotGraph()}, \code{plotNetwork()}, \code{plotManhattan()}, \code{plotMCMCdiag()} or \code{plotCPO()}
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
#' \donttest{
#' # Show the interactive plots. Note that it needs at least 2000*(nbloc+1) iterations 
#' # for the diagnostic plots where nbloc=3 by default 
#' plot(fit)
#' }
#' 
#' ## plot heatmaps of the estimated beta, gamma and Gy
#' plot(fit, estimator=c("beta", "gamma", "Gy"), type="heatmap")
#' 
#' ## plot estimated graph of responses Gy
#' plot(fit, estimator="Gy", type="graph")
#' 
#' ## plot network between response variables and associated predictors
#' plot(fit, estimator=c("gamma", "Gy"), type="network")
#' 
#' ## print Manhattan-like plots
#' plot(fit, estimator="gamma", type="Manhattan")
#' 
#' ## print MCMC diagnostic plots
#' plot(fit, estimator="logP", type="diagnostics")
#' 
#' @export
plot.BayesSUR <- function(x, estimator = NA, type = NA, ...){
  
  if(!inherits(x, "BayesSUR")) 
    stop("Use only with \"BayesSUR\" objects")
  if(sum(!estimator %in% c(NA, 'beta', 'gamma', 'Gy', 'logP', 'CPO'))) 
    stop("'estimator' must be in c(NA, 'beta', 'gamma', 'Gy', 'logP', 'CPO')!")
  if(!sum(type %in% c(NA, "heatmap", "graph", "Manhattan", "network", "diagnostics")))
    stop("Please specify correct type!")
  
  if(!is.na(type)){
    
    if(!( ((sum(estimator %in% c("beta", "gamma", "Gy")) > 0) & (type=="heatmap")) |
          ((length(estimator)==1) & (estimator[1] == "Gy") & (type=="graph")) |
          ((length(estimator)==2) & (sum(estimator %in% c("gamma", "Gy"))==2) & (type=="network")) |
          ((length(estimator)==1) & (estimator[1] == "gamma") & (type == "Manhattan")) |
          ((length(estimator)==1) & (estimator[1] == "logP") & (type == "diagnostics")) |
          ((length(estimator)==1) & (estimator[1] == "CPO") & (type == "diagnostics")) )
       )
    stop("Please specify correct argument!")
    
    ## refer to function plotEstimator()
    if( (sum(estimator %in% c("beta", "gamma", "Gy")) > 0) & (type=="heatmap") )
      plotEstimator(x, estimator, ...)
    
    ## refer to function plotGraph()
    if( (length(estimator)==1) & (estimator[1] == "Gy") & (type=="graph") )
      plotGraph(x, ...)
    
    ## refer to function plotNetwork()
    if( (length(estimator)==2) & (sum(estimator %in% c("gamma", "Gy"))==2) & (type=="network") )
      plotNetwork(x, ...)
    
    ## refer to function plotManhattan()
    if( (length(estimator)==1) & (estimator[1] == "gamma") & (type == "Manhattan") )
      plotManhattan(x, ...)
    
    ## refer to function plotMCMCdiag()
    if( (length(estimator)==1) & (estimator[1] == "logP")  & (type == "diagnostics") )
      plotMCMCdiag(x, ...)
    
    ## refer to function plotCPO()
    if( (length(estimator)==1) & (estimator[1] == "CPO")  & (type == "diagnostics") )
      plotCPO(x, ...)
    
  }else{
    ## print plots interactively
    if( is.na(estimator[1]) ){
      show <- rep(FALSE, 5)
      show[1:4] <- TRUE
      
      if(x$input$covariancePrior=="IG")
        show[ 2:3 ] <- FALSE
      
      devAskNewPage(TRUE)
      
      if (show[1L]) {
        dev.hold()
        plotEstimator(x, estimator=c("beta", "gamma", "Gy"), header="\nEstimators", ...)
        dev.flush()
      }
      if (show[2L]) {
        dev.hold()
        plotGraph(x, ...)
        dev.flush()
      }
      if (show[3L]) {
        dev.hold()
        plotNetwork(x, header="\n\nNetwork respresentation", ...)
        dev.flush()
      }
      if (show[4L]) {
        dev.hold()
        plotManhattan(x, header="\n\nManhattan-like plots", ...)
        dev.flush()
      }
      if (show[5L]) {
        dev.hold()
        plotMCMCdiag(x, header="\nMCMC diagnostic plots", ...)
        dev.flush()
      }
      
      devAskNewPage(options("device.ask.default")[[1]])
    }else{
      stop("Please specify correct argument!")
    }
  }
  
}