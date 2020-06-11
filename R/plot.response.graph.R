#' @title plot the estimated graph for multiple response variables
#' @description
#' Plot the estimated graph for multiple response variables from a "BayesSUR" class object.
#' @importFrom igraph V E plot.igraph graph_from_adjacency_matrix V<-
#' @importFrom graphics par
#' @name plot.response.graph
#' @param x an object of class \code{get.estimator} with \code{estimator="Gy"}
#' @param PmaxResponse cutpoint for thresholding the learning structure matrix of multiple response variables. Default is 0.5
#' @param PtrueResponse true adjacency matrix for the structure of multiple response variables
#' @param name.responses A vector for the node names. The default is "NA" only to show the locations. Value "auto" show the response names from the orginal data. 
#' @param edge.weight draw weighted edges after thresholding at 0.5. The defaul value "FALSE" is not to draw weigthed edges
#' @param label.color label color. Default is "black"
#' @param node.size node size. Default is 30
#' @param node.color node color. Default is "dodgerblue
#' @param ... other arguments
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
#' ## check output
#' # show the graph relationship between responses
#' Gy <- get.estimator(fit, estimator = "Gy")
#' plot(Gy)
#' 
#' @export
plot.response.graph <- function(x, PmaxResponse=0.5, PtrueResponse=NULL, name.responses=NA, edge.weight=FALSE, label.color="black", node.size=30, node.color="dodgerblue", ...){
  
  Gy_hat <- x
  if(!is.na(name.responses)){
    rownames(Gy_hat) <- colnames(Gy_hat) <- name.responses
  }
  
  if(edge.weight){
    Gy_thresh <- Gy_hat
    Gy_thresh[Gy_hat<=PmaxResponse] <- 0
  }else{
    Gy_thresh <- as.matrix( Gy_hat > PmaxResponse )
  }
  
  net <- graph_from_adjacency_matrix(  Gy_thresh, weighted=T, mode="undirected", diag=F)
  V(net)$size <- node.size
  V(net)$label.color <- label.color
  V(net)$color <- node.color
  
  if( !is.null(PtrueResponse) ){
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))    
    par(mfrow=c(1,2))
    netTRUE <- graph_from_adjacency_matrix(  PtrueResponse, weighted=T, mode="undirected", diag=F)
    V(netTRUE)$label.color <- label.color
    V(netTRUE)$color <- node.color
    V(netTRUE)$size <- node.size
    plot.igraph(netTRUE, main = "True graph of responses", edge.width=E(netTRUE)$weight*ifelse(edge.weight,2,1), vertex.frame.color=NA,cex.main=1.5, ...)
  }
  plot.igraph(net, main = "Estimated graph of responses", edge.width=E(net)$weight*ifelse(edge.weight,2,1), vertex.frame.color=NA,cex.main=1.5, ...)
  
}


