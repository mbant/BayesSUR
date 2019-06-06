#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title plotResponseGraph
#' @description
#' Show the relationship between responses
#' @name plotResponseGraph
#' @param object fitted "runSUR" model
#' @param PmaxResponse cutpoint for thresholding the learning structure matrix of multiple response variables. Default is 0.5
#' @param PtrueResponse true adjacency matrix for the structure of multiple response variables
#' @export
plotResponseGraph <- function(object, PmaxResponse=0.5, PtrueResponse=NULL, response.name=NULL){
  #library(igraph)
  G0_hat <- as.matrix( read.table(object$output$G) )
  
  if(!is.null(response.name)){
    rownames(G0_hat) <- colnames(G0_hat) <- response.name
  }else{
    rownames(G0_hat) <- colnames(G0_hat) <- names(read.table(object$output$Y,header=T))
  }
   
  G0_thresh <- as.matrix( G0_hat > PmaxResponse )
  net <- graph_from_adjacency_matrix(  G0_thresh, weighted=T, mode="undirected", diag=F)
  V(net)$size = 30
  
  if( !is.null(PtrueResponse) ){
    par(mfrow=c(1,2))
    netTRUE <- graph_from_adjacency_matrix(  PtrueResponse, weighted=T, mode="undirected", diag=F)
    V(netTRUE)$size <- 30
    plot.igraph(netTRUE, main = "True graph of responses")
  }
  plot.igraph(net, main = "Estimated graph of responses")
  par(mfrow=c(1,1))
  
}


