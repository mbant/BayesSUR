#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title plotResponseGraph
#' @description
#' Show the relationship between responses
#' @name plotResponseGraph
#' @param object fitted "runSUR" model
#' @param PmaxResponse cutpoint for thresholding the learning structure matrix of multiple response variables. Default is 0.5
#' @param PtrueResponse true adjacency matrix for the structure of multiple response variables
#' @param response.name A vector for the node names
#' @param edge.weight draw weighted edges after thresholding at 0.5. The defaul value "FALSE" is not to draw weigthed edges
#' @param label.color label color. Default is "black"
#' @param node.size node size. Default is 30
#' @param node.color node color. Default is "dodgerblue
#' @export
plotResponseGraph <- function(object, PmaxResponse=0.5, PtrueResponse=NULL, response.name=NULL, edge.weight=FALSE, label.color="black", node.size=30, node.color="dodgerblue"){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  Gy_hat <- as.matrix( read.table(object$output$G) )
  
  if(!is.null(response.name)){
    rownames(Gy_hat) <- colnames(Gy_hat) <- response.name
  }else{
    rownames(Gy_hat) <- colnames(Gy_hat) <- names(read.table(object$output$Y,header=T))
  }
  
  if(edge.weight){
    Gy_thresh <- Gy_hat
    Gy_thresh[Gy_hat<=PmaxResponse] <- 0
  }else{
    Gy_thresh <- as.matrix( Gy_hat > PmaxResponse )
  }
  
  net <- graph_from_adjacency_matrix(  Gy_thresh, weighted=T, mode="undirected", diag=F)
  V(net)$size = node.size
  V(net)$label.color <- label.color
  V(net)$color <- node.color
  
  if( !is.null(PtrueResponse) ){
    par(mfrow=c(1,2))
    netTRUE <- graph_from_adjacency_matrix(  PtrueResponse, weighted=T, mode="undirected", diag=F)
    V(netTRUE)$label.color <- label.color
    V(netTRUE)$color <- node.color
    V(netTRUE)$size <- node.size
    plot.igraph(netTRUE, main = "True graph of responses", edge.width=E(netTRUE)$weight*ifelse(edge.weight,2,1))
  }
  plot.igraph(net, main = "Estimated graph of responses", edge.width=E(net)$weight*ifelse(edge.weight,2,1))
  par(mfrow=c(1,1))
  
}


