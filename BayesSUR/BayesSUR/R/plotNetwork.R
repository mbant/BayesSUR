#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title plotNetwork
#' @description
#' Network representation of the associations between responses and features
#' @name plotNetwork
#' @param object fitted "runSUR" model
#' @param includeResponse A vector of the response names which are shown in the network. 
#' @param excludeResponse A vector of the response names which are not shown in the network. 
#' @param includePredictor A vector of the predictor names which are shown in the network. 
#' @param excludePredictor A vector of the predictor names which are not shown in the network. 
#' @param MatrixGamma A matrix or dataframe of the latent indicator variable. Default is "NULL" and to extrate it from object of class inheriting from "runSUR"
#' @param PmaxPredictor cutpoint for thresholding the estimated latent indicator variable. Default is 0.5
#' @param PmaxResponse cutpoint for thresholding the learning structure matrix of multiple response variables. Default is 0.5
#' @param nodesizePredictor node size of Predictors in the output graph. Default is 15
#' @param nodesizePredictor node size of response variables in the output graph. Default is 25
#' @param no.isolates remove isolated nodes from responses graph and Full graph, may get problem if there are also isolated Predictors
#' @param lineup A ratio of the heights between responses' area and Predictors'
#' @param edgewith.response the edge width betwen response nodes
#' @param edgewith.predictor the edge width betwen the predictor and response node
#' @export 
plotNetwork <- function(object, includeResponse=NULL, excludeResponse=NULL, includePredictor=NULL, excludePredictor=NULL, 
                        MatrixGamma=NULL, PmaxPredictor=0.5, PmaxResponse=0.5, nodesizePredictor=15, nodesizeResponse=25, no.isolates=FALSE,
                        lineup=0.8, gray.alpha=0.6, edgewith.response=5, edgewith.predictor=2){
  devAskNewPage(FALSE)
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  
  gamma_hat <- as.matrix( read.table(object$output$gamma) )
  colnames(gamma_hat) <- names(read.table(object$output$Y,header=T))
  rownames(gamma_hat) <- names(read.table(object$output$X,header=T))
  
  # select the required resposes and predictors to plot the network
  excludeResponse.idx <- rep(FALSE, ncol(gamma_hat))
  excludePredictor.idx <- rep(FALSE, nrow(gamma_hat))
  if(!is.null(includeResponse)) 
    excludeResponse.idx <- c(!(colnames(gamma_hat) %in% includeResponse))
  if(!is.null(excludeResponse)) 
    excludeResponse.idx <- c(excludeResponse.idx | c(colnames(gamma_hat) %in% excludeResponse))
  if(!is.null(includePredictor)) 
    excludePredictor.idx <- c(!(rownames(gamma_hat) %in% includePredictor))
  if(!is.null(excludePredictor))
    excludePredictor.idx <- c(excludePredictor.idx | c(rownames(gamma_hat) %in% excludePredictor))
  
  gamma_hat <- gamma_hat[!excludePredictor.idx,!excludeResponse.idx]
  gamma_thresh <- as.data.frame(as.matrix(gamma_hat>PmaxPredictor)+0)
  gamma_thresh <- gamma_thresh[rowSums(gamma_thresh)!=0,]
  
  G0_hat <- as.matrix( read.table(object$output$G) )
  G0_hat <- G0_hat[!excludeResponse.idx,!excludeResponse.idx]
  rownames(G0_hat) <- colnames(G0_hat) <-  colnames(gamma_hat)
  G0_thresh <- as.data.frame(as.matrix(G0_hat>PmaxResponse)+0)
  
  plotSEMgraph(G0_thresh, t(gamma_thresh), nodesizeSNP=nodesizePredictor, nodesizeMET=nodesizeResponse, no.isolates=no.isolates, 
               lineup=lineup, gray.alpha=gray.alpha, edgewith.response=edgewith.response, edgewith.predictor=edgewith.predictor)

}
plotSEMgraph <- function(ADJmatrix,GAMmatrix,nodesizeSNP=15,nodesizeMET=25,no.isolates=FALSE,
                         lineup=0.8,gray.alpha=0.6,edgewith.response=5,edgewith.predictor=2){
  
  # ADJmatrix must be a square qxq adjacency matrix (or data frame)
  qq <- dim(ADJmatrix)[1]
  if(dim(ADJmatrix)[2] != qq) stop("adjacency matrix not square")
  
  # GAMmatrix must be a qxp binary matrix (or data frame)
  pp <- dim(GAMmatrix)[2]
  if(dim(GAMmatrix)[1] != qq) stop("Gamma and Adjacency have different no. q")
  
  print(paste("p = ",pp))
  print(paste("q = ",qq))
  print(paste("some snp names",names(GAMmatrix)[1:4]))
  
  # join mets block (adjency) and lower triangle (gamma)
  semgraph <- rbind(ADJmatrix,t(GAMmatrix))
  
  # add zero blocks for lower triangle and snp block
  zeroblock <- matrix(rep(0,pp*(pp+qq)),nrow=qq+pp,ncol=pp)
  zeroblock <- data.frame(zeroblock)
  colnames(zeroblock) <- colnames(GAMmatrix)
  rownames(zeroblock) <- rownames(semgraph)
  
  semgraph <- cbind(semgraph,zeroblock)
  #print(semgraph[201:206,1:6])
  
  # igraph objects
  graphADJ <- graph.adjacency(as.matrix(ADJmatrix),weight=T,diag=FALSE,mode="undirected")
  graphSEM <- graph.adjacency(as.matrix(semgraph),weight=T,diag=FALSE,mode="directed")
  
  # don't plot isolated nodes?
  if(no.isolates){
    graphADJ <- delete.vertices(graphADJ,degree(graphADJ)==0) 
    graphSEM <- delete.vertices(graphSEM,degree(graphSEM)==0) 
    print("******* Removing isolated nodes from Adjacency and Full SEM, may get problem if there are also isolated SNPs.")
  } 
  
  # get co-ords for undirected edges using layout function (scaled)
  lladj <- layout.fruchterman.reingold(graphADJ)
  lmax <- max(lladj[,1])
  lmin <- min(lladj[,1])
  lladj[,1] <- (lladj[,1]-lmin)/(lmax-lmin)
  lmax <- max(lladj[,2])
  lmin <- min(lladj[,2])
  lladj[,2] <- (lladj[,2]-lmin)/(lmax-lmin)
  
  # plot adjacency only
  #plot(graphADJ,vertex.size=15,edge.width=2,edge.color="black",layout=lladj)
  
  # line up snps
  lymax <- max(lladj[,2]) + (max(lladj[,2])-min(lladj[,2]))*(lineup-1)/2
  lymin <- min(lladj[,2]) + (max(lladj[,2])-min(lladj[,2]))*(1-lineup)/2
  llsnps <- matrix( c(rep(-0.5,pp),lymin+(1:pp)*1.0*(lymax-lymin)/pp), nrow=pp, ncol=2)
  
  llsem <- rbind(lladj,llsnps)
  
  ### plot SEM  
  
  # plot snps and mets nodes differently
  # set node sizes directly in graph object
  V(graphSEM)$size <- c( rep(nodesizeMET,qq), rep(nodesizeSNP,pp) )
  
  n.edgeADJ <- gsize(graphADJ)
  n.edgeGAM <- gsize(graphSEM) - n.edgeADJ
   
  plot(graphSEM,edge.arrow.size=0.5, edge.width=c(rep(edgewith.response,2*n.edgeADJ),rep(edgewith.predictor,2*n.edgeGAM)),
       edge.color=c(rep(gray(0),2*n.edgeADJ),rep(gray(0.7, alpha=gray.alpha),2*n.edgeGAM)),layout=llsem)
  
}
