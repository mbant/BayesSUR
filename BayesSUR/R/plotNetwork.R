#' @title plot the network representation of the associations between responses and predictors
#' @description
#' Plot the network representation of the associations between responses and predictors, based on the estimated gamma matrix and graph of responses from a "BayesSUR" class object.
#' @importFrom graphics text 
#' @importFrom grDevices gray 
#' @importFrom igraph V E gsize layout_in_circle plot.igraph degree layout.fruchterman.reingold delete.vertices graph.adjacency
#' @name plotNetwork
#' @param object an object of class "BayesSUR"
#' @param includeResponse A vector of the response names which are shown in the network
#' @param excludeResponse A vector of the response names which are not shown in the network
#' @param includePredictor A vector of the predictor names which are shown in the network
#' @param excludePredictor A vector of the predictor names which are not shown in the network
#' @param MatrixGamma A matrix or dataframe of the latent indicator variable. Default is \code{NULL} and to extrate it from object of class inheriting from an object of class "BayesSUR"
#' @param PmaxPredictor cutpoint for thresholding the estimated latent indicator variable. Default is 0.5
#' @param PmaxResponse cutpoint for thresholding the learning structure matrix of multiple response variables. Default is 0.5
#' @param nodesizePredictor node size of Predictors in the output graph. Default is 15
#' @param nodesizeResponse node size of response variables in the output graph. Default is 25
#' @param no.isolates remove isolated nodes from responses graph and Full graph, may get problem if there are also isolated Predictors
#' @param lineup A ratio of the heights between responses' area and Predictors'
#' @param gray.alpha the opacity. The default is 0.6
#' @param edgewith.response the edge width betwen response nodes
#' @param edgewith.predictor the edge width betwen the predictor and response node
#' @param edge.weight draw weighted edges after thresholding at 0.5. The defaul value \code{FALSE} is not to draw weigthed edges
#' @param label.predictor A vector of the names of predictors
#' @param label.response A vector of the names of response variables
#' @param color.predictor color of the predictor nodes
#' @param color.response color of the reponse nodes
#' @param name.predictors a subtitle for the predictors
#' @param name.responses a subtitle for the responses
#' @param vertex.frame.color The color of the frame of the vertices. If you don't want vertices to have a frame, supply NA as the color name
#' @param layoutInCircle place vertices on a circle, in the order of their vertex ids. The default is \code{FALSE}
#' @param header the main title
#' @param ... other arguments
#' 
#' @examples
#' data("example_eQTL", package = "BayesSUR")
#' hyperpar <- list( a_w = 2 , b_w = 5 )
#' 
#' fit <- BayesSUR(Y = example_eQTL[["blockList"]][[1]], 
#'                 X = example_eQTL[["blockList"]][[2]],
#'                 data = example_eQTL[["data"]], outFilePath = tempdir(),
#'                 nIter = 100, burnin = 50, nChains = 2, gammaPrior = "hotspot",
#'                 hyperpar = hyperpar, tmpFolder = "tmp/" )
#' 
#' ## check output
#' # show the Network representation of the associations between responses and features
#' plotNetwork(fit)
#' 
#' @export 
plotNetwork <- function(object, includeResponse=NULL, excludeResponse=NULL, includePredictor=NULL, excludePredictor=NULL, 
                        MatrixGamma=NULL, PmaxPredictor=0.5, PmaxResponse=0.5, nodesizePredictor=5, nodesizeResponse=25, no.isolates=FALSE,
                        lineup=1, gray.alpha=0.6, edgewith.response=5, edgewith.predictor=2, edge.weight=FALSE, label.predictor=NULL,
                        label.response=NULL, color.predictor=NULL,color.response=NULL, name.predictors=NULL,name.responses=NULL, 
                        vertex.frame.color=NA,layoutInCircle=FALSE, header="", ...){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  
  gamma_hat <- as.matrix( read.table(object$output$gamma) )
  colnames(gamma_hat) <- names(read.table(object$output$Y,header=T))
  rownames(gamma_hat) <- names(read.table(object$output$X,header=T))
  
  if(sum(colnames(gamma_hat)==paste("V",1:ncol(gamma_hat),sep="")) == ncol(gamma_hat))
    colnames(gamma_hat) <- paste("Y",1:ncol(gamma_hat),sep="")
  if(sum(rownames(gamma_hat)==paste("V",1:nrow(gamma_hat),sep="")) == nrow(gamma_hat))
    rownames(gamma_hat) <- paste("X",1:nrow(gamma_hat),sep="")
  
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
  
  Gy_hat <- as.matrix( read.table(object$output$G) )
  Gy_hat <- Gy_hat[!excludeResponse.idx,!excludeResponse.idx]
  
  if(edge.weight){
    Gy_thresh <- Gy_hat
    Gy_thresh[Gy_hat<=PmaxResponse] <- 0
    
    gamma_thresh <- gamma_hat
    gamma_thresh[gamma_hat<=PmaxPredictor] <- 0
  }else{
    Gy_thresh <- as.matrix( Gy_hat > PmaxResponse )
    gamma_thresh <- as.matrix(gamma_hat>PmaxPredictor)
  }
  
  gamma_thresh <- matrix(gamma_thresh[rowSums(gamma_thresh)!=0,], ncol=ncol(gamma_hat))
  colnames(gamma_thresh) <- colnames(gamma_hat)
  rownames(gamma_thresh) <- rownames(gamma_hat)[rowSums(gamma_hat>PmaxPredictor)!=0]
  rownames(Gy_thresh) <- colnames(Gy_thresh) <-  colnames(gamma_hat) 
  
  plotSEMgraph(Gy_thresh, t(gamma_thresh), nodesizeSNP=nodesizePredictor, nodesizeMET=nodesizeResponse, no.isolates=no.isolates, 
               lineup=lineup, gray.alpha=gray.alpha, edgewith.response=edgewith.response, edgewith.predictor=edgewith.predictor,edge.weight=edge.weight,
               label.predictor=label.predictor,label.response=label.response, color.predictor=color.predictor,color.response=color.response, 
               name.predictors=name.predictors,name.responses=name.responses, vertex.frame.color=vertex.frame.color,layoutInCircle=layoutInCircle,...)
  title(paste("\n\n",header,sep=""), outer=T)
  
}
plotSEMgraph <- function(ADJmatrix,GAMmatrix,nodesizeSNP=2,nodesizeMET=25,no.isolates=FALSE,
                         lineup=1,gray.alpha=0.6,edgewith.response=5,edgewith.predictor=2,
                         label.predictor=NULL,label.response=NULL, color.predictor=NULL,color.response=NULL, 
                         name.predictors=NULL,name.responses=NULL,edge.weight=FALSE, vertex.frame.color=NA,layoutInCircle=FALSE,...){
  
  # ADJmatrix must be a square qxq adjacency matrix (or data frame)
  qq <- dim(ADJmatrix)[1]
  if(dim(ADJmatrix)[2] != qq) stop("adjacency matrix not square")
  
  # GAMmatrix must be a qxp binary matrix (or data frame)
  pp <- dim(GAMmatrix)[2]
  if(dim(GAMmatrix)[1] != qq) stop("Gamma and Adjacency have different no. q")
  
  # print(paste("p = ",pp))
  # print(paste("q = ",qq))
  # print(paste("some snp names",names(GAMmatrix)[1:4]))
  
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
  graphADJ <- graph.adjacency(as.matrix(ADJmatrix),weighted=TRUE,diag=FALSE,mode="undirected")
  graphSEM <- graph.adjacency(as.matrix(semgraph),weighted=TRUE,diag=FALSE,mode="directed")
  
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
  
  V(graphSEM)$label.color <- "black"
  
  V(graphSEM)$color <- c(rep("dodgerblue", nrow(GAMmatrix)), rep("red", ncol(GAMmatrix)))
  if(!is.null(color.predictor)) V(graphSEM)$color[-c(1:nrow(GAMmatrix))] <- color.predictor
  if(!is.null(color.response)) V(graphSEM)$color[1:nrow(GAMmatrix)] <- color.response
  
  V(graphSEM)$label <- c(rownames(GAMmatrix), colnames(GAMmatrix))
  if(!is.null(label.predictor)) V(graphSEM)$label[-c(1:nrow(GAMmatrix))] <- label.predictor
  if(!is.null(label.response)) V(graphSEM)$label[1:nrow(GAMmatrix)] <- label.response
  
  if(edge.weight){
    #edge.width=edge.betweenness(graphSEM)
    edge.width=E(graphSEM)$weight*ifelse(edge.weight,5,1)
  }else{
    edge.width=c(rep(edgewith.response,2*n.edgeADJ),rep(edgewith.predictor,2*n.edgeGAM))
  }
  
  if(!layoutInCircle){
    layoutSEM <- llsem
  }else{
    layoutSEM <- layout_in_circle(graphSEM)
  }
  
  plot.igraph(graphSEM,edge.arrow.size=0.5, edge.width=edge.width, vertex.frame.color=vertex.frame.color,
       edge.color=c(rep(gray(0),2*n.edgeADJ),rep(gray(0.7, alpha=gray.alpha),2*n.edgeGAM)),layout=layoutSEM,...)
  
  if(!is.null(name.predictors)) text(-1,-1.3,name.predictors,cex=1.2)
  if(!is.null(name.responses)) text(0.4,-1.3,name.responses,cex=1.2)
}
