#' @title plot network representation of the associations between responses and predictors
#' @description
#' Plot the network representation of the associations between responses and 
#' predictors, based on the estimated gamma matrix and graph of responses 
#' from a "BayesSUR" class object.
#' 
#' @importFrom graphics text
#' @importFrom grDevices gray
#' @importFrom igraph V E gsize layout_in_circle plot.igraph degree 
#' @importFrom igraph layout.fruchterman.reingold delete.vertices 
#' @importFrom igraph graph.adjacency delete.edges ecount V<-
#' @name plotNetwork
#' @param x an object of class \code{BayesSUR}
#' @param includeResponse A vector of the response names which are shown in the network
#' @param excludeResponse A vector of the response names which are not shown in the network
#' @param includePredictor A vector of the predictor names which are shown in the network
#' @param excludePredictor A vector of the predictor names which are not shown in the network
#' @param MatrixGamma A matrix or dataframe of the latent indicator variable. 
#' Default is \code{NULL} and to extrate it from object of class inheriting 
#' from an object of class \code{BayesSUR}
#' @param PmaxPredictor cutpoint for thresholding the estimated latent 
#' indicator variable. Default is 0.5
#' @param PmaxResponse cutpoint for thresholding the learning structure matrix 
#' of multiple response variables. Default is 0.5
#' @param nodesizePredictor node size of Predictors in the output graph. 
#' Default is 15
#' @param nodesizeResponse node size of response variables in the output graph. 
#' Default is 25
#' @param no.isolates remove isolated nodes from responses graph and full 
#' graph, may get problem if there are also isolated Predictors
#' @param lineup A ratio of the heights between responses' area and predictors'
#' @param gray.alpha the opacity. The default is 0.6
#' @param edgewith.response the edge width between response nodes
#' @param edgewith.predictor the edge width between the predictor and response node
#' @param edge.weight draw weighted edges after thresholding at 0.5. The 
#' default value \code{FALSE} is not to draw weighted edges
#' @param label.predictor A vector of the names of predictors
#' @param label.response A vector of the names of response variables
#' @param color.predictor color of the predictor nodes
#' @param color.response color of the response nodes
#' @param name.predictors A subtitle for the predictors
#' @param name.responses A subtitle for the responses
#' @param vertex.frame.color color of the frame of the vertices. If you don't 
#' want vertices to have a frame, supply NA as the color name
#' @param layoutInCircle place vertices on a circle, in the order of their 
#' vertex ids. The default is \code{FALSE}
#' @param header the main title
#' @param ... other arguments
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
#' ## check output
#' # draw network representation of the associations between responses and covariates
#' plotNetwork(fit)
#'
#' @export
plotNetwork <- function(x, includeResponse = NULL, excludeResponse = NULL, 
                        includePredictor = NULL, excludePredictor = NULL,
                        MatrixGamma = NULL, PmaxPredictor = 0.5, 
                        PmaxResponse = 0.5, nodesizePredictor = 2, 
                        nodesizeResponse = 15, no.isolates = FALSE,
                        lineup = 1.2, gray.alpha = 0.6, edgewith.response = 5, 
                        edgewith.predictor = 2, edge.weight = FALSE, 
                        label.predictor = NULL, label.response = NULL, 
                        color.predictor = NULL, color.response = NULL, 
                        name.predictors = NULL, name.responses = NULL,
                        vertex.frame.color = NA, layoutInCircle = FALSE, 
                        header = "", ...) {
  if (!inherits(x, "BayesSUR")) {
    stop("Use only with a \"BayesSUR\" object")
  }

  if (PmaxPredictor < 0 || PmaxPredictor > 1) {
    stop("Please specify correct argument 'PmaxPredictor' in [0,1]!")
  }
  if (PmaxResponse < 0 || PmaxResponse > 1) {
    stop("Please specify correct argument 'PmaxResponse' in [0,1]!")
  }

  x$output[-1] <- paste(x$output$outFilePath, x$output[-1], sep = "")
  covariancePrior <- x$input$covariancePrior
  if (covariancePrior == "HIW") {
    Gy_hat <- as.matrix(read.table(x$output$Gy))
  } else {
    stop("Gy is only estimated with hyper-inverse Wishart prior for the 
         covariance matrix of responses!")
  }
  gamma_hat <- as.matrix(read.table(x$output$gamma))
  colnames(gamma_hat) <- names(read.table(x$output$Y, header = TRUE))
  rownames(gamma_hat) <- colnames(read.table(x$output$X, header = TRUE))

  if (sum(colnames(gamma_hat) == 
          paste("V", seq_len(ncol(gamma_hat)), sep = "")) == ncol(gamma_hat)) {
    colnames(gamma_hat) <- paste("Y", seq_len(ncol(gamma_hat)), sep = "")
  }
  if (sum(rownames(gamma_hat) == 
          paste("V", seq_len(nrow(gamma_hat)), sep = "")) == nrow(gamma_hat)) {
    rownames(gamma_hat) <- paste("X", seq_len(nrow(gamma_hat)), sep = "")
  }

  # select the required resposes and predictors to plot the network
  excludeResponse.idx <- rep(FALSE, ncol(gamma_hat))
  excludePredictor.idx <- rep(FALSE, nrow(gamma_hat))
  if (!is.null(includeResponse)) {
    excludeResponse.idx <- c(!(colnames(gamma_hat) %in% includeResponse))
  }
  if (!is.null(excludeResponse)) {
    excludeResponse.idx <- c(excludeResponse.idx | 
                               c(colnames(gamma_hat) %in% excludeResponse))
  }
  if (!is.null(includePredictor)) {
    excludePredictor.idx <- c(!(rownames(gamma_hat) %in% includePredictor))
  }
  if (!is.null(excludePredictor)) {
    excludePredictor.idx <- c(excludePredictor.idx | 
                                c(rownames(gamma_hat) %in% excludePredictor))
  }

  gamma_hat <- gamma_hat[!excludePredictor.idx, !excludeResponse.idx]

  Gy_hat <- Gy_hat[!excludeResponse.idx, !excludeResponse.idx]

  if (edge.weight) {
    Gy_thresh <- Gy_hat
    Gy_thresh[Gy_hat <= PmaxResponse] <- 0

    gamma_thresh <- gamma_hat
    gamma_thresh[gamma_hat <= PmaxPredictor] <- 0
  } else {
    Gy_thresh <- as.matrix(Gy_hat > PmaxResponse)
    gamma_thresh <- as.matrix(gamma_hat > PmaxPredictor)
  }

  if (sum(rowSums(gamma_thresh) != 0) == 0) {
    stop(paste("There were no predictors with mPIP gamma > ", PmaxPredictor, 
               ". Not able to draw a network!", sep = ""))
  }

  gamma_thresh <- 
    matrix(gamma_thresh[rowSums(gamma_thresh) != 0, ], ncol = ncol(gamma_hat))
  colnames(gamma_thresh) <- colnames(gamma_hat)
  rownames(gamma_thresh) <- 
    rownames(gamma_hat)[rowSums(gamma_hat > PmaxPredictor) != 0]
  rownames(Gy_thresh) <- colnames(Gy_thresh) <- colnames(gamma_hat)

  plotSEMgraph(Gy_thresh, 
               t(gamma_thresh), 
               nodesizeSNP = nodesizePredictor, 
               nodesizeMET = nodesizeResponse, 
               no.isolates = no.isolates,  
               edgewith.response = edgewith.response, 
               edgewith.predictor = edgewith.predictor, 
               edge.weight = edge.weight, 
               label.predictor = label.predictor, 
               label.response = label.response, 
               color.predictor = color.predictor, 
               color.response = color.response,
               name.predictors = name.predictors, 
               name.responses = name.responses, 
               vertex.frame.color = vertex.frame.color, 
               layoutInCircle = layoutInCircle, ...)
  title(paste("\n\n", header, sep = ""), outer = TRUE)
}
plotSEMgraph <- function(ADJmatrix, 
                         GAMmatrix, 
                         nodesizeSNP = 2, 
                         nodesizeMET = 25, 
                         no.isolates = FALSE,
                         lineup = 1, 
                         gray.alpha = 0.6, 
                         edgewith.response = 5, 
                         edgewith.predictor = 2,
                         label.predictor = NULL, 
                         label.response = NULL, 
                         color.predictor = NULL, 
                         color.response = NULL,
                         name.predictors = NULL, 
                         name.responses = NULL, 
                         edge.weight = FALSE, 
                         vertex.frame.color = NA, 
                         layoutInCircle = FALSE, ...) {
  ## give warnings for re-defined arguments
  if (exists("edge.width")) 
    warning("Argument 'edge.width' was re-defined into new argments 
            'edgewith.response' and 'edgewith.predictor' in this function!")
  if (exists("edge.color")) 
    warning("Argument 'edge.color' cannot be changed  in this function!")
  if (exists("edge.arrow.size")) 
    warning("Argument 'edge.arrow.size' cannot be changed in this function!")

  # ADJmatrix must be a square qxq adjacency matrix (or data frame)
  qq <- dim(ADJmatrix)[1]
  if (dim(ADJmatrix)[2] != qq) stop("adjacency matrix not square")

  # GAMmatrix must be a qxp binary matrix (or data frame)
  pp <- dim(GAMmatrix)[2]
  if (dim(GAMmatrix)[1] != qq) stop("Gamma and Adjacency have different no. q")

  # join mets block (adjency) and lower triangle (gamma)
  semgraph <- rbind(ADJmatrix, t(GAMmatrix))

  # add zero blocks for lower triangle and snp block
  zeroblock <- matrix(rep(0, pp * (pp + qq)), nrow = qq + pp, ncol = pp)
  zeroblock <- data.frame(zeroblock)
  colnames(zeroblock) <- colnames(GAMmatrix)
  rownames(zeroblock) <- rownames(semgraph)

  semgraph <- cbind(semgraph, zeroblock)

  # igraph objects
  graphADJ <- graph.adjacency(as.matrix(ADJmatrix), weighted = TRUE, 
                              diag = FALSE, mode = "undirected")
  graphSEM <- graph.adjacency(as.matrix(semgraph), weighted = TRUE, 
                              diag = FALSE, mode = "directed")

  # don't plot isolated nodes?
  if (no.isolates) {
    graphADJ <- delete.vertices(graphADJ, degree(graphADJ) == 0)
    graphSEM <- delete.vertices(graphSEM, degree(graphSEM) == 0)
    message("Removing isolated nodes from Adjacency and Full SEM, may get 
            problem if there are also isolated SNPs.")
  }

  # get co-ords for undirected edges using layout function (scaled)
  lladj <- layout.fruchterman.reingold(graphADJ)
  lmax <- max(lladj[, 1])
  lmin <- min(lladj[, 1])
  lladj[, 1] <- (lladj[, 1] - lmin) / (lmax - lmin)
  lmax <- max(lladj[, 2])
  lmin <- min(lladj[, 2])
  lladj[, 2] <- (lladj[, 2] - lmin) / (lmax - lmin)

  # plot adjacency only
  # plot(graphADJ,vertex.size=15,edge.width=2,edge.color="black",layout=lladj)

  # line up snps
  lymax <- max(lladj[, 2]) + (max(lladj[, 2]) - 
                                min(lladj[, 2])) * (lineup - 1) / 2
  lymin <- min(lladj[, 2]) + (max(lladj[, 2]) - 
                                min(lladj[, 2])) * (1 - lineup) / 2
  llsnps <- matrix(c(rep(-0.5, pp), lymin + (1:pp) * 1.0 * 
                       (lymax - lymin) / pp), nrow = pp, ncol = 2)

  llsem <- rbind(lladj, llsnps)

  ### plot SEM

  # plot snps and mets nodes differently
  # set node sizes directly in graph object
  V(graphSEM)$size <- c(rep(nodesizeMET, qq), rep(nodesizeSNP, pp))

  n.edgeADJ <- gsize(graphADJ)
  n.edgeGAM <- gsize(graphSEM) - n.edgeADJ

  V(graphSEM)$label.color <- "black"

  V(graphSEM)$color <- c(rep("dodgerblue", nrow(GAMmatrix)), 
                         rep("red", ncol(GAMmatrix)))
  if (!is.null(color.predictor)) 
    V(graphSEM)$color[-c(seq_len(nrow(GAMmatrix)))] <- color.predictor
  if (!is.null(color.response)) 
    V(graphSEM)$color[seq_len(nrow(GAMmatrix))] <- color.response

  V(graphSEM)$label <- c(rownames(GAMmatrix), colnames(GAMmatrix))
  if (!is.null(label.predictor)) 
    V(graphSEM)$label[-c(seq_len(nrow(GAMmatrix)))] <- label.predictor
  if (!is.null(label.response)) 
    V(graphSEM)$label[seq_len(nrow(GAMmatrix))] <- label.response

  if (edge.weight) {
    edge.width <- E(graphSEM)$weight * ifelse(edge.weight, 5, 1)
  } else {
    edge.width <- c(rep(edgewith.response, 2 * n.edgeADJ), 
                    rep(edgewith.predictor, 2 * n.edgeGAM))
  }

  if (!layoutInCircle) {
    layoutSEM <- llsem
  } else {
    layoutSEM <- layout_in_circle(graphSEM)
  }

  # plot undirected graph between response variables
  graphSEMresponses <- delete.edges(
    graphSEM, E(graphSEM)[(1:ecount(graphSEM))[-c(1:(2 * n.edgeADJ))]])
  plot.igraph(graphSEMresponses, 
              edge.arrow.size = 0, 
              edge.width = edge.width[1:(2 * n.edgeADJ)], 
              vertex.frame.color = vertex.frame.color, 
              edge.color = rep(gray(0), 2 * n.edgeADJ), 
              layout = layoutSEM, ...)
  
  # plot directed graph between predictors and response variables
  graphSEMpredictor2responses <- 
    delete.edges(graphSEM, 
                 E(graphSEM)[(1:ecount(graphSEM))[c(1:(2 * n.edgeADJ))]])
  plot.igraph(graphSEMpredictor2responses, 
              edge.arrow.size = 0.5, 
              edge.width = edge.width[-c(1:(2 * n.edgeADJ))], 
              vertex.frame.color = vertex.frame.color, 
              edge.color = rep(gray(0.7, alpha = gray.alpha), 2 * n.edgeGAM), 
              layout = layoutSEM, add = TRUE, ...)

  if (!is.null(name.predictors)) text(-1, -1.3, name.predictors, cex = 1.2)
  if (!is.null(name.responses)) text(0.4, -1.3, name.responses, cex = 1.2)
}
