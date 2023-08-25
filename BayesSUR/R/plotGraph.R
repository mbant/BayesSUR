#' @title plot graph for response variables
#' @description
#' Plot the estimated graph for multiple response variables from a 
#' \code{BayesSUR} class object.
#' @importFrom igraph E plot.igraph graph_from_adjacency_matrix
#' @importFrom graphics par
#' @name plotGraph
#' 
#' @param x either an object of class \code{BayesSUR} (default) or a symmetric 
#' numeric matrix representing an adjacency matrix for a given graph structure.
#' If x is an adjacency matrix, argument \code{main="Given graph of responses"} 
#' by default.
#' @param Pmax a value for thresholding the learning structure matrix of 
#' multiple response variables. Default is 0.5
#' @param main an overall title for the plot
#' @param edge.width edge width. Default is 2
#' @param edge.weight draw weighted edges after thresholding at 0.5. The 
#' default value \code{FALSE} is not to draw weighted edges
#' @param vertex.label character vector used to label the nodes
#' @param vertex.label.color label color. Default is \code{"black"}
#' @param vertex.size node size. Default is 30
#' @param vertex.color node color. Default is \code{"dodgerblue"}
#' @param vertex.frame.color node color. Default is \code{"NA"}
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
#' # show the graph relationship between responses
#' plotGraph(fit, estimator = "Gy")
#'
#' @export
plotGraph <- function(x, Pmax = 0.5, main = "Estimated graph of responses", 
                      edge.width = 2, edge.weight = FALSE, vertex.label = NULL, 
                      vertex.label.color = "black", vertex.size = 30, 
                      vertex.color = "dodgerblue", vertex.frame.color = NA, ...) {
  if (!inherits(x, "BayesSUR")) {
    if (is.matrix(x) && is.numeric(x)) {
      if (!((dim(x)[1] == dim(x)[2]) && (sum(dim(x)) > 2))) {
        stop("Use only with a \"BayesSUR\" object or numeric square matrix")
      }
      Gy_hat <- x
      if (!is.null(vertex.label)) {
        rownames(Gy_hat) <- colnames(Gy_hat) <- vertex.label
      }
      if (main == "Estimated graph of responses") {
        main <- "Given graph of responses"
      }
    } else {
      stop("Use only with a \"BayesSUR\" object or numeric square matrix")
    }
  } else {
    x$output[-1] <- paste0(x$output$outFilePath, x$output[-1])
    covariancePrior <- x$input$covariancePrior
    if (covariancePrior == "HIW") {
      Gy_hat <- as.matrix(read.table(x$output$Gy))
    } else {
      stop("Gy is only estimated with hyper-inverse Wishart prior for the 
           covariance matrix of responses!")
    }

    if (!is.null(vertex.label)) {
      rownames(Gy_hat) <- colnames(Gy_hat) <- vertex.label
    } else {
      rownames(Gy_hat) <- colnames(Gy_hat) <- 
        names(read.table(x$output$Y, header = TRUE))
    }
  }

  if (Pmax < 0 || Pmax > 1) {
    stop("Please specify correct argument 'Pmax' in [0,1]!")
  }

  if (edge.weight) {
    Gy_thresh <- Gy_hat
    Gy_thresh[Gy_hat <= Pmax] <- 0
  } else {
    Gy_thresh <- as.matrix(Gy_hat > Pmax)
  }

  net <- graph_from_adjacency_matrix(Gy_thresh, weighted = TRUE, 
                                     mode = "undirected", diag = FALSE)

  if (edge.weight) {
    plot.igraph(net, 
                main = main, 
                edge.width = E(net)$weight * 2, 
                vertex.label = vertex.label, 
                vertex.color = vertex.color, 
                vertex.frame.color = vertex.frame.color, ...)
  } else {
    plot.igraph(net, 
                main = main, 
                edge.width = edge.width, 
                vertex.label = vertex.label, 
                vertex.color = vertex.color, 
                vertex.frame.color = vertex.frame.color, 
                vertex.label.color = vertex.label.color, 
                vertex.size = vertex.size, ...)
  }
}
