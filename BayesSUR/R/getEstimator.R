#' @title extract the posterior mean of parameters
#' @description
#' Extract the posterior mean of the parameters of a \code{BayesSUR} class object.
#' @name getEstimator
#' 
#' @param object an object of class \code{BayesSUR}
#' @param estimator the name of one estimator. Default is the latent indicator 
#' estimator "\code{gamma}". Other options "\code{beta}", "\code{Gy}", 
#' "\code{CPO}" and "\code{logP}"
#' correspond the marginal (conditional) coefficient matrix if 
#' \code{beta.type="marginal"}(\code{"conditional"}), response graph and 
#' conditional predictive ordinate (CPO) respectively
#' @param Pmax threshold that truncate the estimator "\code{gamma}" or 
#' "\code{Gy}". Default is \code{0}. If \code{Pmax=0.5} and 
#' \code{beta.type="conditional"}, it gives median probability model betas
#' @param beta.type the type of output beta. Default is \code{marginal}, giving 
#' marginal beta estimation. If \code{beta.type="conditional"}, it gives beta 
#' estimation conditional on gamma=1
#'
#' @return Return the estimator from an object of class \code{BayesSUR}. It is 
#' a matrix if the length of argument \code{marginal} is greater than 1. 
#' Otherwise, it is a list
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
#' # extract the posterior mean of the coefficients matrix
#' beta_hat <- getEstimator(fit, estimator = "beta")
#'
#' @export
getEstimator <- function(object, estimator = "gamma", Pmax = 0, 
                         beta.type = "marginal") {
  object$output[-1] <- 
    paste(object$output$outFilePath, object$output[-1], sep = "")
  
  if (sum(!estimator %in% c("gamma", "beta", "Gy", "CPO", "logP")) > 0) {
    stop("Please specify correct 'estimator'!")
  } else {
    ret <- rep(list(NULL), length(estimator))
    names(ret) <- estimator
  }

  if (Pmax < 0 || Pmax > 1) {
    stop("Please specify correct argument 'Pmax' in [0,1]!")
  }

  if ("gamma" %in% estimator) {
    ret$gamma <- as.matrix(read.table(object$output$gamma))
    if (Pmax > 0) {
      ret$gamma[ret$gamma <= Pmax] <- 0
    }
    rownames(ret$gamma) <- colnames(read.table(object$output$X, header = TRUE))
    colnames(ret$gamma) <- colnames(read.table(object$output$Y, header = TRUE))
  }

  if ("beta" %in% estimator) {
    ret$beta <- as.matrix(read.table(object$output$beta))

    if (sum(beta.type %in% c("marginal", "conditional")) > 0) {
      if (beta.type == "conditional") {
        gammas <- as.matrix(read.table(object$output$gamma))

        if ("X0" %in% names(object$output)) {
          X0 <- as.matrix(read.table(object$output$X0))
          ret$beta[-seq_len(ncol(X0)), ] <- 
            (gammas >= Pmax) * ret$beta[-seq_len(ncol(X0)), ] / gammas
        } else {
          ret$beta <- (gammas >= Pmax) * ret$beta / gammas
        }
        ret$beta[is.na(ret$beta)] <- 0
      }
    } else {
      stop("Please specify correct beta.type!")
    }

    colnames(ret$beta) <- colnames(read.table(object$output$Y, header = TRUE))
    if ("X0" %in% names(object$output)) {
      rownames(ret$beta) <- 
        c(colnames(read.table(object$output$X0, header = TRUE)), 
          colnames(read.table(object$output$X, header = TRUE)))
    } else {
      rownames(ret$beta) <- 
        colnames(read.table(object$output$X, header = TRUE))
    }
  }

  if ("Gy" %in% estimator) {
    covariancePrior <- object$input$covariancePrior
    if (covariancePrior == "HIW") {
      ret$Gy <- as.matrix(read.table(object$output$Gy))
    } else {
      stop("Gy is only estimated with hyper-inverse Wishart prior for the 
           covariance matrix of responses!")
    }
    if (Pmax > 0) {
      ret$Gy[ret$Gy <= Pmax] <- 0
    }
    rownames(ret$Gy) <- colnames(ret$Gy) <- 
      names(read.table(object$output$Y, header = TRUE))
  }

  if ("CPO" %in% estimator) {
    if (is.null(object$output$CPO)) {
      stop("Please specify argument output_CPO in BayesSUR()!")
    }

    ret$CPO <- as.matrix(read.table(object$output$CPO))

    rownames(ret$CPO) <- 
      rownames(as.matrix(read.table(object$output$Y, header = TRUE)))
    colnames(ret$CPO) <- 
      colnames(as.matrix(read.table(object$output$Y, header = TRUE)))
  }

  if ("logP" %in% estimator) {
    ret$logP <- t(as.matrix(read.table(object$output$logP)))
  }

  if (length(estimator) > 1) {
    return(ret)
  } else {
    return(ret[[1]])
  }
}
