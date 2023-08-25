#' @title plot MCMC diagnostic plots
#' @description
#' Show trace plots and diagnostic density plots of a fitted model object of class \code{BayesSUR}.
#' @importFrom graphics par plot.default legend title matplot
#' @importFrom stats density
#' @importFrom grDevices hcl.colors
#' @name plotMCMCdiag
#' @param x an object of class \code{BayesSUR}
#' @param nbloc number of splits for the last half iterations after substracting 
#' burn-in length
#' @param HIWg diagnostic plot of the response graph. Default is \code{NULL}. 
#' \code{HIW="degree"} prints the diagnostic of the degrees of response nodes. 
#' \code{HIW="edges"} prints the diagnostic of every edge between two responses. 
#' \code{HIW="lik"} prints the diagnostic of the posterior likelihoods of the 
#' hyperparameters related to the response relationships
#' @param header the main title
#' @param ... other arguments for the plots of the log-likelihood and model size
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
#' plotMCMCdiag(fit)
#'
#' @export
plotMCMCdiag <- function(x, nbloc = 3, HIWg = NULL, header = "", ...) {
  if (!inherits(x, "BayesSUR")) {
    stop("Use only with a \"BayesSUR\" object")
  }

  x$output[-1] <- paste0(x$output$outFilePath, x$output[-1])
  logP <- t(as.matrix(read.table(x$output$logP)))
  model_size <- as.matrix(read.table(x$output$model_size))
  ncol_Y <- ncol(read.table(x$output$gamma))
  nIter <- x$input$nIter

  covariancePrior <- x$input$covariancePrior
  if (covariancePrior == "HIW" && is.null(x$output$Gvisit)) {
    Gvisit <- as.matrix(read.table(x$output$Gvisit))
  }

  if (nIter <= 1) {
    stop("The diagosis only shows results from more than one MCMC iteration!")
  }
  if (nIter < 4000) {
    message("NOTE: The diagosis only shows results of two iteration points due 
            to less than 4000 MCMC iterations!")
  }

  if (nIter >= 4000) {
    logP <- logP[, ncol(logP) - floor(nIter / 1000) - 1 + 1:floor(nIter / 1000)]
  } else {
    logP <- logP[, c(1, ncol(logP))]
  }
  if (is.null(HIWg)) {
    Ptau.indx <- ifelse(covariancePrior != "IG", 7, 3)
    Plik.indx <- ifelse(covariancePrior != "IG", 10, 5)
    
    model_size <- model_size
    if (nIter >= 4000) {
      model_size <- rowSums(model_size[nrow(model_size) - floor(nIter / 1000) - 
                                         1 + 1:floor(nIter / 1000), ])
    } else {
      model_size <- rowSums(model_size[c(1, nrow(model_size)), ])
    }

    dens.all <- density(logP[Ptau.indx, ]) 
    if (nIter >= 4000) {
      dens.first <- density(logP[Ptau.indx, 1:floor(ncol(logP) / 2)])
      dens.last <- 
        density(logP[Ptau.indx, (1 + floor(ncol(logP) / 2)):ncol(logP)])

      ymin <- min(dens.all$y, dens.first$y, dens.last$y)
      ymax <- max(dens.all$y, dens.first$y, dens.last$y)
      xmin <- min(dens.all$x, dens.first$x, dens.last$x)
      xmax <- max(dens.all$x, dens.first$x, dens.last$x)
    } else {
      ymin <- min(dens.all$y)
      ymax <- max(dens.all$y)
      xmin <- min(dens.all$x)
      xmax <- max(dens.all$x)
      nbloc <- 1
    }

    ### nsplit number of split of the sweep
    mid <- floor(floor(ncol(logP) / 2) / nbloc)
    ymax2 <- xmin2 <- xmax2 <- list.dens <- NULL
    for (i in 1:nbloc) {
      dens <- density(logP[Ptau.indx, 
                           (ifelse(nbloc == 1, 0, floor(ncol(logP) / 2)) +
                              1 + mid * (i - 1)):ncol(logP)])
      ymax2 <- max(ymax2, dens$y)
      xmin2 <- min(xmin2, dens$x)
      xmax2 <- max(xmax2, dens$x)
      list.dens <- c(list.dens, list(dens))
    }

    ### plot the figures
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mfrow = c(2, 2))

    if (nbloc > 1) {
      plot.default(logP[Plik.indx, ], xlab = "Iterations (*1000)", 
                   ylab = "Log likelihood (posterior)", type = "l", lty = 1, ...)
    } else {
      plot.default(logP[Plik.indx, ] ~ c(1, nIter), xlab = "Iterations", 
                   ylab = "Log likelihood (posterior)", type = "l", lty = 1, ...)
    }

    if (nbloc > 1) {
      plot.default(model_size, xlab = "Iterations (*1000)", 
                   ylab = "Model size", type = "l", lty = 1, ...)
    } else {
      plot.default(model_size ~ c(1, nIter), xlab = "Iterations", 
                   ylab = "Model size", type = "l", lty = 1, ...)
    }

    title0 <- expression(paste("Log Posterior Distribution: log ", 
                               P(gamma ~ group("|", list(Y, .), ""))))
    plot.default(dens.all, main = "", col = "black", xlim = c(xmin, xmax), 
                 ylim = c(ymin, ymax), xlab = title0, ylab = "", 
                 type = "l", lty = 1)
    if (nIter >= 4000) {
      par(new = TRUE)
      plot.default(dens.first, main = "", col = "red", xlim = c(xmin, xmax), 
                   ylim = c(ymin, ymax), xlab = "", ylab = "", 
                   type = "l", lty = 1)
      par(new = TRUE)
      plot.default(dens.last, main = "", col = "green", xlim = c(xmin, xmax), 
                   ylim = c(ymin, ymax), xlab = "", , ylab = "Density", 
                   type = "l", lty = 1)
      if (nbloc > 1) {
        legend("topleft", title = "iteration", 
               legend = paste0(c("ALL", "First half", "Last half"), " = [", 
                              c(1, 1, floor(ncol(logP) / 2) * 1000 + 1), ":", 
                              c(ncol(logP), floor((ncol(logP)) / 2), 
                                ncol(logP)) * 1000, "]"), 
               col = 1:3, lty = 1, text.col = 1:3, cex = 0.8)
      }
    }
    for (i in 1:nbloc) {
      plot.default(list.dens[[i]], col = i, xlim = c(xmin2, xmax2), 
                   ylim = c(ymin, ymax2), xlab = title0, ylab = "", 
                   type = "l", lty = 1, main = "")
      if (nbloc > 1) par(new = TRUE)
    }
    title(ylab = "Density")
    if (nbloc > 1) {
      legend("topleft", title = "moving window", 
             legend = paste0("set ", 1:nbloc, " = [", 
                             (floor((ncol(logP)) / 2) + mid * (nbloc:1 - 1)) * 
                               1000 + 1, ":", (ncol(logP)) * 1000, "]"), 
             col = 1:nbloc, lty = 1, text.col = 1:nbloc, cex = 0.8)
    }
  } else {
    if (covariancePrior != "HIW") {
      stop("The argument HIWg only works for the model with hyper-inverse 
           Wishart prior on the covariance!")
    }

    if (HIWg == "degree") {
      m <- ncol_Y
      node1 <- node2 <- NULL
      for (i in 1:(m - 1)) {
        node1 <- c(node1, rep(i, m - i))
        node2 <- c(node2, (i + 1):m)
      }
      nodes <- cbind(node1, node2)
      node.degree <- matrix(0, nrow = nrow(Gvisit), ncol = m)
      for (i in 1:m) {
        node.degree[, i] <- 
          rowSums(Gvisit[, which(nodes == i, arr.ind = TRUE)[, 1]])
      }

      matplot(node.degree, type = "l", lty = 1, col = hcl.colors(m), 
              xlab = "Iterations (*1000)", ylab = "degree", 
              main = "Response degrees", xlim = c(1, nrow(Gvisit) * 1.1))
      legend("topright", legend = 1:m, col = hcl.colors(m), lty = 1, 
             text.col = hcl.colors(m), cex = 1 / m * 4)
    }

    if (substr(HIWg, 1, 4) == "edge") {
      m <- ncol_Y
      node1 <- node2 <- NULL
      for (i in 1:(m - 1)) {
        node1 <- c(node1, rep(i, m - i))
        node2 <- c(node2, (i + 1):m)
      }

      if (HIWg == "edge") {
        matplot(Gvisit, type = "l", lty = 1, col = hcl.colors(ncol(Gvisit)), 
                xlab = "Iterations (*1000)", ylab = "", 
                main = "Edges selection", xlim = c(1, nrow(Gvisit) * 1.1))
        legend("topright", legend = paste0(node1, "-", node2), 
               col = hcl.colors(ncol(Gvisit)), lty = 1,
               text.col = hcl.colors(ncol(Gvisit)), cex = 1 / m * 2)
      } else {
        plot.default(
          Gvisit[, which(paste0(node1, node2) == substr(HIWg, 5, nchar(HIWg)))], 
          type = "l", lty = 1, xlab = "Iterations (*1000)", ylab = "", 
          main = paste0("Edge-", substr(HIWg, 5, nchar(HIWg)), " selection"))
      }
    }

    if (HIWg == "lik") {
      Gvisit <- t(logP[1:4, ])
      matplot(Gvisit, type = "l", lty = 1, col = seq_len(ncol(Gvisit)), 
              xlab = "Iterations (*1000)", ylab = "Log likelihood (posterior)", 
              main = "Likelihoods of graph learning")
      legend("topright", legend = c("tau", "eta", "JT", "SigmaRho"), 
             col = seq_len(ncol(Gvisit)), lty = 1, 
             text.col = seq_len(ncol(Gvisit)), cex = 0.8)
    }
  }
  title(paste0("\n", header), outer = TRUE)
}
