#' @title plot heatmap of estimators
#' @description
#' Plot the posterior mean estimators from a \code{BayesSUR} class object, 
#' including the coefficients beta, latent indicator variable gamma and 
#' graph of responses.
#' @importFrom graphics axis box text mtext par image
#' @importFrom grDevices colorRampPalette dev.off grey
#' @importFrom tikzDevice tikz
#' @name plotEstimator
#' 
#' @param x an object of class \code{BayesSUR}
#' @param estimator print the heatmap of estimators. The value "beta" is for 
#' the estimated coefficients matrix, "gamma" for the latent indicator matrix 
#' and "Gy" for the graph of responses
#' @param colorScale.gamma value palette for gamma
#' @param colorScale.beta a vector of three colors for diverging color schemes
#' @param legend.cex.axis magnification of axis annotation relative to cex
#' @param name.responses a vector of the response names. The default is 
#' \code{NA} only to show the locations. The value \code{"auto"} show the 
#' response names from the orginal data.
#' @param name.predictors a vector of the predictor names. The default is 
#' \code{NA} only to show the locations. The value \code{"auto"} show the 
#' predictor names from the orginal data.
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param fig.tex print the figure through LaTex. Default is \code{FALSE}
#' @param output the file name of printed figure
#' @param header the main title
#' @param header.cex size of the main title for all estimators
#' @param cex.main size of the title for each estimator
#' @param mgp the margin line (in mex units) for the axis title, axis labels 
#' and axis line
#' @param title.beta a title for the printed "beta"
#' @param title.gamma a title for the printed "gamma"
#' @param title.Gy a title for the printed "Gy"
#' @param tick a logical value specifying whether tickmarks and an axis line 
#' should be drawn. Default is \code{FALSE}
#' @param beta.type the type of output beta. Default is \code{marginal}, giving 
#' marginal beta estimation. If \code{beta.type="conditional"}, it gives beta 
#' estimation conditional on gamma=1
#' @param Pmax threshold that truncate the estimator "\code{gamma}" or 
#' "\code{Gy}". Default is \code{0}. If \code{Pmax=0.5} and 
#' \code{beta.type="conditional"}, it gives median probability model betas
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
#' # Plot the estimators from the fitted object
#' plotEstimator(fit, estimator = c("beta", "gamma", "Gy"))
#'
#' \dontrun{
#' ## Set up temporary work directory for saving a pdf figure
#' # td <- tempdir()
#' # oldwd <- getwd()
#' # setwd(td)
#'
#' ## Produce authentic math formulas in the graph
#' # plotEstimator(fit, estimator = c("beta", "gamma", "Gy"), fig.tex = TRUE)
#' # system(paste(getOption("pdfviewer"), "ParamEstimator.pdf"))
#' # setwd(oldwd)
#' }
#'
#' @export
plotEstimator <- function(x, estimator = NULL, 
                          colorScale.gamma = grey((100:0) / 100), 
                          colorScale.beta = c("blue", "white", "red"), 
                          legend.cex.axis = 1, name.responses = NA,
                          name.predictors = NA, xlab = "", ylab = "", 
                          fig.tex = FALSE, output = "ParamEstimator", 
                          header = "", header.cex = 2, tick = FALSE, 
                          mgp = c(2.5, 1, 0), cex.main = 1.5, 
                          title.beta = NA, title.gamma = NA, title.Gy = NA, 
                          beta.type = "marginal", Pmax = 0, ...) {
  if (!inherits(x, "BayesSUR")) {
    stop("Use only with a \"BayesSUR\" object")
  }
  if (sum(!estimator %in% c("beta", "gamma", "Gy")) > 0) {
    stop("Please specify correct argument estimator!")
  }
  
  beta_hat <- getEstimator(x, estimator = "beta", 
                           beta.type = beta.type, Pmax = Pmax)
  x$output[-1] <- paste0(x$output$outFilePath, x$output[-1])
  gamma_hat <- as.matrix(read.table(x$output$gamma))
  response_idx <- seq_len(ncol(gamma_hat))
  predictor_idx <- seq_len(nrow(gamma_hat))
  nonpen <- nrow(beta_hat) - nrow(gamma_hat)
  if (nonpen > 0) {
    rownames(beta_hat) <- c(colnames(read.table(x$output$X0, header = TRUE)), 
                            colnames(read.table(x$output$X, header = TRUE)))
  } else {
    rownames(beta_hat) <- colnames(read.table(x$output$X, header = TRUE))
  }
  colnames(beta_hat) <- colnames(read.table(x$output$Y, header = TRUE))

  covariancePrior <- x$input$covariancePrior
  if ((covariancePrior == "HIW") && ("Gy" %in% estimator)) {
    Gy_hat <- as.matrix(read.table(x$output$Gy))
    colnames(Gy_hat) <- rownames(Gy_hat) <- 
      colnames(read.table(x$output$Y, header = TRUE))
  }

  ## BUG TO BE FIXED!!!
  # specify the labels of axes
  if (is.na(name.responses)[1]) name.responses <- response_idx
  if (name.responses[1] == "auto") name.responses <- colnames(beta_hat)
  if (is.character(name.responses)) {
    if (length(name.responses) != ncol(beta_hat)) {
      stop("The length of the given response names are not consistent with the data!")
    }
  }

  if (is.na(name.predictors)) name.predictors <- predictor_idx
  if (name.predictors[1] == "auto") name.predictors <- rownames(beta_hat)
  if (is.character(name.predictors)) {
    if (length(name.predictors) != nrow(beta_hat)) {
      stop("The length of the given predictor names are not consistent with the data!")
    }
  }

  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if (!fig.tex) {
    par(mar = c(6, 6, 5.1, 4.1))
    par(mfrow = c(1, sum(estimator %in% c("beta", "gamma", "Gy"))))

    if ("beta" %in% estimator) {
      # floor(100*constant)+100-1 colours that your want in the legend bar which has the white middle color
      if (sd(beta_hat) == 0) {
        colorbar <- colorRampPalette(c(colorScale.beta[2], colorScale.beta[2]))(1000)
      } else {
        colorbar <- c(colorRampPalette(c(colorScale.beta[1], colorScale.beta[2]))(
          floor(1000 / (-(max(beta_hat) - min(beta_hat)) / min(beta_hat) - 1))), 
          colorRampPalette(c(colorScale.beta[2], colorScale.beta[3]))(1000)[-1])
      }
      if (is.na(title.beta)) title.beta <- expression(hat(bold(B)))

      image(
        z = beta_hat, x = predictor_idx, y = response_idx, 
        col = colorbar, mgp = mgp,
        axes = ifelse(is.na(name.responses)[1], TRUE, FALSE), 
        xlab = xlab, ylab = ylab, main = title.beta, 
        cex.main = cex.main, cex.lab = 1.5, ...)
      box()
      vertical.image.legend(color = colorbar, 
                            zlim = c(min(beta_hat), max(beta_hat)), 
                            legend.cex.axis = legend.cex.axis)
      if (!is.na(name.responses)[1]) {
        par(las = 2, cex.axis = 1)
        axis(2, at = response_idx, labels = name.responses, tick = tick)
        axis(1, at = predictor_idx, labels = name.predictors, tick = tick)
      }
    }
    if ("gamma" %in% estimator) {
      if (is.na(title.gamma)) title.gamma <- expression(hat(bold(Gamma)))
      image(
        z = gamma_hat, x = predictor_idx, y = response_idx, 
        col = colorScale.gamma, mgp = mgp,
        axes = ifelse(is.na(name.responses)[1], TRUE, FALSE), 
        xlab = xlab, ylab = ylab, main = title.gamma, 
        cex.main = cex.main, cex.lab = 1.5, ...)
      box()
      vertical.image.legend(color = colorScale.gamma, zlim = c(0, 1), 
                            legend.cex.axis = legend.cex.axis)
      if (!is.na(name.responses)[1]) {
        par(las = 2, cex.axis = 1)
        axis(2, at = response_idx, labels = name.responses, tick = tick)
        if (nonpen > 0) {
          name.predictors <- name.predictors[-c(1:nonpen)]
        }
        axis(1, at = predictor_idx, labels = name.predictors, tick = tick)
      }
    }

    if ("Gy" %in% estimator) {
      if (is.na(title.Gy)) title.Gy <- "Estimated graph of responses"
      image(
        z = Gy_hat + diag(ncol(Gy_hat)), 
        x = response_idx, y = response_idx, 
        col = colorScale.gamma, mgp = mgp,
        axes = ifelse(is.na(name.responses)[1], TRUE, FALSE), 
        xlab = ylab, ylab = ylab, main = title.Gy, 
        cex.main = cex.main, cex.lab = 1.5, ...)
      box()
      vertical.image.legend(color = colorScale.gamma, 
                            zlim = c(min(Gy_hat), max(Gy_hat)), 
                            legend.cex.axis = legend.cex.axis)
      if (!is.na(name.responses)[1]) {
        par(las = 2, cex.axis = 1)
        axis(2, at = response_idx, labels = name.responses, tick = tick)
        axis(1, at = response_idx, labels = name.responses, tick = tick)
      }
    }
    title(paste0("\n", header), cex.main = header.cex, outer = TRUE)
  } else {
    options(tikzMetricPackages = c("\\usepackage{amsmath}", "\\usepackage{bm}", 
                                   "\\usetikzlibrary{calc}"))
    tikz(paste0(output, ".tex"),
      width = 3.6 * sum(estimator %in% c("beta", "gamma", "Gy")), 
      height = 4, 
      standAlone = TRUE,
      packages = c("\\usepackage{tikz}", 
                   "\\usepackage{amsmath}", 
                   "\\usepackage{bm}", 
                   "\\usepackage[active,tightpage,psfixbb]{preview}",
                   "\\PreviewEnvironment{pgfpicture}"))

    par(mfrow = c(1, sum(estimator %in% c("beta", "gamma", "Gy"))))
    par(mar = c(6, 6, 4, 4) + 0.1)
    if ("beta" %in% estimator) {
      # floor(100*constant)+100-1 colours that your want in the legend bar which has the white middle color
      colorbar <- 
        c(colorRampPalette(c(colorScale.beta[1], colorScale.beta[2]))(
          floor(1000 / (-(max(beta_hat) - min(beta_hat)) / min(beta_hat) - 1))), 
          colorRampPalette(c(colorScale.beta[2], colorScale.beta[3]))(1000)[-1])
      if (is.na(title.beta)) 
        title.beta <- paste("Estimator", "$\\hat{\\bm{B}}$")

      image(
        z = beta_hat, x = predictor_idx, y = response_idx, col = colorbar, 
        axes = ifelse(is.na(name.responses)[1], TRUE, FALSE), mgp = mgp,
        xlab = xlab, ylab = ylab, main = title.beta, 
        cex.main = cex.main, cex.lab = 1.5, ...)
      box()
      vertical.image.legend(color = colorbar, 
                            zlim = c(min(beta_hat), max(beta_hat)), 
                            legend.cex.axis = legend.cex.axis)
      if (!is.na(name.responses)[1]) {
        par(las = 2, cex.axis = 1)
        axis(2, at = response_idx, labels = name.responses, tick = tick)
        # opar <- par(cex.axis=1)
        axis(1, at = predictor_idx, labels = name.predictors, tick = tick)
      }
    }
    if ("gamma" %in% estimator) {
      if (is.na(title.gamma)) 
        title.gamma <- paste("Estimator", "$\\hat{\\mathbf{\\Gamma}}$")
      image(
        z = gamma_hat, x = predictor_idx, y = seq_len(ncol(gamma_hat)), 
        col = colorScale.gamma, 
        axes = ifelse(is.na(name.responses)[1], TRUE, FALSE), 
        xlab = xlab, ylab = ylab, main = title.gamma, mgp = mgp,
        cex.main = cex.main, cex.lab = 1.5, ...)
      box()
      vertical.image.legend(color = colorScale.gamma, 
                            zlim = c(0, 1), legend.cex.axis = legend.cex.axis)
      if (!is.na(name.responses)[1]) {
        par(las = 2, cex.axis = 1)
        axis(2, at = seq_len(ncol(gamma_hat)), 
             labels = name.responses, tick = tick)
        if (nonpen > 0) {
          name.predictors <- name.predictors[-c(1:nonpen)]
        }
        axis(1, at = predictor_idx, labels = name.predictors, tick = tick)
      }
    }

    if ("Gy" %in% estimator) {
      if (is.na(title.Gy)) 
        title.Gy <- paste("Estimator", "$\\hat{\\mathcal{G}}$")
      image(
        z = Gy_hat + diag(ncol(Gy_hat)), 
        x = seq_len(nrow(Gy_hat)), 
        y = seq_len(nrow(Gy_hat)), 
        col = colorScale.gamma, 
        axes = ifelse(is.na(name.responses)[1], TRUE, FALSE), 
        mgp = mgp,
        xlab = ylab, 
        ylab = ylab, 
        main = title.Gy, 
        cex.main = cex.main, 
        cex.lab = 1.5, ...)
      box()
      vertical.image.legend(color = colorScale.gamma, 
                            zlim = c(min(Gy_hat), max(Gy_hat)), 
                            legend.cex.axis = legend.cex.axis)

      if (!is.na(name.responses)[1]) {
        par(las = 2, cex.axis = 1)
        axis(2, at = seq_len(ncol(Gy_hat)), 
             labels = name.responses, tick = tick)
        axis(1, at = seq_len(nrow(Gy_hat)), 
             labels = name.responses, tick = tick)
      }
    }
    title(paste0("\n", header), cex.main = header.cex, outer = TRUE)
    dev.off()
    tools::texi2pdf(paste0(output, ".tex"))
  }
}
# the function vertical.image.legend() is adapted from the R package "aqfig"
vertical.image.legend <- function(zlim, color, legend.cex.axis = 1) {
  starting.par.settings <- par(no.readonly = TRUE)
  mai <- par("mai")
  fin <- par("fin")
  x.legend.fig <- c(1 - (mai[4] / fin[1]), 1)
  y.legend.fig <- c(mai[1] / fin[2], 1 - (mai[3] / fin[2]))
  x.legend.plt <- c(x.legend.fig[1] + (0.18 * (x.legend.fig[2] -
    x.legend.fig[1])), x.legend.fig[2] - (0.6 * (x.legend.fig[2] -
    x.legend.fig[1])))
  y.legend.plt <- y.legend.fig
  cut.pts <- seq(zlim[1], zlim[2], length = length(color) + 1)
  z <- (cut.pts[seq_along(color)] + cut.pts[2:(length(color) + 1)]) / 2
  par(new = TRUE, pty = "m", plt = c(x.legend.plt, y.legend.plt))
  # If z is not increasing, only two values
  if (all(diff(z) > 0)) {
    image(
      x = 1.5, y = z, z = matrix(z, nrow = 1, ncol = length(color)),
      col = color, xlab = "", ylab = "", xaxt = "n", yaxt = "n"
    )
    axis(4, mgp = c(3, 0.2, 0), las = 2, cex.axis = legend.cex.axis, tcl = -0.1)
    box()
  }
  mfg.settings <- par()$mfg
  par(starting.par.settings)
  par(mfg = mfg.settings, new = FALSE)
}
