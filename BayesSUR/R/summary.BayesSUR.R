#' @title summary method for class \code{BayesSUR}
#' @description
#' Summary method for class \code{BayesSUR}. It includes the argument matching 
#' information, Top predictors/responses on average mPIP across all 
#' responses/predictors, elpd estimates, MCMC specification, model 
#' specification and hyper-parameters. The summarized number of the selected 
#' variable corresponds to the posterior mean of the latent indicator variable 
#' thresholding at 0.5 by default.
#'
#' @importFrom Matrix Matrix
#' @name summary.BayesSUR
#' @param object an object of class \code{BayesSUR}
#' @param Pmax threshold that truncates the estimated coefficients based on 
#' thresholding the estimated latent indicator variable. Default is 0.5
#' @param ... other arguments
#'
#' @return Return a result summary from an object of class \code{BayesSUR}, 
#' including the CPOs, number of selected predictors with mPIP>\code{Pmax}, 
#' top 10 predictors on average mPIP across all responses, top 10 responses on 
#' average mPIP across all predictors, Expected log pointwise predictive 
#' density (elpd) estimates, MCMC specification, model specification (i.e., 
#' covariance prior and gamma prior) and hyper-parameters.
#'
#' @examples
#' data(exampleEQTL, package = "BayesSUR")
#' hyperpar <- list(a_w = 2, b_w = 5)
#'
#' set.seed(9173)
#' fit <- BayesSUR(
#'   Y = exampleEQTL[["blockList"]][[1]],
#'   X = exampleEQTL[["blockList"]][[2]],
#'   data = exampleEQTL[["data"]], outFilePath = tempdir(),
#'   nIter = 10, burnin = 0, nChains = 1, gammaPrior = "hotspot",
#'   hyperpar = hyperpar, tmpFolder = "tmp/", output_CPO = TRUE
#' )
#'
#' ## check output
#' # show the summary information
#' summary(fit)
#'
#' @export
summary.BayesSUR <- function(object, Pmax = 0.5, ...) {
  if (Pmax < 0 || Pmax > 1) {
    stop("Please specify correct argument 'Pmax' in (0,1)!")
  }

  ans <- list(status = object$status)
  if (is.null(object$output$CPO)) {
    ans$elpd <- NA
  } else {
    ans$elpd <- c(elpd(object, method = "LOO"), elpd(object, method = "WAIC"))
  }

  object$output[-1] <- 
    paste(object$output$outFilePath, object$output[-1], sep = "")
  if (is.null(object$output$CPO)) {
    ans$CPO <- NA
  } else {
    ans$CPO <- 
      summary.default(as.vector(as.matrix(read.table(object$output$CPO))))[-4]
  }

  gamma <- as.matrix(read.table(object$output$gamma))
  ans$df <- sum(gamma > Pmax)

  # extract top 10 covariates based on average mPIP across all responses
  mean.predictors <- rowMeans(gamma)
  top10.predictors <- mean.predictors[
    sort.list(mean.predictors, decreasing = TRUE)[1:min(10, nrow(gamma))]]
  names(top10.predictors) <- names(read.table(object$output$X, header = TRUE))[
    sort.list(mean.predictors, decreasing = TRUE)[1:min(10, nrow(gamma))]]

  # extract top 10 response variables based on average mPIP across responses
  mean.responses <- colMeans(gamma)
  top10.responses <- mean.responses[
    sort.list(mean.responses, decreasing = TRUE)[1:min(10, ncol(gamma))]]
  names(top10.responses) <- names(read.table(object$output$Y, header = TRUE))[
    sort.list(mean.responses, decreasing = TRUE)[1:min(10, ncol(gamma))]]

  ans$chainParameters <- object$input[1:3]
  ans$modelParameters <- object$input[4:9]
  ans$hyperParameters <- object$input$hyperparameters
  ans$outputFiles <- object$output
  ans$outputFiles["outFilePath"] <- NULL

  cat("\nCall:\n  ", 
      paste(unlist(strsplit(deparse(object$call), ","))[1:3], 
            c(",", ",", ", ...)"), sep = "", collapse = ""), "\n", sep = "")
  cat("\nCPOs:\n")
  print(ans$CPO)
  cat("\nNumber of selected predictors (mPIP > ", Pmax, "): ", 
      sum(gamma > Pmax), " of ", ncol(gamma), "x", nrow(gamma), 
      "\n", sep = "")
  cat("\nTop", min(10, nrow(gamma)), 
      "predictors on average mPIP across all responses:\n")
  print(top10.predictors)

  cat("\nTop", min(10, ncol(gamma)), 
      "responses on average mPIP across all predictors:\n")
  print(top10.responses)

  cat("\nExpected log pointwise predictive density (elpd) estimates:\n", 
      "  elpd.LOO = ", ans$elpd[1], ",  elpd.WAIC = ", ans$elpd[2], "\n", 
      sep = "")
  cat("\nMCMC specification:\n", "  iterations = ", ans$chainParameters$nIter, 
      ",  burn-in = ", ans$chainParameters$burnin, 
      ",  chains = ", ans$chainParameters$nChains,
      "\n  gamma local move sampler: ", ans$modelParameters$gammaSampler, 
      "\n  gamma initialisation: ", ans$modelParameters$gammaInit, "\n",
      sep = "")
  cat("\nModel specification:\n", "  covariance prior: ", 
      ans$modelParameters$covariancePrior, "\n  gamma prior: ", 
      ans$modelParameters$gammaPrior, "\n", sep = "")

  if (is.null(ans$hyperParameters)) {
    cat("\nHyper-parameters:\n")
    print(unlist(object$input$hyperParameters))
  }
  cat("\n")

  invisible(ans)
}
