#' @title print method for class \code{BayesSUR}
#' @description
#' Print a short summary of a \code{BayesSUR} class object. It includes the 
#' argument matching information, number of selected predictors based on 
#' thresholding the posterior mean of the latent indicator variable at 0.5 
#' by default.
#' 
#' @name print.BayesSUR
#' @param x an object of class \code{BayesSUR}
#' @param Pmax threshold that truncates the estimated coefficients based on 
#' thresholding the estimated latent indicator variable. Default is 0.5
#' @param ... other arguments
#'
#' @return Return a short summary from an object of class \code{BayesSUR}, 
#' including the number of selected predictors with mPIP>\code{Pmax} and the 
#' expected log pointwise predictive density estimates (i.e., elpd.LOO and 
#' elpd.WAIC).
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
#'   hyperpar = hyperpar, tmpFolder = "tmp/", output_CPO = TRUE
#' )
#'
#' ## check output
#' # show the print information
#' print(fit)
#'
#' @export
print.BayesSUR <- function(x, Pmax = 0.5, ...) {
  gamma <- as.matrix(read.table(paste(x$output$outFilePath, 
                                      x$output$gamma, sep = "")))

  call.string <- unlist(strsplit(deparse(x$call), ","))
  call.string[which(call.string == " ") + 1] <- 
    substr(call.string[which(call.string == " ") + 1], 4, 
           nchar(call.string[which(call.string == " ") + 1]))
  call.string <- call.string[call.string != " "]

  if (length(call.string) <= 3) {
    cat("\nCall:\n ", paste(call.string, 
                            c(rep(",", length(call.string) - 1), ""), 
                            sep = "", collapse = ""), "\n", sep = "")
  } else {
    cat("\nCall:\n ", paste(call.string[1:3], c(",", ",", ", ...)"), sep = "", 
                            collapse = ""), "\n", sep = "")
  }

  cat("\nNumber of selected predictors (mPIP > ", Pmax, "): ", 
      sum(gamma > Pmax), " of ", ncol(gamma), "x", nrow(gamma), 
      "\n", sep = "")
  cat("\nExpected log pointwise predictive density (elpd):\n", 
      " elpd.LOO = ", elpd(x, method = "LOO"), 
      ",  elpd.WAIC = ", elpd(x, method = "WAIC"), "\n\n", sep = "")
}
