#' @title Fitting BayesSUR models
#' @description
#' Main function of the package. Fits a range of models introduced in the 
#' package vignette \code{BayesSUR.pdf}. Returns an object of S3 class 
#' \code{BayesSUR}. There are three options for the prior on the residual 
#' covariance matrix (i.e., independent inverse-Gamma, inverse-Wishart and 
#' hyper-inverse Wishart) and three options for the prior on the latent 
#' indicator variable (i.e., independent Bernoulli, hotspot and Markov random 
#' field). So there are nine models in total. See details for their combinations.
#'
#' @useDynLib BayesSUR
#' @aliases BayesSUR-package
#' @importFrom utils head tail read.table write.table
#' @importFrom Rcpp sourceCpp
#' @importFrom xml2 as_xml_document write_xml
#' @importFrom parallel detectCores
#'
#' @name BayesSUR
#' @param data a numeric matrix with variables on the columns and observations 
#' on the rows, if arguments \code{Y} and \code{X} (and possibly \code{X_0}) 
#' are vectors. Can be \code{NULL} if arguments \code{Y} and \code{X} (and 
#' possibly \code{X_0}) are numeric matrices
#' @param Y,X vectors of indices (with respect to the data matrix) for the 
#' outcomes (\code{Y}) and the predictors to select (\code{X}) respectively; 
#' if the \code{data} argument is \code{NULL}, these needs to be numeric 
#' matrices containing the data instead, with variables on the columns and 
#' observations on the rows
#' @param X_0 vectors of indices (with respect to the data matrix) for the 
#' fixed predictors that are not selected, i.e. always included in the model; 
#' if the data argument is not provided, this needs to be a numeric matrix 
#' containing the data instead, with variables on the columns and observations 
#' on the rows
#' @param betaPrior string indicating the prior for regression coefficients; it 
#' has to be either \code{independent} for independent spike-and-slab priors 
#' (only slab part for \code{X_0} if specified), or \code{reGroup} for weakly 
#' normal priors for mandatory variables (random effects) and spike-and-slab 
#' priors for other variables of Zhao (2023)
#' @param covariancePrior string indicating the prior for the covariance $C$; 
#' it has to be either \code{HIW} for the hyper-inverse-Wishar (which will 
#' result in a sparse covariance matrix), \code{IW} for the inverse-Wishart 
#' prior (dense covariance) or \code{IG} for independent inverse-Gamma on all 
#' the diagonal elements and 0 otherwise. See the details for the model 
#' specification
#' @param gammaPrior string indicating the gamma prior to use, either 
#' \code{hotspot} (default) for the Hotspot prior of Bottolo (2011), \code{MRF} 
#' for the Markov Random Field prior or \code{hierarchical} for a simpler 
#' hierarchical prior. See the details for the model specification
#' @param nIter number of iterations for the MCMC procedure. Default 10000
#' @param burnin number of iterations to discard at the start of the chain. 
#' Default is 5000
#' @param nChains number of parallel tempered chains to run (default 2). The 
#' temperature is adapted during the burnin phase
#' @param outFilePath path to where the output files are to be written
#' @param gammaSampler string indicating the type of sampler for gamma, either 
#' \code{bandit} for the Thompson sampling inspired samper or \code{MC3} for 
#' the usual MC^3 sampler.  See Russo et al.(2018) or Madigan and York (1995) 
#' for details
#' @param gammaInit gamma initialisation to either all-zeros (\code{0}), all 
#' ones (\code{1}), MLE-informed (\code{MLE}) or (default) randomly (\code{R})
#' @param mrfG either a matrix or a path to the file containing (the edge list 
#' of) the G matrix for the MRF prior on gamma (if necessary)
#' @param standardize logical flag for X variable standardization. Default is 
#' \code{standardize=TRUE}. Coefficients are returned on the standardized scale
#' @param standardize.response logical flag for Y standardization. Default is 
#' \code{standardize.response=TRUE}
#' @param hyperpar a list of named hypeparameters to use instead of the default 
#' values. Valid names are mrf_d, mrf_e, a_sigma, b_sigma, a_tau, b_tau, nu, 
#' a_eta, b_eta, a_o, b_o, a_pi, b_pi, a_w and b_w. Their default values are 
#' a_w=2, b_w=5, a_omega=2, b_omega=1, a_o=2, b_o=p-2, a_pi=2, b_pi=1, nu=s+2, 
#' a_tau=0.1, b_tau=10, a_eta=0.1, b_eta=1, a_sigma=1, b_sigma=1, mrf_d=-3 and 
#' mrf_e=0.03. See the vignette for more information
#' @param maxThreads maximum threads used for parallelization. Default is 1. 
#' Reproducibility of results with \code{set.seed()} is only guaranteed if 
#' \code{maxThreads=1}
#' @param tick an integer used for printing the iteration index and some updated 
#' parameters every tick-th iteration. Default is 1000 
#' @param output_gamma allow (\code{TRUE}) or suppress (\code{FALSE}) the 
#' output for  gamma. See the return value below for more information
#' @param output_beta allow (\code{TRUE}) or suppress (\code{FALSE}) the output 
#' for beta. See the return value below for more information
#' @param output_Gy allow (\code{TRUE}) or suppress (\code{FALSE}) the output 
#' for Gy. See the return value below for more information
#' @param output_sigmaRho allow (\code{TRUE}) or suppress (\code{FALSE}) the 
#' output for sigmaRho. See the return value below for more information
#' @param output_pi allow (\code{TRUE}) or suppress (\code{FALSE}) the output 
#' for pi. See the return value below for more information
#' @param output_tail allow (\code{TRUE}) or suppress (\code{FALSE}) the output 
#' for tail (hotspot tail probability). See the return value below for more 
#' information
#' @param output_model_size allow (\code{TRUE}) or suppress (\code{FALSE}) the 
#' output for model_size. See the return value below for more information
#' @param output_model_visit allow (\code{TRUE}) or suppress (\code{FALSE}) the 
#' output for all visited models over the MCMC iterations. Default is 
#' \code{FALSE}. See the return value below for more information
#' @param output_CPO allow (\code{TRUE}) or suppress (\code{FALSE}) the output 
#' for (scaled) conditional predictive ordinates (\code{*_CPO_out.txt}),
#' CPO with joint posterior predictive of the response variables 
#' (\code{*_CPOsumy_out.txt}) and widely applicable information criterion 
#' (\code{*_WAIC_out.txt}). See the return value below for more information
#' @param output_Y allow (\code{TRUE}) or suppress (\code{FALSE}) the output 
#' for responses dataset Y
#' @param output_X allow (\code{TRUE}) or suppress (\code{FALSE}) the output 
#' for predictors dataset X
#' @param tmpFolder the path to a temporary folder where intermediate data 
#' files are stored (will be erased at the end of the chain). It is specified 
#' relative to \code{outFilePath}
#'
#' @details The arguments \code{covariancePrior} and \code{gammaPrior} specify 
#' the model HRR, dSUR or SSUR with different gamma prior. Let 
#' \eqn{\gamma_{jk}} be latent indicator variable of each coefficient and 
#' \eqn{C} be covariance matrix of response variables. The nine models 
#' specified through the arguments \code{covariancePrior} and 
#' \code{gammaPrior} are as follows.
#' \tabular{cccc}{
#'                 \tab \eqn{\gamma_{jk}}~Bernoulli \tab \eqn{\gamma_{jk}}~hotspot \tab \eqn{\gamma}~MRF \cr
#'   \eqn{C}~indep \tab HRR-B                       \tab HRR-H                     \tab HRR-M           \cr
#'   \eqn{C}~IW   \tab dSUR-B                      \tab dSUR-H                    \tab dSUR-M          \cr
#'   \eqn{C}~HIW   \tab SSUR-B                      \tab SSUR-H                    \tab SSUR-M
#' }
#'
#' @return An object of class \code{BayesSUR} is saved as 
#' \code{obj_BayesSUR.RData} in the output file, including the following 
#' components:
#' \itemize{
#' \item status - the running status
#' \item input - a list of all input parameters by the user
#' \item output - a list of the all output filenames:
#' \itemize{
#' \item "\code{*_logP_out.txt}" - contains each row for the \eqn{1000t}-th iteration's log-likelihoods of parameters, i.e., Tau, Eta, JunctionTree, SigmaRho, O, Pi, Gamma, W, Beta and data conditional log-likelihood depending on the models.
#' \item "\code{*_gamma_out.txt}" - posterior mean of the latent indicator matrix.
#' \item "\code{*_pi_out.txt}" - posterior mean of the predictor effects (prospensity) by decomposing the probability of the latent indicator.
#' \item "\code{*_hotspot_tail_p_out.txt}" - posterior mean of the hotspot tail probability. Only available for the hotspot prior on the gamma.
#' \item "\code{*_beta_out.txt}" - posterior mean of the coefficients matrix.
#' \item "\code{*_Gy_out.txt}" - posterior mean of the response graph. Only available for the HIW prior on the covariance.
#' \item "\code{*_sigmaRho_out.txt}" - posterior mean of the transformed parameters. Not available for the IG prior on the covariance.
#' \item "\code{*_model_size_out.txt}" - contains each row for the\eqn{1000t}-th iteration's model sizes of the multiple response variables.
#' \item "\code{*_model_visit_gy_out.txt}" - contains each row for the nonzero indices of the vectorized estimated graph matrix for each iteration.
#' \item "\code{*_model_visit_gamma_out.txt}" - contains each row for the nonzero indices of the vectorized estimated gamma matrix for each iteration.
#' \item "\code{*_CPO_out.txt}" - the (scaled) conditional predictive ordinates (CPO).
#' \item "\code{*_CPOsumy_out.txt}" - the (scaled) conditional predictive ordinates (CPO) with joint posterior predictive of the response variables.
#' \item "\code{*_WAIC_out.txt}" - the widely applicable information criterion (WAIC).
#' \item "\code{*_Y.txt}" - responses dataset.
#' \item "\code{*_X.txt}" - predictors dataset.
#' \item "\code{*_X0.txt}" - fixed predictors dataset.
#' }
#' \item call - the matched call.
#' }
#'
#' @references Russo D, Van Roy B, Kazerouni A, Osband I, Wen Z (2018). \emph{A tutorial on Thompson sampling.} Foundations and Trends in Machine Learning, 11: 1-96.
#' @references Madigan D, York J (1995). \emph{Bayesian graphical models for discrete data.} International Statistical Review, 63: 215–232.
#' @references Bottolo L, Banterle M, Richardson S, Ala-Korpela M, Jarvelin MR, Lewin A (2020). \emph{A computationally efficient Bayesian seemingly unrelated regressions model for high-dimensional quantitative trait loci discovery.} Journal of Royal Statistical Society: Series C, 70: 886-908.
#' @references Zhao Z, Banterle M, Bottolo L, Richardson S, Lewin A, Zucknick M (2021). \emph{BayesSUR: An R package for high-dimensional multivariate Bayesian variable and covariance selection in linear regression.} Journal of Statistical Software, 100: 1–32.
#' @references Zhao Z, Banterle M, Lewin A, Zucknick M (2023). \emph{Multivariate Bayesian structured variable selection for pharmacogenomic studies.} Journal of the Royal Statistical Society: Series C (Applied Statistics), qlad102.
#'
#' @examples
#' data("exampleEQTL", package = "BayesSUR")
#' hyperpar <- list(a_w = 2, b_w = 5)
#' set.seed(9173)
#' fit <- BayesSUR(
#'   Y = exampleEQTL[["blockList"]][[1]],
#'   X = exampleEQTL[["blockList"]][[2]],
#'   data = exampleEQTL[["data"]], outFilePath = tempdir(),
#'   nIter = 5, burnin = 0, nChains = 1, gammaPrior = "hotspot",
#'   hyperpar = hyperpar, tmpFolder = "tmp/", output_CPO = TRUE
#' )
#'
#' ## check output
#' # show the summary information
#' summary(fit)
#'
#' # show the estimated beta, gamma and graph of responses Gy
#' plot(fit, estimator = c("beta", "gamma", "Gy"), type = "heatmap")
#'
#' \dontrun{
#' ## Set up temporary work directory for saving a pdf figure
#' # td <- tempdir()
#' # oldwd <- getwd()
#' # setwd(td)
#'
#' ## Produce authentic math formulas in the graph
#' # plot(fit, estimator = c("beta", "gamma", "Gy"), type = "heatmap", fig.tex = TRUE)
#' # system(paste(getOption("pdfviewer"), "ParamEstimator.pdf"))
#' # setwd(oldwd)
#' }
#'
#' @export
BayesSUR <- function(data = NULL, Y, X, X_0 = NULL,
                     covariancePrior = "HIW", gammaPrior = "hotspot", 
                     betaPrior = "independent",
                     nIter = 10000, burnin = 5000, nChains = 2,
                     outFilePath = "", 
                     gammaSampler = "bandit", gammaInit = "R", mrfG = NULL,
                     standardize = TRUE, standardize.response = TRUE, 
                     maxThreads = 1, tick = 1000, 
                     output_gamma = TRUE, output_beta = TRUE, output_Gy = TRUE, 
                     output_sigmaRho = TRUE, output_pi = TRUE, 
                     output_tail = TRUE, output_model_size = TRUE, 
                     output_model_visit = FALSE, output_CPO = FALSE, 
                     output_Y = TRUE, output_X = TRUE, 
                     hyperpar = list(), tmpFolder = "tmp/") {
  # Check the directory for the output files
  if (outFilePath == "") {
    stop("Please specify a directory to save all output files!")
  }
  
  outFilePathLength <- nchar(outFilePath)
  if (substr(outFilePath, outFilePathLength, outFilePathLength) != "/") {
    outFilePath <- paste(outFilePath, "/", sep = "")
  }
  if (!file.exists(outFilePath)) {
    dir.create(outFilePath)
  }
  
  # Create temporary directory
  tmpFolderLength <- nchar(tmpFolder)
  if (substr(tmpFolder, tmpFolderLength, tmpFolderLength) != "/") {
    tmpFolder <- paste(tmpFolder, "/", sep = "")
  }
  tmpFolder <- paste(outFilePath, tmpFolder, sep = "")
  if (!file.exists(tmpFolder)) {
    dir.create(tmpFolder)
  }
  
  ## Check the input: reasoning is that the user provides either
  # a data matrix or a data path-to-file
  #     - in this case Y, X (and X_0) need to be provided as vectors of indexes
  # if the data matrix is not provided, the user will give 2/3 matrices for Y, X (and X_0)
  #     - in which case we write those in order into a new joint file
  # everything else throws an error
  
  # check the formula
  cl <- match.call()
  
  # we'll check in reverse order, is data NULL?
  if (is.null(data)) {
    # Y,X (and if there X_0) need to be valid numeric matrices then
    # check Y and X have comfortable number of observations
    npY <- dim(Y)
    if ((!is.numeric(Y)) || is.null(npY)) {
      my_stop("If 'data' is NULL, Y should be a numeric matrix", tmpFolder)
    }
    
    npX <- dim(X)
    if ((!is.numeric(X)) || is.null(npX) || (npX[1] != npY[1])) {
      my_stop("If 'data' is NULL, X should be a numeric matrix and the same 
              number of rows of Y", tmpFolder)
    }
    
    if (is.null(X_0)) {
      X_0 <- matrix(NA, nrow = npY[1], ncol = 0)
    } else {
      npX0 <- dim(X_0)
      if ((!is.numeric(X_0)) || is.null(npX0) || (npX0[1] != npY[1])) {
        my_stop("If 'data' is NULL and X_0 is provided, X_0 should be a numeric 
                matrix and the same number of rows of Y", tmpFolder)
      }
    }
    
    # Standarize the data
    if (standardize) {
      X <- scale(X)
      X_0 <- scale(X_0)
    }
    if (standardize.response) Y <- scale(Y)
    
    # Write the three down in a single data file
    write.table(cbind(Y, X, X_0), paste(sep = "", tmpFolder, "data.txt"), 
                row.names = FALSE, col.names = FALSE)
    data <- paste(sep = "", tmpFolder, "data.txt")
    
    blockLabels <- c(rep(0, ncol(Y)), rep(1, ncol(X)), rep(2, ncol(X_0)))
    
    # Write the data in a output file
    write.table(Y, paste(sep = "", outFilePath, "data_Y.txt"), 
                row.names = FALSE, col.names = TRUE)
    write.table(X, paste(sep = "", outFilePath, "data_X.txt"), 
                row.names = FALSE, col.names = TRUE)
    write.table(X_0, paste(sep = "", outFilePath, "data_X0.txt"), 
                row.names = FALSE, col.names = TRUE)
  } else { # data is not null, so the user wants to use it to input the data
    
    # is the data given as matrix?
    ## If it's valid matrix, simply write it and re-assign the variable data to hold its path
    
    npData <- dim(data)
    if (is.numeric(data) || (!is.null(npData)) || (npData[2] >= 2)) {
      # Standarize the data
      if (standardize) {
        data[, X] <- scale(data[, X])
        data[, X_0] <- scale(data[, X_0])
      }
      if (standardize.response) data[, Y] <- scale(data[, Y])
      
      # Write the Y and X data in a output file
      write.table(data[, Y], paste(sep = "", outFilePath, "data_Y.txt"), 
                  row.names = FALSE, col.names = TRUE)
      write.table(data[, X], paste(sep = "", outFilePath, "data_X.txt"), 
                  row.names = FALSE, col.names = TRUE)
      write.table(data[, X_0], paste(sep = "", outFilePath, "data_X0.txt"), 
                  row.names = FALSE, col.names = TRUE)
      
      write.table(data, paste(sep = "", tmpFolder, "data.txt"), 
                  row.names = FALSE, col.names = FALSE)
      data <- paste(sep = "", tmpFolder, "data.txt")
    } else {
      my_stop("Y should be NULL or a numeric matrix with 2 or more columns!")
    }
    
    # is the data given as a string?
    if (is.character(data) && length(data) == 1) {
      if (substr(data, 1, 1) == "~") {
        data <- path.expand(data)
      }
    }
    
    ## at this point data contains the path to a file that exists
    # try and read one line to check dimensions
    dataHeader <- read.table(data, header = FALSE, nrows = 1)
    nVariables <- ncol(dataHeader)
    
    ## Y, X (and X_0) should be some fixed variables that needs to be included in the model
    if (is.null(X_0)) {
      X_0 <- as.numeric(c())
    }
    
    # be sure that they are vectors
    if (!(is.vector(Y, "numeric") && is.vector(X, "numeric") && 
          is.vector(X_0, "numeric"))) {
      my_stop("When the 'data' argument is set, Y,X and X_0 need to be 
              corresponding index vectors!", tmpFolder)
    }
    
    # check thay do not overlap
    if (length(c(intersect(Y, X), intersect(Y, X_0), intersect(X_0, X))) != 0) {
      my_stop("Y, X and X_0 need to be distinct index vectors!", tmpFolder)
    }
    
    # check if dimensions correspond -- higher dimensions gets an error
    if (length(c(Y, X, X_0)) > nVariables) {
      my_stop("When the 'data' argument is set, Y,X and X_0 need to be 
              corresponding index vectors!", tmpFolder)
    }
    # equal dimensions are ok, but lower dimensions means some columns of the data will be disregarded ( set to -1  )
    
    # We can now init the blockList
    blockLabels <- rep(NA, nVariables)
    blockLabels[Y] <- 0
    blockLabels[X] <- 1
    if (length(X_0) > 0) {
      blockLabels[X_0] <- 2
    }
    
    blockLabels[is.na(blockLabels)] <- -1
  }
  
  # cleanup file PATHS
  dataLength <- nchar(data)
  if (dataLength == 0) {
    my_stop("Please provide a correct path to a plain-text (.txt) file", 
            tmpFolder)
  }
  
  # magicly strip '/' from the start and '.txt' from the end of the data file name
  dataString <- 
    head(strsplit(tail(strsplit(data, split = c("/"))[[1]], 1), ".txt")[[1]], 1)
  
  ## Then init the structure graph
  # Consider that the indexes are written so that Y is 0 , X is 1 and (if there) X_0 is 2
  if (length(X_0) > 0) {
    structureGraph <- structureGraph <- matrix(c(0, 0, 0, 1, 0, 0, 2, 0, 0), 
                                               3, 3, byrow = TRUE)
  } else {
    structureGraph <- structureGraph <- matrix(c(0, 0, 1, 0), 2, 2, 
                                               byrow = TRUE)
  }
  
  ## Finally write blockLabels and structureGraph to a file
  write.table(blockLabels, paste(sep = "", tmpFolder, "blockLabels.txt"), 
              row.names = FALSE, col.names = FALSE)
  blockList <- paste(sep = "", tmpFolder, "blockLabels.txt")
  write.table(structureGraph, paste(sep = "", tmpFolder, "structureGraph.txt"), 
              row.names = FALSE, col.names = FALSE)
  structureGraph <- paste(sep = "", tmpFolder, "structureGraph.txt")
  
  # check how burnin was given
  if (burnin < 0) {
    my_stop("Burnin must be positive or 0", tmpFolder)
  } else {
    if (burnin > nIter) {
      my_stop("Burnin might not be greater than nIter", tmpFolder)
    } else {
      if (burnin < 1) { # given as a fraction
        burnin <- ceiling(nIter * burnin) # the zero case is taken into account here as well
      }
    }
  } # else assume is given as an absolute number
  
  ###############################
  ## prepare the print of hyperparameters corresponding the specified model
  hyperpar.all <- list(a_w = 2, b_w = 5, 
                       a_o = 2, b_o = sum(blockLabels == 1) - 2, 
                       a_pi = NA, b_pi = NA, 
                       nu = sum(blockLabels == 0) + 2, 
                       a_tau = 0.1, b_tau = 10, 
                       a_eta = 0.1, b_eta = 1, 
                       a_sigma = 1, b_sigma = 1, 
                       mrf_d = -3, mrf_e = 0.03)
  if (toupper(gammaPrior) %in% c("HOTSPOT", "HOTSPOTS", "HS")) {
    hyperpar.all$a_pi <- 2
    hyperpar.all$b_pi <- 1
    
    if (toupper(covariancePrior) %in% c("INDEPENDENT", "INDEP", "IG")) {
      hyperpar.all <- hyperpar.all[-c(7:11, 14:15)]
    }
    if (toupper(covariancePrior) %in% c("DENSE", "IW")) {
      hyperpar.all <- hyperpar.all[-c(10:11, 12:15)]
    }
    if (toupper(covariancePrior) %in% c("SPARSE", "HIW")) {
      hyperpar.all <- hyperpar.all[-c(12:15)]
    }
  }
  if (toupper(gammaPrior) %in% c("HIERARCHICAL", "H")) {
    hyperpar.all$a_pi <- 1
    hyperpar.all$b_pi <- sum(blockLabels == 0) - 1
    
    if (toupper(covariancePrior) %in% c("INDEPENDENT", "INDEP", "IG")) {
      hyperpar.all <- hyperpar.all[-c(3:6, 7:11, 14:15)]
    }
    if (toupper(covariancePrior) %in% c("DENSE", "IW")) {
      hyperpar.all <- hyperpar.all[-c(3:6, 10:11, 12:15)]
    }
    if (toupper(covariancePrior) %in% c("SPARSE", "HIW")) {
      hyperpar.all <- hyperpar.all[-c(3:6, 12:15)]
    }
  }
  if (toupper(gammaPrior) %in% c("MRF", "MARKOV RANDOM FIELD")) {
    if (is.null(mrfG)) {
      my_stop("Argument 'mrfG' was specified!", tmpFolder)
    }
    
    if (toupper(covariancePrior) %in% c("INDEPENDENT", "INDEP", "IG")) {
      hyperpar.all <- hyperpar.all[-c(3:6, 7:13)]
    }
    if (toupper(covariancePrior) %in% c("DENSE", "IW")) {
      hyperpar.all <- hyperpar.all[-c(3:6, 10:13)]
    }
    if (toupper(covariancePrior) %in% c("SPARSE", "HIW")) {
      hyperpar.all <- hyperpar.all[-c(3:6, 12:13)]
    }
  }
  if (toupper(betaPrior) == "REGROUP") {
    hyperpar.all$a_w0 <- hyperpar.all$a_w
    hyperpar.all$b_w0 <- hyperpar.all$b_w
  }
  
  if (length(hyperpar) > 0) {
    for (i in seq_along(hyperpar)) {
      if (names(hyperpar)[[i]] %in% names(hyperpar.all)) {
        hyperpar.all[[which(names(hyperpar.all) == names(hyperpar)[[i]])]] <- 
          hyperpar[[i]]
      }
      #if (!is.null(hyperpar$a_omega)) hyperpar.all$a_pi <- hyperpar$a_omega
      #if (!is.null(hyperpar$b_omega)) hyperpar.all$b_pi <- hyperpar$b_omega
    }
  }
  
  if (toupper(gammaPrior) %in% c("HIERARCHICAL", "H")) {
    ## re-use the hyperparameter names a_pi and b_pi for 
    ## Bernoulli prior's hyperparameters a_omega and b_omega
    #names(hyperpar.all)[names(hyperpar.all) == "a_pi"] <- "a_omega"
    #names(hyperpar.all)[names(hyperpar.all) == "b_pi"] <- "b_omega"
    if (!is.null(hyperpar$a_omega)) hyperpar.all$a_pi <- hyperpar$a_omega
    if (!is.null(hyperpar$b_omega)) hyperpar.all$b_pi <- hyperpar$b_omega
  }
  
  # method to use
  if (toupper(covariancePrior) %in% c("SPARSE", "HIW")) {
    covariancePrior <- "HIW"
  } else if (toupper(covariancePrior) %in% c("DENSE", "IW")) {
    covariancePrior <- "IW"
  } else if (toupper(covariancePrior) %in% c("INDEPENDENT", "INDEP", "IG")) {
    covariancePrior <- "IG"
  } else {
    my_stop("Unknown covariancePrior argument: only sparse (HIW), dense(IW) or 
            independent (IG) are available", tmpFolder)
  }
  
  # mrfG and gammaPrior
  if (gammaPrior == "") {
    if (is.null(mrfG)) {
      message("Using default prior for Gamma - hotspot prior\n")
      # mrfG=""
      gammaPrior <- "hotspot"
    } else {
      message("No value for gammaPrior was specified, but mrfG was given - 
              choosing MRF prior\n")
      gammaPrior <- "MRF"
    }
  } else {
    if (toupper(gammaPrior) %in% c("HOTSPOT", "HOTSPOTS", "HS")) {
      gammaPrior <- "hotspot"
    } else if (toupper(gammaPrior) %in% c("MRF", "MARKOV RANDOM FIELD")) {
      gammaPrior <- "MRF"
    } else if (toupper(gammaPrior) %in% c("HIERARCHICAL", "H")) {
      gammaPrior <- "hierarchical"
    } else {
      my_stop("Unknown gammaPrior argument: only hotspot, MRF or hierarchical 
              are available", tmpFolder)
    }
  }
  
  # if mrfG is not a string
  if (!(is.character(mrfG) && length(mrfG) == 1)) {
    # if it's a matrix
    if ((is.numeric(mrfG) || is.data.frame(mrfG)) && !is.null(dim(mrfG))) {
      if (ncol(mrfG) == 2) {
        mrfG <- cbind(mrfG, rep(1, nrow(mrfG)))
      }
      write.table(mrfG, paste(sep = "", outFilePath, "mrfG.txt"), 
                  row.names = FALSE, col.names = FALSE)
      mrfG <- paste(sep = "", outFilePath, "mrfG.txt")
    } else if (is.null(mrfG)) {
      # save a meaningless mrfG.txt file to pass the parameter to C++
      mrfG <- matrix(c(0, 0, 0), ncol = 3)
      write.table(mrfG, paste(sep = "", outFilePath, "mrfG.txt"), 
                  row.names = FALSE, col.names = FALSE)
      mrfG <- paste(sep = "", outFilePath, "mrfG.txt")
    } else {
      my_stop("Unknown mrfG argument: check the help function for possibile 
              values", tmpFolder)
    }
  }
  
  ## Set up the XML file for hyperparameters
  xml <- as_xml_document(
    list(hyperparameters = list(
      lapply(hyperpar.all, function(x) list(x)) # every element in the list should be a list
    ))
  )
  hyperParFile <- paste(sep = "", tmpFolder, "hyperpar.xml")
  write_xml(xml, file = hyperParFile)
  
  ## Create the return object
  ret <- list(status = 1, input = list(), output = list())
  class(ret) <- "BayesSUR"
  
  # Copy the inputs
  ret$input["nIter"] <- nIter
  ret$input["burnin"] <- burnin
  ret$input["nChains"] <- nChains
  ret$input["covariancePrior"] <- covariancePrior
  ret$input["gammaPrior"] <- gammaPrior
  ret$input["gammaSampler"] <- gammaSampler
  ret$input["gammaInit"] <- gammaInit
  ret$input["mrfG"] <- mrfG
  
  ret$input$hyperParameters <- hyperpar.all
  
  methodString <-
    switch(covariancePrior,
           "HIW" = "SSUR",
           "IW"  = "dSUR",
           "IG"  = "HRR"
    )
  
  ret$call <- cl
  
  # Prepare path to outputs
  ret$output["outFilePath"] <- outFilePath
  
  ret$output["logP"] <- 
    paste(sep = "", dataString, "_", methodString, "_logP_out.txt")
  
  if (output_gamma) {
    ret$output["gamma"] <- 
      paste(sep = "", dataString, "_", methodString, "_gamma_out.txt")
  }
  
  if (gammaPrior %in% c("hierarchical", "hotspot") && output_pi) {
    ret$output["pi"] <- 
      paste(sep = "", dataString, "_", methodString, "_pi_out.txt")
  }
  
  if (gammaPrior == "hotspot" && output_tail) {
    ret$output["tail"] <- 
      paste(sep = "", dataString, "_", methodString, "_hotspot_tail_p_out.txt")
  }
  
  if (output_beta) {
    ret$output["beta"] <- 
      paste(sep = "", dataString, "_", methodString, "_beta_out.txt")
  }
  
  if (covariancePrior == "HIW" && output_Gy) {
    ret$output["Gy"] <- 
      paste(sep = "", dataString, "_", methodString, "_Gy_out.txt")
    ret$output["Gvisit"] <- 
      paste(sep = "", dataString, "_", methodString, "_Gy_visit.txt")
  }
  
  if (covariancePrior %in% c("HIW", "IW") && output_sigmaRho) {
    ret$output["sigmaRho"] <- 
      paste(sep = "", dataString, "_", methodString, "_sigmaRho_out.txt")
  }
  
  if (output_model_size) {
    ret$output["model_size"] <- 
      paste(sep = "", dataString, "_", methodString, "_model_size_out.txt")
  }
  
  if (output_CPO) {
    ret$output["CPO"] <- 
      paste(sep = "", dataString, "_", methodString, "_CPO_out.txt")
    ret$output["CPOsumy"] <- 
      paste(sep = "", dataString, "_", methodString, "_CPOsumy_out.txt")
    ret$output["WAIC"] <- 
      paste(sep = "", dataString, "_", methodString, "_WAIC_out.txt")
  }
  
  if (output_Y) {
    ret$output["Y"] <- paste(sep = "", "data_Y.txt")
  }
  
  if (output_X) {
    ret$output["X"] <- paste(sep = "", "data_X.txt")
    if (length(X_0) > 0) {
      ret$output["X0"] <- paste(sep = "", "data_X0.txt")
    }
  }
  
  # set.seed(seed)
  # betaPrior="independent"
  
  # set number of threads
  maxThreads <- min(maxThreads, detectCores())
  
  ret$status <- BayesSUR_internal(data, mrfG, blockList, structureGraph, 
                                  hyperParFile, outFilePath, nIter, burnin, 
                                  nChains, covariancePrior, gammaPrior, 
                                  gammaSampler, gammaInit, betaPrior, 
                                  maxThreads, tick, output_gamma, output_beta, 
                                  output_Gy, output_sigmaRho, output_pi, 
                                  output_tail, output_model_size, output_CPO, 
                                  output_model_visit)
  
  if (toupper(gammaPrior) %in% c("HIERARCHICAL", "H")) {
    ## bring back the re-used hyperparameter names a_pi and b_pi for 
    ## Bernoulli prior's hyperparameters a_omega and b_omega
    names(ret$input$hyperParameters)[names(ret$input$hyperParameters) == "a_pi"] <- "a_omega"
    names(ret$input$hyperParameters)[names(ret$input$hyperParameters) == "b_pi"] <- "b_omega"
  }
  ## save fitted object
  obj_BayesSUR <- list(status = ret$status, input = ret$input, 
                       output = ret$output, call = ret$call)
  save(obj_BayesSUR, file = paste(sep = "", outFilePath, "obj_BayesSUR.RData"))

  if (outFilePath != tmpFolder) {
    unlink(tmpFolder, recursive = TRUE)
  }

  return(ret)
}


my_stop <- function(msg, tmpFolder) {
  unlink(tmpFolder, recursive = TRUE)
  stop(msg)
}
