### Build a new version of the package
#remove.packages("BayesSUR")
#Rcpp::compileAttributes(pkgdir="BayesSUR/")
#tools::resaveRdaFiles("BayesSUR/data/example_GDSC.rda", compress="xz")
#devtools::document("BayesSUR")
#devtools::build("BayesSUR", vignettes=TRUE, args="--compact-vignettes=both")

### CRAN check
#devtools::check("BayesSUR", build_args="--compact-vignettes=both")


## Install the package
library(devtools)
install_github("mbant/BayesSUR/BayesSUR")
#install.packages("BayesSUR_0.1.24.tar.gz",repos = NULL,type = "source")


#####################################################################################################
## load the package
library(BayesSUR)
data(example_eQTL, package = "BayesSUR")
str(example_eQTL)

# fit a SSUR model with hotspot prior
attach(example_eQTL)
fit <- BayesSUR(data = data, Y = blockList[[1]],
              X = blockList[[2]], outFilePath = "results/", 
              nIter = 2000, nChains = 5, covariancePrior = "HIW", burnin=1000,
              gammaPrior = "hotspot")

summary(fit)

# show the interaction of plots
plot(fit)

# show the estimated beta, gamma and G_y
plotEstimator(fit, fig.tex=TRUE)
system(paste(getOption("pdfviewer"), "ParamEstimator.pdf"))

# show the relationship of responses
plotResponseGraph(fit)

# show the network representation of the associations between responses and features
plotNetwork(fit)

# show the manhattan plot
plotManhattan(fit)

# check the convergence of the algorithm
plotMCMCdiag(fit)
