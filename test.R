### Build a new version of the package
#remove.packages("BayesSUR")
#Rcpp::compileAttributes(pkgdir="BayesSUR/")
#tools::resaveRdaFiles("BayesSUR/data/example_GDSC.rda", compress="xz")
#devtools::document("BayesSUR")
#devtools::build("BayesSUR", vignettes=TRUE, args="--compact-vignettes=both")

# $ R CMD check ./BayesSUR_1.0-2.tar.gz --as-cran
# devtools::check("BayesSUR", cran=TRUE, build_args="--compact-vignettes=both")
# $ ./BVS_Reg --covariancePrior HIW --gammaPrior MRF --dataFile data_test_GDSC.txt --blockFile blockLabels.txt
# --structureGraphFile structureGraph.txt --outFilePath results/ --nChains 2 --nIter 50 --burnin 0 --mrfGFile mrfG.txt


## Install the package
library(devtools)
install_github("mbant/BayesSUR/BayesSUR")
#install.packages("BayesSUR_1.0-2.tar.gz",repos = NULL,type = "source",build_vignettes = TRUE)


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
