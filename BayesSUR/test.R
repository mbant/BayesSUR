## Build a new version of the package
remove.packages("BayesSUR")
Rcpp::compileAttributes(pkgdir = "/Users/zhiz/Downloads/BayesSUR/BayesSUR/"); devtools::document("/Users/zhiz/Downloads/BayesSUR/BayesSUR")
devtools::build("/Users/zhiz/Downloads/BayesSUR/BayesSUR")#,vignettes=TRUE)

## Install the package
install.packages("/Users/zhiz/Downloads/BayesSUR/BayesSUR_0.1.1.tar.gz",repos = NULL,type = "source")


#####################################################################################################
## Test the installation
library(BayesSUR)
data(example_data, package = "BayesSUR")

hyperpar = list(mrf_e=-3, mrf_d=3/10, b_pi = 0.2 , a_pi = 0.1 , b_o = 5 , a_o = 1 )

fit = runSUR(data = example_data[["data"]],
                Y = example_data[["blockList"]][[1]],
                X = example_data[["blockList"]][[2]][11:150],
                outFilePath = "results/",hyperpar=hyperpar,
                nIter = 2000, nChains = 2, covariancePrior = "HIW", gammaPrior = "hotspot")

## check output
est_gamma = getEstimator( fit, "gamma")
est_G0 = getEstimator( fit, "G0")
est_beta = getEstimator( fit, "beta")
model_size = getEstimator( fit, "model_size")

sumint.fit = summary(fit)

plotEstimator(fit)

plotResponseGraph(fit, PtrueResponse=example_data$G0)

plotNetwork(fit)

mcmcDiag(fit)


greyscale = grey((100:0)/100)
data(example_ground_truth, package = "BayesSUR")
image(example_ground_truth[["gamma"]],col=greyscale)
plot.default( model_size[,1] , type="l",ylim = c(0,150))
points(1:ncol(example_ground_truth[["gamma"]]),colSums(example_ground_truth[["gamma"]]),col="red",pch=20)
