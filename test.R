## Build a new version of the package
remove.packages("R2SSUR")
Rcpp::compileAttributes(pkgdir = "/Users/zhiz/Downloads/R2SSUR/R2SSUR/"); devtools::document("/Users/zhiz/Downloads/R2SSUR/R2SSUR")
devtools::build("/Users/zhiz/Downloads/R2SSUR/R2SSUR")#,vignettes=TRUE)

## Install the package
install.packages("/Users/zhiz/Downloads/R2SSUR/R2SSUR_0.1.7.tar.gz",repos = NULL,type = "source")


#####################################################################################################
## Test the installation
data(example_data, package = "R2SSUR")

mrfGFile = example_data[["mrfG"]]
hyperpar = list(mrf_e=-3, mrf_d=3/10, b_pi = 0.2 , a_pi = 0.1 , b_o = 5 , a_o = 1 )

fit = R2SSUR::runSSUR(data = example_data[["data"]],
                Y = example_data[["blockList"]][[1]],
                X = example_data[["blockList"]][[2]][11:150],
                X_0 = example_data[["blockList"]][[2]][1:10],
                outFilePath = "results/",hyperpar=hyperpar,
                nIter = 2000, nChains = 2, covariancePrior = "IW", gammaPrior = "hotspot")

## check output
est_gamma = getEst( fit, "gamma")
est_G0 = getEst( fit, "G0")
est_beta = getEst( fit, "beta")
model_size = getEst( fit, "model_size")

sumint.fit = R2SSUR::summary(fit)

R2SSUR::plot(fit)

R2SSUR::mcmcDiag(fit)


greyscale = grey((100:0)/100)
data(example_ground_truth, package = "R2SSUR")
image(example_ground_truth[["gamma"]],col=greyscale)
plot.default( model_size[,1] , type="l",ylim = c(0,150))
points(1:ncol(example_ground_truth[["gamma"]]),colSums(example_ground_truth[["gamma"]]),col="red",pch=20)
