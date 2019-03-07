## Build a new version of the package
remove.packages("R2SSUR")
Rcpp::compileAttributes(pkgdir = "R2SSUR/"); devtools::document("R2SSUR")
devtools::build("R2SSUR")#,vignettes=TRUE)

## Install the package
install.packages("R2SSUR_0.1.5.tar.gz",repos = NULL,type = "source")


#####################################################################################################
## Test the installation
data(example_data, package = "R2SSUR")

# mrfGFile = as.matrix( read.table("mrf.txt") )

fit = R2SSUR::runSSUR(data = example_data[["data"]],
                Y = example_data[["blockList"]][[1]],
                X = example_data[["blockList"]][[2]][11:150],
                X_0 = example_data[["blockList"]][[2]][1:10],
                outFilePath = "results/",
                nIter = 3000, nChains = 2, covariancePrior = "IW", gammaPrior = "hotspot" )

## check output
greyscale = grey((1000:0)/1000)
data(example_ground_truth, package = "R2SSUR")

est_gamma = as.matrix( read.table(fit$output$gamma) )
est_G = as.matrix( read.table(fit$output$G) )
est_beta = as.matrix( read.table(fit$output$beta) )
s = ncol(est_G)

par(mfrow=c(1,3))
image(est_gamma,col=greyscale); image(example_ground_truth[["gamma"]],col=greyscale)
image(est_beta,col=greyscale)
#image((est_G+diag(s))[s:1,],col=greyscale); image(example_ground_truth[["G"]][s:1,],col=greyscale)
par(mfrow=c(1,1))


plot( read.table(fit$output$model_size)[,1] , type="l",ylim = c(0,150))
points(1:ncol(example_ground_truth[["gamma"]]),colSums(example_ground_truth[["gamma"]]),col="red",pch=20)
