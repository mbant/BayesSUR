## Build a new version of the package
remove.packages("R2SSUR")
Rcpp::compileAttributes(pkgdir = "R2SSUR/"); devtools::document("R2SSUR")
devtools::build("R2SSUR")#,vignettes=TRUE)

## Install the package
install.packages("R2SSUR_0.1.3.tar.gz",repos = NULL,type = "source")


#####################################################################################################
## Test the installation
data(example_data, package = "R2SSUR")

R2SSUR::runSSUR(example_data[["data"]],outFilePath = "results/",
                blockList = example_data[["blockList"]], structureGraph = example_data[["structureGraph"]],
                nIter = 100, nChains = 2, method = "SUR", sparse = TRUE )

## check output
greyscale = grey((1000:0)/1000)
data(example_ground_truth, package = "R2SSUR")

est_gamma = as.matrix( read.table("results/data_SSUR_gamma_out.txt") )
est_G = as.matrix( read.table("results/data_SSUR_G_out.txt") )
s = ncol(est_G)

par(mfrow=c(2,2))
image(est_gamma,col=greyscale); image(example_ground_truth[["gamma"]],col=greyscale)
image((est_G+diag(s))[s:1,],col=greyscale); image(example_ground_truth[["G"]][s:1,],col=greyscale)
par(mfrow=c(1,1))

