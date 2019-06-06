## Build a new version of the package
remove.packages("BayesSUR")
Rcpp::compileAttributes(pkgdir = "/Users/zhiz/Downloads/BayesSUR/BayesSUR/"); devtools::document("/Users/zhiz/Downloads/BayesSUR/BayesSUR")
devtools::build("/Users/zhiz/Downloads/BayesSUR/BayesSUR")#,vignettes=TRUE)

## Install the package
install.packages("/Users/zhiz/Downloads/BayesSUR/BayesSUR_0.1.2.tar.gz",repos = NULL,type = "source")


#####################################################################################################
## Test the installation
library(BayesSUR)
data(example_data, package = "BayesSUR")
str(example_data)

# show the simulated gamma matrix and G_0
layout(matrix(1:2, ncol=2))
image(z=example_data$gamma, x=1:150, y=1:10, col=grey(1:0), xlab="SNPs Index", 
      ylab="Responses", main=mtext(bquote(True~" "~gamma)));box()
image(z=t(example_data$G), x=1:10, y=1:10, col=grey(1:0), xlab="Responses", 
      ylab="Responses", main="True graph of responses");box()

fit <- runSUR(data = example_data$data, Y = example_data$blockList[[1]],
              X = example_data$blockList[[2]], outFilePath = "results/", 
              nIter = 100000, nChains = 4, covariancePrior = "HIW", 
              gammaPrior = "hotspot")

str(fit)
#
# show the estimated beta, gamma and G_0
plotEstimator(fit);hat_gamma=getEstimator(fit);sum(hat_gamma>.5);max(hat_gamma)

# show the relationship of responses
plotResponseGraph(fit, PtrueResponse=example_data$G0, 
                  response.name=paste("GEX",1:ncol(example_data$G0),sep=""))

# show the network representation of the associations between responses and features
plotNetwork(fit, PmaxCovariate=0.5, lineup=1.2)
# show the manhattan plot
manhattan(fit)

# check the convergence of the algorithm
mcmcDiag(fit)
