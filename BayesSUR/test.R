# ## Build a new version of the package
# remove.packages("BayesSUR")
# Rcpp::compileAttributes(pkgdir = "/Users/zhiz/Downloads/BayesSUR/BayesSUR/"); devtools::document("/Users/zhiz/Downloads/BayesSUR/BayesSUR")
# devtools::build("/Users/zhiz/Downloads/BayesSUR/BayesSUR")#,vignettes=TRUE)

## Install the package
install.packages("/Users/zhiz/Downloads/BayesSUR/BayesSUR_0.1.6.tar.gz",repos = NULL,type = "source")


#####################################################################################################
## load the package
library(BayesSUR)
data(example_eQTL, package = "BayesSUR")
str(example_eQTL)

# show the simulated gamma matrix and G_0
attach(example_eQTL)
layout(matrix(1:2, ncol=2))
image(z=gamma, x=1:150, y=1:10, col=grey(1:0), xlab="SNPs Index", 
      ylab="Responses", main=mtext(bquote(True~" "~gamma)));box()
image(z=t(G0), x=1:10, y=1:10, col=grey(1:0), xlab="Responses", 
      ylab="Responses", main="True graph of responses");box()

# fit a SSUR model with hotspot prior
fit <- runSUR(data = data, Y = blockList[[1]],
              X = blockList[[2]], outFilePath = "results/", 
              nIter = 2000, nChains = 5, covariancePrior = "HIW", burnin=1000,
              gammaPrior = "hotspot")

str(summary(fit))

# show the interaction of plots
plot(fit)

# show the estimated beta, gamma and G_0
plotEstimator(fit);

# show the relationship of responses
plotResponseGraph(fit, PtrueResponse=G0, response.name=paste("GEX",1:ncol(G0),sep=""))

# show the network representation of the associations between responses and features
plotNetwork(fit,label.predictor = NA,lineup=1.5,nodesizePredictor=2,nodesizeResponse=15,
            name.predictors="SNPs", name.responses="Gene expression",edge.weight=TRUE)

# show the manhattan plot
plotManhattan(fit)

# check the convergence of the algorithm
plotMCMCdiag(fit)
