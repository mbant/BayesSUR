## Build a new version of the package
remove.packages("BayesSUR")
Rcpp::compileAttributes(pkgdir = "/Users/zhiz/Downloads/BayesSUR/BayesSUR/"); devtools::document("/Users/zhiz/Downloads/BayesSUR/BayesSUR")
devtools::build("/Users/zhiz/Downloads/BayesSUR/BayesSUR")#,vignettes=TRUE)

## Install the package
install.packages("/Users/zhiz/Downloads/BayesSUR/BayesSUR_0.1.6.tar.gz",repos = NULL,type = "source")


#####################################################################################################
## Test the installation
library(BayesSUR)
data(example_data, package = "BayesSUR")
str(example_data)

# show the simulated gamma matrix and the response graph \mathcal{G}
attach(example_data)
layout(matrix(1:2, ncol=2))
image(z=gamma, x=1:150, y=1:10, col=grey(1:0), xlab="SNPs Index", 
      ylab="Responses", main=mtext(bquote(True~" "~gamma)));box()
image(z=t(G), x=1:10, y=1:10, col=grey(1:0), xlab="Responses", 
      ylab="Responses", main="True graph of responses");box()

fit <- runSUR(data = data, Y = blockList[[1]],
              X = blockList[[2]], outFilePath = "results/", 
              nIter = 3000, nChains = 5, covariancePrior = "HIW", 
              gammaPrior = "hotspot")

str(summary(fit))
# show the interaction of plots
plot(fit)
# show the estimated beta, gamma and G_0
plotEstimator(fit);hat_gamma=getEstimator(fit);sum(hat_gamma>.5);max(hat_gamma)
# show the relationship of responses
#pdf("ResponseGraph.pdf", height=4, width=7)
plotResponseGraph(fit, PtrueResponse=G0, response.name=paste("GEX",1:ncol(G0),sep=""))
#dev.off()
# show the network representation of the associations between responses and features
#pdf("ResponseNetwork.pdf", height=10, width=10)
plotNetwork(fit,label.predictor = NA,lineup=1.5,nodesizePredictor=2,nodesizeResponse=15,
            name.predictors="SNPs", name.responses="Gene expression",edge.weight=TRUE)
#dev.off()
# show the manhattan plot
plotManhattan(fit)

# check the convergence of the algorithm
plotMCMCdiag(fit)
