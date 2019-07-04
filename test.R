### Build a new version of the package
#remove.packages("BayesSUR")
#Rcpp::compileAttributes(pkgdir = "/Users/zhiz/Downloads/BayesSUR/BayesSUR/"); devtools::document("/Users/zhiz/Downloads/BayesSUR/BayesSUR")
#devtools::build("/Users/zhiz/Downloads/BayesSUR/BayesSUR")#,vignettes=TRUE)

## Install the package
#library(devtools)
#install_github("mbant/BayesSUR/BayesSUR")
install.packages("BayesSUR_0.1.14.tar.gz",repos = NULL,type = "source")


#####################################################################################################
## load the package
library(BayesSUR)
data(example_eQTL, package = "BayesSUR")
str(example_eQTL)

# show the simulated gamma matrix and G_y
attach(example_eQTL)
options(tikzMetricPackages = c("\\usepackage{amsmath}","\\usepackage{bm}", "\\usetikzlibrary{calc}"))
tikz('ParamTrue.tex',width=5.5,height=3, standAlone = TRUE,packages = c("\\usepackage{tikz}","\\usepackage{amsmath}","\\usepackage{bm}"))
layout(matrix(1:2, ncol=2))
image(z=gamma, x=1:150, y=1:10, col=grey(1:0), xlab="SNPs Index", 
      ylab="Responses", main=paste("True","$\\bm{\\gamma}$"));box()
image(z=t(Gy), x=1:10, y=1:10, col=grey(1:0), xlab="Responses", 
      ylab="Responses", main=paste("True","$\\mathcal{G}$"));box()
dev.off()
tools::texi2pdf("ParamTrue.tex")
system(paste(getOption("pdfviewer"), "ParamTrue.pdf"))

# fit a SSUR model with hotspot prior
fit <- runSUR(data = data, Y = blockList[[1]],
              X = blockList[[2]], outFilePath = "results/", 
              nIter = 2000, nChains = 5, covariancePrior = "HIW", burnin=1000,
              gammaPrior = "hotspot")

str(summary(fit))

# show the interaction of plots
plot(fit)

# show the estimated beta, gamma and G_y
plotEstimator(fit, fig.tex=TRUE)
system(paste(getOption("pdfviewer"), "ParamEstimator.pdf"))

# show the relationship of responses
plotResponseGraph(fit, PtrueResponse=Gy, response.name=paste("GEX",1:ncol(Gy),sep=""))

# show the network representation of the associations between responses and features
plotNetwork(fit,label.predictor = NA,lineup=1.5,nodesizePredictor=2,nodesizeResponse=15,
            name.predictors="SNPs", name.responses="Gene expression",edge.weight=TRUE)

# show the manhattan plot
plotManhattan(fit)

# check the convergence of the algorithm
plotMCMCdiag(fit)
