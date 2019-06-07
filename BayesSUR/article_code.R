library(BayesSUR)
library(igraph)
data(example_data, package = "BayesSUR")
str(example_data)
# show the simulated gamma matrix and G_0
layout(matrix(1:2, ncol=2))
image(z=example_data$gamma, x=1:150, y=1:10, col=grey(1:0), xlab="SNPs Index", 
      ylab="Responses", main=mtext(bquote(True~" "~gamma)));box()
image(z=t(example_data$G), x=1:10, y=1:10, col=grey(1:0), xlab="Responses", 
      ylab="Responses", main="True graph of responses");box()
 
#hyperpar = list(mrf_e=-3, mrf_d=3/10, b_pi = 0.2 , a_pi = 0.1 , b_o = 5 , a_o = 1 )
#hyperpar = list(mrf_e=-3, mrf_d=3/10, b_pi = 0.02, a_pi = 0.01, b_o = 5 , a_o = 1 )
#hyperpar = list(mrf_e=-3, mrf_d=3/10, b_pi = 100, a_pi = 200 , b_o = 8 , a_o = 7 )
#(mean_pi=hyperpar$a_pi/hyperpar$b_pi); hyperpar$a_pi/hyperpar$b_pi^2
#(mean_o=hyperpar$a_o/(hyperpar$a_o+hyperpar$b_o)); hyperpar$a_o*hyperpar$b_o/(hyperpar$a_o+hyperpar$b_o)^2/(hyperpar$a_o+hyperpar$b_o+1)
#mean_pi * mean_o
#sum(example_data$gamma)/prod(dim(example_data$gamma)); sum(example_data$gamma)

fit <- runSUR(data = example_data$data, Y = example_data$blockList[[1]],
              X = example_data$blockList[[2]], outFilePath = "results/", 
              nIter = 1000, nChains = 4, covariancePrior = "HIW", 
              gammaPrior = "hotspot")
str(fit)
# show the estimated beta, gamma and G_0
plotEstimator(fit);hat_gamma=getEstimator(fit);sum(hat_gamma>.5);max(hat_gamma)
# show the relationship of responses
#pdf("ResponseGraph.pdf", height=4, width=7)
plotResponseGraph(fit, PtrueResponse=example_data$G0, 
                  response.name=paste("GEX",1:ncol(example_data$G0),sep=""))
#dev.off()
# show the network representation of the associations between responses and features
plotNetwork(fit, PmaxPredictor=0.8, lineup=1.2)
# show the manhattan plot
manhattan(fit)

# check the convergence of the algorithm
mcmcDiag(fit)


#=============
# Analyze the GDSC database
#=============
load("example_GDSC.rda")
fit <- runSUR(data = example_GDSC$data,
              Y = example_GDSC$blockList[[1]],
              X = example_GDSC$blockList[[3]],
              X_0 = example_GDSC$blockList[[2]],
              outFilePath = "results/", nIter = 2000, nChains = 2, 
              covariancePrior = "HIW", gammaPrior = "MRF", mrfG=example_data$mrfG)
# show the estimated beta, gamma and G_0
hat.gamma <- getEstimator(fit, "gamma")
plotEstimator(fit)
plotResponseGraph(fit)
plotNetwork(fit)
mcmcDiag(fit)
manhattan(fit)

pdf("ResponseEstimator.pdf", height=8, width=8);plotEstimator(fit);dev.off()
# show the relationship of responses
pdf("ResponseNetwork0.pdf", height=8, width=8);plotResponseGraph(fit);dev.off()
# show the network representation of the associations between responses and features
pdf("ResponseNetwork1.pdf", height=8, width=10);plotNetwork(fit, PmaxCovariate=0.5);dev.off()
pdf("ResponseNetwork2.pdf", height=8, width=10);plotNetwork(fit, PmaxCovariate=0.6);dev.off()
pdf("ResponseNetwork3.pdf", height=8, width=10);plotNetwork(fit, PmaxCovariate=0.7);dev.off()
# show the manhattan plot
#manhattan(fit)

# check the convergence of the algorithm
pdf("ResponseMCMC.pdf", height=8, width=8);mcmcDiag(fit);dev.off()
