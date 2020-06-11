#' @title show trace plots and diagnostic density plots
#' @description
#' Show trace plots and diagnostic density plots of a fitted model object of class "BayesSUR".
#' @importFrom graphics par plot.default legend title matplot
#' @importFrom stats density
#' @importFrom grDevices hcl.colors
#' @name plot.MCMCdiag
#' @param x an object of class \code{get.estimator} with \code{estimator="logP"}
#' @param nbloc number of splits for the last half iterations after substracting burn-in length
#' @param HIWg diagnostic plot of the response graph. Default is \code{NULL}. \code{HIW="degree"} prints the diagnostic of the degrees of response nodes. \code{HIW="edges"} prints the diagnostic 
#' of every edge between two responses. \code{HIW="lik"} prints the diagnostic of the posterior likelihoods of the hyperparameters related to the response relationships
#' @param header the main title
#' @param ... other arguments for the plots of the log-likelihood and model size
#' 
#' @examples
#' data("example_eQTL", package = "BayesSUR")
#' hyperpar <- list( a_w = 2 , b_w = 5 )
#' 
#' set.seed(9173)
#' fit <- BayesSUR(Y = example_eQTL[["blockList"]][[1]], 
#'                 X = example_eQTL[["blockList"]][[2]],
#'                 data = example_eQTL[["data"]], outFilePath = tempdir(),
#'                 nIter = 100, burnin = 50, nChains = 2, gammaPrior = "hotspot",
#'                 hyperpar = hyperpar, tmpFolder = "tmp/" )
#' 
#' ## check output
#' MCMCdiag <- get.estimator(fit, estimator = "logP")
#' plot(MCMCdiag)
#' 
#' @export
plot.MCMCdiag <- function(x, nbloc=3, HIWg=NULL, header="", ...){
  
  if(x$nIter <= 1)
    stop("The diagosis only shows results from more than one MCMC iteration!")
  if(x$nIter < 4000)
    message("NOTE: The diagosis only shows results of two iteration points due to less than 4000 MCMC iterations!")
  
  nIter <- x$nIter
  logP <- x$logP
  ncol_Y <- x$ncol_Y
  
  if(nIter >= 4000){
    logP <- logP[,ncol(logP)-floor(nIter/1000)-1+1:floor(nIter/1000)]
  }else{
    logP <- logP[,c(1,ncol(logP))]
  }
  if(is.null(HIWg)){
    Ptau.indx <- ifelse(x$covariancePrior!="IG", 7, 3)
    Plik.indx <- ifelse(x$covariancePrior!="IG", 10, 5)
    #nChain <- x$input$nChains
    model_size <-x$model_size
    if(nIter >= 4000){
      model_size <- rowSums(  model_size[nrow(model_size)-floor(nIter/1000)-1+1:floor(nIter/1000),] ) 
    }else{
      model_size <- rowSums(  model_size[c(1,nrow(model_size)),] ) 
    }
    
    dens.all <- density(logP[Ptau.indx,])     #   ,main="log-marginal",col="red")
    if(nIter >= 4000){
      dens.first <- density(logP[Ptau.indx, 1:floor(ncol(logP)/2)])#  ,main="",col="blue",add=TRUE)
      dens.last <- density(logP[Ptau.indx, (1+floor(ncol(logP)/2)):ncol(logP)])#  ,main="",col="black",add=TRUE)
      
      ymin <- min(dens.all$y,dens.first$y,dens.last$y)
      ymax <- max(dens.all$y,dens.first$y,dens.last$y)
      xmin <- min(dens.all$x,dens.first$x,dens.last$x)
      xmax <- max(dens.all$x,dens.first$x,dens.last$x)
    }else{
      ymin <- min(dens.all$y)
      ymax <- max(dens.all$y)
      xmin <- min(dens.all$x)
      xmax <- max(dens.all$x)
      nbloc <- 1
    }
    
    ###nsplit number of split of the sweep
    mid <- floor(floor(ncol(logP)/2)/nbloc)
    ymax2 <- xmin2 <- xmax2 <- list.dens <- NULL
    for (i in 1:nbloc){
      dens <- density(logP[Ptau.indx, (ifelse(nbloc==1,0,floor(ncol(logP)/2))+1+mid*(i-1)):ncol(logP)])
      ymax2 <- max(ymax2,dens$y)
      xmin2 <- min(xmin2,dens$x)
      xmax2 <- max(xmax2,dens$x)
      list.dens <- c(list.dens,list(dens))
    }
    
    ###plot the figures
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))    
    par(mfrow=c(2,2))
    
    if(nbloc>1){
      plot.default(logP[Plik.indx,], xlab="Iterations (*1000)", ylab="Log likelihood (posterior)", type="l", lty=1, ...)
    }else{
      plot.default(logP[Plik.indx,]~c(1,nIter), xlab="Iterations", ylab="Log likelihood (posterior)", type="l", lty=1, ...)
    }
    #legend("topleft", legend=paste("Chain",1:nChain), lty=1, cex=0.5)
    
    if(nbloc>1){
      plot.default(model_size, xlab="Iterations (*1000)", ylab="Model size", type="l", lty=1, ...)
    }else{
      plot.default(model_size~c(1,nIter), xlab="Iterations", ylab="Model size", type="l", lty=1, ...)
    }
    #legend("topleft", legend=paste("Chain",1:nChain), lty=1, cex=0.5)
    
    title.0 <- expression(paste("Log Posterior Distribution: log ",P(gamma~group("|",list(Y,.),""))))
    #plot.default(dens.all,main=title.0,col="black",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="", type="l", lty=1)
    plot.default(dens.all,main="",col="black",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=title.0,ylab="", type="l", lty=1)
    if(nIter >= 4000){ 
      par(new=TRUE) 
      plot.default(dens.first,main="",col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="", type="l", lty=1)
      par(new=TRUE) 
      plot.default(dens.last,main="",col="green",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",,ylab="Density", type="l", lty=1)
      if(nbloc>1) 
        legend("topleft",title="iteration",legend=paste(c("ALL","First half","Last half")," = [",c(1,1,floor(ncol(logP)/2)*1000+1),":",c(ncol(logP),floor((ncol(logP))/2),ncol(logP))*1000,"]",sep=""),col=1:3,lty=1,text.col=1:3, cex=0.8)
    }
    for (i in 1:nbloc){
      #plot.default(list.dens[[i]],col=i,xlim=c(xmin2,xmax2),ylim=c(ymin,ymax2),xlab="",ylab="", type="l", lty=1,main=title.0)
      plot.default(list.dens[[i]],col=i,xlim=c(xmin2,xmax2),ylim=c(ymin,ymax2),xlab=title.0,ylab="", type="l", lty=1,main="")
      if(nbloc>1) par(new=TRUE) 
    }
    title(ylab="Density")
    if(nbloc>1) 
      legend("topleft",title="moving window",legend=paste("set ",1:nbloc," = [",(floor((ncol(logP))/2)+mid*(nbloc:1-1))*1000+1,":",(ncol(logP))*1000,"]",sep=""),col=1:nbloc,lty=1,text.col=1:nbloc, cex=0.8)
  }else{
    if(x$covariancePrior != "HIW")
      stop("The argument HIWg only works for the model with hyper-inverse Wishart prior on the covariance!")
    
    if(HIWg == "degree"){
      Gvisit <- x$Gvisit
      
      m <- ncol_Y
      node1 <- node2 <- NULL
      for(i in 1:(m-1)){
        node1 <- c(node1, rep(i, m-i))
        node2 <- c(node2, (i+1):m)
      }
      nodes <- cbind(node1, node2)
      node.degree <- matrix(0, nrow=nrow(Gvisit), ncol=m)
      for(i in 1:m)
        node.degree[,i] <- rowSums(Gvisit[,which(nodes==i, arr.ind=TRUE)[,1]])
      
      matplot(node.degree, type="l", lty=1, col=hcl.colors(m), xlab="Iterations (*1000)", ylab="degree", main="Response degrees", xlim=c(1,nrow(Gvisit)*1.1))
      legend("topright",legend=1:m,col=hcl.colors(m),lty=1,text.col=hcl.colors(m), cex=1/m*4)
    }
    
    if(substr(HIWg, 1, 4) == "edge"){
      Gvisit <- x$Gvisit
      
      m <- ncol_Y
      node1 <- node2 <- NULL
      for(i in 1:(m-1)){
        node1 <- c(node1, rep(i, m-i))
        node2 <- c(node2, (i+1):m)
      }
      
      if(HIWg == "edge"){
        matplot(Gvisit, type="l", lty=1, col=hcl.colors(ncol(Gvisit)), xlab="Iterations (*1000)", ylab="", main="Edges selection", xlim=c(1,nrow(Gvisit)*1.1))
        legend("topright",legend=paste(node1,"-",node2,sep=""),col=hcl.colors(ncol(Gvisit)),lty=1,text.col=hcl.colors(ncol(Gvisit)), cex=1/m*2)
      }else{
        plot.default(Gvisit[,which(paste(node1,node2,sep="")==substr(HIWg,5,nchar(HIWg)))], type="l", lty=1, xlab="Iterations (*1000)", ylab="", main=paste("Edge-",substr(HIWg,5,nchar(HIWg))," selection",sep=""))
      }
    }
    
    if(HIWg == "lik"){
      Gvisit <- t(logP[1:4,])
      matplot(Gvisit, type="l", lty=1, col=1:ncol(Gvisit), xlab="Iterations (*1000)", ylab="Log likelihood (posterior)", main="Likelihoods of graph learning")
      legend("topright",legend=c("tau", "eta", "JT", "SigmaRho"),col=1:ncol(Gvisit),lty=1,text.col=1:ncol(Gvisit), cex=0.8)
    }
  }
  title(paste("\n",header,sep=""), outer=T)
  
}