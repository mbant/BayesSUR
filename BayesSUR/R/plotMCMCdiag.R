#' @title show trace plots and diagnostic density plots
#' @description
#' Show trace plots and diagnostic density plots of a fitted model object of class "BayesSUR".
#' @importFrom graphics par plot.default legend title
#' @importFrom stats density
#' @name plotMCMCdiag
#' @param object an object of class "BayesSUR"
#' @param nbloc number of splits for the last half iterations after substracting burn-in length
#' @param header the main title
#' @param ... other arguments for the plots of the log-likelihood and model size
#' 
#' @examples
#' data("example_eQTL", package = "BayesSUR")
#' hyperpar <- list( a_w = 2 , b_w = 5 )
#' 
#' fit <- BayesSUR(Y = example_eQTL[["blockList"]][[1]], 
#'                 X = example_eQTL[["blockList"]][[2]],
#'                 data = example_eQTL[["data"]], outFilePath = tempdir(),
#'                 nIter = 100, burnin = 0, nChains = 2, gammaPrior = "hotspot",
#'                 hyperpar = hyperpar, tmpFolder = "tmp/" )
#' 
#' ## check output
#' # show the diagnosis plots with at least 4000 iterations
#' \donttest{
#' plotMCMCdiag(fit)
#' }
#' 
#' @export
plotMCMCdiag <- function(object, nbloc=3, header="", ...){
  
  if(object$input$nIter < 4000)
    stop("The argument nIter should be at least 4000 for the diagnostic plots!")
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  
  nIter <- object$input$nIter
  logP <- t( as.matrix( read.table(object$output$logP) ) )
  logP <- logP[,ncol(logP)-floor(nIter/1000)-1+1:floor(nIter/1000)]
  Ptau.indx <- ifelse(object$input$covariancePrior!="IG", 7, 3)
  Plik.indx <- ifelse(object$input$covariancePrior!="IG", 10, 5)
  #nChain <- object$input$nChains
  model_size <- as.matrix( read.table(object$output$model_size) )
  model_size <- rowSums(  model_size[nrow(model_size)-floor(nIter/1000)-1+1:floor(nIter/1000),] ) 
  
  dens.all <- density(logP[Ptau.indx,])     #   ,main="log-marginal",col="red")
  dens.first <- density(logP[Ptau.indx, 1:floor(ncol(logP)/2)])#  ,main="",col="blue",add=TRUE)
  dens.last <- density(logP[Ptau.indx, (1+floor(ncol(logP)/2)):ncol(logP)])#  ,main="",col="black",add=TRUE)
  ymin <- min(dens.all$y,dens.first$y,dens.last$y)
  ymax <- max(dens.all$y,dens.first$y,dens.last$y)
  xmin <- min(dens.all$x,dens.first$x,dens.last$x)
  xmax <- max(dens.all$x,dens.first$x,dens.last$x)
  
  ###nsplit number of split of the sweep
  mid <- floor(floor(ncol(logP)/2)/nbloc)
  ymax2 <- NULL
  xmin2 <- NULL
  xmax2 <-NULL
  list.dens <- NULL
  for (i in 1:nbloc){
    dens <- density(logP[Ptau.indx, (floor(ncol(logP)/2)+1+mid*(i-1)):ncol(logP)])
    ymax2 <- max(ymax2,dens$y)
    xmin2 <- min(xmin2,dens$x)
    xmax2 <- max(xmax2,dens$x)
    list.dens <- c(list.dens,list(dens))
  }
  
  ###plot the figures
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))    
  par(mfrow=c(2,2))
  
  plot.default(logP[Plik.indx,], xlab="Iterations (*1000)", ylab="Log likelihood (posterior)", type="l", lty=1, ...)
  #legend("topleft", legend=paste("Chain",1:nChain), lty=1, cex=0.5)
  
  plot.default(model_size, xlab="Iterations (*1000)", ylab="Model size", type="l", lty=1, ...)
  #legend("topleft", legend=paste("Chain",1:nChain), lty=1, cex=0.5)
  
  title.0 <- expression(paste("Log Posterior Distribution: log ",P(gamma~group("|",list(Y,.),""))))
  plot.default(dens.all,main=title.0,col="black",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="", type="l", lty=1)
  par(new=TRUE) 
  plot.default(dens.first,main="",col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="", type="l", lty=1)
  par(new=TRUE) 
  plot.default(dens.last,main="",col="green",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",,ylab="Density", type="l", lty=1)
  legend("topleft",title="iteration",legend=paste(c("ALL","First half","Last half")," = [",c(1,1,floor(ncol(logP)/2)*1000+1),":",c(ncol(logP),floor((ncol(logP))/2),ncol(logP))*1000,"]",sep=""),col=1:3,lty=1,text.col=1:3, cex=0.8)
  
  for (i in 1:nbloc){
    plot.default(list.dens[[i]],col=i,xlim=c(xmin2,xmax2),ylim=c(ymin,ymax2),xlab="",ylab="", type="l", lty=1,main=title.0)
    par(new=TRUE) 
  }
  title(ylab="Density")
  legend("topleft",title="moving window",legend=paste("set ",1:nbloc," = [",(floor((ncol(logP))/2)+mid*(nbloc:1-1))*1000+1,":",(ncol(logP))*1000,"]",sep=""),col=1:nbloc,lty=1,text.col=1:nbloc, cex=0.8)
  
  title(paste("\n",header,sep=""), outer=T)
  
}