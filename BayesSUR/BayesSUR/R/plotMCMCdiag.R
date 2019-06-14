#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title getEstimator
#' @description
#' Diagnose the convergence of the MCMC iterations
#' @name plotMCMCdiag
#' @param object fitted "runSUR" model
#' @param nbloc number of splits for the last half iterations after substracting burn-in length
#' @export
plotMCMCdiag <- function(object, nbloc=4){
  
  devAskNewPage(FALSE)
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  
  logP <- t( as.matrix( read.table(object$output$logP) ) )
  Ptau.indx <- ifelse(fit$input$covariancePrior!="IG", 7, 3)
  Plik.indx <- ifelse(fit$input$covariancePrior!="IG", 10, 5)
  #nChain <- object$input$nChains
  model_size <- rowSums(as.matrix( read.table(object$output$model_size) )) 
  
  dens.all <- density(logP[Ptau.indx, -c(1:2)])     #   ,main="log-marginal",col="red")
  dens.first <- density(logP[Ptau.indx, 2+1:floor(ncol(logP)/2)])#  ,main="",col="blue",add=TRUE)
  dens.last <- density(logP[Ptau.indx, -c(1:2,1:floor(ncol(logP)/2))])#  ,main="",col="black",add=TRUE)
  ymin <- min(dens.all$y,dens.first$y,dens.last$y)
  ymax <- max(dens.all$y,dens.first$y,dens.last$y)
  xmin <- min(dens.all$x,dens.first$x,dens.last$x)
  xmax <- max(dens.all$x,dens.first$x,dens.last$x)
  
  ###nsplit number of split of the sweep
  mid <- floor(floor((ncol(logP)-2)/2)/nbloc)
  ymax2 <- NULL
  xmin2 <- NULL
  xmax2 <-NULL
  list.dens <- NULL
  for (i in 1:nbloc){
    dens <- density(logP[Ptau.indx, 2+floor((ncol(logP)-2)/2)+(1+mid*(i-1)):(mid*nbloc)])
    ymax2 <- max(ymax2,dens$y)
    xmin2 <- min(xmin2,dens$x)
    xmax2 <- max(xmax2,dens$x)
    list.dens <- c(list.dens,list(dens))
  }
  
  ###plot the figures
  par(mfrow=c(2,2))
  
  plot(logP[Plik.indx,-c(1:2)], xlab="Iterations (*1000)", ylab="Log likelihood (posterior)", type="l", lty=1)
  #legend("topleft", legend=paste("Chain",1:nChain), lty=1, cex=0.5)
  
  plot(model_size[-c(1:2)], xlab="Iterations (*1000)", ylab="Model size", type="l", lty=1)
  #legend("topleft", legend=paste("Chain",1:nChain), lty=1, cex=0.5)
  
  title.0 <- expression(paste("Log Posterior Distribution: log ",P(gamma~group("|",list(Y,.),""))))
  plot(dens.all,main=title.0,col="black",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="")
  par(new=TRUE)
  plot(dens.first,main="",col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="")
  par(new=TRUE)
  plot(dens.last,main="",col="green",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="")
  legend("topleft",title="sweep",legend=paste(c("ALL","First half","Last half")," = [",c(1,1,(ncol(logP)-2)/2*1000+1),":",c(ncol(logP)-2,floor((ncol(logP)-2)/2),ncol(logP)-2)*1000,"]",sep=""),col=1:3,lty=1,text.col=1:3, cex=0.5)
  
  for (i in 1:nbloc){
    plot(list.dens[[i]],col=i,xlim=c(xmin2,xmax2),ylim=c(ymin,ymax2),xlab="",main=title.0)
    par(new=TRUE)  
  }
  legend("topleft",title="moving window",legend=paste("set ",1:nbloc," = [",(floor((ncol(logP)-2)/2)+mid*(nbloc:1-1))*1000+1,":",(ncol(logP)-2)*1000,"]",sep=""),col=1:nbloc,lty=1,text.col=1:nbloc, cex=0.5)
  
  par(mfrow=c(1,1))
  
  # we might also need other diagnostic plots here together, like temperature, etc.
  
}