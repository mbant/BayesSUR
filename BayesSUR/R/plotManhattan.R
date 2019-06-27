
#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title plotManhattan
#' @description
#' Manhattan plot
#' @param object the object from the runSUR
#' @param x.loc a vector of features distance
#' @param show.marker a vector of covariate names which are shown in the Manhattan plot
#' @param xlab1 a title for the x axis of Manhattan plot for marginal posterior inclusion probabilities
#' @param ylab1 a title for the y axis of Manhattan plot for marginal posterior inclusion probabilities
#' @param xlab2 a title for the x axis of Manhattan plot for numbers of responses 
#' @param ylab2 a title for the y axis of Manhattan plot for numbers of responses 
#' @param threshold threshold for showing number of response variables significantly associated with each feature
#' @param show.all.xlab logical value for showing all labels on the x axis. The defaul value is to show 5 labels on the x axis
#' @param las graphical parameter of plot.default
#' @param cex.axis graphical parameter of plot.default
#' @export
plotManhattan <- function(object, x.loc=FALSE, show.marker=NULL, xlab1="", ylab1="mPIP", xlab2="features", ylab2="No. of responses",threshold=0.5,
                           show.all.xlab=FALSE, las=0, cex.axis=1){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  gamma <- as.matrix( read.table(object$output$gamma) )
  
  if(x.loc){
    x.loc <- 1:dim(gamma)[1]
    names(x.loc) <- 1:dim(gamma)[1]
  }else{
    x.loc <- 1:dim(gamma)[1]
    if(is.null(names(x.loc))) names(x.loc) <- colnames(read.table(object$output$X,header=T))
  }
  if(show.all.xlab){
    n.xlab <- length(x.loc)
  }else{
    n.xlab <- 5
  } 
    
  par(mfrow=c(2,1))
  # Manhattan plot for marginal posterior inclusion probabilities (mPIP) 
  par(mar=c(1,4,6.5,2))
  for(i in 1:(dim(gamma)[2]-1)){
    plot(gamma[,i]~x.loc, xlim=c(min(x.loc),max(x.loc)), ylim=c(0,1), xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', pch=19)
    par(new=T)
  }
  plot(gamma[,dim(gamma)[2]]~x.loc, xlim=c(min(x.loc),max(x.loc)), ylim=c(0,1), xaxt = 'n', ylab=ylab1, xlab=xlab1, main="\n\nManhattan Plot", pch=19)
  axis(1, at=x.loc[seq(1,max(x.loc),length=n.xlab)], labels=names(x.loc)[seq(1,max(x.loc),length=n.xlab)], las=las, cex.axis=cex.axis)
  
  # Manhattan plot for numbers of responses 
  par(mar=c(5,4,3,2))
  no.gamma <- rowSums(gamma>=threshold)
  plot(no.gamma~x.loc, xlim=c(min(x.loc),max(x.loc)), ylim=c(0,max(no.gamma)+0.3), type='n', xaxt = 'n', ylab=ylab2, xlab=xlab2, main="")
  axis(1, at=x.loc[seq(1,max(x.loc),length=n.xlab)], labels=names(x.loc)[seq(1,max(x.loc),length=n.xlab)], las=las, cex.axis=cex.axis)
  segments(x.loc, 0, x.loc, no.gamma)
  if(!is.null(show.marker)) text(no.gamma[names(x.loc) %in% show.marker]+0.2~x.loc[names(x.loc) %in% show.marker], labels=names(x.loc)[names(x.loc) %in% show.marker])
  
  par(mfrow=c(1,1))
  
}
