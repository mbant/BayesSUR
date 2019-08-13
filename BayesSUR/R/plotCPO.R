
#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title plotCPO
#' @description
#' Plot the conditional predictive ordinate (CPO) which is the leave-one-out cross-validation predictive density. The CPO is a handy posterior predictive check because it may be used to identify outliers, influential observations, and for hypothesis testing across different non-nested models.
#' @param object the object from the runSUR
#' @param xlab a title for the x axis 
#' @param ylab a title for the y axis 
#' @param outlier.thresh threshold for the CPOs. The default is 0.01.
#' @param outlier.mark mark the outliers with the response names. The default is \code{FALSE}
#' @param scale.CPO scaled CPOs which is divided by their maximum. The default is \code{TRUE}
#' @param x.loc a vector of features distance
#' @param axis.label a vector of predictor names which are shown in CPO plot. The default is \code{NULL} only showing the indices. The value \code{"auto"} show the predictor names from the orginal data.
#' @param mark.pos the location of the marked text relative to the point
#' @param las graphical parameter of plot.default
#' @param cex.axis graphical parameter of plot.default
#' @param mark.color the color of the marked text. The default color is red.
#' @param mark.cex the fontsize of the marked text. The default fontsize is 0.8.
#' @export
plotCPO <- function(object, outlier.mark=TRUE, outlier.thresh=0.01, scale.CPO=TRUE, x.loc=FALSE, axis.label=NULL, las=0, cex.axis=1, mark.pos=c(0,-.01), mark.color=2, mark.cex=0.8){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  CPO <- as.matrix( read.table(object$output$CPO) )
  if(scale.CPO) CPO <- CPO/max(CPO)
    
  if(is.null(axis.label)){
    x.loc <- 1:nrow(CPO)
    names(x.loc) <- 1:nrow(CPO)
  }else{
    if(axis.label[1] == "auto"){
      x.loc <- 1:nrow(CPO)
      names(x.loc) <- rownames(read.table(object$output$Y,header=T))
    }else{
      if(!x.loc[1]){
        x.loc <- 1:length(axis.label)
      }else{
        if( length(axis.label) != length(x.loc) )
          stop("The given predictor names are not consistent with the data")
      }
      names(x.loc) <- axis.label
    }
  }
  
    plot(as.vector(CPO) ~ rep(1:nrow(CPO), times=ncol(CPO)), xlim=c(1,nrow(CPO)), ylim=c(0,max(CPO)), xaxt = 'n',bty = "n", ylab = ifelse(scale.CPO,"scaled CPOs","CPOs"), xlab = "", main="Conditional predictive ordinate", pch=19)
    axis(1, at=x.loc, labels=names(x.loc), las=las, cex.axis=cex.axis); box()
    
    # mark the names of the specified response variables corresponding to the given responses
    if(outlier.mark){
      if(min(CPO) > outlier.thresh)
        stop("The minimum CPO is larger than the outlier.thresh!")
      name.responses <- colnames(read.table(object$output$Y,header=T))
      text(rep(1:nrow(CPO), times=ncol(CPO))[which(as.vector(CPO) <= outlier.thresh)]+mark.pos[1], as.vector(CPO[CPO<outlier.thresh])+mark.pos[2], labels=rep(name.responses, each=nrow(CPO))[as.vector(CPO) < outlier.thresh], col=mark.color, cex=mark.cex)
      abline(h=outlier.thresh, lty=2, col=mark.color)
    }
}
