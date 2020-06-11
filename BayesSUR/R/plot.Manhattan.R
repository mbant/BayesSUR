#' @title plot Manhattan-like plots for marginal posterior inclusion probabilities (mPIP) and numbers of responses of association for predictors
#' @description
#' Plot Manhattan-like plots for marginal posterior inclusion probabilities (mPIP) and numbers of responses of association for predictors of a "BayesSUR" class object.
#' @importFrom graphics axis box text par plot.default segments
#' @name plot.Manhattan
#' @param x an object of class \code{get.estimator} with \code{estimator="gamma"}
#' @param which if it's value "1" showing the Manhattan-like plot of the marginal posterior inclusion probabilities (mPIP). If it's value "2" showing the Manhattan-like plot of the number of responses. The default is to show both figures.
#' @param x.loc a vector of features distance
#' @param axis.label a vector of predictor names which are shown in the Manhattan-like plot. The value "NULL" only showing the indices. The default "auto" show the predictor names from the orginal data.
#' @param mark.responses a vector of response names which are shown in the Manhattan-like plot for the mPIP
#' @param mark.pos the location of the marked text relative to the point
#' @param xlab1 a title for the x axis of Manhattan-like plot for the mPIP
#' @param ylab1 a title for the y axis of Manhattan-like plot for the mPIP
#' @param xlab2 a title for the x axis of Manhattan-like plot for the numbers of responses 
#' @param ylab2 a title for the y axis of Manhattan-like plot for the numbers of responses 
#' @param threshold threshold for showing number of response variables significantly associated with each feature
#' @param las graphical parameter of plot.default
#' @param cex.axis graphical parameter of plot.default
#' @param mark.color the color of the marked text. The default color is red.
#' @param mark.cex the fontsize of the marked text. The default fontsize is 0.8.
#' @param header the main title
#' @param ... other arguments
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
#' # show the Manhattan-like plots
#' gamma <- get.estimator(fit, estimator="gamma")
#' plot(gamma)
#' 
#' @export
plot.Manhattan <- function(x, which=c(1,2), x.loc=FALSE, axis.label="auto", mark.responses=NULL, xlab1="Predictors", ylab1="mPIP", xlab2="Predictors", ylab2="No. of responses",
                          threshold=0.5,las=0, cex.axis=1, mark.pos=c(0,0), mark.color=2, mark.cex=0.8, header="", ...){
  
  gamma <- x
  if(is.null(axis.label)){
    x.loc <- 1:nrow(gamma)
    names(x.loc) <- 1:nrow(gamma)
  }else{
    #name.predictors <- colnames(read.table(x$output$X,header=T))
    name.predictors <- rownames(gamma)
    if(axis.label[1] == "auto"){
      x.loc <- 1:nrow(gamma)
      names(x.loc) <- name.predictors
    }else{
      if( (!match(axis.label, name.predictors)[1]) & (!x.loc[1]) )
        stop("The given predictor names are not consistent with the data")
      
      if(!x.loc[1]) x.loc <- match( axis.label, name.predictors )
      names(x.loc) <- axis.label
    }
  }
  
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))    
  if(sum(which==c(1,2)) == 2){
    par(mfrow=c(2,1))
  }else{
    par(mfrow=c(1,1))
  }
    
  # Manhattan plot for marginal posterior inclusion probabilities (mPIP) 
  if(1 %in% which){
    par(mar=c(4,4,4,2)) 
  
  plot.default(as.vector(gamma) ~ rep(1:nrow(gamma), times=ncol(gamma)), xlim=c(1,nrow(gamma)), ylim=c(0,max(gamma)), xaxt = 'n',bty = "n", ylab=ylab1, xlab=xlab1, main="", pch=19, ...)
  axis(1, at=x.loc, labels=names(x.loc), las=las, cex.axis=cex.axis); box()
  
  # mark the names of the specified response variables corresponding to the given predictors
  if(!is.null(mark.responses)){
    name.responses <- colnames(read.table(x$output$Y,header=T))
    if(!is.na(match(mark.responses, name.responses)[1])){
      text(rep(x.loc,times=length(mark.responses))+mark.pos[1], as.vector(gamma[x.loc,name.responses %in% mark.responses])+mark.pos[2], labels=rep(mark.responses, each=length(x.loc)), col=mark.color, cex=mark.cex)
    }else{
      stop("The given response names are not consistent with the data")
    }
  }
  }
  
  # Manhattan plot for numbers of responses 
  if(2 %in% which){
    par(mar=c(5,4,3,2))
  no.gamma <- rowSums(gamma>=threshold)
  plot.default(no.gamma ~ c(1:nrow(gamma)), xlim=c(1,nrow(gamma)), ylim=c(0,max(no.gamma)+0.3), type='n', xaxt = 'n', ylab=ylab2, xlab=xlab2, main="", ...)
  segments(1:nrow(gamma), 0, 1:nrow(gamma), no.gamma)
  axis(1, at=x.loc, labels=names(x.loc), las=las, cex.axis=cex.axis)
  }

  title(paste("\n\n",header,sep=""), outer=T)
  
}
