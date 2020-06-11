#' @title plot the posterior mean estimators
#' @description
#' Plot the posterior mean estimators from a "BayesSUR" class object, including the coefficients beta, latent indicator variable gamma and graph of responses.
#' @importFrom graphics axis box text mtext par image
#' @importFrom grDevices colorRampPalette dev.off grey
#' @importFrom tikzDevice tikz
#' @name plot.estimator
#' @param x an object of class \code{get.estimator} with \code{estimator=c("beta","gamma","Gy")}
#' @param estimator print the heatmap of estimators. The value "beta" is for the estimated coefficients matrix, "gamma" for the latent indicator matrix and "Gy" for the graph of responses
#' @param colorScale.gamma value palette for gamma
#' @param colorScale.beta a vector of three colors for diverging color schemes
#' @param legend.cex.axis magnification of axis annotation relative to cex
#' @param name.responses a vector of the response names. The default is "NA" only to show the locations. The value "auto" show the response names from the orginal data. 
#' @param name.predictors a vector of the predictor names. The default is "NA" only to show the locations. The value "auto" show the predictor names from the orginal data.
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param fig.tex print the figure through LaTex. Default is "FALSE"
#' @param output the file name of printed figure
#' @param header the main title
#' @param header.cex size of the main title for all estimators
#' @param cex.main size of the title for each estimator
#' @param title.beta a title for the printed "beta" if \code{fig.tex=TRUE}
#' @param title.gamma a title for the printed "gamma" if \code{fig.tex=TRUE}
#' @param title.Gy a title for the printed "Gy" if \code{fig.tex=TRUE}
#' @param tick a logical value specifying whether tickmarks and an axis line should be drawn. Default is "FALSE"
#' @param mgp the margin line (in mex units) for the axis title, axis labels and axis line
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
#' # Plot the estimators from the fitted object
#' \dontrun{
#' estimators <- get.estimator(fit, estimator = c("beta","gamma","Gy"))
#' plot(estimators)
#' 
#' #Set up temporary work directory for saving a pdf figure
#' td <- tempdir()
#' oldwd <- getwd()
#' setwd(td)
#' 
#' # Produce authentic math formulas in the graph
#' plot(estimators, fig.tex = TRUE)
#' system(paste(getOption("pdfviewer"), "ParamEstimator.pdf"))
#' setwd(oldwd)
#' }
#' 
#' @export
plot.estimator <- function(x, estimator=NULL, colorScale.gamma=grey((100:0)/100), colorScale.beta=c("blue","white","red"), legend.cex.axis=1, name.responses=NA, 
                          name.predictors=NA, xlab="", ylab="", fig.tex=FALSE, output="ParamEstimator", header="", header.cex=2, tick=FALSE, mgp=c(2.5,1,0),
                          title.beta=paste("Estimator","$\\hat{\\bm{B}}$"), title.gamma=paste("Estimator","$\\hat{\\mathbf{\\Gamma}}$"),
                          title.Gy=paste("Estimator","$\\hat{\\mathcal{G}}$"), cex.main=1.5,...){
  
  beta_hat <- x$beta
  gamma_hat <- x$gamma
  nonpen <- nrow(beta_hat) - nrow(gamma_hat)
  
  if(is.null(estimator)){
    estimator <- names(x)
  }else{
    if(sum(! estimator %in% names(x))>0)
      stop("Please specify correct argument estimator!")
  }
  
  # specify the labels of axes
  if(is.na(name.responses)[1]) name.responses <- 1:ncol(beta_hat)
  if(name.responses[1] == "auto") name.responses <- colnames(beta_hat)
  if(is.character(name.responses)){
    if(length(name.responses) != ncol(beta_hat)){
      stop("The length of the given response names are not consistent with the data!")
    }
  }
  
  if(is.na(name.predictors)) name.predictors <- 1:nrow(beta_hat)
  if(name.predictors[1] == "auto") name.predictors <- rownames(beta_hat)
  if(is.character(name.predictors)){
    if(length(name.predictors) != nrow(beta_hat)){
      stop("The length of the given predictor names are not consistent with the data!")
    }
  }
  
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))    
  if(!fig.tex){
    
    par(mar=c(6,6,5.1,4.1))
    par(mfrow=c(1,sum(estimator %in% c("beta", "gamma", "Gy"))))
    
    if("beta" %in% estimator){
      # floor(100*constant)+100-1 colours that your want in the legend bar which has the white middle colour
      colorbar <- c(colorRampPalette(c(colorScale.beta[1], colorScale.beta[2]))(floor(1000/(-(max(beta_hat)-min(beta_hat))/min(beta_hat)-1))), colorRampPalette(c(colorScale.beta[2],colorScale.beta[3]))(1000)[-1])
      
      image(z=beta_hat, x=1:nrow(beta_hat), y=1:ncol(beta_hat), col=colorbar, mgp=mgp,  
            axes=ifelse(is.na(name.responses)[1],TRUE,FALSE), xlab=xlab, ylab=ylab,main=expression(hat(bold(B))),cex.main=cex.main,cex.lab=1.5,...);box()
      vertical.image.legend(col=colorbar, zlim=c(min(beta_hat),max(beta_hat)), legend.cex.axis=legend.cex.axis)
      if(!is.na(name.responses)[1]){
        par(las=2, cex.axis=1)
        axis(2, at = 1:ncol(beta_hat), labels=name.responses, tick=tick)
        #opar <- par(cex.axis=1)
        axis(1, at = 1:nrow(beta_hat), labels=name.predictors, tick=tick)
      }
    }
    if("gamma" %in% estimator){
      image(z=gamma_hat, x=1:nrow(gamma_hat), y=1:ncol(gamma_hat), col=colorScale.gamma, mgp=mgp, 
            axes=ifelse(is.na(name.responses)[1],TRUE,FALSE), xlab=xlab, ylab=ylab,main=expression(hat(bold(Gamma))),cex.main=cex.main,cex.lab=1.5,...);box()
      vertical.image.legend(col=colorScale.gamma, zlim=c(0,1), legend.cex.axis=legend.cex.axis)
      if(!is.na(name.responses)[1]){
        par(las=2, cex.axis=1)
        axis(2, at = 1:ncol(gamma_hat), labels=name.responses, tick=tick)
        #opar <- par(cex.axis=1)
        if(nonpen > 0)
          name.predictors <- name.predictors[-c(1:nonpen)]
        axis(1, at = 1:nrow(gamma_hat), labels=name.predictors, tick=tick)
      }
    }
    
    if("Gy" %in% estimator){
      Gy_hat <- x$Gy
      
      image(z=Gy_hat+diag(ncol(Gy_hat)), x=1:nrow(Gy_hat), y=1:nrow(Gy_hat), col=colorScale.gamma, mgp=mgp, 
            axes=ifelse(is.na(name.responses)[1],TRUE,FALSE), xlab=ylab, ylab=ylab,main="Estimated graph of responses",cex.main=cex.main,cex.lab=1.5,...);box()
      vertical.image.legend(col=colorScale.gamma, zlim=c(min(Gy_hat),max(Gy_hat)), legend.cex.axis=legend.cex.axis)
      if(!is.na(name.responses)[1]){
        par(las=2, cex.axis=1)
        axis(2, at = 1:ncol(Gy_hat), labels=name.responses, tick=tick)
        #opar <- par(cex.axis=1)
        axis(1, at = 1:nrow(Gy_hat), labels=name.responses, tick=tick)
      }
    }
    title(paste("\n",header,sep=""), cex.main=header.cex, outer=T)
    
  }else{
    
    options(tikzMetricPackages = c("\\usepackage{amsmath}","\\usepackage{bm}","\\usetikzlibrary{calc}"))
    tikz(paste(output,".tex",sep=""),width=3.6*sum(estimator %in% c("beta", "gamma", "Gy")),height=4,standAlone=TRUE,
         packages=c("\\usepackage{tikz}","\\usepackage{amsmath}","\\usepackage{bm}","\\usepackage[active,tightpage,psfixbb]{preview}","\\PreviewEnvironment{pgfpicture}"))
      
    par(mfrow=c(1,sum(estimator %in% c("beta", "gamma", "Gy"))))
    par(mar = c(6, 6, 4, 4) + 0.1)
    if("beta"%in% estimator){
      # floor(100*constant)+100-1 colours that your want in the legend bar which has the white middle colour
      colorbar <- c(colorRampPalette(c(colorScale.beta[1], colorScale.beta[2]))(floor(1000/(-(max(beta_hat)-min(beta_hat))/min(beta_hat)-1))), colorRampPalette(c(colorScale.beta[2],colorScale.beta[3]))(1000)[-1])
      
      image(z=beta_hat, x=1:nrow(beta_hat), y=1:ncol(beta_hat), col=colorbar, axes=ifelse(is.na(name.responses)[1],TRUE,FALSE), mgp=mgp,  
            xlab=xlab, ylab=ylab, main=title.beta,cex.main=cex.main,cex.lab=1.5,...);box()
      vertical.image.legend(col=colorbar, zlim=c(min(beta_hat),max(beta_hat)), legend.cex.axis=legend.cex.axis)
      if(!is.na(name.responses)[1]){
        par(las=2, cex.axis=1)
        axis(2, at = 1:ncol(beta_hat), labels=name.responses, tick=tick)
        #opar <- par(cex.axis=1)
        axis(1, at = 1:nrow(beta_hat), labels=name.predictors, tick=tick)
      }
    }
    if("gamma" %in% estimator){
      image(z=gamma_hat, x=1:nrow(gamma_hat), y=1:ncol(gamma_hat), col=colorScale.gamma, axes=ifelse(is.na(name.responses)[1],TRUE,FALSE), mgp=mgp,  
            xlab=xlab, ylab=ylab, main=title.gamma,cex.main=cex.main,cex.lab=1.5,...);box()
      vertical.image.legend(col=colorScale.gamma, zlim=c(0,1), legend.cex.axis=legend.cex.axis)
      if(!is.na(name.responses)[1]){
        par(las=2, cex.axis=1)
        axis(2, at = 1:ncol(gamma_hat), labels=name.responses, tick=tick)
        #opar <- par(cex.axis=1)
        if(nonpen > 0)
          name.predictors <- name.predictors[-c(1:nonpen)]
        axis(1, at = 1:nrow(gamma_hat), labels=name.predictors, tick=tick)
      }
    }
    
    if("Gy" %in% estimator){
      Gy_hat <- x$Gy
      
      image(z=Gy_hat+diag(ncol(Gy_hat)), x=1:nrow(Gy_hat), y=1:nrow(Gy_hat), col=colorScale.gamma, axes=ifelse(is.na(name.responses)[1],TRUE,FALSE), mgp=mgp, 
            xlab=ylab, ylab=ylab, main=title.Gy,cex.main=cex.main,cex.lab=1.5,...);box()
      vertical.image.legend(col=colorScale.gamma, zlim=c(min(Gy_hat),max(Gy_hat)), legend.cex.axis=legend.cex.axis)
      
      if(!is.na(name.responses)[1]){
        par(las=2, cex.axis=1)
        axis(2, at = 1:ncol(Gy_hat), labels=name.responses, tick=tick)
        #opar <- par(cex.axis=1)
        axis(1, at = 1:nrow(Gy_hat), labels=name.responses, tick=tick)
      }
    }
    title(paste("\n",header,sep=""), cex.main=header.cex, outer=T)
    dev.off()
    tools::texi2pdf(paste(output,".tex",sep=""))
  }
  
}
# the function vertical.image.legend is orginally from the R package "aqfig"
vertical.image.legend <- function (zlim, col, legend.cex.axis=1) 
{
  starting.par.settings <- par(no.readonly = TRUE)
  #on.exit(par(starting.par.settings))
  mai <- par("mai")
  fin <- par("fin")
  x.legend.fig <- c(1 - (mai[4]/fin[1]), 1)
  y.legend.fig <- c(mai[1]/fin[2], 1 - (mai[3]/fin[2]))
  x.legend.plt <- c(x.legend.fig[1] + (0.18 * (x.legend.fig[2] - 
                                                 x.legend.fig[1])), x.legend.fig[2] - (0.6 * (x.legend.fig[2] - 
                                                                                                x.legend.fig[1])))
  y.legend.plt <- y.legend.fig
  cut.pts <- seq(zlim[1], zlim[2], length = length(col) + 1)
  z <- (cut.pts[1:length(col)] + cut.pts[2:(length(col) + 1)])/2
  par(new = TRUE, pty = "m", plt = c(x.legend.plt, y.legend.plt))
  # If z is not inincreasing, only two values
  if(all(diff(z) > 0)){
    image(x = 1.5, y = z, z = matrix(z, nrow = 1, ncol = length(col)), 
        col = col, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    axis(4, mgp = c(3, 0.2, 0), las = 2, cex.axis = legend.cex.axis, tcl = -0.1)
    box()
  }
  mfg.settings <- par()$mfg
  par(starting.par.settings)
  par(mfg = mfg.settings, new = FALSE)
 
}
