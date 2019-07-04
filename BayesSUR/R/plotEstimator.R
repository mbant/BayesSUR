#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title plotEstimator
#' @description
#' Plot the estimators from the object of fitted Bayesian Seemingly Unrelated Regression
#' @name plotEstimator
#' @param object fitted "runSUR" model
#' @param estimator print the heatmap of estimators. Default "all" is to print all estimators. The value "beta" is for the estimated coefficients matrix, "gamma" for the latent indicator matrix and "Gy" for the graph of responses
#' @param colorScale.gamma value palette for gamma
#' @param colorScale.beta a vector of three colors for diverging color schemes
#' @param legend.cex.axis magnification of axis annotation relative to cex
#' @param name.responses a vector of the response names
#' @param fig.tex print the figure through LaTex. Default is "FALSE"
#' @param output the file name of printed figure
#' @export
plotEstimator <- function(object, estimator="all", colorScale.gamma=grey((100:0)/100), colorScale.beta=c("blue","white","red"), legend.cex.axis=1, name.responses=NA, fig.tex=FALSE, output="ParamEstimator"){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  beta_hat <- as.matrix( read.table(object$output$beta) )
  gamma_hat <- as.matrix( read.table(object$output$gamma) )
  
  if(!fig.tex){
    par(mar=c(6,6,3.1,2.5))
    if(estimator=="all"){
      par(mfrow=c(1,ifelse(toupper(object$input$covariancePrior)=="HIW",3,2))) 
    }
    
    if(estimator=="all" | estimator=="beta"){
      # floor(100*constant)+100-1 colours that your want in the legend bar which has the white middle colour
      colorbar <- c(colorRampPalette(c(colorScale.beta[1], colorScale.beta[2]))(floor(1000/(-(max(beta_hat)-min(beta_hat))/min(beta_hat)-1))), colorRampPalette(c(colorScale.beta[2],colorScale.beta[3]))(1000)[-1])
      
      image(z=beta_hat, x=1:nrow(beta_hat), y=1:ncol(beta_hat), col=colorbar, xlab="", ylab="",main=mtext(bquote(hat(bold(beta)))),cex.main=1.5);box()
      vertical.image.legend(col=colorbar, zlim=c(min(beta_hat),max(beta_hat)), legend.cex.axis=legend.cex.axis)
    }
    if(estimator=="all" | estimator=="gamma"){
      image(z=gamma_hat, x=1:nrow(gamma_hat), y=1:ncol(gamma_hat), col=colorScale.gamma, xlab="", ylab="",main=mtext(bquote(hat(gamma))),cex.main=1.5);box()
      vertical.image.legend(col=colorScale.gamma, zlim=c(min(gamma_hat),max(gamma_hat)), legend.cex.axis=legend.cex.axis)
    }
    
    if(toupper(object$input$covariancePrior) == "HIW"){
      if(estimator=="all" | estimator=="Gy"){
        Gy_hat <- as.matrix( read.table(object$output$G) )
        
        image(z=Gy_hat+diag(ncol(Gy_hat)), x=1:nrow(Gy_hat), y=1:nrow(Gy_hat), col=colorScale.gamma, 
              axes=ifelse(is.na(name.responses)[1],TRUE,FALSE), xlab="", ylab="",main="Estimated graph of responses",cex.main=1.5);box()
        vertical.image.legend(col=colorScale.gamma, zlim=c(min(Gy_hat),max(Gy_hat)), legend.cex.axis=legend.cex.axis)
        if(!is.na(name.responses)[1]){
          par(las=2)
          par(cex.axis=1)
          axis(2, at = 1:dim(Gy_hat)[2], label=name.responses)
          par(cex.axis=1)
          axis(1, at = 1:dim(Gy_hat)[2], label=name.responses)
        }
      }
    }
    #par(mfrow=c(1,1))
  }else{
    #library(tikzDevice)
    options(tikzMetricPackages = c("\\usepackage{amsmath}","\\usepackage{bm}","\\usetikzlibrary{calc}"))
    if(estimator=="all"){
      tikz(paste(output,".tex",sep=""),width=7,height=2.3,standAlone=TRUE,packages=c("\\usepackage{tikz}","\\usepackage{amsmath}",
            "\\usepackage{bm}","\\usepackage[active,tightpage,psfixbb]{preview}","\\PreviewEnvironment{pgfpicture}"))
      par(mfrow=c(1,ifelse(toupper(object$input$covariancePrior)=="HIW",3,2))) 
    }else{
      tikz(paste(output,".tex",sep=""),width=4,height=4,standAlone=TRUE,packages=c("\\usepackage{tikz}","\\usepackage{amsmath}",
            "\\usepackage{bm}","\\usepackage[active,tightpage,psfixbb]{preview}","\\PreviewEnvironment{pgfpicture}"))
    }
    par(mar=c(6,6,3.1,2.5))
    if(estimator=="all" | estimator=="beta"){
      # floor(100*constant)+100-1 colours that your want in the legend bar which has the white middle colour
      colorbar <- c(colorRampPalette(c(colorScale.beta[1], colorScale.beta[2]))(floor(1000/(-(max(beta_hat)-min(beta_hat))/min(beta_hat)-1))), colorRampPalette(c(colorScale.beta[2],colorScale.beta[3]))(1000)[-1])
      
      image(z=beta_hat, x=1:nrow(beta_hat), y=1:ncol(beta_hat), col=colorbar, xlab="", ylab="", main=paste("Estimator","$\\hat{\\bm{B}}$"),cex.main=1.5);box()
      vertical.image.legend(col=colorbar, zlim=c(min(beta_hat),max(beta_hat)), legend.cex.axis=legend.cex.axis)
    }
    if(estimator=="all" | estimator=="gamma"){
      image(z=gamma_hat, x=1:nrow(gamma_hat), y=1:ncol(gamma_hat), col=colorScale.gamma, xlab="", ylab="", main=paste("Estimator","$\\hat{\\mathbf{\\Gamma}}$"),cex.main=1.5);box()
      vertical.image.legend(col=colorScale.gamma, zlim=c(min(gamma_hat),max(gamma_hat)), legend.cex.axis=legend.cex.axis)
    }
    
    if(toupper(object$input$covariancePrior) == "HIW"){
      if(estimator=="all" | estimator=="Gy"){
        Gy_hat <- as.matrix( read.table(object$output$G) )
        
        image(z=Gy_hat+diag(ncol(Gy_hat)), x=1:nrow(Gy_hat), y=1:nrow(Gy_hat), col=colorScale.gamma, 
              axes=ifelse(is.na(name.responses)[1],TRUE,FALSE), xlab="", ylab="", main=paste("Estimator","$\\hat{\\mathcal{G}}$"),cex.main=1.5);box()
        vertical.image.legend(col=colorScale.gamma, zlim=c(min(Gy_hat),max(Gy_hat)), legend.cex.axis=legend.cex.axis)
        
        if(!is.na(name.responses)[1]){
          par(las=2)
          par(cex.axis=1)
          axis(2, at = 1:dim(Gy_hat)[2], label=name.responses)
          par(cex.axis=1)
          axis(1, at = 1:dim(Gy_hat)[2], label=name.responses)
        }
      }
    }
    dev.off()
    tools::texi2pdf(paste(output,".tex",sep=""))
  }
  par(mfrow=c(1,1))
}

# the function vertical.image.legend is orginally from the R package "aqfig"
vertical.image.legend <- function (zlim, col, legend.cex.axis=1) 
{
  starting.par.settings <- par(no.readonly = TRUE)
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
  image(x = 1.5, y = z, z = matrix(z, nrow = 1, ncol = length(col)), 
        col = col, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  axis(4, mgp = c(3, 0.2, 0), las = 2, cex.axis = legend.cex.axis, tcl = -0.1)
  box()
  mfg.settings <- par()$mfg
  par(starting.par.settings)
  par(mfg = mfg.settings, new = FALSE)
}
