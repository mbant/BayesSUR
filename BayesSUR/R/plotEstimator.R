#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title plotEstimator
#' @description
#' Plot the estimators from the object of fitted Bayesian Seemingly Unrelated Regression
#' @name plotEstimator
#' @param object fitted "runSUR" model
#' @param colorScale.gamma value palette for gamma
#' @param colorScale.beta a vector of three colors for diverging color schemes
#' @param legend.cex.axis magnification of axis annotation relative to cex
#' @param fig.tex print the figure through LaTex. Default is "FALSE"
#' @export
plotEstimator <- function(object, colorScale.gamma=grey((100:0)/100), colorScale.beta=c("blue","white","red"), legend.cex.axis=1, fig.tex=FALSE){
  
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  beta_hat <- as.matrix( read.table(object$output$beta) )
  gamma_hat <- as.matrix( read.table(object$output$gamma) )
  
  
  # floor(100*constant)+100-1 colours that your want in the legend bar which has the white middle colour
  colorbar <- c(colorRampPalette(c(colorScale.beta[1], colorScale.beta[2]))(floor(1000/(-(max(beta_hat)-min(beta_hat))/min(beta_hat)-1))), colorRampPalette(c(colorScale.beta[2],colorScale.beta[3]))(1000)[-1])
   
  if(!fig.tex){
    par(mfrow=c(1,ifelse(toupper(object$input$covariancePrior)=="HIW",3,2))) 
    
    image(z=beta_hat, x=1:nrow(beta_hat), y=1:ncol(beta_hat), col=colorbar, xlab="", ylab="",main=mtext(bquote(hat(bold(beta)))));box()
    vertical.image.legend(col=colorbar, zlim=c(min(beta_hat),max(beta_hat)), legend.cex.axis=legend.cex.axis)
    image(z=gamma_hat, x=1:nrow(gamma_hat), y=1:ncol(gamma_hat), col=colorScale.gamma, xlab="", ylab="",main=mtext(bquote(hat(gamma))));box()
     vertical.image.legend(col=colorScale.gamma, zlim=c(min(gamma_hat),max(gamma_hat)), legend.cex.axis=legend.cex.axis)
  
    if(toupper(object$input$covariancePrior) == "HIW"){
       G0_hat <- as.matrix( read.table(object$output$G) )
       image(z=G0_hat+diag(ncol(G0_hat)), x=1:nrow(G0_hat+diag(ncol(G0_hat))), y=1:ncol(G0_hat+diag(ncol(G0_hat))), col=colorScale.gamma, xlab="", ylab="",main="Estimated graph of responses");box()
       vertical.image.legend(col=colorScale.gamma, zlim=c(min(G0_hat),max(G0_hat)), legend.cex.axis=legend.cex.axis)
    }
    par(mfrow=c(1,1))
  }else{
    #library(tikzDevice)
    options(tikzMetricPackages = c("\\usepackage{amsmath}","\\usepackage{bm}","\\usetikzlibrary{calc}"))
    tikz('ParamEstimator.tex',width=6,height=2.5,standAlone=TRUE,packages=c("\\usepackage{tikz}","\\usepackage{amsmath}","\\usepackage{bm}"))
    par(mfrow=c(1,ifelse(toupper(object$input$covariancePrior)=="HIW",3,2))) 
    
    image(z=beta_hat, x=1:nrow(beta_hat), y=1:ncol(beta_hat), col=colorbar, xlab="", ylab="", main=paste("Estimator","$\\hat{\\mathbf{B}}$"));box()
    vertical.image.legend(col=colorbar, zlim=c(min(beta_hat),max(beta_hat)), legend.cex.axis=legend.cex.axis)
    image(z=gamma_hat, x=1:nrow(gamma_hat), y=1:ncol(gamma_hat), col=colorScale.gamma, xlab="", ylab="", main=paste("Estimator","$\\hat{\\Gamma}$"));box()
    vertical.image.legend(col=colorScale.gamma, zlim=c(min(gamma_hat),max(gamma_hat)), legend.cex.axis=legend.cex.axis)
    
    if(toupper(object$input$covariancePrior) == "HIW"){
      G0_hat <- as.matrix( read.table(object$output$G) )
      image(z=G0_hat+diag(ncol(G0_hat)), x=1:nrow(G0_hat+diag(ncol(G0_hat))), y=1:ncol(G0_hat+diag(ncol(G0_hat))), col=colorScale.gamma, xlab="", ylab="", main=paste("Estimator","$\\hat{\\mathcal{G}}$"));box()
      vertical.image.legend(col=colorScale.gamma, zlim=c(min(G0_hat),max(G0_hat)), legend.cex.axis=legend.cex.axis)
    }
    dev.off()
    tools::texi2pdf("ParamEstimator.tex")
  }

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
