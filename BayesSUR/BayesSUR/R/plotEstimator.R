#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title plotEstimator
#' @description
#' Plot the estimators from the object of fitted Bayesian Seemingly Unrelated Regression
#' @name plotEstimator
#' @param object fitted "runSUR" model
#' @param colorScale.betay value palette for beta
#' @param colorScale.gamma value palette for gamma
#' @export
plotEstimator <- function(object, colorScale.gamma=grey((100:0)/100)){
  
  devAskNewPage(FALSE)
  object$output[-1] <- paste(object$output$outFilePath,object$output[-1],sep="")
  beta_hat <- as.matrix( read.table(object$output$beta) )
  gamma_hat <- as.matrix( read.table(object$output$gamma) )
  
  par(mfrow=c(1,ifelse(toupper(object$input$covariancePrior)=="HIW",3,2)))
  #num.blue <- sum(beta_hat<(-0.0001)) 
  #num.gray <- sum(beta_hat>0.0001) 
  
  # blue.r <- blue.g <- blue.b <- floor(seq(153,204,length=num.blue))
  # blue.g <- floor(seq(204,229,length=num.blue))
  # blue.b<- rep(255, num.blue)
  # white.r <- white.g <- white.b <- rep(255, sum(beta_hat>=-0.001 & beta_hat<=0.001))
  # gray.r <- gray.g <- gray.b <- floor(seq(255,0,length=num.gray))
  # color.r <- c(blue.r, white.r, gray.r)
  # color.g <- c(blue.g, white.g, gray.g)
  # color.b <- c(blue.b, white.b, gray.b)
  # colorbar <- rep(NA, prod(dim(beta_hat)))
  # for(i in 1:prod(dim(beta_hat))) colorbar[i] <- rgb(color.r[i],color.g[i],color.b[i], maxColorValue = 255)
  # colorbar <- c(colorRampPalette(c("blue", "white"))(100/(1+min(beta_hat)/(max(beta_hat)-min(beta_hat)))), grey((100:0)/100)[-1])
  
  # floor(100*constant)+100-1 colours that your want in the legend bar which has the white middle colour
  colorbar <- c(colorRampPalette(c("blue", "white"))(floor(100/(-(max(beta_hat)-min(beta_hat))/min(beta_hat)-1))), colorRampPalette(c("white","red"))(100)[-1])
  
  #colorbar <- c(rgb(floor(seq(153,204,length=num.blue)) ,floor(seq(204,229,length=num.blue)) ,255 , maxColorValue = 255), rep(rgb(255,255,255,maxColorValue = 255), sum(beta_hat>=-0.0001 & beta_hat<=0.0001)),
  #              rgb(floor(seq(255,0,length=num.gray)) ,floor(seq(255,0,length=num.gray)) ,floor(seq(255,0,length=num.gray)), maxColorValue = 255))
  #colorbar <- colorbar[sort.list(beta_hat)]
  image(beta_hat, col=colorbar, axes = FALSE, main=mtext(bquote(hat(bold(beta)))));box()
  vertical.image.legend(col=colorbar, zlim=c(min(beta_hat),max(beta_hat)))
  image(gamma_hat, col=colorScale.gamma, axes = FALSE, main=mtext(bquote(hat(gamma))));box()
  vertical.image.legend(col=colorScale.gamma, zlim=c(0,1))
  
  if(toupper(object$input$covariancePrior) == "HIW"){
    G0_hat <- as.matrix( read.table(object$output$G) )
    image((G0_hat+diag(ncol(G0_hat))), col=colorScale.gamma, axes = FALSE, main="Estimated graph of responses");box()
    vertical.image.legend(col=colorScale.gamma, zlim=c(0,1))
  }
  
  par(mfrow=c(1,1))

}

# the function vertical.image.legend is orginally from the R package "aqfig"
vertical.image.legend <- function (zlim, col) 
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
  axis(4, mgp = c(3, 0.2, 0), las = 2, cex.axis = 0.5, tcl = -0.1)
  box()
  mfg.settings <- par()$mfg
  par(starting.par.settings)
  par(mfg = mfg.settings, new = FALSE)
}