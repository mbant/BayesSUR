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
  num.blue <- sum(beta_hat<(-0.0001)) 
  num.gray <- sum(beta_hat>0.0001) 
  
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
  #colorbar <- c(colorRampPalette(c("blue", "white"))(100/(1+min(beta_hat)/(max(beta_hat)-min(beta_hat)))), grey((100:0)/100)[-1])
  
  colorbar <- c(colorRampPalette(c("blue", "white"))(floor(100/(-(max(beta_hat)-min(beta_hat))/min(beta_hat)-1))), grey(seq(1,min(max(beta_hat),1),length=100))[-1])
  
  #colorbar <- c(rgb(floor(seq(153,204,length=num.blue)) ,floor(seq(204,229,length=num.blue)) ,255 , maxColorValue = 255), rep(rgb(255,255,255,maxColorValue = 255), sum(beta_hat>=-0.0001 & beta_hat<=0.0001)),
  #              rgb(floor(seq(255,0,length=num.gray)) ,floor(seq(255,0,length=num.gray)) ,floor(seq(255,0,length=num.gray)), maxColorValue = 255))
  #colorbar <- colorbar[sort.list(beta_hat)]
  image(beta_hat, col=colorbar, axes = FALSE, main=mtext(bquote(hat(bold(beta)))));box()
  image(gamma_hat, col=colorScale.gamma, axes = FALSE, main=mtext(bquote(hat(gamma))));box()
  
  if(toupper(object$input$covariancePrior) == "HIW"){
    G0_hat <- as.matrix( read.table(object$output$G) )
    image((G0_hat+diag(ncol(G0_hat))), col=colorScale.gamma, axes = FALSE, main="Estimated graph of responses");box()
  }
  
  par(mfrow=c(1,1))

}