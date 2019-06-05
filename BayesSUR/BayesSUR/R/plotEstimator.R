#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title plotEstimator
#' @description
#' Plot the estimators from the object of fitted Bayesian Seemingly Unrelated Regression
#' @name plotEstimator
#' @param object fitted "runSUR" model
#' @param greyscale grey value palette
#' @export
plotEstimator <- function(object, greyscale=grey((0:100)/100)){
  
  beta_hat <- as.matrix( read.table(object$output$beta) )
  gamma_hat <- as.matrix( read.table(object$output$gamma) )
  
  par(mfrow=c(1,ifelse(toupper(object$input$covariancePrior)=="HIW",3,2)))
  image(beta_hat, col=greyscale, axes = FALSE, main=mtext(bquote(hat(bold(beta)))));box()
  image(gamma_hat, col=greyscale, axes = FALSE, main=mtext(bquote(hat(gamma))));box()
  
  if(toupper(object$input$covariancePrior) == "HIW"){
    G0_hat <- as.matrix( read.table(object$output$G) )
    image((G0_hat+diag(ncol(G0_hat))), col=greyscale, axes = FALSE, main=mtext(bquote(hat(G)[0])));box()
  }
  
  par(mfrow=c(1,1))

}