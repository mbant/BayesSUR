#' R2SSUR -- Bayesian Sparse Seemingly Unrelated Regression
#' Plot the results from the object of fitted Bayesian Sparse Seemingly Unrelated Regression
#' @export
plot <- function(object, y=NULL, greyscale=grey((1000:0)/1000)){
  library(igraph)
  
  if(toupper(object$input$covariancePrior) == "HIW"){
    G0_hat = as.matrix( read.table(object$output$G) )
  }
  
  gamma_hat = as.matrix( read.table(object$output$gamma) )
  beta_hat = as.matrix( read.table(object$output$beta) )
  
  threshold = 0.2
  sigmaRho = as.matrix( read.table(object$output$sigmaRho) > threshold )
  net = graph_from_adjacency_matrix(  sigmaRho, weighted=T, mode="undirected", diag=F); V(net)$size = 30
  
  par(mfrow=c(ifelse(is.null(y),2,3), 2))
  plot.igraph(net, main = "Network of the response variables")
  if(!is.null(y)){
    plot.default(1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
    image(x=1:dim(y)[2],y=1:dim(y)[2],cor(y)[,ncol(y):1],axes=FALSE,xlab="",ylab="",main="Correlation between the response variables",col=greyscale);box()
    #image(x=1:dim(y)[2],y=1:dim(y)[2],cor(y)[,ncol(y):1],axes=FALSE,xlab="",ylab="",zlim=c(-1,1),main="Correlation between the response variables",col=greyscale);box()
  }
  if(toupper(object$input$covariancePrior) == "HIW"){
    image((G0_hat+diag(ncol(G0_hat)))[ncol(G0_hat):1,], col=greyscale, axes = FALSE, main=mtext(bquote(hat(G)[0])));box()
  }else{
    plot.default(1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  }
  image(gamma_hat, col=greyscale, axes = FALSE, main=mtext(bquote(hat(Gamma))));box()
  image(beta_hat, col=greyscale, axes = FALSE, main=mtext(bquote(hat(bold(Beta)))));box()
  par(mfrow=c(1,1))

  return(0)
}