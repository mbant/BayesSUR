% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotNetwork.R
\name{plotNetwork}
\alias{plotNetwork}
\title{plotNetwork}
\usage{
plotNetwork(object, includeResponse = NULL, excludeResponse = NULL,
  includePredictor = NULL, excludePredictor = NULL,
  MatrixGamma = NULL, PmaxPredictor = 0.5, PmaxResponse = 0.5,
  nodesizePredictor = 15, nodesizeResponse = 25, no.isolates = FALSE,
  lineup = 0.8, gray.alpha = 0.6)
}
\arguments{
\item{object}{fitted "runSUR" model}

\item{includeResponse}{A vector of the response names which are shown in the network.}

\item{excludeResponse}{A vector of the response names which are not shown in the network.}

\item{includePredictor}{A vector of the predictor names which are shown in the network.}

\item{excludePredictor}{A vector of the predictor names which are not shown in the network.}

\item{MatrixGamma}{A matrix or dataframe of the latent indicator variable. Default is "NULL" and to extrate it from object of class inheriting from "runSUR"}

\item{PmaxPredictor}{cutpoint for thresholding the estimated latent indicator variable. Default is 0.5}

\item{PmaxResponse}{cutpoint for thresholding the learning structure matrix of multiple response variables. Default is 0.5}

\item{nodesizePredictor}{node size of Predictors in the output graph. Default is 15}

\item{no.isolates}{remove isolated nodes from responses graph and Full graph, may get problem if there are also isolated Predictors}

\item{lineup}{A ratio of the heights between responses' area and Predictors'}

\item{nodesizePredictor}{node size of response variables in the output graph. Default is 25}
}
\description{
Network representation of the associations between responses and features
}
\details{
BayesSUR -- Bayesian Seemingly Unrelated Regression
}