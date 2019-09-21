// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(igraph)]]
// [[Rcpp::depends(Matrix)]]
//' @title BayesSUR_internal
//' @description
//' Run a SUR Bayesian sampler -- internal function
//' @name BayesSUR_internal
//' @param dataFile path to data file
//' @param outFilePath path to where the output is to be written
//' @param nIter number of iterations
//' @param nChains number of parallel chains to run
//
// NOTE THAT THIS IS BASICALLY JUST A WRAPPER

#include "drive.h"
#include <RcppArmadillo.h>

using Rcpp::Rcerr;


// [[Rcpp::export(rng=false)]]
int BayesSUR_internal(const std::string& dataFile, const std::string& mrfGFile, const std::string& blockFile, const std::string& structureGraphFile, const std::string& hyperParFile, const std::string& outFilePath,
                    unsigned int nIter=10, unsigned int burnin=0, unsigned int nChains=1,
                    const std::string& covariancePrior="HIW", 
                    const std::string& gammaPrior="hotspot", const std::string& gammaSampler="bandit", 
                    const std::string& gammaInit = "MLE",
                    const std::string& betaPrior="independent",
                    bool output_gamma = true, bool output_beta = true, bool output_G = true, bool output_sigmaRho = true, 
                    bool output_pi = true, bool output_tail = true, bool output_model_size = true, bool output_CPO = true )
{
  int status {1};
  
  try
  {
    status =  drive(dataFile,mrfGFile,blockFile,structureGraphFile,hyperParFile,outFilePath,nIter,burnin,nChains,
                    covariancePrior,gammaPrior,gammaSampler,gammaInit,betaPrior, output_gamma, output_beta,
                    output_G, output_sigmaRho, output_pi, output_tail, output_model_size, output_CPO);
  }
  catch(const std::exception& e)
  {
    Rcerr << e.what() << '\n'; // we can use Rcerr here because we're reaching here from R for sure
  }
  
  return status;
  
}
