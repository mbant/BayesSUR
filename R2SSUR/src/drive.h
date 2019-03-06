#ifndef DRIVESUR
#define DRIVESUR

#include <vector>
#include <iostream>
#include <string>
#include <armadillo>
// #include <tgmath.h>
// #include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "global.h"
#include "utils.h"
#include "distr.h"

#include "ESS_Sampler.h"
#include "HESS_Chain.h"
#include "SUR_Chain.h"


struct Chain_Data
{
	// Data
	Utils::SUR_Data surData;

	// Misc MCMC quantities
	unsigned int nChains = 1 , nIter = 10 , burnin = 0;
	
	// Parameter and sampler types
	Covariance_Type covariance_type;
	Gamma_Type gamma_type;
	Beta_Type beta_type;
	Gamma_Sampler_Type gamma_sampler_type;

	arma::mat mrfG;

	// init for some variables
	arma::mat betaInit;
	arma::umat gammaInit;

	// file names and paths
	std::string filePrefix , outFilePath;

};
	

int drive_SUR( Chain_Data& chainData );

int drive_HESS( Chain_Data& chainData );

int drive( const std::string& dataFile, const std::string& blockFile, const std::string& structureGraphFile, const std::string& outFilePath,  
			unsigned int nIter, unsigned int burnin, unsigned int nChains,
			const std::string& covariancePrior, 
			const std::string& gammaPrior, const std::string& gammaSampler, const std::string& gammaInit, const std::string& mrfGFile ,
			const std::string& betaPrior );

#endif