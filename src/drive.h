#include <vector>
#include <iostream>
#include <string>
#include <armadillo>
#include <tgmath.h>
#include <limits>
#include <omp.h>

#include "global.h"
#include "utils.h"
#include "distr.h"

#include "ESS_Sampler.h"
#include "HESS_Chain.h"
#include "SSUR_Chain.h"
#include "dSUR_Chain.h"

#ifndef DRIVESUR
#define DRIVESUR

struct Chain_Data
{
	Utils::SUR_Data surData;
	unsigned int nChains , nIter ;
	std::string gammaSampler;
	
	bool usingGPrior;

	arma::mat betaInit;
	arma::umat gammaInit;

	std::string filePrefix , outFilePath;
};
	

int drive_SSUR( Chain_Data& chainData );

int drive_dSUR( Chain_Data& chainData );

int drive_HESS( Chain_Data& chainData );

int drive( const std::string& dataFile, const std::string& blockFile, const std::string& structureGraphFile, const std::string& outFilePath,  
			unsigned int nIter, unsigned int nChains,
			const std::string& method, const std::string& gammaSampler, const std::string& gammaInit, bool usingGPrior );

#endif