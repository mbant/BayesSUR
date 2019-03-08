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
	
using Utils::Chain_Data;

template<typename T>
void setHyperParameters( ESS_Sampler<T>& chain, const Chain_Data& chainData );

int drive_SUR( Chain_Data& chainData );

int drive_HESS( Chain_Data& chainData );

int drive( const std::string& dataFile, const std::string& blockFile, const std::string& structureGraphFile, const std::string& hyperParFile, const std::string& outFilePath,  
			unsigned int nIter, unsigned int burnin, unsigned int nChains,
			const std::string& covariancePrior, 
			const std::string& gammaPrior, const std::string& gammaSampler, const std::string& gammaInit, const std::string& mrfGFile ,
			const std::string& betaPrior,
			bool output_gamma, bool output_beta, bool output_G, bool output_sigmaRho, bool output_pi, bool output_tail, bool output_model_size );

#endif