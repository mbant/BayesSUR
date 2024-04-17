#ifndef DRIVESUR
#define DRIVESUR

#include <vector>
#include <string>

#ifndef CCODE
	#include <RcppArmadillo.h>
	#include <fstream>
#else
	#include <armadillo>
	#include <iostream>
#endif

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "global.h"
#include "utils.h"
#include "distr.h"

#include "ESS_Sampler.h"
#include "HRR_Chain.h"
#include "SUR_Chain.h"
	
using Utils::Chain_Data;

template<typename T>
void setHyperParameters( ESS_Sampler<T>& chain, const Chain_Data& chainData );

int drive_SUR( Chain_Data& chainData );

int drive_HRR( Chain_Data& chainData );

int drive( const std::string& dataFile, const std::string& mrfGFile, const std::string& blockFile, const std::string& structureGraphFile, const std::string& hyperParFile, const std::string& outFilePath,  
			unsigned int nIter, unsigned int burnin, unsigned int nChains,
			const std::string& covariancePrior, 
			const std::string& gammaPrior, const std::string& gammaSampler, const std::string& gammaInit,
			const std::string& betaPrior, const int maxThreads, const int tick,
			bool output_gamma, bool output_beta, bool output_Gy, bool output_sigmaRho, bool output_pi, bool output_tail, bool output_model_size,
            bool output_CPO, bool output_model_visit );

#endif
