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

int drive_SSUR( arma::mat& Y , arma::mat& X , unsigned int& nChains , unsigned int& nIter , 
				std::string& inFile , std::string& outFilePath , std::string& gammaSampler , bool usingGPrior );

int drive_dSUR( arma::mat& Y , arma::mat& X , unsigned int& nChains , unsigned int& nIter , 
				std::string& inFile , std::string& outFilePath , std::string& gammaSampler , bool usingGPrior );


int drive_HESS( arma::mat& Y , arma::mat& X , unsigned int& nChains , unsigned int& nIter , 
				std::string& inFile , std::string& outFilePath , std::string& gammaSampler , bool usingGPrior );

int drive( unsigned int nIter, unsigned int s, unsigned int p, unsigned int nChains, std::string inFile,
			std::string outFilePath, std::string method, std::string gammaSampler, bool usingGPrior );

#endif