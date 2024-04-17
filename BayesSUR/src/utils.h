#ifndef UTILS
#define UTILS

#ifdef CCODE
	#include <iostream>
	#include <armadillo>
#else
	#include <RcppArmadillo.h>
#endif

#include <memory>
#include <string>
#include <cmath>
#include <limits>

#include "global.h"
#include "Parameter_types.h"
#include "pugixml.hpp"

namespace Utils{

	struct SUR_Data
	{	
		std::shared_ptr<arma::mat> data;
		std::shared_ptr<arma::mat> mrfG;
		unsigned int nObservations, nOutcomes, nPredictors, nVSPredictors, nFixedPredictors;
		std::shared_ptr<arma::uvec> outcomesIdx, VSPredictorsIdx, fixedPredictorsIdx;

		arma::ivec blockLabels;
		arma::umat structureGraph;

		std::shared_ptr<arma::umat> missingDataArrayIdx;
		std::shared_ptr<arma::uvec> completeCases;

		SUR_Data() // use this constructor to instanciate all the object at creation (to be sure pointers point to *something*)
		{
			data = std::make_shared<arma::mat>();
			mrfG = std::make_shared<arma::mat>();
			outcomesIdx = std::make_shared<arma::uvec>();
			VSPredictorsIdx = std::make_shared<arma::uvec>();
			fixedPredictorsIdx = std::make_shared<arma::uvec>();

			missingDataArrayIdx = std::make_shared<arma::umat>();
			completeCases = std::make_shared<arma::uvec>();
		}
	};

	struct Chain_Data
	{
		// Data
		Utils::SUR_Data surData;

		// Misc MCMC quantities
		unsigned int nChains = 1 , nIter = 10 , burnin = 0;
        
        int maxThreads = 1 ;
        int tick = 1000 ;
		
		// Parameter and sampler types
		Covariance_Type covariance_type;
		Gamma_Type gamma_type;
		Beta_Type beta_type;
		Gamma_Sampler_Type gamma_sampler_type;

		// hyperparameters
		arma::mat mrfG;

		double mrfD		= std::nan("0"), mrfE 	= std::nan("0") ;
		double sigmaA 	= std::nan("0"), sigmaB = std::nan("0") ;
		double tauA 	= std::nan("0"), tauB 	= std::nan("0") ;
		double nu		= std::nan("0");
		double etaA 	= std::nan("0"), etaB 	= std::nan("0") ;
		double oA 		= std::nan("0"), oB 	= std::nan("0") ;
		double piA 		= std::nan("0"), piB 	= std::nan("0") ;
		double wA 		= std::nan("0"), wB 	= std::nan("0") ;
        double w0A      = std::nan("0"), w0B    = std::nan("0") ;

		// init for some variables
		arma::mat betaInit;
		arma::umat gammaInit;

		// file names and paths
		std::string filePrefix , outFilePath;

		// outputs
		bool output_gamma, output_beta, output_sigmaRho,
			output_Gy, output_pi, output_tail, output_model_size, output_CPO, output_model_visit;
        
	};

	class badFile : public std::exception
	{
		const char * what () const throw ()
		{
			return "The file is either missing or in a wrong format, make sure you're feeding plaintext files.";
		}
	};

	class badBlocks : public std::exception
	{
		const char * what () const throw ()
		{
			return "The program needs at least two blocks, one of responses and one of covariates, to form a regression model.";
		}
	};

	class badSURGraph : public std::exception
	{
		const char * what () const throw ()
		{
			return "Bad structure graph, only one outcome is permitted, check your graph structure or use the SEM version.";
		}
	};

	class badGraph : public std::exception
	{
		const char * what () const throw ()
		{
			return "Bad structure Graph.";
		}
	};

	class badRead : public std::exception
	{
		const char * what () const throw ()
		{
			return "Unknown error: Something went wrong while reading the files. Check your input.";
		}
	};
	

	bool readData(const std::string& dataFileName, std::shared_ptr<arma::mat> data);
    
	bool readGmrf(const std::string& mrfGFileName, std::shared_ptr<arma::mat> mrfG);

	bool readGraph(const std::string& graphFileName, arma::umat& graph);

	bool readBlocks(const std::string& blocksFileName, arma::ivec& blockLabels);

	void removeDisposable(std::shared_ptr<arma::mat> data, arma::ivec& blockLabels);

	void getBlockDimensions(const arma::ivec& blockLabels, const arma::umat& structureGraph,
							const std::shared_ptr<arma::mat>& data, const std::shared_ptr<arma::mat>& mrfG, 
							unsigned int& nObservations,
							unsigned int& nOutcomes, std::shared_ptr<arma::uvec> outcomesIndexes, 
							unsigned int& nPredictors, unsigned int& nVSPredictors, unsigned int& nFixedPredictors,
							std::shared_ptr<arma::uvec> VSPredictorsIndexes, std::shared_ptr<arma::uvec> fixedPredictorsIndexes);

	/* Computes the set-difference from two vectors of indexes */
	arma::uvec arma_setdiff_idx(const arma::uvec& x, const arma::uvec& y);

	void initMissingData(std::shared_ptr<arma::mat> data, std::shared_ptr<arma::umat> missingDataArrayIndexes, std::shared_ptr<arma::uvec> completeCases, bool print=false );

	void formatData(const std::string& dataFileName, const std::string& mrfGFileName, const std::string& blockFileName, const std::string& structureGraphFileName, 
					SUR_Data& surData );

	void readHyperPar(const std::string& hyperParFile, Chain_Data& chainData );

	template <typename T> int sgn(T val)
	{
		return (T(0) < val) - (val < T(0));
	}

	inline double internalRound( double val )
	{
		if( val < 0 ) return ceil(val - 0.5);
		return floor(val + 0.5);
	}

	inline double round( double val , unsigned int decimals )
	{
		return internalRound( val * std::pow(10.,decimals) ) / std::pow(10.,decimals) ;
	}

	double logspace_add(const arma::vec& logv);
	double logspace_add(double a,double b);

	arma::uvec nonZeroLocations_row( arma::sp_umat X);  // if you pass a row subview
	arma::uvec nonZeroLocations_col( arma::sp_umat X); // if you pass a col subview

	// arma::uvec nonZeroLocations_col( arma::sp_umat X, unsigned int column);
	// arma::uvec nonZeroLocations_row( arma::sp_umat X, unsigned int row);
		
}

#endif
