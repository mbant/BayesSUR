#ifndef UTILS
#define UTILS

#include <iostream>
#include <string>
#include <memory>
#include <armadillo>


namespace Utils{

	struct SUR_Data
	{	
		std::shared_ptr<arma::mat> data;
		unsigned int nObservations, nOutcomes, nPredictors, nVSPredictors, nFixedPredictors;
		std::shared_ptr<arma::uvec> outcomesIdx, VSPredictorsIdx, fixedPredictorsIdx;

		arma::ivec blockLabels;
		arma::umat structureGraph;

		std::shared_ptr<arma::umat> missingDataArrayIdx;
		std::shared_ptr<arma::uvec> completeCases;

		SUR_Data() // use this constructor to instanciate all the object at creation (to be sure pointers point to *something*)
		{
			data = std::make_shared<arma::mat>();
			outcomesIdx = std::make_shared<arma::uvec>();
			VSPredictorsIdx = std::make_shared<arma::uvec>();
			fixedPredictorsIdx = std::make_shared<arma::uvec>();

			missingDataArrayIdx = std::make_shared<arma::umat>();
			completeCases = std::make_shared<arma::uvec>();
		}
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
    
	bool readGmrf(const std::string& mrfGFileName, arma::mat& mrfG);

	bool readGraph(const std::string& graphFileName, arma::umat& graph);

	bool readBlocks(const std::string& blocksFileName, arma::ivec& blockLabels);

	void removeDisposable(std::shared_ptr<arma::mat> data, arma::ivec& blockLabels);

	void getBlockDimensions(const arma::ivec& blockLabels, const arma::umat& structureGraph,
							const std::shared_ptr<arma::mat>& data, unsigned int& nObservations,
							unsigned int& nOutcomes, std::shared_ptr<arma::uvec> outcomesIndexes, 
							unsigned int& nPredictors, unsigned int& nVSPredictors, unsigned int& nFixedPredictors,
							std::shared_ptr<arma::uvec> VSPredictorsIndexes, std::shared_ptr<arma::uvec> fixedPredictorsIndexes);

	/* Computes the set-difference from two vectors of indexes */
	arma::uvec arma_setdiff_idx(const arma::uvec& x, const arma::uvec& y);

	void initMissingData(std::shared_ptr<arma::mat> data, std::shared_ptr<arma::umat> missingDataArrayIndexes, std::shared_ptr<arma::uvec> completeCases, bool print=false );

	void formatData(const std::string& dataFileName, const std::string& blockFileName, const std::string& structureGraphFileName, 
					SUR_Data& surData );

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
