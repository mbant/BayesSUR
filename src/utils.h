#ifndef UTILS
#define UTILS

#include <iostream>
#include <string>
#include <armadillo>

namespace Utils{

	bool readData(std::string fileName, unsigned int &nOutcomes, unsigned int &nPredictors, unsigned int &nObservations, arma::mat &Y, arma::mat& X);

	bool readDataSEM(std::string fileName, arma::mat &data, arma::ivec &blockIndexes, arma::ivec &variableType,
			arma::uvec &missingDataIndexes, arma::uvec &nOutcomes, unsigned int &nPredictors, unsigned int &nObservations);

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