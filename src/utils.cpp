#include "utils.h"

#include <iostream>
#include <string>
#include <armadillo>
#include <cmath>

#include <limits>

namespace Utils{


	bool readData(std::string fileName, unsigned int &nOutcomes, unsigned int &nPredictors, unsigned int &nObservations, arma::mat &Y, arma::mat& X)
	{

		bool status = X.load(fileName,arma::raw_ascii);

		if(!status){ 
			std::cout<< "Somethign went wrong while reading "<<fileName<<std::endl;
			return false;
		}else{

			if(X.n_cols < (nOutcomes+nPredictors) ){
				std::cout<< "Less columns than the sum of the specified outcomes and predictors. Specify the correct numbers. "<<std::endl;
				return false;
			}

			Y = X.cols(0,nOutcomes-1);
			X.shed_cols(0,nOutcomes-1);

			nObservations = Y.n_rows;

			if( nPredictors < X.n_cols){
				X.shed_cols(nPredictors,X.n_cols-1);
				std::cout<< "More predictors read then what specified in p=" << nPredictors << " -- I kept the firsts and discarded the others." << std::endl;
				nPredictors = X.n_cols;
			}

			if( nOutcomes < Y.n_cols){
				Y.shed_cols(nOutcomes,Y.n_cols-1);
				std::cout<< "More outcomes read then what specified in q=" << nOutcomes << " -- I kept the firsts and discarded the others." << std::endl;
				nOutcomes = Y.n_cols;
			}
		}

		return true;
	}


	bool readDataSEM(std::string fileName, arma::mat &data, arma::ivec &blockIndexes, arma::ivec &variableType,
		arma::uvec &missingDataIndexes, arma::uvec &nOutcomes, unsigned int &nPredictors, unsigned int &nObservations)
	{

		bool status = data.load(fileName,arma::raw_ascii);
		arma::uvec jnk;
		arma::ivec ijnk;
		
		if(!status){ 
			std::cout<< "Somethign went wrong while reading "<<fileName<<std::endl;
			return false;
		}else{

			// the first row is a row of blockIndexes
			blockIndexes = arma::conv_to<arma::ivec>::from( data.row(0)/*.t()*/ );
			data.shed_row(0);
			
			// checks on the blockIndexes
			// index 0 stands for the Xs, predictors
			// index 1+ are the upper-level outcomes
			// so we always need at least some zeros and some ones
			ijnk = arma::unique(blockIndexes);
			if( arma::max( blockIndexes ) < 1 || jnk.n_elem < 2 ) // more indepth check would be length of positive indexes..
			{
				std::cout<< "You need to define at least two blocks -- Xs (block 0) and Ys (block 1)"<<std::endl;
				return false;
			}


			// the second row is a row of variableType
			variableType = arma::conv_to<arma::ivec>::from( data.row(0)/*.t()*/ );
			data.shed_row(0);
			

			// miscellanea variables
			for( int i=0, n=arma::max(blockIndexes) ; i<n ; ++i )
			{
				jnk = arma::find( blockIndexes == i );   // meh, hate this temporary but is needed
				nOutcomes(i) = jnk.n_elem;
			}
			// first index is the number of predictors
			nPredictors = nOutcomes(0);
			nOutcomes.shed_row(0);
			
			nObservations = data.n_rows;

			// Now deal with NANs
			if( data.has_nan() )
			{
				missingDataIndexes = arma::find_nonfinite(data);
				data(missingDataIndexes).fill( arma::datum::nan );  // This makes all the ind values into valid armadillo NANs (should be ok even without, but..)
			}
			else
			{
				missingDataIndexes.set_size(0);
			}
			
		}

		return true;
	}


	// sgn is defined in the header in order for it to be visible

	double logspace_add(const arma::vec& logv)
	{

		if( logv.is_empty() )
			return std::numeric_limits<double>::lowest();
		double m;
		if( logv.has_inf() || logv.has_nan() ) // || logv.is_empty()
		{
			return logspace_add(logv.elem(arma::find_finite(logv)));
		}else{ 
			m = arma::max(logv);
			return m + std::log( (double)arma::as_scalar( arma::sum( arma::exp( - (m - logv) ) ) ) );
		}
	}

	double logspace_add(double a, double b)
	{

		if(a <= std::numeric_limits<float>::lowest())
      		return b;
    	if(b <= std::numeric_limits<float>::lowest())
      		return a;
    	return std::max(a, b) + std::log( (double)(1. + std::exp( (double)-std::abs((double)(a - b)) )));
	}

	arma::uvec nonZeroLocations_col( arma::sp_umat X)
	{
		std::vector<arma::uword> locations;

		for( arma::sp_umat::const_iterator it = X.begin(); it != X.end(); ++it)
		{
			locations.push_back(it.row());
		}

		return arma::uvec(locations);
	}

	arma::uvec nonZeroLocations_row( arma::sp_umat X)
	{
		std::vector<arma::uword> locations;

		for( arma::sp_umat::const_iterator it = X.begin(); it != X.end(); ++it)
		{
			locations.push_back(it.col());
		}

		return arma::uvec(locations);
	}

}