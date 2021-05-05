#include "utils.h"

#ifndef CCODE
	using Rcpp::Rcout;
	using Rcpp::Rcerr;
#else
	#define Rcout std::cout
	#define Rcerr std::cerr
#endif

namespace Utils{

	bool readData(const std::string& dataFileName, std::shared_ptr<arma::mat> data)
	{

		bool status = data->load(dataFileName,arma::raw_ascii);
		if( !status )
			throw badFile();

		return status;
	}
    
    bool readGmrf(const std::string& mrfGFileName, std::shared_ptr<arma::mat> mrfG)
    {
        
        bool status = mrfG->load(mrfGFileName,arma::raw_ascii);
        if( !status )
            throw badFile();

        return status;
        
    }

	bool readGraph(const std::string& graphFileName, arma::umat& graph)
	{

		bool status = graph.load(graphFileName,arma::raw_ascii);
		if( !status )
			throw badFile();

		return status;
	}

	bool readBlocks(const std::string& blocksFileName, arma::ivec& blockLabels)
	{

		bool status = blockLabels.load(blocksFileName,arma::raw_ascii);
		if( !status )
			throw badFile();

		// checks on the blockLabels
		// index 0 stands for the Xs, predictors
		// index 1+ are the upper-level outcomes
		// -1 are for variable to be excluded from analysis
		// so we always need at least some zeros and some ones
		arma::ivec uniqueblockLabels = arma::unique(blockLabels);

		if( arma::max( blockLabels ) < 1 || uniqueblockLabels.n_elem < 2 ) // more indepth check would be length of positive indexes..
			throw badBlocks();

		// remember R produces these labels and thus checks the dimensions (need to correspond to data dimensions) before entering the C++ code
		return status;
	}

	void removeDisposable(std::shared_ptr<arma::mat> data, arma::ivec& blockLabels)
	{
		arma::uword shedIdx;
		while( arma::any( arma::find(blockLabels < 0)) )
		{
			shedIdx = arma::as_scalar(arma::find(blockLabels < 0 , 1 , "first"));
			
			// then shed the rest		
			data->shed_col( shedIdx );
			blockLabels.shed_row( shedIdx ); //shed the blockIdx as well!
		}
		return;
	}


	void getBlockDimensions(const arma::ivec& blockLabels, const arma::umat& structureGraph,
							const std::shared_ptr<arma::mat>& data, const std::shared_ptr<arma::mat>& mrfG, unsigned int& nObservations,
							unsigned int& nOutcomes, std::shared_ptr<arma::uvec> outcomeIndexes, 
							unsigned int& nPredictors, unsigned int& nVSPredictors, unsigned int& nFixedPredictors,
							std::shared_ptr<arma::uvec> VSPredictorsIndexes, std::shared_ptr<arma::uvec> fixedPredictorsIndexes)
	{
		// dimensions variables
		nObservations = data->n_rows;

		// define the structure
		// structureGraph(i,j) != 0 means an arrow j->i
		arma::uvec allOutcomeLabels = arma::find( arma::sum(structureGraph,0)!=0 );  // so search for sum over columns different from zero
		unsigned int nEquations = allOutcomeLabels.n_elem;
		if ( nEquations > 1 || nEquations == 0 )
			throw badSURGraph();


		// outcomes
		arma::uword outcomeLabel{ allOutcomeLabels(0) };
		(*outcomeIndexes) = arma::find( blockLabels == outcomeLabel );   // groups in the graph are ordered by their position in the blockList
		nOutcomes = outcomeIndexes->n_elem;
		
		// Predictors
		arma::uvec VSPredictorsLabels = arma::find( structureGraph.col(outcomeLabel) == 1 );
		arma::uvec fixedPredictorsLabels = arma::find( structureGraph.col(outcomeLabel) == 2 );


		// reset
		VSPredictorsIndexes->clear();
		fixedPredictorsIndexes->clear();

		for( auto label : VSPredictorsLabels )
		{
			VSPredictorsIndexes->insert_rows( VSPredictorsIndexes->n_elem , arma::find( blockLabels == label ) );
		}

		for( auto label : fixedPredictorsLabels )
		{
			fixedPredictorsIndexes->insert_rows( fixedPredictorsIndexes->n_elem , arma::find( blockLabels == label ) );
		}

		nVSPredictors = VSPredictorsIndexes->n_elem;
		nFixedPredictors = fixedPredictorsIndexes->n_elem;
		nPredictors = nVSPredictors + nFixedPredictors;

		return;
	}


	/* Computes the set-difference from two vectors of indexes */
	arma::uvec arma_setdiff_idx(const arma::uvec& x, const arma::uvec& y)
	{

		arma::uvec ux = arma::unique(x);
		arma::uvec uy = arma::unique(y);

		for (size_t j = 0; j < uy.n_elem; j++) {
			arma::uvec q1 = arma::find(ux == uy[j]);
			if (!q1.empty()) {
				ux.shed_row(q1(0));
			}
		}

		return ux;
	}	

	void initMissingData(std::shared_ptr<arma::mat> data, std::shared_ptr<arma::umat> missingDataArrayIndexes, std::shared_ptr<arma::uvec> completeCases, bool print )
	{

		const unsigned int nObservations = data->n_rows;
		arma::uvec missingDataVecIdx;

		// Now deal with NANs
		if( data->has_nan() )
		{
			missingDataVecIdx = arma::find_nonfinite(*data);
			(*data)(missingDataVecIdx).fill( arma::datum::nan );  // This makes all the ind values into valid armadillo NANs (should be ok even without, but..)
		}

		// Init the missing data array in a more readable way, i.e. in a row,column format
		if( missingDataVecIdx.n_elem > 0 )
		{
			// create an array of indexes with rows and columns
			(*missingDataArrayIndexes) = arma::umat(missingDataVecIdx.n_elem,2);
			for( unsigned int j=0, n=missingDataVecIdx.n_elem; j<n; ++j)
			{
				(*missingDataArrayIndexes)(j,1) = std::floor( missingDataVecIdx(j) / nObservations ); // this is the corresponding column
				(*missingDataArrayIndexes)(j,0) = missingDataVecIdx(j) - (*missingDataArrayIndexes)(j,1) * nObservations; // this is the row
			}
			(*completeCases) = Utils::arma_setdiff_idx( arma::regspace<arma::uvec>(0, nObservations-1)   , missingDataArrayIndexes->col(0) );
		
		}else{

			(*missingDataArrayIndexes) = arma::umat(0,2);
			(*completeCases) = arma::regspace<arma::uvec>(0, nObservations-1);
		} 

		if( print )
		{
			Rcout << 100. * missingDataVecIdx.n_elem/(double)(nObservations*data->n_cols) <<"% of missing data.." << '\n';
			Rcout << 100. * completeCases->n_elem/(double)(nObservations) <<"% of Complete cases" << '\n';
		}

		return;
	}
/*
	void standardiseData( std::shared_ptr<arma::mat> data, const std::shared_ptr<arma::uvec>& outcomeIndexes, 
						const std::shared_ptr<arma::uvec>& VSPredictorsIndexes, const std::shared_ptr<arma::uvec>& fixedPredictorsIndexes)
	{
		// standardise data
		data->each_col( [](arma::vec& a){ a = ( a - arma::mean(a) ) / arma::stddev(a); } ); // no need to separate between Ys and Xs sincw=e we stdaze each column independently
		// i'll keep the extra function arguments in case something changes

	}
*/

	void formatData( const std::string& dataFileName, const std::string& mrfGFileName, const std::string& blockFileName, const std::string& structureGraphFileName,  SUR_Data& surData )
	{
		
		bool status = readData(dataFileName,surData.data);
		status = status && readGmrf(mrfGFileName,surData.mrfG);
		status = status && readBlocks(blockFileName,surData.blockLabels);
		status = status && readGraph(structureGraphFileName,surData.structureGraph);

		if( !status )
			throw badRead();
		
		removeDisposable(surData.data,surData.blockLabels); // variables indexed by negative indexes are deemed unnecessary by the user, remove them


		getBlockDimensions( surData.blockLabels, surData.structureGraph, surData.data, surData.mrfG, surData.nObservations,
							surData.nOutcomes, surData.outcomesIdx, surData.nPredictors, surData.nVSPredictors, surData.nFixedPredictors,
							surData.VSPredictorsIdx, surData.fixedPredictorsIdx);

		initMissingData( surData.data, surData.missingDataArrayIdx, surData.completeCases, false );

		//standardiseData( surData.data, surData.outcomesIdx, surData.VSPredictorsIdx, surData.fixedPredictorsIdx );

	}

	void readHyperPar(const std::string& hyperParFile, Chain_Data& chainData )
	{
		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_file( hyperParFile.c_str() );

		if ( result )
		{
			// We'll use Xpaths to get variables
			// the results of empty queries will be a nan after the evaluate_number() and thus we'll preserve NANs for non specified parameters
			pugi::xpath_query query_mrf_d("/hyperparameters/mrf_d");
			pugi::xpath_query query_mrf_e("/hyperparameters/mrf_e");

			pugi::xpath_query query_a_sigma("/hyperparameters/a_sigma");
			pugi::xpath_query query_b_sigma("/hyperparameters/b_sigma");

			pugi::xpath_query query_a_tau("/hyperparameters/a_tau");
			pugi::xpath_query query_b_tau("/hyperparameters/b_tau");
			pugi::xpath_query query_nu("/hyperparameters/nu");

			pugi::xpath_query query_a_eta("/hyperparameters/a_eta");
			pugi::xpath_query query_b_eta("/hyperparameters/b_eta");

			pugi::xpath_query query_a_o("/hyperparameters/a_o");
			pugi::xpath_query query_b_o("/hyperparameters/b_o");

			pugi::xpath_query query_a_pi("/hyperparameters/a_pi");
			pugi::xpath_query query_b_pi("/hyperparameters/b_pi");

			pugi::xpath_query query_a_w("/hyperparameters/a_w");
			pugi::xpath_query query_b_w("/hyperparameters/b_w");

            pugi::xpath_query query_a_w0("/hyperparameters/a_w0");
            pugi::xpath_query query_b_w0("/hyperparameters/b_w0");

			// get results
			chainData.mrfD = query_mrf_d.evaluate_number(doc);
			chainData.mrfE = query_mrf_e.evaluate_number(doc);

			chainData.sigmaA = query_a_sigma.evaluate_number(doc);
			chainData.sigmaB = query_b_sigma.evaluate_number(doc);

			chainData.tauA = query_a_tau.evaluate_number(doc);
			chainData.tauB = query_b_tau.evaluate_number(doc);
			chainData.nu = query_nu.evaluate_number(doc);

			chainData.etaA = query_a_eta.evaluate_number(doc);
			chainData.etaB = query_b_eta.evaluate_number(doc);

			chainData.oA  = query_a_o.evaluate_number(doc);
			chainData.oB = query_b_o.evaluate_number(doc);

			chainData.piA  = query_a_pi.evaluate_number(doc);
			chainData.piB = query_b_pi.evaluate_number(doc);

			chainData.wA  = query_a_w.evaluate_number(doc);
			chainData.wB = query_b_w.evaluate_number(doc);

            chainData.w0A  = query_a_w0.evaluate_number(doc);
            chainData.w0B = query_b_w0.evaluate_number(doc);

			std::vector<std::string> valid_top_level = {"hyperparameters"}; // ,"model","chain"}; ?
			std::vector<std::string> valid_hyperpar = {"mrf_d","mrf_e","a_sigma","b_sigma","a_tau","b_tau","nu",
													"a_eta","b_eta","a_o","b_o","a_pi","b_pi","a_w","b_w","a_w0","b_w0"};

			for (pugi::xml_node node = doc.first_child(); node; node = node.next_sibling())
			{
				if ( std::find(valid_top_level.begin(), valid_top_level.end(), node.name() ) == valid_top_level.end() )

				Rcout << "\n\n\tWarning: " << node.name() << " not recognised as a valid top level node - only 'hyperparameters' is valid" << '\n'; // ,model and chain
			}

			for (pugi::xml_node node = doc.child("hyperparameters").first_child(); node; node = node.next_sibling())
			{
				if ( std::find(valid_hyperpar.begin(), valid_hyperpar.end(), node.name() ) == valid_hyperpar.end() )
				{
					Rcout << "\n\n\tWARNING: " << node.name() << " was not recognised as a valid hyperaparameter" << '\n';
					Rcout << "\t"<<node.name() << ": " << node.child_value() << " disregarded .. " << '\n';
					Rcout << "\tValid hyperparameters are: \n\t";
					for( auto& v : valid_hyperpar ) 
						Rcout << v << ", ";
					Rcout << " --- see the documentation for more details " << '\n' << '\n';
					
				}
			}

		}
		else{
			Rcout << '\n' << "No hyperparameter input file was given (or wrong format detected), so default values will be used." << '\n';
		}

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
