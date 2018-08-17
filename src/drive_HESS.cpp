#include <vector>
#include <iostream>
#include <string>
#include <armadillo>
#include <cmath>
#include <limits>
#include <omp.h>

#include "global.h"
#include "utils.h"
#include "distr.h"
#include "HESS.h"

extern omp_lock_t RNGlock; //defined in global.h
extern std::vector<std::mt19937_64> rng;


int main(int argc, char *  argv[])
{

	omp_init_lock(&RNGlock);  // RNG lock for the parallel part

	omp_set_lock(&RNGlock);
    // one thread at a time stuff
    omp_unset_lock(&RNGlock);

	unsigned int nIter = 10; // default number of iterations
	unsigned int s=1,p=1;      // might read them from a meta-data file, but for the moment is easier like this..
	unsigned int nChains = 1;
	double deltaTempRatio = 1.;

	std::string inFile = "data.txt";
	std::string outFilePath = "";
	std::string omegaInitPath = "";

	unsigned int method = 1; // Defaul is our novel "bandit" method

	/*
	0: MC^3 -- BASE algorithm, simple randow walker with add-delete and swap move
	1: Bandit -- Novel method
	*/

    // ### Read and interpret command line (to put in a separate file / function?)
    int na = 1;
    while(na < argc)
    {
		if ( 0 == strcmp(argv[na],"--method") || 0 == strcmp(argv[na],"--algorithm") || 0 == strcmp(argv[na],"--algo")  )
		{
			method = std::stoi(argv[++na]);
			if(method > 1 || method < 0 )
			{
				std::cout << "Invalid method argument ("<<method<<"), see README.md\nDefaulting to bandit sampler\n"<<std::flush;
				method = 1;
			}
			if (na+1==argc) break;
			++na;
		}
		else if ( 0 == strcmp(argv[na],"--nIter") )
		{
			nIter = atoi(argv[++na]);
			if (na+1==argc) break;
			++na;
		}
		else if ( 0 == strcmp(argv[na],"--nOutcomes") )
		{
			s = std::stoi(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == strcmp(argv[na],"--nPredictors") )
		{
			p = std::stoi(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == strcmp(argv[na],"--nChains") )
		{
			nChains = std::stoi(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		/* DEPRECATED, ALWAYS INITIALISE WITH A SEQUENCE AT RANDOM, MOREOVER I NEED MULTIPLE ENGINES NOW.. */
		// else if ( 0 == strcmp(argv[na],"--seed") )
		// {
		// 	seed = atoi(argv[++na]); // use the next
		// 	if( seed < 0 )
		// 	{
		// 		seed = std::time(NULL);
		// 	}
		// 	if (na+1==argc) break; // in case it's last, break
		// 	++na; // otherwise augment counter
		// }
		else if ( 0 == strcmp(argv[na],"--inFile") )
		{
			inFile = ""+std::string(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == strcmp(argv[na],"--outFilePath") )
		{
			outFilePath = std::string(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == strcmp(argv[na],"--omegaInitPath") )
		{
			omegaInitPath = ""+std::string(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == strcmp(argv[na],"--deltaTempRatio") )
		{
			deltaTempRatio = std::stof(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else
        {
          	    std::cout << "Unknown option: " << argv[na] << std::endl;
          	    return(1); //this is exit if I'm in a function elsewhere
        }
    }//end reading from command line


	// ### Read the data
	unsigned int n;
	arma::mat Y, X;

	if( Utils::readData(inFile, s, p, n, Y, X) ){
		std::cout << "Reading successfull!" << std::endl;
	}else{
		std::cout << "OUCH! EXITING --- " << std::endl;
		return 1;
	}

	// if( Utils::readDataSEM(inFile, data, blockIdx, varType, NAIdx,
	// 	s, p, n) ){
	// 	std::cout << "Reading successfull!" << std::endl;
	// }else{
	// 	std::cout << "OUCH! EXITING --- " << std::endl;
	// 	return 1;
	// }


	// ############# Init the RNG generator/engine
	std::random_device r;
	unsigned int nThreads = omp_get_max_threads();
	if(nChains < nThreads)
		nThreads = nChains;
	omp_set_num_threads(nThreads);

	rng.reserve(nThreads);  // reserve the correct space for the vector of rng engines
	std::seed_seq seedSeq;	// and declare the seedSequence
	std::vector<unsigned int> seedInit(8);
	long long int seed = r();

	// seed all the engines
	for(unsigned int i=0; i<nThreads; ++i)
	{
		// seedInit = {r(), r(), r(), r(), r(), r(), r(), r()};
		// seedSeq.generate(seedInit.begin(), seedInit.end()); // init with a sequence of 8 TRUE random values
		rng[i].seed(seed + i*(1000*(p*s*3+s*s)*nIter)); // 1000 rng per (p*s*3+s*s) variables each loop .. is this...ok? is random better?
		// rng[i].seed(r());
	}


	// ############################
	// HESS is the model where we have a (possibly) multivariate response with diagonal covariance matrix but we do VS separately for each outcome
	// Beta is a vector/matrix of (1+p)xq coefficients, R is a diagonal qxq matrix (or a qx1 vector of sigma coeffs) and gamma is a pxq matrix
	// We sample using a MH step for gamma while Beta and R are integrated out
	// ############################

	// ### Set hyper-parameters of distributions involved

	// ## Set Prior Distributions  -- (and init values)

	// # R -- diagonal matrix with elements coming from IG(a_r_0,b_r_0)
	double a_r_0 = 0.1;
	double b_r_0 = 0.01;

	// # Beta s --- (Alpha s, intercepts, are included here even though there's no variable selection on them)
	// Beta comes from a Matrix Normal distribution with parameters MN(mean_beta_0,W_0,R)

	double m_beta_0 = 0;  // this is common for all the coefficients
	double sigma2_beta_0 = 10;

	arma::mat mean_beta_0(p+1, s); mean_beta_0.fill(m_beta_0);
	arma::mat W_0 = arma::eye( p+1, p+1 ) * sigma2_beta_0;

	// Add to X the intercept
	X.insert_cols(0, arma::ones<arma::vec>(n) );

	// XtX covariates only
	arma::mat XtX = X.t() * X; // this is needed for crossover for example
	arma::mat covariatesCorrelation = arma::inv( arma::diagmat( arma::sqrt(XtX.diag()) ) ) * XtX * arma::inv( arma::diagmat( arma::sqrt(XtX.diag()) ) );

	// # gamma  (p+1xs elements, but the first row doesn't move, always 1, cause are the intercepts)
	// gamma_jk (j=1,..p  i.e. not 0 - k=0,s-1) comes from a Bernoulli distribution of parameters omega_j, who are all coming from a Beta(a_0,b_0)
	arma::vec a_0(p); a_0.fill( 5. ); 		// the average of a beta is E[p]=a/(a+b), this way E[p] <= 1.% FOR EVERY OUTCOME
	arma::vec b_0(p); b_0.fill( std::max( (double)p , 500.) - 5. );

	// ### Initialise and Start the chain

	// ## Initialise chain traces -- INIT HERE
	arma::mat omega_init(p,s);
	arma::umat gamma_init(p,s);

	if( omegaInitPath == "" )
	{
		omega_init = arma::ones<arma::mat>(p,s)/(double)p; // This might be read from file
		gamma_init = arma::zeros<arma::umat>(p,s);
	}else{
		omega_init.load(omegaInitPath,arma::raw_ascii);
		for(unsigned int j=0; j<p; ++j)
			for(unsigned int k=0; k<s; ++k)
				gamma_init(j,k) = Distributions::randBernoulli(omega_init(j,k));
	}

	// arma::cube mcmc_omega(p,s,nIter); mcmc_omega.slice(0) = omega_init; //inits
	// arma::ucube mcmc_gamma(p,s,nIter); mcmc_gamma.slice(0) = gamma_init;

	arma::vec mcmc_logPrior(nIter); mcmc_logPrior(0) = Model::logPrior(omega_init, gamma_init, a_0, b_0);
	arma::vec mcmc_logLik(nIter); mcmc_logLik(0) = Model::logLikelihood(Y, X, gamma_init, a_r_0, b_r_0, W_0, 1.);

	// ## Defines proposal parameters and temporary variables for the MCMC
	double accCount = 0.; 
	arma::vec accCount_tmp = arma::zeros<arma::vec>(nChains);

	// arma::mat omega_curr;				//current state of the main chain not needed
	// arma::umat gamma_curr;
	// double logLik_curr,logPrior_curr;

	arma::cube omega_state(p,s,nChains); 	// Current state of ALL the chains
	arma::ucube gamma_state(p,s,nChains);	// no need to keep track of anything else for all the non-main chains
	arma::vec logPrior_state(nChains);
	arma::vec logLik_state(nChains);

	for(unsigned int m=0; m<nChains ; ++m)
	{
		omega_state.slice(m) = omega_init;
		gamma_state.slice(m) = gamma_init;
		logPrior_state(m) = mcmc_logPrior(0);
		logLik_state(m) = mcmc_logLik(0);
	}

	// std::cout << std::endl<< std::endl<< Model::logLikelihood(Y, X, gamma_init, a_r_0, b_r_0, W_0, 1.)<< std::endl;


	// # Define temperture ladder
	double maxTemperature = arma::min(a_0) - 2.; // due to the tempered gibbs moves
	double temperatureRatio = 2.; // std::max( 100., (double)n );

	arma::vec temperature(nChains);
	temperature(0) = 1.;

	for(unsigned int m=1; m<nChains; ++m)
		temperature(m) = std::min( maxTemperature, temperature(m-1)*temperatureRatio );

	for(unsigned int m=0; m<nChains ; ++m)
	{
		logLik_state(m) = mcmc_logLik(0)/temperature(m);
	}

	unsigned int nGlobalUpdates = floor( nChains/2 );

	unsigned int countGlobalUpdates = 0;    // count the global updates that happen on the first chain
	double accCountGlobalUpdates = 0.;		// count the accepted ones

	unsigned int globalType;


	// OUTPUT STATE


	// Re-define inFile so that I can use it in the output
	std::size_t slash = inFile.find("/");  // remove the path from inFile
	while( slash != std::string::npos )
	{
		inFile.erase(inFile.begin(),inFile.begin()+slash+1);
		slash = inFile.find("/");
		// std::cout << inFile << std::endl;
	}
	inFile.erase(inFile.end()-4,inFile.end());  // remomve the .txt from inFile


	// open new files in append mode
	// std::ofstream omegaOutFile; omegaOutFile.open( outFilePath+inFile+"_SSUR_omega_out.txt" , std::ios_base::trunc); omegaOutFile.close();
	std::ofstream gammaOutFile; gammaOutFile.open( outFilePath+inFile+"_HESS_gamma_out.txt" , std::ios_base::trunc); gammaOutFile.close();

	// Output to file the current state
	arma::umat gamma_out = gamma_state.slice(0); // out var for the gammas
	arma::mat omega_out(omega_state.slice(0)); // out var for the omegas


	// could use save btI need to sum and then normalise so I'd need to store another matrix for each...

	gammaOutFile.open( outFilePath+inFile+"_HESS_gamma_out.txt" , std::ios_base::trunc);
	gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out)) << std::flush;
	gammaOutFile.close();

	// All-exchange operator hyper pars

	// Crossover operator hyper pars (see http://www3.stat.sinica.edu.tw/statistica/oldpdf/A10n21.pdf
	double pXO_0 = 0.1, pXO_1 = 0.2 , pXO_2 = 0.2;
	double p11 = pXO_0*pXO_0 + (1.-pXO_0)*(1.-pXO_0) ,p12 = 2.*pXO_0*(1.-pXO_0) ,p21= pXO_1*(1.-pXO_2) + pXO_2*(1.-pXO_1) ,p22 = pXO_1*pXO_2 + (1.-pXO_1)*(1.-pXO_2);


	// NON-BANDIT related things
	unsigned int nUpdates = p/10; //arbitrary nunmber, should I use something different?

	// BANDIT ONLY SECTION
	arma::cube alpha_z;	arma::cube beta_z;	arma::cube zeta; std::vector<arma::vec> mismatch;
	std::vector<arma::vec> normalised_mismatch; std::vector<arma::vec> normalised_mismatch_backwards;

	if( method == 1)
	{
		nUpdates = 4; // for Bandit this must be SMALL (as it scales with nUpdates! and has a different meaning anyway)

		// Starting vaues for the Bandit tuning parameters
		// stating proposals are beta( 0.5 , 0.5 ) so that they're centered in 0.5 with spikes at the extremes
		// ALL OF THESE NEED TO HAVE A COPY FOR EACH CHAIN!!
		alpha_z = arma::cube(p,s,nChains); alpha_z.fill(0.5);
		beta_z = arma::cube (p,s,nChains); beta_z.fill(0.5);

		// these need to fill, they'll be overwritten anyway
		zeta = arma::cube(p,s,nChains);
		mismatch = std::vector<arma::vec>(nChains);
		normalised_mismatch = std::vector<arma::vec>(nChains);
		normalised_mismatch_backwards = std::vector<arma::vec>(nChains);
		for(unsigned int i = 0; i<nChains; ++i)
		{
			mismatch[i] = arma::vec(p*s);
			normalised_mismatch[i] = arma::vec(p*s);
			normalised_mismatch_backwards[i] = arma::vec(p*s);
		}
		// I still wonder why .slice() is an instance of arma::mat and return a reference to it
		// while .col() is a subview and has hald the method associated with a Col object ...
	}
	// END Bandit only section

	// ###########################################################
	// ###########################################################
	// ## Start the MCMC
	// ###########################################################
	// ###########################################################

	std::cout << "Starting "<< nChains <<" (parallel) chain(s) for " << nIter << " iterations:" << std::endl;

	for(unsigned int i=1; i < nIter ; ++i)
	{

		switch(method){

			case 0:
					#pragma omp parallel for num_threads(nThreads)
					for(unsigned int m=0; m<nChains ; ++m)
					{
							Model::MC3_SUR_MCMC_step(Y,X, omega_state.slice(m),gamma_state.slice(m),logPrior_state(m),logLik_state(m),
											a_r_0, b_r_0, W_0, a_0, b_0, accCount_tmp(m), nUpdates, temperature(m)); // in ESS accCount could be a vec, one for each chain
					}// end parallel updates
					break;

			case 1:
					#pragma omp parallel for num_threads(nThreads)
					for(unsigned int m=0; m<nChains ; ++m)
					{
							Model::bandit_SUR_MCMC_step(Y,X, omega_state.slice(m),gamma_state.slice(m),logPrior_state(m),logLik_state(m),
											a_r_0, b_r_0, W_0, a_0, b_0, accCount_tmp(m), nUpdates, temperature(m) ,
											zeta.slice(m), alpha_z.slice(m), beta_z.slice(m),
											mismatch[m], normalised_mismatch[m], normalised_mismatch_backwards[m]); // in ESS accCount could be a vec, one for each chain
					}// end parallel updates
					break;

					// there is no default since by default method = 2 and hence this last one is the default anyway
		}


		// UPDATE OUTPUT STATE
		gamma_out += gamma_state.slice(0); // the result of the whole procedure is now my new mcmc point, so add that up
		omega_out += omega_state.slice(0);


		// mcmc_omega.slice(i) = omega_state.slice(0); mcmc_gamma.slice(i) = gamma_state.slice(0); 
		mcmc_logPrior(i) = logPrior_state(0); mcmc_logLik(i) = logLik_state(0);

		// ####################
		// ## Global moves
		if( nChains > 1 )
		{
			for(unsigned int k=0; k < nGlobalUpdates ; ++k)  // repeat global updates
			{

				// # Global move
				// Select the type of exchange/crossOver to apply

				globalType = Distributions::randIntUniform(0,6);   // (nChains-1)*(nChains-2)/2 is the number of possible chain combinations with nChains

				switch(globalType){

					case 0: break;

					// -- Exchange
					case 1: Model::exchange_step(gamma_state, omega_state, logPrior_state, logLik_state,
							temperature, nChains, nGlobalUpdates, accCountGlobalUpdates, countGlobalUpdates);
							break;

					case 2: Model::nearExchange_step(gamma_state, omega_state, logPrior_state, logLik_state,
							temperature, nChains, nGlobalUpdates, accCountGlobalUpdates, countGlobalUpdates);
							break;

					case 3: Model::allExchange_step(gamma_state, omega_state, logPrior_state, logLik_state,
							temperature, nChains, nGlobalUpdates, accCountGlobalUpdates, countGlobalUpdates);
							break;

					// -- CrossOver
					case 4: Model::adapCrossOver_step(gamma_state, omega_state, logPrior_state, logLik_state,
								a_r_0, b_r_0, W_0, a_0, b_0, Y, X, temperature, 
								pXO_0, pXO_1, pXO_2, p11, p12, p21, p22,
								nChains, nGlobalUpdates, accCountGlobalUpdates, countGlobalUpdates);
							break;

					case 5: Model::uniformCrossOver_step(gamma_state, omega_state, logPrior_state, logLik_state,
								a_r_0, b_r_0, W_0, a_0, b_0, Y, X, temperature,
								nChains, nGlobalUpdates, accCountGlobalUpdates, countGlobalUpdates);
							break;

					case 6: Model::blockCrossOver_step(gamma_state, omega_state, logPrior_state, logLik_state,
								a_r_0, b_r_0, W_0, a_0, b_0, Y, X, covariatesCorrelation, temperature, 0.25,
								nChains, nGlobalUpdates, accCountGlobalUpdates, countGlobalUpdates);
							break;
				}


			} // end "K" Global Moves

			// Update the output vars
			// mcmc_omega.slice(i) = omega_state.slice(0);
			// mcmc_gamma.slice(i) = gamma_state.slice(0);
			mcmc_logPrior(i) = logPrior_state(0);
			mcmc_logLik(i) = logLik_state(0);

			// ## Update temperature ladder
			if ( i % 10 == 0 )
			{
				if( (accCountGlobalUpdates / (double)countGlobalUpdates) > 0.35 )
				{
					temperatureRatio = std::max( 2. , temperatureRatio-deltaTempRatio );
				}else{
					temperatureRatio += deltaTempRatio;
				}

				temperatureRatio = std::min( temperatureRatio , pow( maxTemperature, 1./( (double)nChains - 1.) ) );

				// std::cout << "---------------------------------- \n"<< temperature << '\n'<< std::flush;
				for(unsigned int m=1; m<nChains; ++m)
				{
					// untempered lik and prior
					logLik_state(m) = logLik_state(m)*temperature(m);
					logPrior_state(m) = logPrior_state(m)*temperature(m);

					// std::cout << '\n'<< temperature(m-1) << " " << temperature(m-1)*temperatureRatio  << '\n' << '\n'<< std::flush;
					temperature(m) = std::min( maxTemperature, temperature(m-1)*temperatureRatio );

					// re-tempered lik and prior
					logLik_state(m) = logLik_state(m)/temperature(m);
					logPrior_state(m) = logPrior_state(m)/temperature(m);

				}
			}


		} //end Global move's section

		// Print something on how the chain is going
		if( (i+1) % 100 == 0 )
		{
			// Update Acc Rate only for main chain
			accCount = accCount_tmp(0)/nUpdates;

			std::cout << " Running iteration " << i+1 << " ... loc.acc.rate ~ " << accCount/(double)i;
			if( nChains > 1 )
				std::cout << " global.acc.rate ~ " << accCountGlobalUpdates / (double)countGlobalUpdates;
			std::cout << /*"  ~~  " << temperature.t() << */ std::endl;
		}

		// Output to files every now and then
		if( (i+1) % (1000) == 0 )
		{
			gammaOutFile.open( outFilePath+inFile+"_HESS_gamma_out.txt" , std::ios_base::trunc);
			gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/((double)i+1.0) << std::flush;
			gammaOutFile.close();
		}

	} // end MCMC


	std::cout << " MCMC ends. Final temperature ratio ~ " << temperatureRatio << "   --- Saving results and exiting" << std::endl;

	// ### Collect results and save them

	// std::size_t slash = inFile.find("/");  // remove the path from inFile
	// while( slash != std::string::npos )
	// {
	// 	inFile.erase(inFile.begin(),inFile.begin()+slash+1);
	// 	slash = inFile.find("/");
	// 	// std::cout << inFile << std::endl;
	// }
	// inFile.erase(inFile.end()-4,inFile.end());  // remomve the .txt from inFile
	// 	// std::cout << inFile << std::endl;

	// mcmc_omega.save(outFilePath+inFile+"_HESS_omega_out.txt",arma::raw_ascii);
	// mcmc_gamma.save(outFilePath+inFile+"_HESS_gamma_out.txt",arma::raw_ascii);
	mcmc_logPrior += mcmc_logLik; mcmc_logPrior.save(outFilePath+inFile+"_HESS_logP_out.txt",arma::raw_ascii);

	omega_out = omega_out/((double)nIter);
	omega_out.save(outFilePath+inFile+"_HESS_omega_out.txt",arma::raw_ascii);

	// Exit
	return 0;
}
