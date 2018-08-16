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

#include "junction_tree.h"
#include "SSUR_Chain.h"

extern omp_lock_t RNGlock; //defined in global.h
extern std::vector<std::mt19937_64> rng;

int main(int argc, char *  argv[])
{
	omp_init_lock(&RNGlock);  // RNG lock for the parallel part

	unsigned int nIter = 10; // default number of iterations
	unsigned int s=1,p=1;      // might read them from a meta-data file, but for the moment is easier like this..
	unsigned int nChains = 1;
	unsigned int M = 20;

	std::string inFile = "data.txt";
	std::string outFilePath = "";
	std::string omegaInitPath = "";

    // ### Read and interpret command line (to put in a separate file / function?)
    int na = 1;
    while(na < argc)
    {
		if ( 0 == strcmp(argv[na],"--nIter") )
		{
			nIter = std::stoi(argv[++na]);
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
		else if ( 0 == strcmp(argv[na],"--M") )
		{
			M = std::stoi(argv[++na]); // use the next
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
		else
    {
	    std::cout << "Unknown option: " << argv[na] << std::endl;
	    return(1); //this is exit if I'm in a function elsewhere
    }
    }//end reading from command line

	std::cout << "Init RNG engine .. " << std::endl;

	// ############# Init the RNG generator/engine
	std::random_device r;
	unsigned int nThreads = omp_get_max_threads();

	rng.reserve(nThreads);  // reserve the correct space for the vector of rng engines
	std::seed_seq seedSeq;	// and declare the seedSequence
	std::vector<unsigned int> seedInit(8);
	long long int seed = r();

	// seed all the engines
	for(unsigned int i=0; i<nThreads; ++i)
	{
		rng[i].seed(seed + i*(1000*(p*s*3+s*s)*M*nIter)); // 1000 rng per (p*s*3+s*s) variables each loop .. is this...ok? is random better?
	}

	// ############

	// ### Read the data
	unsigned int n;
	arma::mat Y, X;
	std::cout << "Trying to read data ...  " << std::flush;

	if( Utils::readData(inFile, s, p, n, Y, X) ){
		std::cout << "Reading successfull!" << std::endl;
	}else{
		std::cout << "OUCH! EXITING --- " << std::endl;
		return 1;
	}

	// The intercept columnto X will be inserted when initialising the chain

	std::cout << "Clearing and initialising output files " << std::endl;
	// Re-define inFile so that I can use it in the output
	std::size_t slash = inFile.find("/");  // remove the path from inFile
	while( slash != std::string::npos )
	{
		inFile.erase(inFile.begin(),inFile.begin()+slash+1);
		slash = inFile.find("/");
	}
	inFile.erase(inFile.end()-4,inFile.end());  // remomve the .txt from inFile


	// clear the content of previous files
	std::ofstream fileClear;

	// open new files in append mode
	fileClear.open(outFilePath+inFile+"_SSUR_logP_out.txt", std::ios::out | std::ios::trunc); fileClear.close(); // clear contents
	std::ofstream logPOutFile; logPOutFile.open( outFilePath+inFile+"_SSUR_logP_out.txt" , std::ios_base::app);

	// open avg files in trunc mode to cut previous content
	std::ofstream gammaOutFile; gammaOutFile.open( outFilePath+inFile+"_SSUR_gamma_out.txt" , std::ios_base::trunc); gammaOutFile.close();
	std::ofstream gOutFile; gOutFile.open( outFilePath+inFile+"_SSUR_G_out.txt" , std::ios_base::trunc); gammaOutFile.close();


	// ****************************************
	// **********  INIT THE CHAIN *************
	// ****************************************
	std::cout << "Initialising the MCMC Chain " << std::endl;

	SSUR_Chain chain( Y, X , 1. );

	// ****************************************

	// INIT THE FILE OUTPUT
	// Output to file the current state
	arma::umat gamma_out = chain.getGamma(); // out var for the gammas
	arma::umat g_out = arma::umat( chain.getGAdjMat() ); // out var for G

	arma::mat beta_out = chain.getBeta(); // out var for the betas
	arma::mat sigmaRho_out = chain.getSigmaRho(); // out var for the sigmas and rhos

	gammaOutFile.open( outFilePath+inFile+"_SSUR_gamma_out.txt" , std::ios_base::trunc);
	gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out)) << std::flush;
	gammaOutFile.close();

	gOutFile.open( outFilePath+inFile+"_SSUR_G_out.txt" , std::ios_base::trunc);
	gOutFile << ( arma::conv_to<arma::mat>::from(g_out) ) << std::flush;   // this might be quite long...
	gOutFile.close();

	logPOutFile << 	chain.getLogPTau() << " " <<
					chain.getLogPEta() <<  " " <<
					chain.getLogPJT() <<  " " <<
					chain.getLogPSigmaRho() <<  " " <<
					chain.getLogPO() <<  " " <<
					chain.getLogPPi() <<  " " <<
					chain.getLogPGamma() <<  " " <<
					chain.getLogPW() <<  " " <<
					chain.getLogPBeta() <<  " " <<
				std::endl << std::flush;
					

	// ########
	// ########
	// ######## Start
	// ########
	// ########

	std::cout << "Starting "<< nChains <<" (parallel) chain(s) for " << nIter << " iterations:" << std::endl << std::flush;

	unsigned int tick = 100; // how may iter for each print?

	for(unsigned int i=1; i < nIter ; ++i)
	{

        // #pragma omp parallel for num_threads(nThreads)
        // for(unsigned int m=0; m<nChains ; ++m)
        // {
				chain.step();

        // } // end parallel updates
		
		// #################### END LOCAL MOVES

		// ## Global moves
		// *** end Global move's section

		// UPDATE OUTPUT STATE
		gamma_out += chain.getGamma(); // the result of the whole procedure is now my new mcmc point, so add that up
		g_out += arma::umat( chain.getGAdjMat() );

		beta_out += chain.getBeta();
		sigmaRho_out += chain.getSigmaRho();

		// Print something on how the chain is going
		if( (i+1) % tick == 0 )
		{

			std::cout << " Running iteration " << i+1 << " ... local Acc Rate: ~ gamma: " << Utils::round( chain.getGammaAccRate() , 3 ) << " -- JT: " << Utils::round( chain.getJTAccRate() , 3 ) << std::endl; 
			// if( nChains > 1 )
				// std::cout << "\033[A" << "\033[105C" << " - global Acc Rate ~ " << 1 << std::endl;;

			// Output to files every now and then
			if( (i+1) % (tick*10) == 0 )
			{

				gammaOutFile.open( outFilePath+inFile+"_SSUR_gamma_out.txt" , std::ios_base::trunc);
				gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/((double)i+1.0) << std::flush;
				gammaOutFile.close();

				gOutFile.open( outFilePath+inFile+"_SSUR_G_out.txt" , std::ios_base::trunc);
				gOutFile << ( arma::conv_to<arma::mat>::from(g_out) )/((double)i+1.0) << std::flush;   // this might be quite long...
				gOutFile.close();

				logPOutFile << 	chain.getLogPTau() << " " <<
								chain.getLogPEta() <<  " " <<
								chain.getLogPJT() <<  " " <<
								chain.getLogPSigmaRho() <<  " " <<
								chain.getLogPO() <<  " " <<
								chain.getLogPPi() <<  " " <<
								chain.getLogPGamma() <<  " " <<
								chain.getLogPW() <<  " " <<
								chain.getLogPBeta() <<  " " <<
							std::endl << std::flush;
			}

		}

	} // end MCMC


	// Print the end
	std::cout << " MCMC ends. " /* << " Final temperature ratio ~ " << temperatureRatio  */<< "  --- Saving results and exiting" << std::endl;

	// ### Collect results and save them
	gammaOutFile.open( outFilePath+inFile+"_SSUR_gamma_out.txt" , std::ios_base::trunc);
	gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/((double)nIter+1.0) << std::flush;
	gammaOutFile.close();

	gOutFile.open( outFilePath+inFile+"_SSUR_G_out.txt" , std::ios_base::trunc);
	gOutFile << ( arma::conv_to<arma::mat>::from(g_out) )/((double)nIter+1.0) << std::flush;   // this might be quite long...
	gOutFile.close();

	logPOutFile << 	chain.getLogPTau() << " " <<
					chain.getLogPEta() <<  " " <<
					chain.getLogPJT() <<  " " <<
					chain.getLogPSigmaRho() <<  " " <<
					chain.getLogPO() <<  " " <<
					chain.getLogPPi() <<  " " <<
					chain.getLogPGamma() <<  " " <<
					chain.getLogPW() <<  " " <<
					chain.getLogPBeta() <<  " " <<
				std::endl << std::flush;

	// ----
	beta_out = beta_out/((double)nIter);
	beta_out.save(outFilePath+inFile+"_SSUR_beta_out.txt",arma::raw_ascii);

	sigmaRho_out = sigmaRho_out/((double)nIter);
	sigmaRho_out.save(outFilePath+inFile+"_SSUR_sigmaRho_out.txt",arma::raw_ascii);
	// -----

	std::cout << "Saved to :   "+outFilePath+inFile+"_SSUR_****_out.txt" << std::endl;
	std::cout << "Final tau : " << chain.getTau() << "    w/ proposal variance: " << chain.getVarTauProposal() << std::endl;
	std::cout << "Final eta : " << chain.getEta() <<  std::endl;
	std::cout << "Final o : " << chain.getO().t() << "       w/ proposal variance: " << chain.getVarOProposal() << std::endl;  
	std::cout << "Final pi : " << chain.getPi().t() << "       w/ proposal variance: " << chain.getVarPiProposal() << std::endl;
	std::cout << "  -- Average Omega : " << arma::accu( chain.getO() * chain.getPi().t() )/((double)p*s) <<  std::endl;
	std::cout << "Final w : " << chain.getW() <<  std::endl << std::endl ;

	// Exit

	std::cout << "DONE, exiting! " << std::endl << std::endl ;
	return 0;
}
