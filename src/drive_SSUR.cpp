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

#include "ESS_Sampler.h"

extern omp_lock_t RNGlock; //defined in global.h
extern std::vector<std::mt19937_64> rng;

int main(int argc, char *  argv[])
{
	omp_init_lock(&RNGlock);  // RNG lock for the parallel part

	unsigned int nIter = 10; // default number of iterations
	unsigned int s=1,p=1;      // might read them from a meta-data file, but for the moment is easier like this..
	unsigned int nChains = 1;

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
	long long int seed = std::chrono::system_clock::now().time_since_epoch().count();

	// seed all the engines
	for(unsigned int i=0; i<nThreads; ++i)
	{
		rng[i] = std::mt19937_64(seed + i*(1000*(p*s*3+s*s)*nIter) );
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
	inFile.erase(inFile.end()-4,inFile.end());  // remomve the .txt from inFile !


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

	ESS_Sampler sampler( Y, X , nChains , "SSUR");

	// ****************************************

	// INIT THE FILE OUTPUT
	// Output to file the current state
	arma::umat gamma_out = sampler[0] -> getGamma(); // out var for the gammas
	arma::umat g_out = arma::umat( sampler[0] -> getGAdjMat() ); // out var for G

	arma::mat beta_out = sampler[0] -> getBeta(); // out var for the betas
	arma::mat sigmaRho_out = sampler[0] -> getSigmaRho(); // out var for the sigmas and rhos

	gammaOutFile.open( outFilePath+inFile+"_SSUR_gamma_out.txt" , std::ios_base::trunc);
	gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out)) << std::flush;
	gammaOutFile.close();

	gOutFile.open( outFilePath+inFile+"_SSUR_G_out.txt" , std::ios_base::trunc);
	gOutFile << ( arma::conv_to<arma::mat>::from(g_out) ) << std::flush;   // this might be quite long...
	gOutFile.close();

	logPOutFile << 	sampler[0] -> getLogPTau() << " " <<
					sampler[0] -> getLogPEta() <<  " " <<
					sampler[0] -> getLogPJT() <<  " " <<
					sampler[0] -> getLogPSigmaRho() <<  " " <<
					sampler[0] -> getLogPO() <<  " " <<
					sampler[0] -> getLogPPi() <<  " " <<
					sampler[0] -> getLogPGamma() <<  " " <<
					sampler[0] -> getLogPW() <<  " " <<
					sampler[0] -> getLogPBeta() <<  " " <<
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

		sampler.step();
		
		// #################### END LOCAL MOVES

		// ## Global moves
		// *** end Global move's section

		// UPDATE OUTPUT STATE
		gamma_out += sampler[0] -> getGamma(); // the result of the whole procedure is now my new mcmc point, so add that up
		g_out += arma::umat( sampler[0] -> getGAdjMat() );

		beta_out += sampler[0] -> getBeta();
		sigmaRho_out += sampler[0] -> getSigmaRho();

		// Print something on how the chain is going
		if( (i+1) % tick == 0 )
		{

			std::cout << " Running iteration " << i+1 << " ... local Acc Rate: ~ gamma: " << Utils::round( sampler[0] -> getGammaAccRate() , 3 ) 
					<< " -- JT: " << Utils::round( sampler[0] -> getJTAccRate() , 3 ) ;
			if( nChains > 1)
				std::cout << " -- Global: " << Utils::round( sampler.getGlobalAccRate() , 3 ) << std::endl; 
			else
				std::cout << std::endl;
				
			// Output to files every now and then
			if( (i+1) % (tick*10) == 0 )
			{

				gammaOutFile.open( outFilePath+inFile+"_SSUR_gamma_out.txt" , std::ios_base::trunc);
				gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/((double)i+1.0) << std::flush;
				gammaOutFile.close();

				gOutFile.open( outFilePath+inFile+"_SSUR_G_out.txt" , std::ios_base::trunc);
				gOutFile << ( arma::conv_to<arma::mat>::from(g_out) )/((double)i+1.0) << std::flush;   // this might be quite long...
				gOutFile.close();

				logPOutFile << 	sampler[0] -> getLogPTau() << " " <<
								sampler[0] -> getLogPEta() <<  " " <<
								sampler[0] -> getLogPJT() <<  " " <<
								sampler[0] -> getLogPSigmaRho() <<  " " <<
								sampler[0] -> getLogPO() <<  " " <<
								sampler[0] -> getLogPPi() <<  " " <<
								sampler[0] -> getLogPGamma() <<  " " <<
								sampler[0] -> getLogPW() <<  " " <<
								sampler[0] -> getLogPBeta() <<  " " <<
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

	logPOutFile << 	sampler[0] -> getLogPTau() << " " <<
					sampler[0] -> getLogPEta() <<  " " <<
					sampler[0] -> getLogPJT() <<  " " <<
					sampler[0] -> getLogPSigmaRho() <<  " " <<
					sampler[0] -> getLogPO() <<  " " <<
					sampler[0] -> getLogPPi() <<  " " <<
					sampler[0] -> getLogPGamma() <<  " " <<
					sampler[0] -> getLogPW() <<  " " <<
					sampler[0] -> getLogPBeta() <<  " " <<
				std::endl << std::flush;

	// ----
	beta_out = beta_out/((double)nIter);
	beta_out.save(outFilePath+inFile+"_SSUR_beta_out.txt",arma::raw_ascii);

	sigmaRho_out = sigmaRho_out/((double)nIter);
	sigmaRho_out.save(outFilePath+inFile+"_SSUR_sigmaRho_out.txt",arma::raw_ascii);
	// -----

	std::cout << "Saved to :   "+outFilePath+inFile+"_SSUR_****_out.txt" << std::endl;
	std::cout << "Final tau : " << sampler[0] -> getTau() << "    w/ proposal variance: " << sampler[0] -> getVarTauProposal() << std::endl;
	std::cout << "Final eta : " << sampler[0] -> getEta() <<  std::endl;
	std::cout << "Final o : " << sampler[0] -> getO().t() << "       w/ proposal variance: " << sampler[0] -> getVarOProposal() << std::endl;  
	std::cout << "Final pi : " << sampler[0] -> getPi().t() << "       w/ proposal variance: " << sampler[0] -> getVarPiProposal() << std::endl;
	std::cout << "  -- Average Omega : " << arma::accu( sampler[0] -> getO() * sampler[0] -> getPi().t() )/((double)p*s) <<  std::endl;
	std::cout << "Final w : " << sampler[0] -> getW() <<  std::endl;
	std::cout << "Final temperature ratio : " << sampler[1]->getTemperature() <<  std::endl << std::endl ;

	// Exit

	std::cout << "DONE, exiting! " << std::endl << std::endl ;
	return 0;
}
