#include "drive.h"

extern omp_lock_t RNGlock; //defined in global.h
extern std::vector<std::mt19937_64> rng;

int drive_SSUR( arma::mat& Y , arma::mat& X , unsigned int& nChains , unsigned int& nIter , 
				std::string& inFile , std::string& outFilePath , std::string& gammaSampler , bool usingGPrior )
{

	// ****************************************
	// **********  INIT THE CHAIN *************
	// ****************************************
	std::cout << "Initialising the MCMC Chain " << std::endl;

	ESS_Sampler<SSUR_Chain> sampler( Y , X , nChains );// this is thus also some sort of default
		// although note that you won't pass the input phase with a different method string

	// *****************************
	// need to use getX because I need the intercept
	arma::mat Q,R; arma::qr(Q,R, *sampler[0]->getX() );
	arma::mat betaInit = arma::solve(R,arma::trans(Q) * Y );
	arma::umat gammaInit = betaInit > 0.5*arma::stddev(arma::vectorise(betaInit));
	gammaInit.shed_row(0);

	sampler[0] -> gammaInit( gammaInit );
    sampler[0] -> updateQuantities();
    sampler[0] -> logLikelihood();
	sampler[0] -> stepSigmaRhoAndBeta();

	// *****************************

	// set when the JT move should start
	unsigned int jtStartIteration = nIter/10;
	for( unsigned int i=0; i<nChains; ++i)
		sampler[i]->setJTStartIteration( jtStartIteration );

	// *****************************

	if( usingGPrior )
		for(unsigned int m=0; m<nChains; ++m)
			sampler[m] -> gPriorInit();

	for(unsigned int m=0; m<nChains; ++m)
		sampler[m] -> setGammaSamplerType(gammaSampler);
	

	// ****************************************

	// INIT THE FILE OUTPUT

	// clear the content of previous files
	std::ofstream logPOutFile; logPOutFile.open(outFilePath+inFile+"logP_out.txt", std::ios::out | std::ios::trunc); logPOutFile.close();
	// openlogP file in append mode
	logPOutFile.open( outFilePath+inFile+"logP_out.txt" , std::ios_base::app); // note we don't close!
	// open avg files in trunc mode to cut previous content
	std::ofstream gammaOutFile; gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc); gammaOutFile.close();
	std::ofstream gOutFile; gOutFile.open( outFilePath+inFile+"G_out.txt" , std::ios_base::trunc); gOutFile.close();
	std::ofstream piOutFile; piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc); piOutFile.close();
	std::ofstream htpOutFile; htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc); htpOutFile.close();

	// Output to file the current state
	arma::umat gamma_out = sampler[0] -> getGamma(); // out var for the gammas
	
	arma::umat g_out = arma::umat( sampler[0] -> getGAdjMat() ); // out var for G
	arma::mat beta_out = sampler[0] -> getBeta(); // out var for the betas
	arma::mat sigmaRho_out  = sampler[0] -> getSigmaRho(); // out var for the sigmas and rhos

	arma::vec tmpVec = sampler[0] -> getPi();
	arma::vec pi_out = tmpVec;
	tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
	arma::vec hotspot_tail_prob_out = tmpVec;


	gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc);
	gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out)) << std::flush;
	gammaOutFile.close();

	gOutFile.open( outFilePath+inFile+"G_out.txt" , std::ios_base::trunc);
	gOutFile << ( arma::conv_to<arma::mat>::from(g_out) ) << std::flush;   // this might be quite long...
	gOutFile.close();

	logPOutFile << 	sampler[0] -> getLogPTau() << " ";
	logPOutFile << 	sampler[0] -> getLogPEta() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPJT() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPSigmaRho() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPO() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPPi() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPGamma() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPW() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPBeta() <<  " ";
	logPOutFile << 	sampler[0] -> getLogLikelihood();
	logPOutFile << 	std::endl << std::flush;

	piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc);
	piOutFile << pi_out << std::flush;
	piOutFile.close();

	htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
	htpOutFile << hotspot_tail_prob_out << std::flush;
	htpOutFile.close();
					

	// ########
	// ########
	// ######## Start
	// ########
	// ########

	std::cout << "Starting "<< nChains <<" (parallel) chain(s) for " << nIter << " iterations:" << std::endl << std::flush;

	unsigned int tick = 1000; // how may iter for each print?

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

		tmpVec = sampler[0] -> getPi();
		pi_out += tmpVec;
		tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
		hotspot_tail_prob_out += tmpVec;

		// Print something on how the chain is going
		if( (i+1) % tick == 0 )
		{

			std::cout << " Running iteration " << i+1 << " ... local Acc Rate: ~ gamma: " << Utils::round( sampler[0] -> getGammaAccRate() , 3 );
			std::cout << " -- JT: " << Utils::round( sampler[0] -> getJTAccRate() , 3 ) ;

			if( nChains > 1)
				std::cout << " -- Global: " << Utils::round( sampler.getGlobalAccRate() , 3 ) << std::endl; 
			else
				std::cout << std::endl;
				
			// Output to files every now and then
			if( (i+1) % (tick*10) == 0 )
			{

				gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc);
				gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/((double)i+1.0) << std::flush;
				gammaOutFile.close();

				gOutFile.open( outFilePath+inFile+"G_out.txt" , std::ios_base::trunc);
				gOutFile << ( arma::conv_to<arma::mat>::from(g_out) )/((double)(i-jtStartIteration)+1.0) << std::flush;   // this might be quite long...
				gOutFile.close();

				logPOutFile << 	sampler[0] -> getLogPTau() << " ";
				logPOutFile << 	sampler[0] -> getLogPEta() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPJT() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPSigmaRho() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPO() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPPi() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPGamma() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPW() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPBeta() <<  " ";
				logPOutFile << 	sampler[0] -> getLogLikelihood();
				logPOutFile << 	std::endl << std::flush;

				piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc);
				piOutFile << pi_out/((double)i+1.0) << std::flush;
				piOutFile.close();

				htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
				htpOutFile << hotspot_tail_prob_out/((double)i+1.0) << std::flush;
				htpOutFile.close();
			}

		}

	} // end MCMC


	// Print the end
	std::cout << " MCMC ends. " /* << " Final temperature ratio ~ " << temperatureRatio  */<< "  --- Saving results and exiting" << std::endl;

	// ### Collect results and save them
	gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc);
	gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/((double)nIter+1.0) << std::flush;
	gammaOutFile.close();

	gOutFile.open( outFilePath+inFile+"G_out.txt" , std::ios_base::trunc);
	gOutFile << ( arma::conv_to<arma::mat>::from(g_out) )/((double)(nIter-jtStartIteration)+1.0) << std::flush;   // this might be quite long...
	gOutFile.close();

	logPOutFile << 	sampler[0] -> getLogPTau() << " ";
	logPOutFile << 	sampler[0] -> getLogPEta() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPJT() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPSigmaRho() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPO() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPPi() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPGamma() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPW() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPBeta() <<  " ";
	logPOutFile << 	sampler[0] -> getLogLikelihood();
	logPOutFile << 	std::endl << std::flush;
	logPOutFile.close();

	// ----
	beta_out = beta_out/((double)nIter);
	beta_out.save(outFilePath+inFile+"beta_out.txt",arma::raw_ascii);

	sigmaRho_out = sigmaRho_out/((double)nIter);
	sigmaRho_out.save(outFilePath+inFile+"sigmaRho_out.txt",arma::raw_ascii);
	// -----

	// -----
	piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc);
	piOutFile << pi_out/((double)nIter) << std::flush;
	piOutFile.close();

	htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
	htpOutFile << hotspot_tail_prob_out/((double)nIter) << std::flush;
	htpOutFile.close();
	// -----


	std::cout << "Saved to :   "+outFilePath+inFile+"****_out.txt" << std::endl;
	std::cout << "Final w : " << sampler[0] -> getW() <<  std::endl;
	std::cout << "Final tau : " << sampler[0] -> getTau() << "    w/ proposal variance: " << sampler[0] -> getVarTauProposal() << std::endl;
	std::cout << "Final eta : " << sampler[0] -> getEta() <<  std::endl;
	// std::cout << "Final o : " << sampler[0] -> getO().t() << "       w/ proposal variance: " << sampler[0] -> getVarOProposal() << std::endl;  
	// std::cout << "Final pi : " << sampler[0] -> getPi().t() << "       w/ proposal variance: " << sampler[0] -> getVarPiProposal() << std::endl;
	std::cout << "  -- Average Omega : " << arma::accu( sampler[0] -> getO() * sampler[0] -> getPi().t() )/((double)(sampler[0]->getP()*sampler[0]->getS())) <<  std::endl;
	if( nChains > 1 ) 
		std::cout << "Final temperature ratio : " << sampler[1]->getTemperature() <<  std::endl << std::endl ;

	// Exit

	std::cout << "DONE, exiting! " << std::endl << std::endl ;
	return 0;
}

int drive_dSUR( arma::mat& Y , arma::mat& X , unsigned int& nChains , unsigned int& nIter , 
				std::string& inFile , std::string& outFilePath , std::string& gammaSampler , bool usingGPrior )
{

	// ****************************************
	// **********  INIT THE CHAIN *************
	// ****************************************
	std::cout << "Initialising the MCMC Chain " << std::endl;

	ESS_Sampler<dSUR_Chain> sampler( Y , X , nChains );// this is thus also some sort of default
		// although note that you won't pass the input phase with a different method string


	// *****************************
	// need to use getX because I need the intercept
	arma::mat Q,R; arma::qr(Q,R, *sampler[0]->getX() );
	arma::mat betaInit = arma::solve(R,arma::trans(Q) * Y );
	arma::umat gammaInit = betaInit > 0.5*arma::stddev(arma::vectorise(betaInit));
	gammaInit.shed_row(0);

	sampler[0] -> gammaInit( gammaInit );
    sampler[0] -> updateQuantities();
    sampler[0] -> logLikelihood();
	sampler[0] -> stepSigmaRhoAndBeta();

	// *****************************

	if( usingGPrior )
		for(unsigned int m=0; m<nChains; ++m)
			sampler[m] -> gPriorInit();

	for(unsigned int m=0; m<nChains; ++m)
		sampler[m] -> setGammaSamplerType(gammaSampler);
		
	// *****************************


	// INIT THE FILE OUTPUT

	// clear the content of previous files
	std::ofstream logPOutFile; logPOutFile.open(outFilePath+inFile+"logP_out.txt", std::ios::out | std::ios::trunc); logPOutFile.close();
	// openlogP file in append mode
	logPOutFile.open( outFilePath+inFile+"logP_out.txt" , std::ios_base::app); // note we don't close!
	// open avg files in trunc mode to cut previous content
	std::ofstream gammaOutFile; gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc); gammaOutFile.close();
	std::ofstream gOutFile; gOutFile.open( outFilePath+inFile+"G_out.txt" , std::ios_base::trunc); gOutFile.close();
	std::ofstream piOutFile; piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc); piOutFile.close();
	std::ofstream htpOutFile; htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc); htpOutFile.close();

	// Output to file the current state
	arma::umat gamma_out = sampler[0] -> getGamma(); // out var for the gammas
	
	arma::mat beta_out = sampler[0] -> getBeta(); // out var for the betas
	arma::mat sigmaRho_out  = sampler[0] -> getSigmaRho(); // out var for the sigmas and rhos

	arma::vec tmpVec = sampler[0] -> getPi();
	arma::vec pi_out = tmpVec;
	tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
	arma::vec hotspot_tail_prob_out = tmpVec;


	gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc);
	gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out)) << std::flush;
	gammaOutFile.close();

	logPOutFile << 	sampler[0] -> getLogPTau() << " ";
	logPOutFile << 	sampler[0] -> getLogPSigmaRho() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPO() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPPi() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPGamma() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPW() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPBeta() <<  " ";
	logPOutFile << 	sampler[0] -> getLogLikelihood();
	logPOutFile << 	std::endl << std::flush;

	piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc);
	piOutFile << pi_out << std::flush;
	piOutFile.close();

	htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
	htpOutFile << hotspot_tail_prob_out << std::flush;
	htpOutFile.close();
					

	// ########
	// ########
	// ######## Start
	// ########
	// ########

	std::cout << "Starting "<< nChains <<" (parallel) chain(s) for " << nIter << " iterations:" << std::endl << std::flush;

	unsigned int tick = 1000; // how may iter for each print?

	for(unsigned int i=1; i < nIter ; ++i)
	{

		sampler.step();
		
		// #################### END LOCAL MOVES

		// ## Global moves
		// *** end Global move's section

		// UPDATE OUTPUT STATE
		gamma_out += sampler[0] -> getGamma(); // the result of the whole procedure is now my new mcmc point, so add that up

		beta_out += sampler[0] -> getBeta();
		sigmaRho_out += sampler[0] -> getSigmaRho();	

		tmpVec = sampler[0] -> getPi();
		pi_out += tmpVec;
		tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
		hotspot_tail_prob_out += tmpVec;

		// Print something on how the chain is going
		if( (i+1) % tick == 0 )
		{

			std::cout << " Running iteration " << i+1 << " ... local Acc Rate: ~ gamma: " << Utils::round( sampler[0] -> getGammaAccRate() , 3 );

			if( nChains > 1)
				std::cout << " -- Global: " << Utils::round( sampler.getGlobalAccRate() , 3 ) << std::endl; 
			else
				std::cout << std::endl;
				
			// Output to files every now and then
			if( (i+1) % (tick*10) == 0 )
			{

				gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc);
				gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/((double)i+1.0) << std::flush;
				gammaOutFile.close();

				logPOutFile << 	sampler[0] -> getLogPTau() << " ";
				logPOutFile << 	sampler[0] -> getLogPSigmaRho() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPO() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPPi() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPGamma() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPW() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPBeta() <<  " ";
				logPOutFile << 	sampler[0] -> getLogLikelihood();
				logPOutFile << 	std::endl << std::flush;

				piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc);
				piOutFile << pi_out/((double)i+1.0) << std::flush;
				piOutFile.close();

				htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
				htpOutFile << hotspot_tail_prob_out/((double)i+1.0) << std::flush;
				htpOutFile.close();
			}

		}

	} // end MCMC


	// Print the end
	std::cout << " MCMC ends. " /* << " Final temperature ratio ~ " << temperatureRatio  */<< "  --- Saving results and exiting" << std::endl;

	// ### Collect results and save them
	gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc);
	gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/((double)nIter+1.0) << std::flush;
	gammaOutFile.close();

	logPOutFile << 	sampler[0] -> getLogPTau() << " ";
	logPOutFile << 	sampler[0] -> getLogPSigmaRho() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPO() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPPi() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPGamma() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPW() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPBeta() <<  " ";
	logPOutFile << 	sampler[0] -> getLogLikelihood();
	logPOutFile << 	std::endl << std::flush;
	logPOutFile.close();

	// ----
	beta_out = beta_out/((double)nIter);
	beta_out.save(outFilePath+inFile+"beta_out.txt",arma::raw_ascii);

	sigmaRho_out = sigmaRho_out/((double)nIter);
	sigmaRho_out.save(outFilePath+inFile+"sigmaRho_out.txt",arma::raw_ascii);
	// -----

	// -----
	piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc);
	piOutFile << pi_out/((double)nIter) << std::flush;
	piOutFile.close();

	htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
	htpOutFile << hotspot_tail_prob_out/((double)nIter) << std::flush;
	htpOutFile.close();
	// -----


	std::cout << "Saved to :   "+outFilePath+inFile+"****_out.txt" << std::endl;
	std::cout << "Final w : " << sampler[0] -> getW() <<  std::endl;
	std::cout << "Final tau : " << sampler[0] -> getTau() << "    w/ proposal variance: " << sampler[0] -> getVarTauProposal() << std::endl;
	// std::cout << "Final o : " << sampler[0] -> getO().t() << "       w/ proposal variance: " << sampler[0] -> getVarOProposal() << std::endl;  
	// std::cout << "Final pi : " << sampler[0] -> getPi().t() << "       w/ proposal variance: " << sampler[0] -> getVarPiProposal() << std::endl;
	std::cout << "  -- Average Omega : " << arma::accu( sampler[0] -> getO() * sampler[0] -> getPi().t() )/((double)(sampler[0]->getP()*sampler[0]->getS())) <<  std::endl;
	if( nChains > 1 ) 
		std::cout << "Final temperature ratio : " << sampler[1]->getTemperature() <<  std::endl << std::endl ;

	// Exit

	std::cout << "DONE, exiting! " << std::endl << std::endl ;
	return 0;
}


int drive_HESS( arma::mat& Y , arma::mat& X , unsigned int& nChains , unsigned int& nIter , 
				std::string& inFile , std::string& outFilePath , std::string& gammaSampler , bool usingGPrior )
{

	// ****************************************
	// **********  INIT THE CHAIN *************
	// ****************************************
	std::cout << "Initialising the MCMC Chain " << std::endl;

	ESS_Sampler<HESS_Chain> sampler( Y , X , nChains );// this is thus also some sort of default
		// although note that you won't pass the input phase with a differetn method string


	// *****************************
	arma::mat Q,R; arma::qr(Q,R, *sampler[0]->getX() );
	arma::mat betaInit = arma::solve(R,arma::trans(Q) * Y );
	arma::umat gammaInit = betaInit > 0.5*arma::stddev(arma::vectorise(betaInit));
	gammaInit.shed_row(0);

	sampler[0] -> gammaInit( gammaInit );
    sampler[0] -> updateGammaMask();
    sampler[0] -> logLikelihood();

	// *****************************

	if( usingGPrior )
		for(unsigned int m=0; m<nChains; ++m)
			sampler[m] -> gPriorInit();

	for(unsigned int m=0; m<nChains; ++m)
		sampler[m] -> setGammaSamplerType(gammaSampler);
		

	// ****************************************

	// INIT THE FILE OUTPUT
	// clear the content of previous files

	std::ofstream logPOutFile; logPOutFile.open(outFilePath+inFile+"logP_out.txt", std::ios::out | std::ios::trunc); logPOutFile.close();
	// openlogP file in append mode
	logPOutFile.open( outFilePath+inFile+"logP_out.txt" , std::ios_base::app); // note we don't close!
	// open avg files in trunc mode to cut previous content
	std::ofstream gammaOutFile; gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc); gammaOutFile.close();
	std::ofstream piOutFile; piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc); piOutFile.close();
	std::ofstream htpOutFile; htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc); htpOutFile.close();

	// Output to file the current state
	arma::umat gamma_out = sampler[0] -> getGamma(); // out var for the gammas

	arma::vec tmpVec = sampler[0] -> getPi();
	arma::vec pi_out = tmpVec;
	tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
	arma::vec hotspot_tail_prob_out = tmpVec;


	gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc);
	gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out)) << std::flush;
	gammaOutFile.close();

	logPOutFile << 	sampler[0] -> getLogPO() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPPi() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPGamma() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPW() <<  " ";
	logPOutFile << 	sampler[0] -> getLogLikelihood();
	logPOutFile << 	std::endl << std::flush;
					
	piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc);
	piOutFile << pi_out << std::flush;
	piOutFile.close();

	htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
	htpOutFile << hotspot_tail_prob_out << std::flush;
	htpOutFile.close();


	// ########
	// ########
	// ######## Start
	// ########
	// ########

	std::cout << "Starting "<< nChains <<" (parallel) chain(s) for " << nIter << " iterations:" << std::endl << std::flush;

	unsigned int tick = 1000; // how may iter for each print?

	for(unsigned int i=1; i < nIter ; ++i)
	{

		sampler.step();
		
		// #################### END LOCAL MOVES

		// ## Global moves
		// *** end Global move's section

		// UPDATE OUTPUT STATE
		gamma_out += sampler[0] -> getGamma(); // the result of the whole procedure is now my new mcmc point, so add that up

		tmpVec = sampler[0] -> getPi();
		pi_out += tmpVec;
		tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
		hotspot_tail_prob_out += tmpVec;

		// Print something on how the chain is going
		if( (i+1) % tick == 0 )
		{

			std::cout << " Running iteration " << i+1 << " ... local Acc Rate: ~ gamma: " << Utils::round( sampler[0] -> getGammaAccRate() , 3 );

			if( nChains > 1)
				std::cout << " -- Global: " << Utils::round( sampler.getGlobalAccRate() , 3 ) << std::endl; 
			else
				std::cout << std::endl;
				
			// Output to files every now and then
			if( (i+1) % (tick*10) == 0 )
			{

				gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc);
				gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/((double)i+1.0) << std::flush;
				gammaOutFile.close();

				logPOutFile << 	sampler[0] -> getLogPO() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPPi() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPGamma() <<  " ";
				logPOutFile << 	sampler[0] -> getLogPW() <<  " ";
				logPOutFile << 	sampler[0] -> getLogLikelihood();
				logPOutFile << 	std::endl << std::flush;

				piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc);
				piOutFile << pi_out/((double)i+1.0) << std::flush;
				piOutFile.close();

				htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
				htpOutFile << hotspot_tail_prob_out/((double)i+1.0) << std::flush;
				htpOutFile.close();
			}

		}

	} // end MCMC


	// Print the end
	std::cout << " MCMC ends. " /* << " Final temperature ratio ~ " << temperatureRatio  */<< "  --- Saving results and exiting" << std::endl;

	// ### Collect results and save them
	gammaOutFile.open( outFilePath+inFile+"gamma_out.txt" , std::ios_base::trunc);
	gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/((double)nIter+1.0) << std::flush;
	gammaOutFile.close();

	logPOutFile << 	sampler[0] -> getLogPO() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPPi() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPGamma() <<  " ";
	logPOutFile << 	sampler[0] -> getLogPW() <<  " ";
	logPOutFile << 	sampler[0] -> getLogLikelihood();
	logPOutFile << 	std::endl << std::flush;
	logPOutFile.close();

	// -----
	piOutFile.open( outFilePath+inFile+"pi_out.txt" , std::ios_base::trunc);
	piOutFile << pi_out/((double)nIter) << std::flush;
	piOutFile.close();

	htpOutFile.open( outFilePath+inFile+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
	htpOutFile << hotspot_tail_prob_out/((double)nIter) << std::flush;
	htpOutFile.close();
	// -----

	std::cout << "Saved to :   "+outFilePath+inFile+"****_out.txt" << std::endl;
	std::cout << "Final w : " << sampler[0] -> getW() << "       w/ proposal variance: " << sampler[0] -> getVarWProposal() << std::endl;  
	// std::cout << "Final o : " << sampler[0] -> getO().t() << "       w/ proposal variance: " << sampler[0] -> getVarOProposal() << std::endl;  
	// std::cout << "Final pi : " << sampler[0] -> getPi().t() << "       w/ proposal variance: " << sampler[0] -> getVarPiProposal() << std::endl;
	std::cout << "  -- Average Omega : " << arma::accu( sampler[0] -> getO() * sampler[0] -> getPi().t() )/((double)(sampler[0]->getP()*sampler[0]->getS())) <<  std::endl;
	if( nChains > 1 ) 
		std::cout << "Final temperature ratio : " << sampler[1]->getTemperature() <<  std::endl << std::endl ;

	// Exit

	std::cout << "DONE, exiting! " << std::endl << std::endl ;
	return 0;

}

int drive( unsigned int nIter, unsigned int s, unsigned int p, unsigned int nChains, std::string inFile,
			std::string outFilePath, std::string method, std::string gammaSampler, bool usingGPrior )
{

	omp_init_lock(&RNGlock);  // RNG lock for the parallel part
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

	// Update the "outFilePath" (inFile variable) with the method's name
	inFile += "_"+method+"_";

	int status;

	// TODO, I hate this, but I can't initialise/instanciate templated classes
	// at runtime so this seems fair (given that the 2 drive functions have their differences in output and stuff...)
	// still if there's a more elegant solution I'd like to find it

	if( method == "SSUR" )
		status = drive_SSUR(Y,X,nChains,nIter,inFile,outFilePath,gammaSampler,usingGPrior);
	else if( method == "dSUR" )
		status = drive_dSUR(Y,X,nChains,nIter,inFile,outFilePath,gammaSampler,usingGPrior);
	else if( method == "HESS" )
		status = drive_HESS(Y,X,nChains,nIter,inFile,outFilePath,gammaSampler,usingGPrior);
	else
		status = drive_SSUR(Y,X,nChains,nIter,inFile,outFilePath,gammaSampler,usingGPrior); // this makes a default, but
			// you shound't reach here if method is wrongly specified

	return status;
}