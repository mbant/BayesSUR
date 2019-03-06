#include <iostream>
#include <string>

// redeclare the drive funciton
int drive( const std::string& dataFile, const std::string& blockFile, const std::string& structureGraphFile, const std::string& outFilePath,  
			unsigned int nIter, unsigned int burnin, unsigned int nChains,
			const std::string& covariancePrior, 
			const std::string& gammaPrior, const std::string& gammaSampler, const std::string& gammaInit, const std::string& mrfGFile ,
			const std::string& betaPrior );

int main(int argc, char *  argv[])
{

	unsigned int nIter = 10; // default number of iterations
	unsigned int burnin = 0;
	unsigned int nChains = 1;

	std::string dataFile = "data.txt";
	std::string blockFile = "blocks.txt";
	std::string structureGraphFile = "structureGraph.txt";

	std::string outFilePath = "";

	std::string covariancePrior = "";
	std::string gammaPrior = "";
	std::string mrfGFile = "";
	std::string gammaSampler = "bandit";
	std::string gammaInit = "MLE";

	std::string betaPrior = "independent";

    // ### Read and interpret command line (to put in a separate file / function?)
    int na = 1;
    while(na < argc)
    {
		if ( 0 == std::string{argv[na]}.compare(std::string{"--covariancePrior"}) )
		{
			covariancePrior = std::string(argv[++na]); // use the next

			if ( covariancePrior == "sparse" || covariancePrior == "Sparse" || covariancePrior == "SPARSE" || covariancePrior == "HIW" || covariancePrior == "hiw" )
				covariancePrior = "HIW";
			else if ( covariancePrior == "dense" || covariancePrior == "Dense" || covariancePrior == "DENSE" || covariancePrior == "IW" || covariancePrior == "iw" ) 
				covariancePrior = "IW";
			else if ( covariancePrior == "independent" || covariancePrior == "Independent" || covariancePrior == "INDEPENDENT" || covariancePrior == "INDEP" || covariancePrior == "indep" || covariancePrior == "IG" || covariancePrior == "ig" ) 
				covariancePrior = "IG";
			else
			{
				std::cout << "Unknown covariancePrior argument: only sparse (HIW), dense(IW) or independent (IG) are available" << std::endl;
			    return(1);
			}		

			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--mrfGFile"}) )
		{
			mrfGFile = ""+std::string(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--gammaPrior"}) )
		{
			gammaPrior = std::string(argv[++na]); // use the next

			if ( gammaPrior == "hotspot" || gammaPrior == "HOTSPOT" || gammaPrior == "hotspots" || gammaPrior == "HOTSPOTS" || gammaPrior == "hs" || gammaPrior == "HS" )
				gammaPrior = "hotspot";
			else if ( gammaPrior == "MRF" || gammaPrior == "mrf" || gammaPrior == "markov random field" || gammaPrior == "Markov Random Field" ) 
				gammaPrior = "MRF";
			else if ( gammaPrior == "hierarchical" || gammaPrior == "h" || gammaPrior == "H" ) 
				gammaPrior = "hierarchical";

			else
			{
				std::cout << "Unknown gammaPrior argument: only hotspot, MRF or hierarchical are available" << std::endl;
			    return(1);
			}

			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--gammaSampler"}) )
		{
			gammaSampler = std::string(argv[++na]); // use the next

			if ( gammaSampler == "MC3" || gammaSampler == "mc3" )
				gammaSampler = "MC3";
			else if ( gammaSampler == "Bandit" || gammaSampler == "bandit" || gammaSampler == "BANDIT" ) 
				gammaSampler = "bandit";
			else
			{
				std::cout << "Unknown gammaSampler method: only Bandit or MC3 are available" << std::endl;
			    return(1);
			}

			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--gammaInit"}) )
		{
			gammaInit = std::string(argv[++na]); // use the next

			if( gammaInit != "R" && gammaInit != "0" && gammaInit != "1" && gammaInit != "MLE")
			{
				std::cout << "Unknown gammaInit method: only allowed:\n\t*\tR: random init (0.5 probability)\n\t*\t0: all elements set to 0\n\t*\t1: all elements set to 1\n\t*\tMLE: computes MLE for beta and init gamma for all significant coeffs" << std::endl;
			    return(1);
			}

			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--betaPrior"}) )
		{
			betaPrior = std::string(argv[++na]); // use the next

			if ( betaPrior == "independent" || betaPrior == "Independent" || betaPrior == "INDEPENDENT" || betaPrior == "indep" || betaPrior == "Indep" || betaPrior == "i" || betaPrior == "I" )
				betaPrior = "independent";
			else if ( betaPrior == "gprior" || betaPrior == "gPrior" || betaPrior == "g-prior" || betaPrior == "G-Prior" || betaPrior == "GPRIOR" ) 
				betaPrior = "g-prior";
			else
			{
				std::cout << "Unknown betaPrior method: only independent is available as of yet" << std::endl;
			    return(1);
			}

			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--nIter"}) )
		{
			nIter = std::stoi(argv[++na]);
			if (na+1==argc) break;
			++na;
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--burnin"}) )
		{
			burnin = std::stoi(argv[++na]);
			if (na+1==argc) break;
			++na;
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--nChains"}) )
		{
			nChains = std::stoi(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--dataFile"}) )
		{
			dataFile = ""+std::string(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--blockFile"}) )
		{
			blockFile = ""+std::string(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--structureGraphFile"}) )
		{
			structureGraphFile = ""+std::string(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--outFilePath"}) )
		{
			outFilePath = std::string(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else
		{
			std::cout << "Unknown option: " << argv[na] << std::endl;
			return(1);
    	}
    }//end reading from command line

	if( gammaPrior == "" )
	{
		if ( mrfGFile == "" )
		{
			std::cout << "Using default prior for Gamma - hotspot prior" << std::endl;
			gammaPrior = "hotspot";
		}
		else
		{
			std::cout << "No value for gammaPrior was specified, but mrfG was given - choosing MRF prior" << std::endl;
			gammaPrior = "MRF";
		}
		
	}

	int status{1};
	
	try
	{
		status =  drive(dataFile,blockFile,structureGraphFile,outFilePath,nIter,burnin,nChains,
			covariancePrior,gammaPrior,gammaSampler,gammaInit,mrfGFile,betaPrior);
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}

	return status;

}