#include <iostream>
#include <string>

// redeclare the drive funciton
int drive( const std::string& dataFile, const std::string& blockFile, const std::string& structureGraphFile, const std::string& outFilePath,  
			unsigned int nIter, unsigned int burnin, unsigned int nChains,
			const std::string& method, const bool sparse, const std::string& gammaSampler, const std::string& gammaInit, bool usingGPrior );

int main(int argc, char *  argv[])
{

	unsigned int nIter = 10; // default number of iterations
	unsigned int burnin = 0;
	unsigned int nChains = 1;

	std::string dataFile = "data.txt";
	std::string blockFile = "blocks.txt";
	std::string structureGraphFile = "structureGraph.txt";

	std::string outFilePath = "";

	std::string method = "SUR";
	bool sparse = false;
	std::string gammaSampler = "Bandit";
	std::string gammaInit = "MLE";
	
	bool usingGPrior = false;

    // ### Read and interpret command line (to put in a separate file / function?)
    int na = 1;
    while(na < argc)
    {
		if ( 0 == std::string{argv[na]}.compare(std::string{"--method"}) )
		{
			method = std::string(argv[++na]); // use the next

			if ( method != "SUR" && method != "HESS" )
			{
				std::cout << "Unknown method: only SUR or HESS are available" << std::endl;
			    return(1); //this is exit if I'm in a function elsewhere
			}

			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--sparse"}) )
		{
			sparse = true;
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--gammaSampler"}) )
		{
			gammaSampler = std::string(argv[++na]); // use the next

			if ( gammaSampler == "MC3" || gammaSampler == "mc3" )
				gammaSampler = "MC3";
			else if ( gammaSampler == "Bandit" || gammaSampler == "bandit" || gammaSampler == "BANDIT" ) 
				gammaSampler = "Bandit";
			else
			{
				std::cout << "Unknown gammaSampler method: only Bandit or MC3 are available" << std::endl;
			    return(1); //this is exit if I'm in a function elsewhere
			}

			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--gammaInit"}) )
		{
			gammaSampler = std::string(argv[++na]); // use the next

			if( gammaInit != "R" && gammaInit != "0" && gammaInit != "1" && gammaInit != "MLE")
			{
				std::cout << "Unknown gammaInit method: only allowed:\n\t*\tR: random init (0.5 probability)\n\t*\t0: all elements set to 0\n\t*\t1: all elements set to 1\n\t*\tMLE: computes MLE for beta and init gamma for all significant coeffs" << std::endl;
			    return(1); //this is exit if I'm in a function elsewhere
			}

			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--gPrior"}) )
		{
			usingGPrior = true;

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
			return(1); //this is exit if I'm in a function elsewhere
    	}
    }//end reading from command line

	return drive(dataFile,blockFile,structureGraphFile,outFilePath,nIter,burnin,nChains,method,sparse,gammaSampler,gammaInit,usingGPrior);

}