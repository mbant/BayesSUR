#include <iostream>
#include <string>

// redeclare the drive funciton
int drive( unsigned int nIter, unsigned int s, unsigned int p, unsigned int nChains, std::string inFile,
			std::string outFilePath, std::string method, std::string gammaSampler, bool usingGPrior );

int main(int argc, char *  argv[])
{

	unsigned int nIter = 10; // default number of iterations
	unsigned int s=1,p=1;      // might read them from a meta-data file, but for the moment is easier like this..
	unsigned int nChains = 1;

	std::string inFile = "data.txt";
	std::string outFilePath = "";

	std::string method = "";
	std::string gammaSampler = "Bandit";
	bool usingGPrior = false;

    // ### Read and interpret command line (to put in a separate file / function?)
    int na = 1;
    while(na < argc)
    {
		if ( 0 == std::string{argv[na]}.compare(std::string{"--method"}) )
		{
			method = std::string(argv[++na]); // use the next

			if( method != "SSUR" && method != "HESS" && method != "dSUR")
			{
				std::cout << "Unknown method: only SSUR, dSUR or HESS are available" << std::endl;
			    return(1); //this is exit if I'm in a function elsewhere
			}

			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--gammaSampler"}) )
		{
			gammaSampler = std::string(argv[++na]); // use the next

			if( gammaSampler != "MC3" && gammaSampler != "mc3" && gammaSampler != "Bandit" && gammaSampler != "bandit")
			{
				std::cout << "Unknown gammaSampler method: only Bandit or MC3 are available" << std::endl;
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
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--nOutcomes"}) )
		{
			s = std::stoi(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--nPredictors"}) )
		{
			p = std::stoi(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--nChains"}) )
		{
			nChains = std::stoi(argv[++na]); // use the next
			if (na+1==argc) break; // in case it's last, break
			++na; // otherwise augment counter
		}
		else if ( 0 == std::string{argv[na]}.compare(std::string{"--inFile"}) )
		{
			inFile = ""+std::string(argv[++na]); // use the next
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


	return drive(nIter,s,p,nChains,inFile,outFilePath,method,gammaSampler,usingGPrior);

}