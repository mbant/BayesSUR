#include <omp.h>
#include <vector>
#include <string>
#include <memory>

#include "SSUR_Chain.h"

#ifndef ESS_SAMPLER_H
#define ESS_SAMPLER_H

class ESS_Sampler{

    public:
        
        // Constructor - nChains and type of MCMC
        ESS_Sampler( arma::mat& , arma::mat& , unsigned int , std::string );
        
        // this gets one of the chains from the vector
        std::shared_ptr<SSUR_Chain>& operator[]( unsigned int );

        // this gets the whole vector if useful for some reason
        std::vector<std::shared_ptr<SSUR_Chain>>& getChains();
        
        // this creates a new chain with the global type
        // void addChain();
        // I mean, more complex things like this would be nice, but let's not for now ... TODO

        void step();

        // Local Operator
        void localStep();
        
        // Global Operators
        // the internal chains class know how to perform them between two chains
        //  this one selectes two chains and ask them to check for global operators to be applied
        std::pair<unsigned int , unsigned int>  randomChainSelect();
        std::pair<unsigned int , unsigned int>  nearChainSelect();

        void globalStep();

        void updateTemperatures(); // we will need extra tracking of the global acceptance rate of the moves

    private:
        unsigned int nChains;

        // Pointer to chains -  
        // we use pointers so that the client can ask for the original object and manipulate them as he wish
        std::vector<std::shared_ptr<SSUR_Chain>> chain;

        // Chain type
        std::string chainType;

        // extra parameters for global moves
        arma::mat corrMatX;
        // should also put here the ones for adapt_XO


};

#endif