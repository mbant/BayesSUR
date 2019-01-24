#include <omp.h>
#include <vector>
#include <string>
#include <memory>

#include "utils.h"
#include "distr.h"

#include "ESS_Atom.h"
#include "HESS_Chain.h"
#include "SSUR_Chain.h"
#include "dSUR_Chain.h"

#ifndef ESS_SAMPLER_H
#define ESS_SAMPLER_H

template<typename T>  //  the template here should be a class derived from ESS_Atom
class ESS_Sampler{

    public:
        
        // Constructor - nChains and type of MCMC
        ESS_Sampler( Utils::SUR_Data& surData , unsigned int nChains_ , double temperatureRatio );
        ESS_Sampler( Utils::SUR_Data& surData , unsigned int nChains_ ) : ESS_Sampler( surData , nChains_ , 1.2 ){}
        
        // this gets one of the chains from the vector
        std::shared_ptr<T>& operator[]( unsigned int );

        // this gets the whole vector if useful for some reason
        std::vector<std::shared_ptr<T>>& getChains();
        
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

        int allExchangeAll_step();
        void swapAll( std::shared_ptr<T>& thisChain , std::shared_ptr<T>& thatChain );


        // Temperature ladder update and getter for the acceptance rate of global updates
        double getGlobalAccRate() const;

        void updateTemperatures(); // we will need extra tracking of the global acceptance rate of the moves

    private:

        unsigned int nChains;

        // Pointer to chains -  
        // we use pointers so that the client can ask for the original object and manipulate them as he wish
        std::vector<std::shared_ptr<T>> chain;

        unsigned int updateCounter; // how often do we update the temperatures?
        unsigned int global_proposal_count, global_acc_count;

};

#include "ESS_Sampler.cpp"

#endif