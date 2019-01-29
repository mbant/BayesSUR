#ifdef _OPENMP
#include <omp.h>
#endif
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
        std::shared_ptr<T> operator[]( unsigned int i ) { return chain[i]; }

        // this gets the size of the chain which should be equal to nChains
        unsigned int size() const { return chain.size(); }
        
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


// Template declarations to force the compiler to build the needed classes
template class ESS_Sampler<SSUR_Chain>;
template class ESS_Sampler<dSUR_Chain>;
template class ESS_Sampler<HESS_Chain>;
// the alternative would be to include implementation details alongside with these declarations in a single file to include in drive.h
// https://stackoverflow.com/questions/8752837/undefined-reference-to-template-class-constructor

#endif