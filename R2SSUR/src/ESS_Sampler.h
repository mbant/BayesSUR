#ifndef ESS_SAMPLER_H
#define ESS_SAMPLER_H

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
#include "SUR_Chain.h"

template<typename T>  //  the template here should be a class derived from ESS_Atom
class ESS_Sampler{

    public:
        
        // Constructor - nChains and type of MCMC
        ESS_Sampler( Utils::SUR_Data& surData , unsigned int nChains_ , double temperatureRatio ,
        	Gamma_Sampler_Type gamma_sampler_type, Gamma_Type gamma_type, Beta_Type beta_type, Covariance_Type covariance_type);

        ESS_Sampler( Utils::SUR_Data& surData , unsigned int nChains_ , double temperatureRatio ) : 
            ESS_Sampler( surData , nChains_ , temperatureRatio , 
               Gamma_Sampler_Type::bandit, Gamma_Type::hotspot , Beta_Type::independent , Covariance_Type::sparse){}

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

// ***********************************
// ***** Implementation
// ***********************************


template<typename T>
ESS_Sampler<T>::ESS_Sampler( Utils::SUR_Data& surData , unsigned int nChains_ , double temperatureRatio ,
    Gamma_Sampler_Type gamma_sampler_type, Gamma_Type gamma_type, Beta_Type beta_type, Covariance_Type covariance_type):
        updateCounter(500), // how often do we update the temperatures?
        global_proposal_count(0),
        global_acc_count(0),
        nChains(nChains_),
        chain(std::vector<std::shared_ptr<T>>(nChains))
{


    // compile-time check that T is one of ESS_Atom's derived classes
    static_assert(std::is_base_of<ESS_Base, T>::value, "type parameter of this class must derive from ESS_Atom");

    for( unsigned int i=0; i<nChains; ++i )
        chain[i] = std::make_shared<T>( surData , 
            gamma_sampler_type, gamma_type, beta_type, covariance_type,
            std::pow( temperatureRatio , (double)i ) );  // default init for now
}

// Example of specialised constructor, might be needed to initialise with more precise arguments depending on the chain type
// template<> ESS_Sampler<SSUR_Chain>::ESS_Sampler( ... ):
//     ...
// {
//     ...
// }

// ********************************
// STEP OPERATORS
// ********************************

template<typename T>
void ESS_Sampler<T>::step()
{ 
    this->localStep();
    this->globalStep();
}

// Local Operator
template<typename T>
void ESS_Sampler<T>::localStep()
{
    // OMP GIVES RUNTIME ERRORS, PROBABLY DUE TO SOME VARIABLE ACCESSED AT THE SAME TIME OR SOMRTHING :/
    // TODO but leave it out for now -- ideally use MPI or similar to distribute the different chains -- shared variables will be a pain though
    // #ifdef _OPENMP    
    // #pragma omp parallel for schedule(static,1)
    // #endif
    for( auto i = chain.begin(); i < chain.end(); ++i )
        (*i) -> step();

    // this sintactic sugar is disabled for omp
    // for( auto i : chain )
    //     i->step();
}

// Global Operators
// the internal chains class know how to perform them between two chains
//  this one selectes two chains and ask them to check for global operators to be applied

template<typename T>
std::pair<unsigned int , unsigned int>  ESS_Sampler<T>::randomChainSelect()
{
    unsigned int chainIdx = 1, firstChain = 0, secondChain = 1;

    // Select the chains to swap
    chainIdx = (nChains>2) ? Distributions::randIntUniform(1, (nChains)*(nChains-1)/2 ) : 1;   // (nChains-1)*(nChains-2)/2 is the number of possible chain combinations with nChains

    for(unsigned int c=1; c<nChains; ++c)
    {
        for(unsigned int r=0; r<c; ++r)
        {
            if( (--chainIdx) == 0 ){
                firstChain = r;
                secondChain = c;
                break;
            }
        }
    }

    return std::pair<unsigned int , unsigned int>( firstChain , secondChain );
}


template<typename T>
std::pair<unsigned int , unsigned int>  ESS_Sampler<T>::nearChainSelect()
{
    unsigned int firstChain = 0, secondChain = 1;

    if( nChains>2 )
    {
        firstChain = Distributions::randIntUniform(1, nChains-2 );  // so not the first (0) or last (nChains-1) indexes
        secondChain = ( Distributions::randU01() < 0.5 ) ? firstChain-1 : firstChain+1 ; // then select a neighbour
    }

    return std::pair<unsigned int , unsigned int>( firstChain , secondChain );
}


template<typename T>
void ESS_Sampler<T>::globalStep()
{
    ++global_proposal_count;
    std::pair<unsigned int , unsigned int> chainIdx = {0,1};

    if( nChains > 1 )
    {
        if( Distributions::randU01() < 0.9 )
        {
            if( Distributions::randU01() < 0.5 )
                chainIdx = randomChainSelect();
            else
                chainIdx = nearChainSelect();

            global_acc_count += chain[chainIdx.first] -> globalStep( chain[chainIdx.second] );
    
        }else
            global_acc_count += allExchangeAll_step();

        if ( (global_proposal_count % updateCounter) == 0 )
            updateTemperatures();
    }
}

template<typename T>
double ESS_Sampler<T>::getGlobalAccRate() const { return ((double)global_acc_count)/((double)global_proposal_count); }


// This below assumes that acceptance ratio is a monotonic function of the temperatre ratio
template<typename T>
void ESS_Sampler<T>::updateTemperatures()
{

    double tempRatio = chain[1]->getTemperature(); // / temperatures(0) = 1

    // check acceptance rate
    // if too high/low , update temperatures

    if( getGlobalAccRate() > 0.3 )
    {
        tempRatio *= 1.1 ;

        for( unsigned int i=1; i < nChains ; ++i )
        {
            chain[i]->setTemperature( tempRatio * chain[i-1]->getTemperature() );
        }

        std::cout << "Temperature ladder updated, new temperature ratio : " << tempRatio << std::endl;

    }else if( getGlobalAccRate() < 0.05 )
    {
        tempRatio = std::max( 1. + 1e-8 , tempRatio * 0.9 ); 

        for( unsigned int i=1; i < nChains ; ++i )
        {
            chain[i]->setTemperature( tempRatio * chain[i-1]->getTemperature() );
        }
        
        std::cout << "Temperature ladder updated, new temperature ratio : " << tempRatio << std::endl;
    }

    // I want to maintain a sort-of-moving avaerage acceptance count for the global moves, so that 
    // when we check if we should update the temperatures we only check in the near-past
    // in order to do this, I simply reset the varaibles each time I update the temperature
    global_proposal_count = 0;
    global_acc_count = 0;

}

template<typename T>
int ESS_Sampler<T>::allExchangeAll_step()
{
    unsigned int nChainCombinations = ((nChains)*(nChains-1)/2);

    arma::vec pExchange( nChainCombinations +1 );
    unsigned int swapIdx, firstChain, secondChain;

    arma::umat indexTable( pExchange.n_elem, 2);
    unsigned int tabIndex = 0;
    indexTable(tabIndex,0) = 0; indexTable(tabIndex,1) = 0;
    tabIndex++;

    for(unsigned int c=1; c<nChains; ++c)
    {
        for(unsigned int r=0; r<c; ++r)
        {
            indexTable(tabIndex,0) = r; indexTable(tabIndex,1) = c;
            tabIndex++;
        }
    }


    pExchange(0) = 0.; // these are log probabilities, remember!
    // #ifdef _OPENMP
    // #pragma omp parallel for private(tabIndex, firstChain, secondChain)
    // #endif
    for(tabIndex = 1; tabIndex <= nChainCombinations; ++tabIndex)
    {

        firstChain = indexTable(tabIndex,0);
        secondChain  = indexTable(tabIndex,1);

        // Swap probability
        pExchange(tabIndex) = ( chain[firstChain]->getLogLikelihood() * chain[firstChain]->getTemperature() -
                            chain[secondChain]->getLogLikelihood() * chain[secondChain]->getTemperature() ) * 
					( 1. / chain[secondChain]->getTemperature() - 1. / chain[firstChain]->getTemperature() );
                    //  no priors because that is not tempered so it cancels out
    }

    // normalise and cumulate the weights
    double logSumWeights = Utils::logspace_add(pExchange); // normaliser
    arma::vec cumulPExchange = arma::cumsum( arma::exp( pExchange - logSumWeights ) ); // this should sum to one

    // Now select which swap happens
    double val = Distributions::randU01();

    swapIdx = 0;
    while( val > cumulPExchange(swapIdx) )
    {
        swapIdx++;
    }

    if( swapIdx != 0 )
    {
        firstChain = indexTable(swapIdx,0);
        secondChain  = indexTable(swapIdx,1);

        swapAll( chain[firstChain] , chain[secondChain] );

        return 1;
    }else
        return 0;

}

template<typename T>
void ESS_Sampler<T>::swapAll( std::shared_ptr<T>& thisChain , std::shared_ptr<T>& thatChain )
{

    // POINTER SWAP
    std::swap ( thisChain , thatChain );

    double swapTemp = thisChain -> getTemperature();
    thisChain -> setTemperature( thatChain -> getTemperature() );
    thatChain -> setTemperature( swapTemp );

}


#endif