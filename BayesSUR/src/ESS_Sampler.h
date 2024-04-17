#ifndef ESS_SAMPLER_H
#define ESS_SAMPLER_H

#ifndef CCODE
#include <RcppArmadillo.h>
using Rcpp::Rcout;
using Rcpp::Rcerr;
#else
#include <armadillo>
#include <iostream>
#define Rcout std::cout
#define Rcerr std::cerr
#endif

#include <vector>
#include <string>
#include <memory>

#include "utils.h"
#include "distr.h"

#include "ESS_Atom.h"
#include "HRR_Chain.h"
#include "SUR_Chain.h"

template<typename T>  //  the template here should be a class derived from ESS_Atom
class ESS_Sampler{
    
public:
    
    // Constructor - nChains and type of MCMC
    ESS_Sampler( Utils::SUR_Data& surData , unsigned int nChains_ , double temperatureRatio ,
                Gamma_Sampler_Type gamma_sampler_type, Gamma_Type gamma_type, Beta_Type beta_type, Covariance_Type covariance_type, 
                bool output_CPO , int maxThreads , int tick , unsigned int burnin_ );
    
    ESS_Sampler( Utils::SUR_Data& surData , unsigned int nChains_ , double temperatureRatio ) :
    ESS_Sampler( surData , nChains_ , temperatureRatio ,
    Gamma_Sampler_Type::bandit, Gamma_Type::hotspot , Beta_Type::independent , Covariance_Type::HIW, false){}
    
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
    
    void setHyperParameters( const Utils::Chain_Data& chainData );
    
private:
    
    unsigned int nChains;
    unsigned int burnin;
    int maxThreads;
    int tick;
    
    // Pointer to chains -
    // we use pointers so that the client can ask for the original object and manipulate them as he wish
    std::vector<std::shared_ptr<T>> chain;
    
    unsigned int updateCounter; // how often do we update the temperatures?
    unsigned int global_proposal_count, global_acc_count, global_count;
    double tmpRand;
    
};

// ***********************************
// ***** Implementation
// ***********************************


template<typename T>
ESS_Sampler<T>::ESS_Sampler( Utils::SUR_Data& surData , unsigned int nChains_ , double temperatureRatio ,
                            Gamma_Sampler_Type gamma_sampler_type, Gamma_Type gamma_type, Beta_Type beta_type, Covariance_Type covariance_type, 
                            bool output_CPO , int maxThreads , int tick, unsigned int burnin_ ):
nChains(nChains_),
burnin(burnin_),
chain(std::vector<std::shared_ptr<T>>(nChains)),
updateCounter(100), // how often do we update the temperatures?
global_proposal_count(0),
global_acc_count(0),
global_count(0)
{
    
    
    // compile-time check that T is one of ESS_Atom's derived classes
    static_assert(std::is_base_of<ESS_Base, T>::value, "type parameter of this class must derive from ESS_Atom");
    
    for( unsigned int i=0; i<nChains; ++i )
        chain[i] = std::make_shared<T>( surData ,
                                       gamma_sampler_type, gamma_type, beta_type, covariance_type, output_CPO, maxThreads, tick,
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
    
    for( unsigned int i=0; i<nChains; ++i )
        chain[i] -> step();
    
    // this sintactic sugar is disabled for omp
    // for( auto i : chain )
    //     i->step();
    
    // and also this for pedantic compilers
    // for( auto i = chain.begin(); i!=chain.end();  ++i)
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
    chainIdx = (nChains>2) ? randIntUniform(1, (nChains)*(nChains-1)/2 ) : 1;   // (nChains-1)*(nChains-2)/2 is the number of possible chain combinations with nChains
    
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
        firstChain = randIntUniform(1, nChains-2 );  // so not the first (0) or last (nChains-1) indexes
        secondChain = ( randU01() < 0.5 ) ? firstChain-1 : firstChain+1 ; // then select a neighbour
    }
    
    return std::pair<unsigned int , unsigned int>( firstChain , secondChain );
}


template<typename T>
void ESS_Sampler<T>::globalStep()
{
    ++global_proposal_count;
    ++global_count;
    std::pair<unsigned int , unsigned int> chainIdx = {0,1};
    
    if( nChains > 1 )
    {
        tmpRand = randU01();
        if( tmpRand < 0.9 )
        {
            if( tmpRand < 0.5 )
                chainIdx = randomChainSelect();
            else
                chainIdx = nearChainSelect();
            
            global_acc_count += chain[chainIdx.first] -> globalStep( chain[chainIdx.second] );
            
        }else
            global_acc_count += allExchangeAll_step();
        
        if ( ((global_proposal_count % updateCounter) == 0) && (global_count <= burnin) )
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
        
        Rcout << "Temperature ladder updated, new temperature ratio : " << tempRatio << std::endl;
        
    }else if( getGlobalAccRate() < 0.05 )
    {
        tempRatio = std::max( 1. + 1e-8 , tempRatio * 0.9 );
        
        for( unsigned int i=1; i < nChains ; ++i )
        {
            chain[i]->setTemperature( tempRatio * chain[i-1]->getTemperature() );
        }
        
        Rcout << "Temperature ladder updated, new temperature ratio : " << tempRatio << std::endl;
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
    
    // build the index table
    arma::vec pExchange( nChainCombinations +1 );
    arma::umat indexTable( pExchange.n_elem, 2);
    
    unsigned int idx = 0;
    indexTable(idx,0) = 0; indexTable(idx,1) = 0;
    idx++;
    
    for(unsigned int c=1; c<nChains; ++c)
    {
        for(unsigned int r=0; r<c; ++r)
        {
            indexTable(idx,0) = r; indexTable(idx,1) = c;
            idx++;
        }
    }
    
    // Compute the swap probabilities
    pExchange(0) = 0.; // these are log probabilities, remember!
   
#ifdef _OPENMP
#pragma omp parallel for
#endif
    
    for(unsigned int tabIndex = 1; tabIndex <= nChainCombinations; ++tabIndex)
    {
        
        unsigned int firstChain = indexTable(tabIndex,0);
        unsigned int secondChain = indexTable(tabIndex,1);
        
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
    double val = randU01();
    
    unsigned int swapIdx = 0;
    while( val > cumulPExchange(swapIdx) )
    {
        swapIdx++;
    }
    
    if( swapIdx != 0 )
    {
        unsigned int firstChain = indexTable(swapIdx,0);
        unsigned int secondChain  = indexTable(swapIdx,1);
        
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

template<typename T>
void ESS_Sampler<T>::setHyperParameters( const Utils::Chain_Data& chainData )
{
    // MRF Prior
    if ( chainData.gamma_type == Gamma_Type::mrf )
    {
        if ( ! std::isnan( chainData.mrfD ) )
        {
            if ( ! std::isnan( chainData.mrfE ) )
            {
                for( auto c : chain )
                    c->setGammaDE( chainData.mrfD, chainData.mrfE );
            }
            else
            {
                for( auto c : chain )
                    c->setGammaD( chainData.mrfD );
            }
        }
        else if ( ! std::isnan( chainData.mrfE ) )
        {
            for( auto c : chain )
                c->setGammaE( chainData.mrfE );
        }
    }
    
    // Hierarchical and Hotspot
    if ( chainData.gamma_type == Gamma_Type::hierarchical || chainData.gamma_type == Gamma_Type::hotspot )
    {
        // a_pi and b_pi
        if ( ! std::isnan( chainData.piA ) )
        {
            if ( ! std::isnan( chainData.piB ) )
            {
                for( auto c : chain )
                    c->setPiAB( chainData.piA, chainData.piB );
            }
            else
            {
                for( auto c : chain )
                    c->setPiA( chainData.piA );
            }
        }
        else if ( ! std::isnan( chainData.piB ) )
        {
            for( auto c : chain )
                c->setPiB( chainData.piB );
        }
        
        // HOTSPOT only
        if ( chainData.gamma_type == Gamma_Type::hotspot )
        {
            // a_o and b_o
            if ( ! std::isnan( chainData.oA ) )
            {
                if ( ! std::isnan( chainData.oB ) )
                {
                    for( auto c : chain )
                        c->setOAB( chainData.oA, chainData.oB );
                }
                else
                {
                    for( auto c : chain )
                        c->setOA( chainData.oA );
                }
            }
            else if ( ! std::isnan( chainData.oB ) )
            {
                for( auto c : chain )
                    c->setOB( chainData.oB );
            }
        }
    }
    
    
    // Covariance Prior
    if ( chainData.covariance_type == Covariance_Type::IG )  // HRR
    {
        // A_Sigma and B_Sigma
        if ( ! std::isnan( chainData.sigmaA ) )
        {
            if ( ! std::isnan( chainData.sigmaB ) )
            {
                for( auto c : chain )
                    c->setSigmaAB( chainData.sigmaA, chainData.sigmaB );
            }
            else
            {
                for( auto c : chain )
                    c->setSigmaA( chainData.sigmaA );
            }
        }
        else if ( ! std::isnan( chainData.sigmaB ) )
        {
            for( auto c : chain )
                c->setSigmaB( chainData.sigmaB );
        }
    }
    else // SUR
    {
        // A_tau and B_tau
        if ( ! std::isnan( chainData.tauA ) )
        {
            if ( ! std::isnan( chainData.tauB ) )
            {
                for( auto c : chain )
                    c->setTauAB( chainData.tauA, chainData.tauB );
            }
            else
            {
                for( auto c : chain )
                    c->setTauA( chainData.tauA );
            }
        }
        else if ( ! std::isnan( chainData.tauB ) )
        {
            for( auto c : chain )
                c->setTauB( chainData.tauB );
        }
        
        // NU
        if ( ! std::isnan( chainData.nu ) )
            for( auto c : chain )
                c->setNu( chainData.nu );
    }
    
    if ( chainData.covariance_type == Covariance_Type::HIW )  // Sparse SUR
    {
        // A_eta and B_eta
        if ( ! std::isnan( chainData.etaA ) )
        {
            if ( ! std::isnan( chainData.etaB ) )
            {
                for( auto c : chain )
                    c->setEtaAB( chainData.etaA, chainData.etaB );
            }
            else
            {
                for( auto c : chain )
                    c->setEtaA( chainData.etaA );
            }
        }
        else if ( ! std::isnan( chainData.etaB ) )
        {
            for( auto c : chain )
                c->setEtaB( chainData.etaB );
        }
    }
    
    // a_w and b_w
    if ( ! std::isnan( chainData.wA ) )
    {
        if ( ! std::isnan( chainData.wB ) )
        {
            for( auto c : chain )
                c->setWAB( chainData.wA, chainData.wB );
        }
        else
        {
            for( auto c : chain )
                c->setWA( chainData.wA );
        }
    }
    else if ( ! std::isnan( chainData.wB ) )
    {
        for( auto c : chain )
            c->setWB( chainData.wB );
    }
    
    // a_w0 and b_w0
    if ( ! std::isnan( chainData.w0A ) )
    {
        if ( ! std::isnan( chainData.w0B ) )
        {
            for( auto c : chain )
                c->setW0AB( chainData.w0A, chainData.w0B );
        }
        else
        {
            for( auto c : chain )
                c->setW0A( chainData.w0A );
        }
    }
    else if ( ! std::isnan( chainData.w0B ) )
    {
        for( auto c : chain )
            c->setW0B( chainData.w0B );
    }
    
}

#endif
