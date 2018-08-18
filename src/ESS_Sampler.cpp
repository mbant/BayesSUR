#include "ESS_Sampler.h"
#include <omp.h>

ESS_Sampler::ESS_Sampler( arma::mat& Y_ , arma::mat& X_ , unsigned int nChains_ , std::string chainType_ )
{
    nChains = nChains_ ;
    chainType = chainType_ ;

    chain = std::vector<std::shared_ptr<SSUR_Chain>>(nChains);

    double tempRatio = 1.2;

    for( unsigned int i=0; i<nChains; ++i )
        chain[i] = std::make_shared<SSUR_Chain>( Y_ , X_ , std::pow( tempRatio , (double)i ) );  // default init for now

    corrMatX = arma::cor( X_ );

    updateCounter = 200; // how often do we update the temperatures?
    global_proposal_count = 0;
    global_acc_count = 0;

}

// this gets one of the chains from the vector
std::shared_ptr<SSUR_Chain>& ESS_Sampler::operator[]( unsigned int i )
{
    return chain[i];
}

// this gets the whole vector if useful for some reason
std::vector<std::shared_ptr<SSUR_Chain>>& ESS_Sampler::getChains()
{
    return chain;
}

// ********************************
// STEP OPERATORS
// ********************************

void ESS_Sampler::step()
{ 
    this->localStep();
    this->globalStep();
}

// Local Operator
void ESS_Sampler::localStep()
{
    // OMP GIVES RUNTIME ERRORS, PROBABLY DUE TO SOME VARIABLE ACCESSED AT THE SAME TIME OR SOMRTHING :/
    // TODO but leave it out for now -- ideally use MPI or similar to distribute the different chains -- shared variables will be a pain though
    // #pragma omp parallel for schedule(static,1)
    for( auto i = chain.begin(); i < chain.end(); ++i )
        (*i) -> step();

    // this sintactic sugar is disabled for omp
    // for( auto i : chain )
    //     i->step();
}

// Global Operators
// the internal chains class know how to perform them between two chains
//  this one selectes two chains and ask them to check for global operators to be applied

std::pair<unsigned int , unsigned int>  ESS_Sampler::randomChainSelect()
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


std::pair<unsigned int , unsigned int>  ESS_Sampler::nearChainSelect()
{
    unsigned int firstChain = 0, secondChain = 1;

    if( nChains>2 )
    {
        firstChain = Distributions::randIntUniform(1, nChains-2 );  // so not the first (0) or last (nChains-1) indexes
        secondChain = ( Distributions::randU01() < 0.5 ) ? firstChain-1 : firstChain+1 ; // then select a neighbour
    }

    return std::pair<unsigned int , unsigned int>( firstChain , secondChain );
}

void ESS_Sampler::globalStep()
{
    unsigned int globalType = 0;
    if( chainType == "SSUR" || chainType == "ssur") 
        globalType = Distributions::randIntUniform(0,10);
    else if( chainType == "HESS" || chainType == "hess") 
        globalType = Distributions::randIntUniform(0,5);

    std::pair<unsigned int , unsigned int> chainIdx;

    if( nChains > 1 )
    {
        switch(globalType){

            case 0: 
                break;

            // -- Exchange and CrossOver
            case 1: 
                ++global_proposal_count;
                chainIdx = randomChainSelect();
                global_acc_count += chain[ chainIdx.first ] -> exchangeGamma_step( chain[ chainIdx.second ] );
                break;
            
            case 2: 
                ++global_proposal_count;
                chainIdx = nearChainSelect();
                global_acc_count += chain[ chainIdx.first ] -> exchangeGamma_step( chain[ chainIdx.second ] );
                break;
            case 3: 
                ++global_proposal_count;
                chainIdx = randomChainSelect();
                global_acc_count += chain[ chainIdx.first ] -> adapt_crossOver_step( chain[ chainIdx.second ] );
                break;
            
            case 4: 
                ++global_proposal_count;
                chainIdx = randomChainSelect();
                global_acc_count += chain[ chainIdx.first ] -> uniform_crossOver_step( chain[ chainIdx.second ] );
                break;
            
            case 5: 
                ++global_proposal_count;
                chainIdx = randomChainSelect();
                global_acc_count += chain[ chainIdx.first ] -> block_crossOver_step( chain[ chainIdx.second ] , corrMatX , 0.25 );
                break;
            
            case 6: 
                ++global_proposal_count;
                chainIdx = randomChainSelect();
                global_acc_count += chain[ chainIdx.first ] -> exchangeJT_step( chain[ chainIdx.second ] );
                break;
            
            case 7: 
                ++global_proposal_count;
                chainIdx = nearChainSelect();
                global_acc_count += chain[ chainIdx.first ] -> exchangeJT_step( chain[ chainIdx.second ] );
                break;
            
            case 8: 
                ++global_proposal_count;
                chainIdx = randomChainSelect();
                global_acc_count += exchangeAll_step( chain[ chainIdx.first ] , chain[ chainIdx.second ] );
                break;
            
            case 9: 
                ++global_proposal_count;
                chainIdx = nearChainSelect();
                global_acc_count += exchangeAll_step( chain[ chainIdx.first ] , chain[ chainIdx.second ] );
                break;

            case 10: 
                ++global_proposal_count;
                global_acc_count += allExchangeAll_step( chain );
                break;

            default: 
                break;
        }

        if ( (global_proposal_count % updateCounter) == 0 )
            updateTemperatures();

    }
}


double ESS_Sampler::getGlobalAccRate() const { return ((double)global_acc_count)/((double)global_proposal_count); }


// This below assumes that acceptance ratio is a monotonic function of the temperatre ratio
void ESS_Sampler::updateTemperatures()
{

    double tempRatio = chain[1]->getTemperature(); // / temperatures(0) = 1

    // check acceptance rate
    // if too high/low , update temperatures

    if( getGlobalAccRate() > 0.35 )
    {
        tempRatio *= 1.1 ;

        for( unsigned int i=1; i < nChains ; ++i )
        {
            chain[i]->setTemperature( tempRatio * chain[i-1]->getTemperature() );
        }

    }else if( getGlobalAccRate() < 0.15 )
    {
        tempRatio *= 0.9 ; 

        for( unsigned int i=1; i < nChains ; ++i )
        {
            chain[i]->setTemperature( tempRatio * chain[i-1]->getTemperature() );
        }
    }

    // I want to maintain a sort-of-moving avaerage acceptance count for the global moves, so that 
    // when we check if we should update the temperatures we only check in the near-past
    // in order to do this, I simply reset the varaibles each time I update the temperature
    global_proposal_count = 0;
    global_acc_count = 0;



}