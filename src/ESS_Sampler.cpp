template<typename T>
ESS_Sampler<T>::ESS_Sampler( Utils::SUR_Data& surData , unsigned int nChains_ , double temperatureRatio ):
    updateCounter(500), // how often do we update the temperatures?
    global_proposal_count(0),
    global_acc_count(0),
    nChains(nChains_),
    chain(std::vector<std::shared_ptr<T>>(nChains_))
{

    // compile-time check that T is one of ESS_Atom's derived classes
    static_assert(std::is_base_of<ESS_Atom<T>, T>::value, "type parameter of this class must derive from ESS_Atom");

    for( unsigned int i=0; i<nChains; ++i )
        chain[i] = std::make_shared<T>( surData , std::pow( temperatureRatio , (double)i ) );  // default init for now
        
}

// Example of specialised constructor, might be needed to initialise with more precise arguments depending on the chain type
// template<> ESS_Sampler<SSUR_Chain>::ESS_Sampler( Utils::SUR_Data& surData , unsigned int nChains_ , double temperatureRatio ):
//     updateCounter(500), // how often do we update the temperatures?
//     global_proposal_count(0),
//     global_acc_count(0),
//     nChains(nChains_),
//     chain(std::vector<std::shared_ptr<SSUR_Chain>>(nChains_))
// {
//     std::string gst = "Bandit";

//     for( unsigned int i=0; i<nChains; ++i )
//         chain[i] = std::make_shared<SSUR_Chain>( surData , gst , false , std::pow( temperatureRatio , (double)i ) );  // default init for now

// }

// this gets one of the chains from the vector
template<typename T>
std::shared_ptr<T>& ESS_Sampler<T>::operator[]( unsigned int i )
{
    return chain[i];
}

// this gets the whole vector if useful for some reason
template<typename T>
std::vector<std::shared_ptr<T>>& ESS_Sampler<T>::getChains()
{
    return chain;
}

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
    // #pragma omp parallel for private(tabIndex, firstChain, secondChain)
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