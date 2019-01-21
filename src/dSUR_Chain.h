#include <iostream>  // for std::cout
#include <string>
#include <vector>
#include <memory>

#include "SSUR_Chain.h"

#ifndef dSUR_CHAIN_H
#define dSUR_CHAIN_H

/************************************
 * dSUR class that works with the ESS_Sampler class (dense SUR)
 * Code is not restricted to dense Gs but rather to fixed Gs (dense is still the main usecase)
 * CRTP used for global exchanges
 ***********************************/

class dSUR_Chain : public SSUR_Chain // derive from SSUR and override only a few methods
{ 

    public:

        // *******************************
        // Constructors
        // *******************************

        // empty, everything is default, except the data
        dSUR_Chain( arma::mat& , arma::mat& , double ); // Y and X  (and temperature that'll have a default value of 1.)
        dSUR_Chain( arma::mat& , arma::mat& , std::string& , bool , double ); // Y and X and gammaSamplerType (and temperature that'll have a default value of 1.)

        // full, every parameter object is initialised from existing objects
        dSUR_Chain( arma::mat& , arma::mat& , // Y and X
                double , JunctionTree& , arma::mat& , // tau, jt, sigmaRho 
                arma::vec& , arma::vec& , arma::umat& , double , arma::mat& , // o, pi, gamma, w, beta
                double ); // temperature

        // full, every parameter object is initialised from existing objects, plus the type of gamma sampler
        dSUR_Chain( arma::mat& , arma::mat& , // Y and X
                double , JunctionTree& jt_init , arma::mat& , // tau, jt, sigmaRho 
                arma::vec& , arma::vec& , arma::umat& , double , arma::mat& , // o, pi, gamma, w, beta
                std::string , bool , double ); // gamma sampler type , temperature

        // ******************************
        // Init Methods
        // ******************************

        void jtInit();
        void jtInit( JunctionTree& );

        // *********************
        // STEP FUNCTION - PERFORM ONE ITERATION FOR THE CHAIN
        // *********************

        // sample sigmaRho given Beta or Beta given sigmaRho
        // return the probability of the move and sample given the provided state rather than given the internal state
        double sampleSigmaRhoGivenBeta( const arma::mat& , arma::mat& ,
                        const arma::umat& , const arma::mat& , const arma::mat& , arma::mat& ); // "quantities"

        // logProbabilities of the above samplers (for the reverse moves)
        double logPSigmaRhoGivenBeta( const arma::mat& , const arma::mat& , 
                        const arma::umat& , const arma::mat& , const arma::mat& , const arma::mat& );

        // it updates all the internal states
        void step();

        // *******************************
        // Global operators between two chains
        // *******************************

        void swapTau( std::shared_ptr<dSUR_Chain>& );
        void swapSigmaRho( std::shared_ptr<dSUR_Chain>& );
        void swapO( std::shared_ptr<dSUR_Chain>& );
        void swapPi( std::shared_ptr<dSUR_Chain>& );
        void swapGamma( std::shared_ptr<dSUR_Chain>& );
        void swapW( std::shared_ptr<dSUR_Chain>& );
        void swapBeta( std::shared_ptr<dSUR_Chain>& );

        int globalStep( std::shared_ptr<dSUR_Chain>& );
        void swapAll( std::shared_ptr<dSUR_Chain>& );

        int exchangeAll_step( std::shared_ptr<dSUR_Chain>& );
        int exchangeGamma_step( std::shared_ptr<dSUR_Chain>& );
        int exchangeJT_step( std::shared_ptr<dSUR_Chain>& );

        int uniform_crossOver_step( std::shared_ptr<dSUR_Chain>& );
        int adapt_crossOver_step( std::shared_ptr<dSUR_Chain>& );
        int block_crossOver_step( std::shared_ptr<dSUR_Chain>& , arma::mat& , double );

        // *******************************
        // Other Methods
        // *******************************

        // update relavant quantities
        arma::mat createRhoU( const arma::mat& , const arma::mat& ); // U , sigmaRho
        void updateRhoU();

        void createQuantities( arma::umat& , arma::mat& , arma::mat& , arma::mat& ,
                const arma::umat& , const arma::mat& , const arma::mat& );
        void updateQuantities();

};

#endif


// technically all the beta functions that use xi are to be re-written, 
// but this mkakes for HUGE code repeats
// I already need to replace most of the global functions so that the global operators behave corerctly
// https://stackoverflow.com/questions/25925641/error-passing-shared-ptrderived-as-shared-ptrbase-without-const