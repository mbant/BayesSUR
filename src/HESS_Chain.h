#include <iostream>  // for std::cout
#include <string>
#include <vector>
#include <memory>

#include "utils.h"
#include "distr.h"
#include "junction_tree.h"

#include "ESS_Atom.h"

#ifndef HESS_CHAIN_H
#define HESS_CHAIN_H

/************************************
 * HESS class that works with the ESS_Sampler class
 * This is the same model from Bottolo et al 2010, so no hyperprior on a_sigma, b_sigma
 * To maintain notation consistent with SSUR tau ( from paper) is names w here
 * CRTP used for global exchanges
 * Note that contrarily from SUR here the model still have sigma in the prior for beta for computaional convenience
 ***********************************/

class HESS_Chain : public ESS_Atom<HESS_Chain>
{ 

    public:

        // *******************************
        // Constructors
        // *******************************

        // empty, everything is default, except the data
        HESS_Chain( arma::mat& , arma::mat& , double ); // Y and X  (and temperature that'll have a default value of 1.)

        // full, every parameter object is initialised from existing objects
        HESS_Chain( arma::mat& , arma::mat& , // Y and X
                arma::vec& , arma::vec& , arma::umat& , double , // o, pi, gamma, w
                double ); // temperature

        // full, every parameter object is initialised from existing objects, plus the type of gamma sampler
        HESS_Chain( arma::mat& , arma::mat& , // Y and X
                arma::vec& , arma::vec& , arma::umat& , double ,  // o, pi, gamma, w
                std::string , bool , double ); // gamma sampler type , gPrior? , temperature

        // *******************************
        // Getters and Setters
        // *******************************
        // getters for big objects (arma::mat, JunctionTree, ...) will return references
        // so the client must be carefull with what it does with it


        // data
        std::shared_ptr<arma::mat> getY() const;
        std::shared_ptr<arma::mat> getX() const;

        // data are best set together as we need to confirm dimension matching
        void setData( std::shared_ptr<arma::mat> , std::shared_ptr<arma::mat> );

        arma::mat& getXtX();
        unsigned int getN() const;
        unsigned int getP() const;
        unsigned int getS() const;
        // no setters as they are linked to the data

        // gPrior
        void gPriorInit(); // g Prior can only be init at the start, so no proper "set" method
        bool getGPrior() const;

        // usefull quantities to keep track of
        arma::umat& getGammaMask();
        void setGammaMask( arma::umat );

        // MCMC related tuning parameters
        double getTemperature() const;
        void setTemperature( double );

        unsigned int getinternalIterationCounter() const;
        //no need to set anything for this, updated by the step function

        std::string getGammaSamplerType();
        void setGammaSamplerType( std::string& );

        // Bandit-sampling related quantities
        unsigned int getNUpdatesBandit() const;
        void setNUpdatesBandit( unsigned int );

        arma::mat& getBanditZeta();
        void setBanditZeta( arma::mat );
        
        arma::mat& getBanditAlpha();
        void setBanditAlpha( arma::mat );

        arma::mat& getBanditBeta();
        void setBanditBeta( arma::mat );
        
        arma::vec& getBanditMismatch();
        void setBanditMismatch( arma::vec );

        arma::vec& getBanditNormalisedMismatch();
        void setBanditNormalisedMismatch( arma::vec );
        
        arma::vec& getBanditNormalisedMismatchBackwards();
        void setBanditNormalisedMismatchBackwards( arma::vec );
        
        // Parameter states etc

        double getSigmaA() const; 
        void setSigmaA( double ); 
        
        double getSigmaB() const; 
        void setSigmaB( double ); 
        
        void sigmaABInit();

        // no setter for this, dedicated setter below

        // o_k
        arma::vec& getO();
        void setO( arma::vec& );
        void setO( arma::vec& , double );

        double getOA() const;
        void setOA( double );

        double getOB() const;
        void setOB( double );

        double getVarOProposal() const;
        void setVarOProposal( double );

        double getOAccRate() const;
        // no setter for this, is updated internally

        double getLogPO() const;
        // no setter for this, dedicated setter below

        // pi_j
        arma::vec& getPi();
        void setPi( arma::vec& );
        void setPi( arma::vec& , double );

        double getPiA() const;
        void setPiA( double );

        double getPiB() const;
        void setPiB( double );

        double getVarPiProposal() const;
        void setVarPiProposal( double );

        double getPiAccRate() const;
        // no setter for this, is updated internally

        double getLogPPi() const;
        // no setter for this, dedicated setter below

        // GAMMA (bandit defined above)
        arma::umat& getGamma();
        void setGamma( arma::umat& );
        void setGamma( arma::umat& , double );
        
        unsigned int getNUpdatesMC3() const;
        void setNUpdatesMC3( unsigned int );
        
        double getGammaAccRate() const;
        // no setter for this, is updated internally
        
        double getLogPGamma() const;
        // no setter for this, dedicated setter below

        // W
        double getW() const;
        void setW( double );
        void setW( double , double );

        double getWA() const;
        void setWA( double );
        
        double getWB() const;
        void setWB( double );

        double getVarWProposal() const;
        void setVarWProposal( double );

        double getWAccRate() const;

        double getLogPW() const;
        // no setter for this, dedicated setter below

        // LOG-LIKELIHOOD FOR THE SSUR MODEL
        double getLogLikelihood() const;
        // setter for this because of exchange operators and stuffs
        void setLogLikelihood( double );

        // log prior and log posterior
        double getJointLogPrior() const;
        double getJointLogPosterior() const;

        // ******************************
        // Init Methods
        // ******************************

        // init for all parameters 
        // (contrarily from the data, these are native inside the class, not pointers)
        // so the init must happen here to avoid copying and memory-waste in general
        // different version help in selecting different values for fixed hyperameters etc..
        void oInit();
        void oInit( arma::vec& );
        void oInit( arma::vec& , double , double , double );

        void piInit();
        void piInit( arma::vec& );
        void piInit( arma::vec& , double , double , double );

        void gammaInit();
        void gammaInit( arma::umat& );

        void wInit();
        void wInit( double );
        void wInit( double , double , double );
        void wInit( double , double , double , double );

        // *****************************
        // Methods for Log Probabilities
        // *****************************

        // LOG PRIORS 
        // logP for all parameters in 3 versions
        // empty for re-computing (and updating) the logP given all current values and hyperparameter values
        // with one argument for computing with a different value given the current hyperapameter values
        // with full arguments for computing the same logP but with different values and different hyperameter values

        // o_k
        double logPO( );
        double logPO( const arma::vec& );
        double logPO( const arma::vec& , double , double );

        // pi_j
        double logPPi( );
        double logPPi( arma::vec& );
        double logPPi( arma::vec& , double , double );

        // GAMMA
        double logPGamma( );
        double logPGamma( const arma::umat& );
        double logPGamma( const arma::umat& , const arma::vec& , const arma::vec& );

        // W
        double logPW( );
        double logPW( double );
        double logPW( double , double , double );

        // LOG LIKELIHOODS
        // logLik for the SSUR model, 2 versions
        // empty for re-computing and updating the logLik given all current values
        double logLikelihood( );  // this is fast and uses (gammaMask, ) XB [contains the betas] (, U) , rhoU and sigmaRho

        // with modified but available gammaMask or with available gamma dna mutantGammaMask
        double logLikelihood( const arma::umat& ); // still fast
        double logLikelihood( arma::umat& , const arma::umat& ); //gammaMask , gamma 

        // with full arguments for computing using different values
        double logLikelihood( const arma::umat& , const double , const double , const double); //gammaMask , w, a_sigma, b_sigma


        // *********************
        // STEP FUNCTION - PERFORM ONE ITERATION FOR THE CHAIN
        // *********************

        // sampler for proposed updates on gamma
        double gammaBanditProposal( arma::umat& , arma::uvec& , unsigned int& ); // steppedGamma , updateIdx, outcomeIdx
        double gammaMC3Proposal( arma::umat& , arma::uvec& , unsigned int& ); // steppedGamma , updateIdx, outcomeIdx

        // update the internal state of each parameter given all the others
        void stepOneO();
        void stepO();
        void stepOnePi();
        void stepPi();
        void stepW();

        void stepGamma();

        // it updates all the internal states
        void step();

        // update all the internal proposal RW variances based on their acceptance rate
        void updateProposalVariances();
        // more complex functions could be defined outide
        // through public methods but this as a baseline is good to have.

        int globalStep( std::shared_ptr<HESS_Chain>& );
        void swapAll( std::shared_ptr<HESS_Chain>& );

        // *******************************
        // Global operators between two chains
        // *******************************
        // assuming nu and other fixed hyperparameters are the same across chains, woudn;t make sense otherwise
        void swapO( std::shared_ptr<HESS_Chain>& );
        void swapPi( std::shared_ptr<HESS_Chain>& );
        void swapGamma( std::shared_ptr<HESS_Chain>& );
        void swapW( std::shared_ptr<HESS_Chain>& );
        
        int exchangeAll_step( std::shared_ptr<HESS_Chain>& );
        int exchangeGamma_step( std::shared_ptr<HESS_Chain>& );
        int exchangeJT_step( std::shared_ptr<HESS_Chain>& );

        int uniform_crossOver_step( std::shared_ptr<HESS_Chain>& );
        int adapt_crossOver_step( std::shared_ptr<HESS_Chain>& );
        int block_crossOver_step( std::shared_ptr<HESS_Chain>& , arma::mat& , double );

        // *******************************
        // Other Methods
        // *******************************

        // update relavant quantities
        arma::umat createGammaMask( const arma::umat& );
        void updateGammaMask();

        // Bandit-sampling related methods
        void banditInit(); // initialise all the private memebers

        // MC3 init
        void MC3Init();

    protected:

        // Data (and related quatities)
        std::shared_ptr<arma::mat> Y;
        std::shared_ptr<arma::mat> X;
        // these are pointers cause they will live on outside the MCMC
        
        bool preComputedXtX;
        arma::mat XtX;

        unsigned int n; // number of samples
        unsigned int p; // number of predictors
        unsigned int s; // number of outcomes

        // usefull quantities to keep track of
        arma::umat gammaMask;
        // these and basically everything else is native of the MCMC so are object defined here, not pointers
        
        // MCMC related tuning parameters
        double temperature;
        unsigned int internalIterationCounter;
        std::string gammaSamplerType;

        // empirical mean/variances to adapt the proposal distribution
        double wEmpiricalMean;
        arma::vec oEmpiricalMean, piEmpiricalMean;
        double wEmpiricalM2;
        arma::vec oEmpiricalM2, piEmpiricalM2; // second moment ( nb, all on the log-scale )
        double var_w_proposal_init, var_o_proposal_init, var_pi_proposal_init ;

        // Bandit-sampling related quantities
        unsigned int n_updates_bandit;
        arma::mat banditZeta;
        arma::mat banditAlpha;
        arma::mat banditBeta;
        arma::vec mismatch;
        arma::vec normalised_mismatch;
        arma::vec normalised_mismatch_backwards;

        double banditLimit;
        double banditIncrement;

        // **************************
        // Parameter states, with their associated parameters from priors and proposal and current logP
        // **************************

        // beta Prior
        bool gPrior;

        // o_k - outcome association propensity 
        // o_k ~ Beta(a_o,b_o) || its proposal is a symmetric normal RWMH in the log-scale
        arma::vec o;
        double a_o, b_o;
        double var_o_proposal;
        double o_acc_count;
        double logP_o;

        // pi_j - predictor association propensity 
        // pi_j ~ Gamma(a_pi,b_pi) || its proposal is a symmetric normal RWMH in the log-scale
        arma::vec pi;
        double a_pi, b_pi;
        double var_pi_proposal;
        double pi_acc_count;
        double logP_pi;

        // GAMMA - variable selection binary indexes
        // gamma_jk ~ Bernulli( omega_jk ), with omega_jk = v_k * u_j
        // its proposal distribution is either classic MC3 or the adaptive Bandit sampler
        arma::umat gamma;
        // prior hyperparameters are all already defined
        // proposal tuning parameters for Bandit are defined in its section
        unsigned int n_updates_MC3;
        double gamma_acc_count;
        double logP_gamma;

        // SIGMAs have a IG prior with parameters
        double a_sigma, b_sigma;

        // W - prior variance for the (slab) Normal component of the beta spike and slab prior
        // BETA slab is N(0,\sigma^2_k w)
        // w ~ InvGamma(a_w,b_w) || its proposal is a symmetric normal RWMH in the log-scale
        // a common value makes sense given common-scale for Xs and Ys
        double w;
        double var_w_proposal; // variance of the normal RW proposal
        double w_acc_count; // acceptance rate of the RW proposal
        double a_w,b_w;
        double logP_w;

        // **************************
        // LOG-LIKELIHOOD FOR THE HESS MODEL
        // **************************
        double log_likelihood;

        // Parameters for the Global moves
        // extra parameters for global moves
        arma::mat corrMatX;
        // should also put here the ones for adapt_XO

};

#endif