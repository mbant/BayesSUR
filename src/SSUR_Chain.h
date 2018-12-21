#include <iostream>  // for std::cout
#include <string>
#include <vector>
#include <memory>

#include "utils.h"
#include "distr.h"
#include "junction_tree.h"

#include "ESS_Atom.h"

#ifndef SSUR_CHAIN_H
#define SSUR_CHAIN_H

/************************************
 * SSUR class that works with the ESS_Sampler class
 * CRTP used for global exchanges
 ***********************************/

class SSUR_Chain : public ESS_Atom<SSUR_Chain>
{ 

    public:

        // *******************************
        // Constructors
        // *******************************

        // empty, everything is default, except the data
        SSUR_Chain( arma::mat& , arma::mat& , double ); // Y and X  (and temperature that'll have a default value of 1.)
        SSUR_Chain( arma::mat& , arma::mat& , std::string& , bool , double ); // Y and X and gammaSamplerType (and temperature that'll have a default value of 1.)

        // full, every parameter object is initialised from existing objects
        SSUR_Chain( arma::mat& , arma::mat& , // Y and X
                double , double , JunctionTree& , arma::mat& , // tau, eta, jt, sigmaRho 
                arma::vec& , arma::vec& , arma::umat& , double , arma::mat& , // o, pi, gamma, w, beta
                double ); // temperature

        // full, every parameter object is initialised from existing objects, plus the type of gamma sampler
        SSUR_Chain( arma::mat& , arma::mat& , // Y and X
                double , double , JunctionTree& , arma::mat& , // tau, eta, jt, sigmaRho 
                arma::vec& , arma::vec& , arma::umat& , double , arma::mat& , // o, pi, gamma, w, beta
                std::string , bool , double ); // gamma sampler type , temperature

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

        arma::mat& getXB();
        void setXB( arma::mat );
        
        arma::mat& getU();
        void setU( arma::mat );

        arma::mat& getRhoU();
        void setRhoU( arma::mat );

        // MCMC related tuning parameters
        double getTemperature() const;
        void setTemperature( double );

        unsigned int getinternalIterationCounter() const;
        //no need to set anything for this, updated by the step function
        unsigned int getJTStartIteration() const;
        void setJTStartIteration( const unsigned int );


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

        // TAU
        double getTau() const;
        void setTau( double );
        void setTau( double , double );
        
        double getTauA() const;
        void setTauA( double );

        double getTauB() const;
        void setTauB( double );

        double getVarTauProposal() const;
        void setVarTauProposal( double );

        double getTauAccRate() const;
        // no setter for this, is updated internally

        double getLogPTau() const;
        // no setter for this, dedicated setter below

        // ETA
        double getEta() const;
        void setEta( double );
        void setEta( double , double );
        
        double getEtaA() const;
        void setEtaA( double );

        double getEtaB() const;
        void setEtaB( double );

        double getLogPEta() const;
        // no setter for this, dedicated setter below

        // JT
        JunctionTree& getJT();
        arma::sp_umat getGAdjMat() const;
        void setJT( JunctionTree& );
        void setJT( JunctionTree& , double );

        unsigned int getNUpdatesJT() const;
        void setNUpdatesJT( unsigned int );

        double getJTAccRate() const;
        // no setter for this, is updated internally

        double getLogPJT() const;
        // no setter for this, dedicated setter below

        // sigmas and rhos
        arma::mat& getSigmaRho();
        void setSigmaRho( arma::mat& );
        void setSigmaRho( arma::mat& , double );

        double getNu() const; 
        void setNu( double ); 
        
        double getLogPSigmaRho();
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

        double getLogPW() const;
        // no setter for this, dedicated setter below

        // BETA
        arma::mat& getBeta();
        void setBeta( arma::mat& );
        void setBeta( arma::mat& , double );

        double getLogPBeta() const;
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
        void tauInit();
        void tauInit( double );
        void tauInit( double , double , double , double );

        void etaInit();
        void etaInit( double );
        void etaInit( double , double , double );

        void jtInit();
        void jtInit( JunctionTree& );

        void sigmaRhoInit();
        void sigmaRhoInit( arma::mat& );
        void sigmaRhoInit( arma::mat& , double );

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

        void betaInit();
        void betaInit( arma::mat& );

        // *****************************
        // Methods for Log Probabilities
        // *****************************

        // LOG PRIORS 
        // logP for all parameters in 3 versions
        // empty for re-computing (and updating) the logP given all current values and hyperparameter values
        // with one argument for computing with a different value given the current hyperapameter values
        // with full arguments for computing the same logP but with different values and different hyperameter values

        // TAU
        double logPTau( );
        double logPTau( double );
        double logPTau( double , double , double );

        // ETA
        double logPEta( );
        double logPEta( double );
        double logPEta( double , double , double );

        // JT
        double logPJT( );
        double logPJT( const JunctionTree& );
        double logPJT( const JunctionTree& , double );

        // sigma + rhos
        double logPSigmaRho( );
        double logPSigmaRho( const arma::mat& );
        double logPSigmaRho( const arma::mat& , double , double , const JunctionTree& );

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

        // BETA
        double logPBeta( );
        double logPBeta( const arma::mat& );
        double logPBeta( const arma::mat& , const arma::umat& , double );
        double logPBetaMask( const arma::mat& , const arma::umat& , double ); // faster version if the gamma mask is available

        // LOG LIKELIHOODS
        // logLik for the SSUR model, 2 versions
        // empty for re-computing and updating the logLik given all current values
        double logLikelihood( );  // this is fast and uses (gammaMask, ) XB [contains the betas] (, U) , rhoU and sigmaRho

        // with modified but available (gammaMask, ) XB [betas] (, U) , rhoU and sigmaRho
        double logLikelihood( const arma::umat& , const arma::mat& , const arma::mat& ,
                                 const arma::mat& , const arma::mat& ); // still fast

        // with full arguments for computing using different values
        // this re-computes everything and update the first four argument passed to new gammaMask,XB,U,rhoU
        double logLikelihood( arma::umat& , arma::mat& , arma::mat& , arma::mat& , //gammaMask,XB,U,rhoU
                              const arma::mat& , const arma::umat& , // beta , gamma 
                              const arma::mat& , const JunctionTree& ); // sigmaRho, jt


        // *********************
        // STEP FUNCTION - PERFORM ONE ITERATION FOR THE CHAIN
        // *********************

        // sample sigmaRho given Beta or Beta given sigmaRho
        // return the probability of the move and sample given the provided state rather than given the internal state
        double sampleSigmaRhoGivenBeta( const arma::mat& , arma::mat& , const JunctionTree& ,
                        const arma::umat& , const arma::mat& , const arma::mat& , arma::mat& ); // "quantities"
        double sampleBetaGivenSigmaRho( arma::mat& , const arma::mat& , const JunctionTree& ,
                        const arma::umat& , arma::mat& , arma::mat& , arma::mat& );
        double sampleBetaKGivenSigmaRho( const unsigned int , arma::mat& , const arma::mat& , const JunctionTree& ,
                        const arma::umat& , arma::mat& , arma::mat& , arma::mat& );
        // this samples only beta_k, so te beta vector for one outcome

        // logProbabilities of the above samplers (for the reverse moves)
        double logPSigmaRhoGivenBeta( const arma::mat& , const arma::mat& , const JunctionTree& ,
                        const arma::umat& , const arma::mat& , const arma::mat& , const arma::mat& );
        double logPBetaGivenSigmaRho( const arma::mat& , const arma::mat& , const JunctionTree& ,
                        const arma::umat& , const arma::mat& , const arma::mat& , const arma::mat& );
        double logPBetaKGivenSigmaRho( const unsigned int , const arma::mat& , const arma::mat& , const JunctionTree& ,
                const arma::umat& , const arma::mat& , const arma::mat& , const arma::mat& );


        // sample sigmaRho given Beta or Beta given sigmaRho
        // simple interface to gibbs sampling that updates the internal states (given the internal states)
        void sampleSigmaRhoGivenBeta();
        void sampleBetaGivenSigmaRho();

        // sampler for proposed updates on gamma
        double gammaBanditProposal( arma::umat& , arma::uvec& , unsigned int& ); // steppedGamma , updateIdx , outcomeIdx
        double gammaMC3Proposal( arma::umat& , arma::uvec& , unsigned int&); // steppedGamma , updateIdx , outcomeIdx


        // update the internal state of each parameter given all the others
        void stepTau();
        void stepEta();
        void stepOneO();
        void stepO();
        void stepOnePi();
        void stepPi();
        void stepW();
        void stepWMH();
        void stepWGibbs();

        void stepJT();
        void stepGamma();
        void stepSigmaRhoAndBeta();

        // it updates all the internal states
        void step();

        // update all the internal proposal RW variances based on their acceptance rate
        void updateProposalVariances();
        // more complex functions could be defined outide
        // through public methods but this as a baseline is good to have.

        int globalStep( std::shared_ptr<SSUR_Chain>& );
        void swapAll( std::shared_ptr<SSUR_Chain>& );

        // *******************************
        // Global operators between two chains
        // *******************************
        // assuming nu and other fixed hyperparameters are the same across chains, woudn;t make sense otherwise
        void swapTau( std::shared_ptr<SSUR_Chain>& );
        void swapEta( std::shared_ptr<SSUR_Chain>& );
        void swapJT( std::shared_ptr<SSUR_Chain>& );
        void swapSigmaRho( std::shared_ptr<SSUR_Chain>& );
        void swapO( std::shared_ptr<SSUR_Chain>& );
        void swapPi( std::shared_ptr<SSUR_Chain>& );
        void swapGamma( std::shared_ptr<SSUR_Chain>& );
        void swapW( std::shared_ptr<SSUR_Chain>& );
        void swapBeta( std::shared_ptr<SSUR_Chain>& );
        
        int exchangeAll_step( std::shared_ptr<SSUR_Chain>& );
        int exchangeGamma_step( std::shared_ptr<SSUR_Chain>& );
        int exchangeJT_step( std::shared_ptr<SSUR_Chain>& );

        int uniform_crossOver_step( std::shared_ptr<SSUR_Chain>& );
        int adapt_crossOver_step( std::shared_ptr<SSUR_Chain>& );
        int block_crossOver_step( std::shared_ptr<SSUR_Chain>& , arma::mat& , double );

        // *******************************
        // Other Methods
        // *******************************

        // update relavant quantities
        arma::umat createGammaMask( const arma::umat& );
        void updateGammaMask();

        arma::mat createXB( const arma::umat& , const arma::mat& ); // gammaMask, beta
        void updateXB(); 

        arma::mat createU( const arma::mat& ); // XB
        void updateU(); 

        arma::mat createRhoU( const arma::mat& , const arma::mat& , const JunctionTree& ); // U , sigmaRho, jt
        void updateRhoU();

        void createQuantities( arma::umat& , arma::mat& , arma::mat& , arma::mat& ,
                const arma::umat& , const arma::mat& , const arma::mat& , const JunctionTree& );
        void updateQuantities();


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
        arma::mat XB;
        arma::mat U;
        arma::mat rhoU;
        // these and basically everything else is native of the MCMC so are object defined here, not pointers
        
        // MCMC related tuning parameters
        double temperature;
        unsigned int internalIterationCounter;
        std::string gammaSamplerType;

        // empirical mean/variances to adapt the proposal distribution
        double tauEmpiricalMean, wEmpiricalMean;
        arma::vec oEmpiricalMean, piEmpiricalMean;
        double tauEmpiricalM2, wEmpiricalM2;
        arma::vec oEmpiricalM2, piEmpiricalM2; // second moment ( nb, all on the log-scale )
        double var_tau_proposal_init, var_o_proposal_init, var_pi_proposal_init, var_w_proposal_init ;

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

        // TAU - prior scale for the covariance matrix C -- a common value makes sense given standardised Ys --
        // tau ~ Gamma(a_tau,b_tau) || its proposal is a symmetric normal RWMH in the log-scale
        double tau;
        double a_tau,b_tau; // fixed prior parameters
        double var_tau_proposal; // variance of the normal RW proposal
        double tau_acc_count; // acceptance rate of the RW proposal
        double logP_tau; // current logP value

        // ETA - edge probability for the graph G
        // eta ~ Beta(a_eta,b_eta) || the update is a Gibbs move given the current graph G (in adjacency matrix form)
        double eta;
        double a_eta,b_eta;
        double logP_eta;

        // JT - the Junction Tree associated with the graph G
        // each edge iid ~ Bernulli(eta) || the update is MH with proposal the JT-sampler of Thomas and Green (2013)
        JunctionTree jt;
        unsigned int n_updates_jt;
        double jt_acc_count;
        double logP_jt;

        // sigmas and rhos - reparametrisation of the residual covariance matrix C
        // C ~ (hyper)InverseWishart( nu , tau * I_s ) || its update is Gibbs from the posterior full-conditional
        // note that this object has sigmas on the diagonal and symmetric rhos as off-diagonal elements
        arma::mat sigmaRho;
        double nu; // prior degree of freedom
        double logP_sigmaRho;

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

        unsigned int jtStartIteration;

        // GAMMA - variable selection binary indexes
        // gamma_jk ~ Bernulli( omega_jk ), with omega_jk = v_k * u_j
        // its proposal distribution is either classic MC3 or the adaptive Bandit sampler
        arma::umat gamma;
        // prior hyperparameters are all already defined
        // proposal tuning parameters for Bandit are defined in its section
        unsigned int n_updates_MC3;
        double gamma_acc_count;
        double logP_gamma;

        // W - prior variance for the Norma component of the beta spike and slab prior
        // w ~ InvGamma(a_w,b_w) || its proposal is a symmetric normal RWMH in the log-scale
        // a common value makes sense given common-scale for Xs and Ys
        double w;
        double a_w,b_w;
        double logP_w;
        double w_acc_count;
        double var_w_proposal;

        // BETA - regression coefficients
        // beta_jk | gamma_jk ~ gamma_jk * Normal( 0 , w ) + (1-gamma_jk) * delta(0) where delta(0) is a Dirac point mass on 0
        // its update is a Gibbs move from the posterior full conditional
        arma::mat beta;
        // prior hyperparameters are all already defined
        double logP_beta;

        bool gPrior;


        // **************************
        // LOG-LIKELIHOOD FOR THE SSUR MODEL
        // **************************
        double log_likelihood;

        // Parameters for the Global moves
        // extra parameters for global moves
        arma::mat corrMatX;
        // should also put here the ones for adapt_XO

};

#endif