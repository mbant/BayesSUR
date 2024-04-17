#ifndef SUR_Chain_H
#define SUR_Chain_H

#include <string>
#include <vector>
#include <memory>

#include "utils.h"
#include "distr.h"
#include "junction_tree.h"

#include "ESS_Atom.h"
#include "Parameter_types.h"

/************************************
 * SSUR class that works with the ESS_Sampler class
 * CRTP used for global exchanges
 ***********************************/


/******************
 * Note that we can go back and forth between parameterisation using (for example)
 * arma::mat originalC =
            arma::inv( arma::eye<arma::mat>(nOutcomes,nOutcomes) - arma::trimatl(sigmaRho,-1) ) *
            arma::diagmat(sigmaRho) *
            arma::inv( arma::eye<arma::mat>(nOutcomes,nOutcomes) - arma::trimatl(sigmaRho,-1) ).t() ;
 *
 * but we don't need to...            
 * ***************/

class SUR_Chain : public ESS_Atom<SUR_Chain>
{ 

    public:

        // *******************************
        // Constructors
        // *******************************

        SUR_Chain( std::shared_ptr<arma::mat> data_, std::shared_ptr<arma::mat> mrfG_, unsigned int nObservations_, 
            unsigned int nOutcomes_, unsigned int nVSPredictors_, unsigned int nFixedPredictors_,
            std::shared_ptr<arma::uvec> outcomesIdx_, std::shared_ptr<arma::uvec> VSPredictorsIdx_,
            std::shared_ptr<arma::uvec> fixedPredictorsIdx_, std::shared_ptr<arma::umat> missingDataArrayIdx_, std::shared_ptr<arma::uvec> completeCases_, 
            Gamma_Sampler_Type gamma_sampler_type_ , Gamma_Type gamma_type_ ,
            Beta_Type beta_type_ , Covariance_Type covariance_type_ , bool output_CPO = false , int maxThreads = 1, int tick = 1000, 
            double externalTemperature = 1. );

        SUR_Chain( Utils::SUR_Data& surData, 
            Gamma_Sampler_Type gamma_sampler_type_ , Gamma_Type gamma_type_ ,
            Beta_Type beta_type_ , Covariance_Type covariance_type_ ,  bool output_CPO = false , int maxThreads = 1, int tick = 1000, 
            double externalTemperature = 1. );

        SUR_Chain( Utils::SUR_Data& surData, double externalTemperature = 1. );
    
        virtual ~SUR_Chain() {};

        // *******************************
        // Getters and Setters
        // *******************************
        // getters for big objects (arma::mat, JunctionTree, ...) will return references
        // so the client must be carefull with what it does with it


        // data
        inline std::shared_ptr<arma::mat> getData() const{ return data ; }
        inline arma::mat& getXtX(){ return XtX ; }

        // mrfG
        inline std::shared_ptr<arma::mat> getMRFG() const{ return mrfG ; } 

        inline unsigned int getN() const{ return nObservations ; }
        inline unsigned int getP() const{ return nFixedPredictors+nVSPredictors ; }
        inline unsigned int getPFixed() const{ return nFixedPredictors ; }
        inline unsigned int getPVS() const{ return nVSPredictors ; }
        inline unsigned int getS() const{ return nOutcomes ; }
        // no setters as they are linked to the data

        // gPrior
        void gPriorInit(); // g Prior can only be init at the start, so no proper "set" method

        // usefull quantities to keep track of
        arma::umat& getGammaMask();
        void setGammaMask( arma::umat );

        arma::mat& getXB();
        void setXB( arma::mat );
        
        arma::mat& getU();
        void setU( arma::mat );

        arma::mat& getRhoU();
        void setRhoU( arma::mat );

        arma::urowvec& getModelSize() const;

        // MCMC related tuning parameters
        double getTemperature() const;
        void setTemperature( double );

        unsigned int getinternalIterationCounter() const;
        //no need to set anything for this, updated by the step function
        unsigned int getJTStartIteration() const;
        void setJTStartIteration( const unsigned int );


        Gamma_Sampler_Type getGammaSamplerType();
        void setGammaSamplerType( Gamma_Sampler_Type );

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

        void setTauAB( double, double );

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

        void setEtaAB( double, double );

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

        void setOAB( double, double );

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

        void setPiAB( double, double );

        double getVarPiProposal() const;
        void setVarPiProposal( double );

        double getPiAccRate() const;
        // no setter for this, is updated internally

        double getLogPPi() const;
        // no setter for this, dedicated setter below

        // MRF
//        inline arma::mat& getMRFG() { return mrf_G; }
//        void setMRFG( arma::mat& mrf_G_ ) { mrf_G = mrf_G_; logPGamma(); }

        // GAMMA (bandit defined above)
        arma::umat& getGamma();
        void setGamma( arma::umat& );
        void setGamma( arma::umat& , double );

        double getGammaD() const;
        void setGammaD( double );

        double getGammaE() const;
        void setGammaE( double );

        void setGammaDE( double, double );        
        
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

        void setWAB( double, double );        

        double getLogPW() const;
        // no setter for this, dedicated setter below

        // W0
        double getW0() const;
        void setW0( double );
        void setW0( double , double );

        double getW0A() const;
        void setW0A( double );
              
        double getW0B() const;
        void setW0B( double );

        void setW0AB( double, double );

        double getLogPW0() const;
        // no setter for this, dedicated setter below
    
        // BETA
        arma::mat& getBeta();
        void setBeta( arma::mat& );
        void setBeta( arma::mat& , double );

        double getLogPBeta() const;
        // no setter for this, dedicated setter below
    
        // PREDICTIV LIKELIHOOD FOR THE SUR MODEL
        arma::mat getPredLikelihood();
        void setPredLikelihood( arma::mat );

        // LOG-LIKELIHOOD FOR THE SUR MODEL
        double getLogLikelihood() const;
        // setter for this because of exchange operators and stuffs
        void setLogLikelihood( double );

        // log prior and log posterior
        double getJointLogPrior() const;
        double getJointLogPosterior() const;

        void setSigmaA( double ){ throw Bad_Covariance_Type( covariance_type ) ; }
        void setSigmaB( double ){ throw Bad_Covariance_Type( covariance_type ) ; }
        void setSigmaAB( double , double ){ throw Bad_Covariance_Type( covariance_type ) ; }

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
        void piInit( arma::vec& , double , double );
        void piInit( arma::vec& , double , double , double );

        void mrfGInit();
//        void mrfGInit( arma::mat& );

        void gammaInit();
        void gammaInit( arma::umat& );

        void wInit();
        void wInit( double );
        void wInit( double , double , double );
        void wInit( double , double , double , double );
    
        void w0Init();
        void w0Init( double );
        void w0Init( double , double , double );
        void w0Init( double , double , double , double );

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
        double logPGamma( const arma::umat& , const arma::vec& );
//        double logPGamma( const arma::umat& , double , double , const arma::mat& );
        double logPGamma( const arma::umat& , double , double );

        // W
        double logPW( );
        double logPW( double );
        double logPW( double , double , double );
    
        // W0
        double logPW0( );
        double logPW0( double );
        double logPW0( double , double , double );

        // BETA
        double logPBeta( );
        double logPBeta( const arma::mat& );
        double logPBeta( const arma::mat& , const arma::umat& , double , double );
        double logPBetaMask( const arma::mat& , const arma::umat& , double , double ); // faster version if the gamma mask is available
    
        // PREDICTIVE LIKELIHOODS
        arma::mat predLikelihood();
        arma::mat predLikelihood(const arma::mat& , const arma::mat& , const arma::mat&);

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
        void stepW0();
        void stepWMH();
        void stepWGibbs();
        void stepW0Gibbs();

        void stepJT();
        void stepGamma();
        void stepSigmaRhoAndBeta();

        // it updates all the internal states
        void step();

        // update all the internal proposal RW variances based on their acceptance rate
        void updateProposalVariances();
        // more complex functions could be defined outide
        // through public methods but this as a baseline is good to have.

        int globalStep( std::shared_ptr<SUR_Chain>& );
        void swapAll( std::shared_ptr<SUR_Chain>& );

        // *******************************
        // Global operators between two chains
        // *******************************
        // assuming nu and other fixed hyperparameters are the same across chains, woudn;t make sense otherwise
        void swapTau( std::shared_ptr<SUR_Chain>& );
        void swapEta( std::shared_ptr<SUR_Chain>& );
        void swapJT( std::shared_ptr<SUR_Chain>& );
        void swapSigmaRho( std::shared_ptr<SUR_Chain>& );
        void swapO( std::shared_ptr<SUR_Chain>& );
        void swapPi( std::shared_ptr<SUR_Chain>& );
        void swapGamma( std::shared_ptr<SUR_Chain>& );
        void swapW( std::shared_ptr<SUR_Chain>& );
        void swapW0( std::shared_ptr<SUR_Chain>& );
        void swapBeta( std::shared_ptr<SUR_Chain>& );
        
        int exchangeAll_step( std::shared_ptr<SUR_Chain>& );
        int exchangeGamma_step( std::shared_ptr<SUR_Chain>& );
        int exchangeJT_step( std::shared_ptr<SUR_Chain>& );

        int uniform_crossOver_step( std::shared_ptr<SUR_Chain>& );
        int adapt_crossOver_step( std::shared_ptr<SUR_Chain>& );
        int block_crossOver_step( std::shared_ptr<SUR_Chain>& , arma::mat& , double );

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

    protected:  // not private, so that they're available to derived classes

        // Data (and related quatities)
        std::shared_ptr<arma::mat> data;
        std::shared_ptr<arma::mat> mrfG;
        std::shared_ptr<arma::uvec> outcomesIdx;

        std::shared_ptr<arma::uvec> predictorsIdx;
        std::shared_ptr<arma::uvec> VSPredictorsIdx;
        std::shared_ptr<arma::uvec> fixedPredictorsIdx;

        std::shared_ptr<arma::umat> missingDataArrayIdx;
        std::shared_ptr<arma::uvec> completeCases;

        // these are pointers cause they will live on outside the MCMC
        
        bool preComputedXtX;
        arma::mat XtX;
        void setXtX();

        unsigned int nObservations; // number of samples
        unsigned int nOutcomes; // number of outcomes
        unsigned int nVSPredictors; // number of predictors to be selected
        unsigned int nFixedPredictors; // number of predictors to be kept no matter what
        bool output_CPO;
        int maxThreads;
        int tick;

        // usefull quantities to keep track of
        arma::umat gammaMask;
        arma::mat XB;
        arma::mat U;
        arma::mat rhoU;
        // these and basically everything else is native of the MCMC so are object defined here, not pointers
        
        // MCMC related tuning parameters
        double temperature;
        unsigned int internalIterationCounter;

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

        // MRF PRIOR
//        arma::mat mrf_G;
        double mrf_d, mrf_e;

        // GAMMA - variable selection binary indexes
        // gamma_jk ~ Bernulli( omega_jk ), with omega_jk = o_k * pi_j
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
    
        double w0;
        double a_w0,b_w0;
        double logP_w0;
        double w0_acc_count;
        double var_w0_proposal;

        // BETA - regression coefficients
        // beta_jk | gamma_jk ~ gamma_jk * Normal( 0 , w ) + (1-gamma_jk) * delta(0) where delta(0) is a Dirac point mass on 0
        // its update is a Gibbs move from the posterior full conditional
        arma::mat beta;
        // prior hyperparameters are all already defined
        double logP_beta;

        // **************************
        // LOG-LIKELIHOOD FOR THE SSUR MODEL
        // **************************
        double log_likelihood;
        // PREDICTIVE LIKELIHOOD
        arma::mat predLik;

        // Parameters for the Global moves
        // extra parameters for global moves
        arma::mat corrMatX;
        // should also put here the ones for adapt_XO

        // Parameter and sampler types
        Covariance_Type covariance_type;
        Gamma_Type gamma_type;
        Beta_Type beta_type;
        Gamma_Sampler_Type gamma_sampler_type;


};

#endif
