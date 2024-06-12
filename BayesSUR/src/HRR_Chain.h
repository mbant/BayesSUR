#ifndef HRR_CHAIN_H
#define HRR_CHAIN_H

#include <string>
#include <vector>
#include <memory>

#include "utils.h"
#include "distr.h"
#include "junction_tree.h"

#include "ESS_Atom.h"
#include "Parameter_types.h"

/************************************
 * HRR class that works with the ESS_Sampler class
 * This is the same model from Bottolo et al 2010, so no hyperprior on a_sigma, b_sigma
 * To maintain notation consistent with SSUR tau ( from paper) is names w here
 * CRTP used for global exchanges
 * Note that contrarily from SUR here the model still have sigma in the prior for beta for computaional convenience
 ***********************************/

class HRR_Chain : public ESS_Atom<HRR_Chain>
{ 

    public:

        // *******************************
        // Constructors
        // *******************************

        HRR_Chain( std::shared_ptr<arma::mat> data_, std::shared_ptr<arma::mat> mrfG_, unsigned int nObservations,
            unsigned int nOutcomes, unsigned int nVSPredictors, unsigned int nFixedPredictors,
            std::shared_ptr<arma::uvec> outcomesIdx_, std::shared_ptr<arma::uvec> VSPredictorsIdx_,
            std::shared_ptr<arma::uvec> fixedPredictorIdx_, std::shared_ptr<arma::umat> missingDataArrayIdx_, std::shared_ptr<arma::uvec> completeCases_, 
            Gamma_Sampler_Type gamma_sampler_type_ , Gamma_Type gamma_type_ ,
            Beta_Type beta_type_ , Covariance_Type covariance_type_ , bool output_CPO = false, int maxThreads = 1, int tick = 1000, 
            double externalTemperature = 1. );

        HRR_Chain( Utils::SUR_Data& surData,
            Gamma_Sampler_Type gamma_sampler_type_ , Gamma_Type gamma_type_ ,
            Beta_Type beta_type_ , Covariance_Type covariance_type_ , bool output_CPO = false, int maxThreads = 1, int tick = 1000, 
            double externalTemperature = 1. );

        HRR_Chain( Utils::SUR_Data& surData, double externalTemperature = 1. );
    
        virtual ~HRR_Chain() {};

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

        arma::urowvec& getModelSize() const;

        // MCMC related tuning parameters
        double getTemperature() const;
        void setTemperature( double );

        unsigned int getinternalIterationCounter() const;
        //no need to set anything for this, updated by the step function

        Gamma_Sampler_Type getGammaSamplerType();
        void setGammaSamplerType( Gamma_Sampler_Type gamma_sampler_type_ );

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
        void setSigmaAB( double, double );

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
        //inline arma::mat& getMRFG() { return mrf_G; }
        //void setMRFG( arma::mat& mrf_G_ ) { mrf_G = mrf_G_; logPGamma(); }

        double getGammaD() const;
        void setGammaD( double );

        double getGammaE() const;
        void setGammaE( double );

        void setGammaDE( double, double );

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

        void setWAB( double, double );

        double getVarWProposal() const;
        void setVarWProposal( double );

        double getWAccRate() const;

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

        double getVarW0Proposal() const;
        void setVarW0Proposal( double );

        double getW0AccRate() const;

        double getLogPW0() const;
        // no setter for this, dedicated setter below
    
        // get Beta, here we get a sample from the posterior for output reasons
        arma::mat& getBeta();// const;
    
        // PREDICTIV LIKELIHOOD FOR THE SUR MODEL
        arma::mat getPredLikelihood();
        void setPredLikelihood( arma::mat );

        // LOG-LIKELIHOOD FOR THE SSUR MODEL
        double getLogLikelihood() const;
        // setter for this because of exchange operators and stuffs
        void setLogLikelihood( double );

        // log prior and log posterior
        double getJointLogPrior() const;
        double getJointLogPosterior() const;

        void setNu( double ){ throw Bad_Covariance_Type( covariance_type ) ; }
        void setTauA( double ){ throw Bad_Covariance_Type( covariance_type ) ; }
        void setTauB( double ){ throw Bad_Covariance_Type( covariance_type ) ; }
        void setTauAB( double , double ){ throw Bad_Covariance_Type( covariance_type ) ; }
        void setEtaA( double ){ throw Bad_Covariance_Type( covariance_type ) ; }
        void setEtaB( double ){ throw Bad_Covariance_Type( covariance_type ) ; }
        void setEtaAB( double , double ){ throw Bad_Covariance_Type( covariance_type ) ; }

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
        void piInit( arma::vec& , double , double );
        void piInit( arma::vec& , double , double , double );

        void mrfGInit();
        //void mrfGInit( arma::mat& );

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
        double logPGamma( const arma::umat& , const arma::vec& );
        double logPGamma( const arma::umat& , double , double );

        // W
        double logPW( );
        double logPW( double );
        double logPW( double , double , double );
    
        // W0
        double logPW0( );
        double logPW0( double );
        double logPW0( double , double , double );
    
        // PREDICTIVE LIKELIHOODS
        arma::mat predLikelihood();
        arma::mat predLikelihood(const arma::mat& , const arma::mat& , const arma::mat&);

        // LOG LIKELIHOODS
        // logLik for the SSUR model, 2 versions
        // empty for re-computing and updating the logLik given all current values
        double logLikelihood( );  // this is fast and uses (gammaMask, ) XB [contains the betas] (, U) , rhoU and sigmaRho

        // with modified but available gammaMask or with available gamma dna mutantGammaMask
        double logLikelihood( const arma::umat& ); // still fast
        double logLikelihood( arma::umat& , const arma::umat& ); //gammaMask , gamma 

        // with full arguments for computing using different values
        double logLikelihood( const arma::umat& , const double , const double , const double , const double); //gammaMask , w, w0, a_sigma, b_sigma


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
        //void stepW0();

        void stepGamma();

        // it updates all the internal states
        void step();

        // update all the internal proposal RW variances based on their acceptance rate
        void updateProposalVariances();
        // more complex functions could be defined outide
        // through public methods but this as a baseline is good to have.

        int globalStep( std::shared_ptr<HRR_Chain>& );
        void swapAll( std::shared_ptr<HRR_Chain>& );

        // *******************************
        // Global operators between two chains
        // *******************************
        // assuming nu and other fixed hyperparameters are the same across chains, woudn;t make sense otherwise
        void swapO( std::shared_ptr<HRR_Chain>& );
        void swapPi( std::shared_ptr<HRR_Chain>& );
        void swapGamma( std::shared_ptr<HRR_Chain>& );
        void swapW( std::shared_ptr<HRR_Chain>& );
        void swapW0( std::shared_ptr<HRR_Chain>& );
        
        int exchangeAll_step( std::shared_ptr<HRR_Chain>& );
        int exchangeGamma_step( std::shared_ptr<HRR_Chain>& );
        int exchangeJT_step( std::shared_ptr<HRR_Chain>& );

        int uniform_crossOver_step( std::shared_ptr<HRR_Chain>& );
        int adapt_crossOver_step( std::shared_ptr<HRR_Chain>& );
        int block_crossOver_step( std::shared_ptr<HRR_Chain>& , arma::mat& , double );

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
        // these and basically everything else is native of the MCMC so are object defined here, not pointers
        
        // MCMC related tuning parameters
        double temperature;
        unsigned int internalIterationCounter;

        // empirical mean/variances to adapt the proposal distribution
        double wEmpiricalMean;
        arma::vec oEmpiricalMean, piEmpiricalMean;
        double wEmpiricalM2;
        arma::vec oEmpiricalM2, piEmpiricalM2; // second moment ( nb, all on the log-scale )
        double var_w_proposal_init, var_o_proposal_init, var_pi_proposal_init ;
    
        double w0EmpiricalMean;
        double w0EmpiricalM2;
        double var_w0_proposal_init;

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

        // MRF PRIOR
        arma::mat mrf_G;
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
    
        double w0;
        double var_w0_proposal; // variance of the normal RW proposal
        double w0_acc_count; // acceptance rate of the RW proposal
        double a_w0,b_w0;
        double logP_w0;

        // **************************
        // LOG-LIKELIHOOD FOR THE HRR MODEL
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
