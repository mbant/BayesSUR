#include "SUR_Chain.h"

// *******************************
// Constructors
// *******************************


SUR_Chain::SUR_Chain( std::shared_ptr<arma::mat> data_, std::shared_ptr<arma::mat> mrfG_, unsigned int nObservations_,
                     unsigned int nOutcomes_, unsigned int nVSPredictors_, unsigned int nFixedPredictors_,
                     std::shared_ptr<arma::uvec> outcomesIdx_, std::shared_ptr<arma::uvec> VSPredictorsIdx_,
                     std::shared_ptr<arma::uvec> fixedPredictorsIdx_, std::shared_ptr<arma::umat> missingDataArrayIdx_, std::shared_ptr<arma::uvec> completeCases_,
                     Gamma_Sampler_Type gamma_sampler_type_ , Gamma_Type gamma_type_ ,
                     Beta_Type beta_type_ , Covariance_Type covariance_type_ , bool output_CPO , int maxThreads , int tick, 
                     double externalTemperature ):
data(data_), mrfG(mrfG_), outcomesIdx(outcomesIdx_), VSPredictorsIdx(VSPredictorsIdx_), fixedPredictorsIdx(fixedPredictorsIdx_),
missingDataArrayIdx(missingDataArrayIdx_), completeCases(completeCases_),
nObservations(nObservations_), nOutcomes(nOutcomes_), nVSPredictors(nVSPredictors_), nFixedPredictors(nFixedPredictors_),
temperature(externalTemperature),internalIterationCounter(0),jtStartIteration(0),
covariance_type(covariance_type_), gamma_type(gamma_type_),beta_type(beta_type_),gamma_sampler_type(gamma_sampler_type_)
{
    
    predictorsIdx = std::make_shared<arma::uvec>(arma::join_vert( *fixedPredictorsIdx, *VSPredictorsIdx ));
    setXtX();
    
    switch ( gamma_sampler_type )
    {
        case Gamma_Sampler_Type::bandit :
            banditInit();
            break;
            
        case Gamma_Sampler_Type::mc3 :
            MC3Init();
            break;
            
        default:
            throw Bad_Gamma_Sampler_Type ( gamma_sampler_type ) ;
    }
    
    tauInit();
    
    etaInit();
    jtInit();
    
    switch ( gamma_type )
    {
        case Gamma_Type::hotspot :
            oInit();
            piInit();
            break;
            
        case Gamma_Type::hierarchical :
            piInit();
            break;
            
        case Gamma_Type::mrf :
            mrfGInit();
            break;
            
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
    
    gammaInit();
    updateGammaMask();
    w0Init();
    wInit();
    
    betaInit();
    sigmaRhoInit();
    
    updateQuantities();
    
    logLikelihood();
    
    // init for sigma rho and beta to reasonable values -- one step of gibbs
    stepSigmaRhoAndBeta();
    if( output_CPO && (temperature = 1.) ){
        predLikelihood();
    }
    
}


SUR_Chain::SUR_Chain( Utils::SUR_Data& surData,
                     Gamma_Sampler_Type gamma_sampler_type_ , Gamma_Type gamma_type_ ,
                     Beta_Type beta_type_ , Covariance_Type covariance_type_  , bool output_CPO , int maxThreads , int tick, 
                     double externalTemperature ):
SUR_Chain(surData.data,surData.mrfG,surData.nObservations,surData.nOutcomes,surData.nVSPredictors,surData.nFixedPredictors,
surData.outcomesIdx,surData.VSPredictorsIdx,surData.fixedPredictorsIdx,surData.missingDataArrayIdx,surData.completeCases,
          gamma_sampler_type_,gamma_type_,beta_type_,covariance_type_,output_CPO,maxThreads,tick,externalTemperature){ }

SUR_Chain::SUR_Chain( Utils::SUR_Data& surData, double externalTemperature ):
SUR_Chain(surData.data,surData.mrfG,surData.nObservations,surData.nOutcomes,surData.nVSPredictors,surData.nFixedPredictors,
surData.outcomesIdx,surData.VSPredictorsIdx,surData.fixedPredictorsIdx,surData.missingDataArrayIdx,surData.completeCases,
          Gamma_Sampler_Type::bandit , Gamma_Type::hotspot , Beta_Type::independent , Covariance_Type::HIW , false,
          externalTemperature){ }


// *******************************
// Getters and Setters
// *******************************

void SUR_Chain::setXtX()
{
    
    // Compute XtX
    if( (nFixedPredictors+nVSPredictors) < 100000 )  // kinda arbitrary value, how can we assess a more sensible one?
    {
        preComputedXtX = true;
        XtX = data->cols( *predictorsIdx ).t() * data->cols( *predictorsIdx );
        //corrMatX = arma::cor( data->submat(arma::regspace<arma::uvec>(0,nObservations-1), *VSPredictorsIdx ) );  // this is only for values to be selected
        corrMatX = arma::cor( data->cols( *VSPredictorsIdx ) );
    }else{
        
        preComputedXtX = false;
        XtX.clear();          // if not precomputed, these two are just reset
        corrMatX.clear();
    }
}

// gPrior
void SUR_Chain::gPriorInit() // g Prior can only be init at the start, so no proper set method
{
    if( internalIterationCounter > 0 )
        throw std::runtime_error(std::string("gPrior can only be initialised at the start of the MCMC"));
    
    
    // unconditional Error, gPrior not implemented yet
    throw std::runtime_error(std::string("gPrior is not fully functional yet, so its use is blocked"));
    
    // set the beta type
    beta_type = Beta_Type::gprior;
    
    // re-initialise the w parameter in line with the new prior, as w now has a different meaning
    wInit( (double)nObservations , 0.5*nOutcomes + nOutcomes -1. , 0.5*nObservations*nOutcomes ); // these values are taken from Lewin 2016
    
    // update internals
    logPW();
    logPBeta();
    
}

// useful quantities to keep track of
arma::umat& SUR_Chain::getGammaMask(){ return gammaMask; }
void SUR_Chain::setGammaMask( arma::umat  externalGammaMask ){ gammaMask =  externalGammaMask ; }

arma::mat& SUR_Chain::getXB(){ return XB; }
void SUR_Chain::setXB( arma::mat externalXB ){ XB = externalXB ; }

arma::mat& SUR_Chain::getU(){ return U ; }
void SUR_Chain::setU( arma::mat externalU ){ U = externalU ; }

arma::mat& SUR_Chain::getRhoU(){ return rhoU ; }
void SUR_Chain::setRhoU( arma::mat externalRhoU ){ rhoU = externalRhoU ; }

arma::urowvec& SUR_Chain::getModelSize() const
{
    static arma::urowvec modelSize;
    modelSize = nFixedPredictors + arma::sum( gamma , 0 ); // 0 is to get the sum of the elements in each column
    return modelSize;
}

// MCMC related tuning parameters
double SUR_Chain::getTemperature() const{ return temperature; }
void SUR_Chain::setTemperature( double temp_ )
{
    log_likelihood = log_likelihood * temperature / temp_ ; // re-temper the log_likelihood first
    temperature = temp_ ;
}

unsigned int SUR_Chain::getinternalIterationCounter() const{ return internalIterationCounter ; }

unsigned int SUR_Chain::getJTStartIteration() const{ return jtStartIteration; }
void SUR_Chain::setJTStartIteration( const unsigned int jts ){ jtStartIteration = jts; }

Gamma_Sampler_Type SUR_Chain::getGammaSamplerType(){ return gamma_sampler_type ; }
void SUR_Chain::setGammaSamplerType( Gamma_Sampler_Type gamma_sampler_type_ )
{
    if( gamma_sampler_type != gamma_sampler_type_ )
    {
        gamma_sampler_type = gamma_sampler_type_ ;
        
        switch ( gamma_sampler_type )
        {
            case Gamma_Sampler_Type::bandit :
                banditInit();
                break;
                
            case Gamma_Sampler_Type::mc3 :
                MC3Init();
                break;
                
            default:
                throw Bad_Gamma_Sampler_Type ( gamma_sampler_type );
        }
    }
}

// Bandit-sampling related quantities
unsigned int SUR_Chain::getNUpdatesBandit() const{ return n_updates_bandit ; }
void SUR_Chain::setNUpdatesBandit( unsigned int n_updates_bandit_ ){ n_updates_bandit = n_updates_bandit_ ; }

arma::mat& SUR_Chain::getBanditZeta(){ return banditZeta; }
void SUR_Chain::setBanditZeta( arma::mat banditZeta_ ){ banditZeta = banditZeta_ ; }

arma::mat& SUR_Chain::getBanditAlpha(){ return banditAlpha ; }
void SUR_Chain::setBanditAlpha( arma::mat banditAlpha_ ){ banditAlpha = banditAlpha_ ; }

arma::mat& SUR_Chain::getBanditBeta(){ return banditBeta ; }
void SUR_Chain::setBanditBeta( arma::mat banditBeta_ ){ banditBeta = banditBeta_ ; }

arma::vec& SUR_Chain::getBanditMismatch(){ return mismatch; }
void SUR_Chain::setBanditMismatch( arma::vec mismatch_ ){ mismatch = mismatch_ ; }

arma::vec& SUR_Chain::getBanditNormalisedMismatch(){ return normalised_mismatch ; }
void SUR_Chain::setBanditNormalisedMismatch( arma::vec normalised_mismatch_ ){ normalised_mismatch = normalised_mismatch_ ; }

arma::vec& SUR_Chain::getBanditNormalisedMismatchBackwards(){ return normalised_mismatch_backwards ; }
void SUR_Chain::setBanditNormalisedMismatchBackwards( arma::vec normalised_mismatch_backwards_ ){ normalised_mismatch_backwards = normalised_mismatch_backwards_ ; }

// Parameter states etc

// TAU
double SUR_Chain::getTau() const{ return tau; }
void SUR_Chain::setTau( double tau_ )
{
    tau = tau_;
    logPTau(); // this updates the internal logP variable
}

void SUR_Chain::setTau( double tau_ , double logP_tau_)
{
    tau = tau_;
    logP_tau = logP_tau_ ;
}

double SUR_Chain::getTauA() const{ return a_tau; }
void SUR_Chain::setTauA( double a_tau_ )
{
    a_tau = a_tau_ ;
    logPTau();
}

double SUR_Chain::getTauB() const{ return b_tau; }
void SUR_Chain::setTauB( double b_tau_ )
{
    b_tau = b_tau_ ;
    logPTau();
}

void SUR_Chain::setTauAB( double a_tau_ , double b_tau_ )
{
    if ( covariance_type != Covariance_Type::HIW && covariance_type != Covariance_Type::IW )
        throw Bad_Covariance_Type( covariance_type );
    
    a_tau = a_tau_ ;
    b_tau = b_tau_ ;
    logPTau();
}

double SUR_Chain::getVarTauProposal() const{ return var_tau_proposal ; }
void SUR_Chain::setVarTauProposal( double var_tau_proposal_ ){ var_tau_proposal = var_tau_proposal_ ; }

double SUR_Chain::getTauAccRate() const{ return tau_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double SUR_Chain::getLogPTau() const{ return logP_tau ; }
// no setter for this, dedicated setter below

// ETA
double SUR_Chain::getEta() const{ return eta ; }
void SUR_Chain::setEta( double eta_ )
{
    eta = eta_ ;
    logPEta(); // this updates the internal logP variable
}

void SUR_Chain::setEta( double eta_ , double logP_eta_ )
{
    eta = eta_ ;
    logP_eta = logP_eta_ ;
}

double SUR_Chain::getEtaA() const{ return a_eta ; }
void SUR_Chain::setEtaA( double a_eta_ )
{
    a_eta = a_eta_ ;
    logPEta();
}

double SUR_Chain::getEtaB() const{ return b_eta ; }
void SUR_Chain::setEtaB( double b_eta_ )
{
    b_eta = b_eta_ ;
    logPEta();
}

void SUR_Chain::setEtaAB( double a_eta_ , double b_eta_ )
{
    if( covariance_type != Covariance_Type::HIW )
        throw Bad_Covariance_Type( covariance_type );
    
    a_eta = a_eta_ ;
    b_eta = b_eta_ ;
    logPEta();
}

double SUR_Chain::getLogPEta() const{ return logP_eta ; }
// no setter for this, dedicated setter below

// JT
JunctionTree& SUR_Chain::getJT(){ return jt; }  // need to check this works correctly TODO (i.e. correct behaviour is returning a copy of jt)
arma::sp_umat SUR_Chain::getGAdjMat() const{ return jt.getAdjMat(); }  // need to check this works correctly TODO (i.e. correct behaviour is returning a copy of jt)
void SUR_Chain::setJT( JunctionTree& externalJT )
{
    jt = externalJT ; // old jt gets destroyed as there's no reference to him anymore, new jt "points" to the foreign object
    if ( covariance_type == Covariance_Type::HIW )
        logPJT();
}

void SUR_Chain::setJT( JunctionTree& externalJT , double logP_jt_ )
{
    jt = externalJT ; // old jt gets destroyed as there's no reference to him anymore, new jt "points" to the foreign object
    logP_jt = logP_jt_ ;
}

unsigned int SUR_Chain::getNUpdatesJT() const{ return n_updates_jt; }
void SUR_Chain::setNUpdatesJT( unsigned int n_updates_jt_ ){ n_updates_jt = n_updates_jt_ ; }

double SUR_Chain::getJTAccRate() const{ return jt_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double SUR_Chain::getLogPJT() const{ return logP_jt ; }
// no setter for this, dedicated setter below

// sigmas and rhos
arma::mat& SUR_Chain::getSigmaRho(){ return sigmaRho ; }
void SUR_Chain::setSigmaRho( arma::mat&  externalSigmaRho )
{
    sigmaRho =  externalSigmaRho ;
    logPSigmaRho();
    // THIS DO NOT ACT ON LOGLIK, SO IT'S IMPORTANT TO REMEMBER TO CHANGE THAT
}

void SUR_Chain::setSigmaRho( arma::mat&  externalSigmaRho , double logP_sigmaRho_ )
{
    sigmaRho =  externalSigmaRho ;
    logP_sigmaRho = logP_sigmaRho_ ;
    // THIS DO NOT ACT ON LOGLIK, SO IT'S IMPORTANT TO REMEMBER TO CHANGE THAT
}

double SUR_Chain::getNu() const{ return nu; }
void SUR_Chain::setNu( double nu_ )
{
    nu = nu_ ;
    logPSigmaRho();
}

double SUR_Chain::getLogPSigmaRho(){ return logP_sigmaRho ; }
// no setter for this, dedicated setter below

// o_k
arma::vec& SUR_Chain::getO(){ return o ; }
void SUR_Chain::setO( arma::vec& o_ )
{
    o = o_ ;
    logPO();
}

void SUR_Chain::setO( arma::vec& o_ , double logP_o_ )
{
    o = o_ ;
    logP_o = logP_o_ ;
}

double SUR_Chain::getOA() const{ return a_o ; }
void SUR_Chain::setOA( double a_o_ )
{
    a_o = a_o_ ;
    logPO();
}

double SUR_Chain::getOB() const{ return b_o ; }
void SUR_Chain::setOB( double b_o_ )
{
    b_o = b_o_ ;
    logPO();
}

void SUR_Chain::setOAB( double a_o_ , double b_o_ )
{
    if ( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type( gamma_type );
    
    a_o = a_o_ ;
    b_o = b_o_ ;
    logPO();
}

double SUR_Chain::getVarOProposal() const{ return var_o_proposal ; }
void SUR_Chain::setVarOProposal( double var_o_proposal_ ){ var_o_proposal = var_o_proposal_ ; }

double SUR_Chain::getOAccRate() const{ return o_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double SUR_Chain::getLogPO() const{ return logP_o ; }
// no setter for this, dedicated setter below

// pi_j
arma::vec& SUR_Chain::getPi(){ return pi ; }
void SUR_Chain::setPi( arma::vec& pi_ )
{
    pi = pi_ ;
    logPPi();
}

void SUR_Chain::setPi( arma::vec& pi_ , double logP_pi_ )
{
    pi = pi_ ;
    logP_pi = logP_pi_ ;
}

double SUR_Chain::getPiA() const{ return a_pi ; }
void SUR_Chain::setPiA( double a_pi_ )
{
    a_pi = a_pi_ ;
    logPPi();
}

double SUR_Chain::getPiB() const{ return b_pi ; }
void SUR_Chain::setPiB( double b_pi_ )
{
    b_pi = b_pi_ ;
    logPPi();
}

void SUR_Chain::setPiAB( double a_pi_ , double b_pi_ )
{
    if ( gamma_type != Gamma_Type::hotspot && gamma_type != Gamma_Type::hierarchical )
        throw Bad_Gamma_Type( gamma_type );
    
    a_pi = a_pi_ ;
    b_pi = b_pi_ ;
    logPPi();
}

double SUR_Chain::getVarPiProposal() const{ return var_pi_proposal ; }
void SUR_Chain::setVarPiProposal( double var_pi_proposal_ ){ var_pi_proposal = var_pi_proposal_ ; }

double SUR_Chain::getPiAccRate() const{ return pi_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double SUR_Chain::getLogPPi() const{ return logP_pi ; }
// no setter for this, dedicated setter below

// GAMMA
arma::umat& SUR_Chain::getGamma(){ return gamma ; }
void SUR_Chain::setGamma( arma::umat& externalGamma )
{
    gamma = externalGamma ;
    logPGamma();
}

void SUR_Chain::setGamma( arma::umat& externalGamma , double logP_gamma_ )
{
    gamma = externalGamma ;
    logP_gamma = logP_gamma_ ;
}

double SUR_Chain::getGammaD() const{ return mrf_d ; }
void SUR_Chain::setGammaD( double mrf_d_ )
{
    if( gamma_type != Gamma_Type::mrf )
        throw Bad_Gamma_Type( gamma_type );
    
    mrf_d = mrf_d_ ;
    logPGamma();
}

double SUR_Chain::getGammaE() const{ return mrf_e ; }
void SUR_Chain::setGammaE( double mrf_e_ )
{
    if( gamma_type != Gamma_Type::mrf )
        throw Bad_Gamma_Type( gamma_type );
    
    mrf_e = mrf_e_ ;
    logPGamma();
}

void SUR_Chain::setGammaDE( double mrf_d_ , double mrf_e_ )
{
    if( gamma_type != Gamma_Type::mrf )
        throw Bad_Gamma_Type( gamma_type );
    
    mrf_d = mrf_d_ ;
    mrf_e = mrf_e_ ;
    logPGamma();
}

unsigned int SUR_Chain::getNUpdatesMC3() const{ return n_updates_MC3 ; }
void SUR_Chain::setNUpdatesMC3( unsigned int n_updates_MC3_ ){ n_updates_MC3 = n_updates_MC3_ ; }

double SUR_Chain::getGammaAccRate() const{ return gamma_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double SUR_Chain::getLogPGamma() const{ return logP_gamma ; }
// no setter for this, dedicated setter below

// W
double SUR_Chain::getW() const{ return w ; }
void SUR_Chain::setW( double w_ )
{
    w = w_ ;
    logPW();
}

void SUR_Chain::setW( double w_ , double logP_w_ )
{
    w = w_ ;
    logP_w = logP_w_ ;
}

double SUR_Chain::getWA() const{ return a_w ; }
void SUR_Chain::setWA( double a_w_ )
{
    a_w = a_w_ ;
    logPW();
}

double SUR_Chain::getWB() const{ return b_w ; }
void SUR_Chain::setWB( double b_w_ )
{
    b_w = b_w_ ;
    logPW();
}

void SUR_Chain::setWAB( double a_w_ , double b_w_ )
{
    a_w = a_w_ ;
    b_w = b_w_ ;
    logPW();
}

double SUR_Chain::getLogPW() const{ return logP_w; }
// no setter for this, dedicated setter below

// W_0
double SUR_Chain::getW0() const{ return w0 ; }
void SUR_Chain::setW0( double w0_ )
{
  w0 = w0_ ;
  logPW0();
}

void SUR_Chain::setW0( double w0_ , double logP_w0_ )
{
  w0 = w0_ ;
  logP_w0 = logP_w0_ ;
}

double SUR_Chain::getW0A() const{ return a_w0 ; }
void SUR_Chain::setW0A( double a_w0_ )
{
  a_w0 = a_w0_ ;
  logPW0();
}

double SUR_Chain::getW0B() const{ return b_w0 ; }
void SUR_Chain::setW0B( double b_w0_ )
{
  b_w0 = b_w0_ ;
  logPW0();
}

void SUR_Chain::setW0AB( double a_w0_ , double b_w0_ )
{
  a_w0 = a_w0_ ;
  b_w0 = b_w0_ ;
  logPW0();
}

double SUR_Chain::getLogPW0() const{ return logP_w0; }
// no setter for this, dedicated setter below

// BETA
arma::mat& SUR_Chain::getBeta(){ return beta ; }
void SUR_Chain::setBeta( arma::mat&  externalBeta )
{
    beta =  externalBeta ;
    logPBeta();
    // THIS DO NOT ACT ON LOGLIK, SO IT'S IMPORTANT TO REMEMBER TO CHANGE THAT
}

void SUR_Chain::setBeta( arma::mat&  externalBeta , double logP_beta_ )
{
    beta =  externalBeta ;
    logP_beta = logP_beta_ ;
    // THIS DO NOT ACT ON LOGLIK, SO IT'S IMPORTANT TO REMEMBER TO CHANGE THAT
}

double SUR_Chain::getLogPBeta() const{ return logP_beta ; }
// no setter for this, dedicated setter below

// PREDICTIVE-LIKELIHOOD OF INDIVIDUALS FOR THE SSUR MODEL
arma::mat SUR_Chain::getPredLikelihood() { return predLik ; }
void SUR_Chain::setPredLikelihood( arma::mat predLik_ ){ predLik = predLik_ ; }

// LOG-LIKELIHOOD FOR THE SSUR MODEL
double SUR_Chain::getLogLikelihood() const{ return log_likelihood ; }
void SUR_Chain::setLogLikelihood( double log_likelihood_ ){ log_likelihood = log_likelihood_ ; }

double SUR_Chain::getJointLogPrior() const
{
    return logP_tau + logP_eta + logP_jt + logP_sigmaRho + logP_o + logP_pi + logP_gamma + logP_w + logP_w0 + logP_beta;
}

double SUR_Chain::getJointLogPosterior() const
{
    return getJointLogPrior() + getLogLikelihood();
}

// ******************************
// Init Methods
// ******************************

// init for all parameters
// (contrarily from the data, these are native inside the class, not pointers)
// so the init must happen here to avoid copying and memory-waste in general
// different version help in selecting different values for fixed hyperameters etc..
void SUR_Chain::tauInit( double tau_init , double a_tau_ , double b_tau_ , double var_tau_proposal_ )
{
    tau = tau_init ;
    
    a_tau = a_tau_;
    b_tau = b_tau_;
    
    var_tau_proposal = var_tau_proposal_;
    
    tau_acc_count = 0.;
    
    logPTau();
}

void SUR_Chain::tauInit()
{
    tauInit( 1. , 0.1 , 10. , 0.02 );
}

void SUR_Chain::tauInit( double tau_init )
{
    tauInit( tau_init , 0.1 , 10. , 0.02 );
}


void SUR_Chain::etaInit( double eta_init , double a_eta_ , double b_eta_ )
{
    eta = eta_init;
    
    a_eta = a_eta_;
    b_eta = b_eta_;
    
    logPEta();
}

void SUR_Chain::etaInit()
{
    etaInit( 0.1 , /*0.5*s*(s-1)*0.1*/ 0.1 , /*0.5*s*(s-1)*(1.-0.1)*/ 1. ); // 0.1  is essentially the expected active proportion
}

void SUR_Chain::etaInit( double eta_init )
{
    etaInit( eta_init , /*0.5*s*(s-1)*0.1*/ 0.1 , /*0.5*s*(s-1)*(1.-0.1)*/ 1. ); // should I use eta_init in the prior parameters?
}


// **** care, due to the need to pass a whole object, for constuctor purposes
// this one function differently than the rest (I need to repeat the code in both function essentally..)
void SUR_Chain::jtInit( JunctionTree& jt_init )
{
    jt = jt_init ;
    jt_acc_count = 0.;
    
    switch ( covariance_type )
    {
        case Covariance_Type::HIW :
            n_updates_jt = 5; // default value, should I pick something different?
            logPJT();
            break;
            
        case Covariance_Type::IW :
            n_updates_jt = 0;
            break;
            
        default:
            throw Bad_Covariance_Type ( covariance_type );
    }
}

void SUR_Chain::jtInit()
{
    jt_acc_count = 0.;
    
    switch ( covariance_type )
    {
        case Covariance_Type::HIW :
            jt = JunctionTree( nOutcomes , "empty" ); //empty constructor, diagonal adj matrix of dimension s
            n_updates_jt = 5; // default value, should I pick something different?
            logPJT();
            break;
            
        case Covariance_Type::IW :
            jt = JunctionTree( nOutcomes , "full" ); // full constructor, dense adj matrix of dimension s
            n_updates_jt = 0;
            break;
            
        default:
            break;
    }
}
// end different from the rest


void SUR_Chain::sigmaRhoInit( arma::mat& sigmaRho_init , double nu_ )
{
    sigmaRho = sigmaRho_init;
    nu = nu_;
    logPSigmaRho();
}

void SUR_Chain::sigmaRhoInit()
{
    arma::mat init = arma::eye<arma::mat>(nOutcomes,nOutcomes);
    sigmaRhoInit( init , nOutcomes+2. );
}

void SUR_Chain::sigmaRhoInit( arma::mat& sigmaRho_init )
{
    sigmaRhoInit( sigmaRho_init , nOutcomes+2. );
}

void SUR_Chain::oInit( arma::vec& o_init , double a_o_ , double b_o_ , double var_o_proposal_ )
{
    
    if( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    o = o_init;
    a_o = a_o_;
    b_o = b_o_;
    var_o_proposal = var_o_proposal_;
    o_acc_count = 0.;
    
    logPO();
}

void SUR_Chain::oInit()
{
    if( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    arma::vec init = arma::ones<arma::vec>(nOutcomes) / std::max( 500. , (double)nVSPredictors ) ;
    oInit( init , 2. , (double)nVSPredictors-2. , 0.005 );
}

void SUR_Chain::oInit( arma::vec& o_init )
{
    if( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    oInit( o_init , 2. , (double)nVSPredictors-2. , 0.005 );
}


// this is only for the hotspot
void SUR_Chain::piInit( arma::vec& pi_init , double a_pi_ , double b_pi_ , double var_pi_proposal_ )
{
    if ( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    pi = pi_init;
    a_pi = a_pi_;
    b_pi = b_pi_;
    var_pi_proposal = var_pi_proposal_;
    pi_acc_count = 0.;
    logPPi();
}

// this is only for the hierarchical
void SUR_Chain::piInit( arma::vec& pi_init , double a_pi_ , double b_pi_ )
{
    if ( gamma_type != Gamma_Type::hierarchical )
        throw Bad_Gamma_Type ( gamma_type );
    
    pi = pi_init;
    a_pi = a_pi_;
    b_pi = b_pi_;
    logPPi();
}

// this is either hotspot or hierarchical
void SUR_Chain::piInit()
{
    arma::vec init = arma::ones<arma::vec>(nVSPredictors) ;
    switch ( gamma_type )
    {
        case Gamma_Type::hotspot :
            piInit( init , 2. , 1. , 0.02 );
            break;
            
        case Gamma_Type::hierarchical :
            piInit( init , 1. , (double)nOutcomes-1. );
            break;
            
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
}

void SUR_Chain::piInit( arma::vec& pi_init )
{
    switch ( gamma_type )
    {
        case Gamma_Type::hotspot :
            piInit( pi_init , 2. , 1. , 0.02 );
            break;
            
        case Gamma_Type::hierarchical :
            piInit( pi_init , 1. , (double)nOutcomes-1. );
            break;
            
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
}

void SUR_Chain::mrfGInit()
{
    if( gamma_type != Gamma_Type::mrf )
        throw Bad_Gamma_Type ( gamma_type );
    
    //    mrf_G = arma::zeros<arma::mat>(0,2);
    mrf_d = -3. ;
    mrf_e = 0.3 ;
}
/*
 void SUR_Chain::mrfGInit( arma::mat& mrf_G_ )
 {
 if( gamma_type != Gamma_Type::mrf )
 throw Bad_Gamma_Type ( gamma_type );
 mrf_G = mrf_G_;
 mrf_d = -3. ;
 mrf_e = 0.2 ;
 }
 */

void SUR_Chain::gammaInit( arma::umat& gamma_init )
{
    gamma = gamma_init;
    gamma_acc_count = 0.;
    logPGamma();
    updateGammaMask();
}

void SUR_Chain::gammaInit()
{
    arma::umat init = arma::zeros<arma::umat>(nVSPredictors,nOutcomes);
    gammaInit( init );
}

void SUR_Chain::wInit( double w_init , double a_w_ , double b_w_ , double var_w_proposal_ )
{
    w = w_init;
    a_w = a_w_;
    b_w = b_w_;
    
    var_w_proposal = var_w_proposal_ ;
    
    w_acc_count = 0.;
    
    logPW();
}

void SUR_Chain::wInit( double w_init , double a_w_ , double b_w_ )
{
    wInit( w_init , a_w_ , b_w_ , 0.02 );
}

void SUR_Chain::wInit( double w_init )
{
    wInit( w_init , 2. , 5. );
}

void SUR_Chain::wInit()
{
    wInit( 1. , 2. , 5. );
}

void SUR_Chain::w0Init( double w0_init , double a_w0_ , double b_w0_ , double var_w0_proposal_ )
{
  w0 = w0_init;
  a_w0 = a_w0_;
  b_w0 = b_w0_;

  var_w0_proposal = var_w0_proposal_ ;

  w0_acc_count = 0.;
    
  logPW0();
}

void SUR_Chain::w0Init( double w0_init , double a_w0_ , double b_w0_ )
{
  w0Init( w0_init , a_w0_ , b_w0_ , 0.02 );
}

void SUR_Chain::w0Init( double w0_init )
{
  w0Init( w0_init , 2. , 5. );
}

void SUR_Chain::w0Init()
{
  w0Init( 1. , 2. , 5. );
}

void SUR_Chain::betaInit( arma::mat& beta_init )
{
    beta = beta_init;
    logPBeta();
}

void SUR_Chain::betaInit()
{
    arma::mat init = arma::zeros<arma::mat>(nFixedPredictors+nVSPredictors,nOutcomes);
    betaInit( init );
}

// *****************************
// Methods for Log Probabilities
// *****************************

// LOG PRIORS
// logP for all parameters in 3 versions
// empty for re-computing (and updating) the logP given all current values and hyperparameter values
// with one argument for computing with a different value given the current hyperapameter values
// with full arguments for computing the same logP but with different values and different hyperameter values

// TAU
double SUR_Chain::logPTau( double tau_ , double a_tau_ , double b_tau_ )
{
    return Distributions::logPDFGamma( tau_ , a_tau_, b_tau_ );
}

double SUR_Chain::logPTau( )
{
    logP_tau = logPTau( tau , a_tau, b_tau );
    return logP_tau;
}

double SUR_Chain::logPTau( double tau_ )
{
    return logPTau( tau_ , a_tau, b_tau );
}

// ETA
double SUR_Chain::logPEta( double eta_ , double a_eta_ , double b_eta_ )
{
    return Distributions::logPDFBeta( eta_ , a_eta_ , b_eta_ );
}

double SUR_Chain::logPEta( )
{
    logP_eta = logPEta( eta , a_eta, b_eta );
    return logP_eta;
}

double SUR_Chain::logPEta( double eta_ )
{
    return logPEta( eta_ , a_eta, b_eta );
}

// JT
double SUR_Chain::logPJT( const JunctionTree& externalJT , double eta_ )
{
    if ( covariance_type != Covariance_Type::HIW )
        throw Bad_Covariance_Type ( covariance_type );
    
    double logP = 0.;
    for(unsigned int k=0; k<(nOutcomes-1); ++k)
    {
        for(unsigned int l=k+1; l<nOutcomes; ++l)
        {
            logP += Distributions::logPDFBernoulli( externalJT.adjacencyMatrix(k,l) , eta_ );
        }
    }
    
    return logP;
}

double SUR_Chain::logPJT( )
{
    if ( covariance_type != Covariance_Type::HIW )
        throw Bad_Covariance_Type ( covariance_type );
    
    logP_jt = logPJT( jt , eta );
    return logP_jt;
}

double SUR_Chain::logPJT( const JunctionTree& externalJT )
{
    if ( covariance_type != Covariance_Type::HIW )
        throw Bad_Covariance_Type ( covariance_type );
    
    return logPJT( externalJT , eta );
}


// sigma + rhos
double SUR_Chain::logPSigmaRho( const arma::mat&  externalSigmaRho , double nu_ , double tau_ , const JunctionTree& externalJT )
{
    double logP = 0.;
    
    double a,b;
    double thisSigmaTT;
    arma::uvec connectedNodes;
    arma::uvec conditioninIndexes;
    arma::uvec singleIdx_l(1);
    unsigned int nConditioninIndexes;
    
    arma::mat rhoVar; // inverse matrix of the residual elements of Sigma in the component
    arma::rowvec rhoMean; // this is the partial Schur complement, needed for the sampler
    
    std::vector<unsigned int> Prime_q,Res_q, Sep_q;
    unsigned int l;
    
    for( unsigned q=0; q < externalJT.perfectCliqueSequence.size(); ++q )
    {
        Sep_q = externalJT.perfectCliqueSequence[q]->getSeparator();
        Prime_q = externalJT.perfectCliqueSequence[q]->getNodes();
        Res_q.clear();
        std::set_difference(Prime_q.begin(), Prime_q.end(),
                            Sep_q.begin(), Sep_q.end(),
                            std::inserter(Res_q, Res_q.begin()));;
        
        for( unsigned int t=0; t<Res_q.size(); ++t )
        {
            l = Res_q[t];
            singleIdx_l(0) = l;
            
            // start computing interesting things
            thisSigmaTT = tau_;
            
            nConditioninIndexes = Sep_q.size() + t;
            conditioninIndexes.zeros( nConditioninIndexes );
            
            if( nConditioninIndexes > 0 )
            {
                if( Sep_q.size() > 0 )
                    conditioninIndexes(arma::span(0,Sep_q.size()-1)) = arma::conv_to<arma::uvec>::from( Sep_q );
                
                for( unsigned int inner_l=0; inner_l<t; ++inner_l )
                {
                    conditioninIndexes( Sep_q.size() + inner_l ) = Res_q[inner_l];
                }
                
                rhoVar = (1./tau_) * arma::eye(nConditioninIndexes,nConditioninIndexes);
                rhoMean = arma::zeros<arma::rowvec>(nConditioninIndexes); // off-diagonals are zeros
            }
            
            // *** Diagonal Element
            
            // Compute parameters
            a = 0.5 * ( nu_ - nOutcomes + nConditioninIndexes + 1. ) ;
            b = 0.5 * thisSigmaTT ;
            
            logP += Distributions::logPDFIGamma(  externalSigmaRho(l,l), a , b );
            
            // *** Off-Diagonal Element(s)
            if( nConditioninIndexes > 0 )
            {
                
                logP += Distributions::logPDFNormal(  externalSigmaRho( conditioninIndexes , singleIdx_l ) ,
                                                    rhoMean.t() ,  externalSigmaRho(l,l) * rhoVar );
            }
            
        } // end loop over elements inside components
    } // end loop over components
    
    return logP;
}

double SUR_Chain::logPSigmaRho( )
{
    logP_sigmaRho = logPSigmaRho( sigmaRho , nu , tau , jt );
    return logP_sigmaRho;
}

double SUR_Chain::logPSigmaRho( const arma::mat&  externalSigmaRho )
{
    return logPSigmaRho(  externalSigmaRho , nu , tau , jt );
}

// o_k
double SUR_Chain::logPO( const arma::vec& o_ , double a_o_ , double b_o_ )
{
    if ( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    double logP = 0.;
    for(unsigned int k=0; k<nOutcomes; ++k)
        logP += Distributions::logPDFBeta( o_(k) , a_o_, b_o_ );
    
    return logP;
}

double SUR_Chain::logPO( )
{
    if ( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    logP_o = logPO( o , a_o , b_o );
    return logP_o;
}

double SUR_Chain::logPO( const arma::vec& o_ )
{
    if ( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    return logPO( o_ , a_o , b_o );
}


// pi_j
double SUR_Chain::logPPi( arma::vec& pi_ , double a_pi_ , double b_pi_ )
{
    if ( gamma_type != Gamma_Type::hotspot && gamma_type != Gamma_Type::hierarchical )
        throw Bad_Gamma_Type ( gamma_type );
    
    double logP = 0.;
    
    switch ( gamma_type )
    {
        case Gamma_Type::hotspot :
            for(unsigned int j=0; j<nVSPredictors; ++j)
                logP += Distributions::logPDFGamma( pi_(j) , a_pi_, b_pi_ );
            break;
            
        case Gamma_Type::hierarchical :
            for(unsigned int j=0; j<nVSPredictors; ++j)
                logP += Distributions::logPDFBeta( pi_(j) , a_pi_, b_pi_ );
            break;
            
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
    return logP;
}

double SUR_Chain::logPPi( )
{
    if ( gamma_type != Gamma_Type::hotspot && gamma_type != Gamma_Type::hierarchical )
        throw Bad_Gamma_Type ( gamma_type );
    
    logP_pi = logPPi( pi , a_pi , b_pi );
    return logP_pi;
}
double SUR_Chain::logPPi( arma::vec& pi_ )
{
    if ( gamma_type != Gamma_Type::hotspot && gamma_type != Gamma_Type::hierarchical )
        throw Bad_Gamma_Type ( gamma_type );
    
    return logPPi( pi_ , a_pi , b_pi );
}

// GAMMA

// this is the hotspot prior
double SUR_Chain::logPGamma( const arma::umat& externalGamma , const arma::vec& o_ , const arma::vec& pi_ )
{
    if( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    double logP = 0.;
    for(unsigned int j=0; j<nVSPredictors; ++j)
    {
        for(unsigned int k=0; k<nOutcomes; ++k)
        {
            if( ( o_(k) * pi_(j) ) > 1 )
                return -std::numeric_limits<double>::infinity();
            
            logP += Distributions::logPDFBernoulli( externalGamma(j,k), o_(k) * pi_(j) );
        }
    }
    return logP;
}

// this is the simpler hierarchical prior
double SUR_Chain::logPGamma( const arma::umat& externalGamma , const arma::vec& pi_ )
{
    if( gamma_type != Gamma_Type::hierarchical )
        throw Bad_Gamma_Type ( gamma_type );
    double logP = 0.;
    for(unsigned int j=0; j<nVSPredictors; ++j)
    {
        logP += Distributions::logPDFBernoulli( externalGamma.row(j).t() , pi_(j) );
        // logP += Distributions::logPDFBinomial( arma::sum( externalGamma.row(j) ) , nOutcomes , pi_(j) ); // do we care about the binomial coeff? I don't think so..
    }
    return logP;
}

// this is the MRF prior
//double SUR_Chain::logPGamma( const arma::umat& externalGamma , double d, double e, const arma::mat& externalMRFG )
double SUR_Chain::logPGamma( const arma::umat& externalGamma , double d, double e )
{
    if( gamma_type != Gamma_Type::mrf )
        throw Bad_Gamma_Type ( gamma_type );
    
    //arma::mat externalMRFG = mrfG->cols( arma::linspace<arma::uvec>(0,2,3) );
    
    double logP = 0.;
    // calculate the linear and quadratic parts in MRF by using all edges of G
    arma::vec gammaVec = arma::conv_to< arma::vec >::from(arma::vectorise(externalGamma));
    double quad_mrf = 0.;
    double linear_mrf = 0.;
    //int count_linear_mrf = 0; // If the MRF graph matrix has diagonals 0, count_linear_mrf is always 0.
    for( unsigned i=0; i < mrfG->n_rows; ++i )
    {
        if( (*mrfG)(i,0) != (*mrfG)(i,1) ){
            quad_mrf += 2.0 * gammaVec( (*mrfG)(i,0) ) * gammaVec( (*mrfG)(i,1) ) * (*mrfG)(i,2);
        }else{
                if( gammaVec( (*mrfG)(i,0) ) == 1 ){
                    linear_mrf += (*mrfG)(i,2); // should this be 'linear_mrf += e * (externalMRFG)(i,2)'?
                    //count_linear_mrf ++;
                }
        }
    }
    //logP = arma::as_scalar( linear_mrf + d * (arma::accu( externalGamma ) - count_linear_mrf) + e * 2.0 * quad_mrf );
    // Should logP be the following?
    logP = arma::as_scalar( d * arma::accu( externalGamma ) + e * (linear_mrf + quad_mrf) );
    
    return logP;
}

// these below are general interfaces
double SUR_Chain::logPGamma( )
{
    switch ( gamma_type )
    {
        case Gamma_Type::hotspot :
            logP_gamma = logPGamma( gamma , o , pi );
            break;
            
        case Gamma_Type::hierarchical :
            
            logP_gamma = logPGamma( gamma , pi );
            break;
            
        case Gamma_Type::mrf :
        {
            //            logP_gamma = logPGamma( gamma , mrf_d , mrf_e , mrf_G );
            logP_gamma = logPGamma( gamma , mrf_d , mrf_e );
            break;
        }
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
    return logP_gamma;
}

double SUR_Chain::logPGamma( const arma::umat& externalGamma )
{
    double logP {0} ;
    
    switch ( gamma_type )
    {
        case Gamma_Type::hotspot :
            logP = logPGamma( externalGamma , o , pi );
            break;
            
        case Gamma_Type::hierarchical :
            logP = logPGamma( externalGamma , pi );
            break;
            
        case Gamma_Type::mrf :
        {
            //            logP = logPGamma( externalGamma , mrf_d , mrf_e , mrf_G );
            logP = logPGamma( externalGamma , mrf_d , mrf_e );
            break;
        }
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
    
    return logP;
}

// W
double SUR_Chain::logPW( double w_ , double a_w_ , double b_w_ )
{
    return Distributions::logPDFIGamma( w_ , a_w_, b_w_ );
}

double SUR_Chain::logPW( )
{
    logP_w = logPW( w , a_w , b_w );
    return logP_w;
}

double SUR_Chain::logPW( double w_ )
{
    return logPW( w_ , a_w , b_w );
}

// W0
double SUR_Chain::logPW0( double w0_ , double a_w0_ , double b_w0_ )
{
  return Distributions::logPDFIGamma( w0_ , a_w0_, b_w0_ );
}

double SUR_Chain::logPW0( )
{
  logP_w0 = logPW0( w0 , a_w0 , b_w0 );
  return logP_w0;
}

double SUR_Chain::logPW0( double w0_ )
{
  return logPW0( w0_ , a_w0 , b_w0 );
}

// BETA
double logPBetaMaskgPriorK( const arma::vec&  externalBeta_k , const double& w_ , const arma::mat& invXtX_k , const double varianceFactor_k )
{
    return Distributions::logPDFNormal( externalBeta_k , (w_/varianceFactor_k) * invXtX_k );
}


double SUR_Chain::logPBetaMask( const arma::mat&  externalBeta , const arma::umat& mask_ , double w_ , double w0_  )
{
    arma::uvec VS_IN_k;
    arma::uvec singleIdx_k(1);
    
    double logP = 0.;
    
    if(mask_.n_rows > 0)
    {
        
        switch ( beta_type )
        {
            case Beta_Type::gprior :
            {
                arma::uvec xi = arma::conv_to<arma::uvec>::from(jt.perfectEliminationOrder);
                arma::vec xtxMultiplier(nOutcomes);
                
                // prepare posterior full conditional's hyperparameters
                xtxMultiplier(xi(nOutcomes-1)) = 0;
                
                for( unsigned int k=0; k < (nOutcomes-1); ++k)
                {
                    xtxMultiplier(xi(k)) = 0;
                    for(unsigned int l=k+1 ; l<nOutcomes ; ++l)
                    {
                        xtxMultiplier(xi(k)) += pow( sigmaRho(xi(l),xi(k)) , 2 ) / sigmaRho(xi(l),xi(l));
                    }
                }
                
                for(unsigned int k=0; k<nOutcomes ; ++k)
                {
                    singleIdx_k(0) = k;
                    VS_IN_k = mask_( arma::find( mask_.col(1) == k ) , arma::zeros<arma::uvec>(1) );
                    
                    if( preComputedXtX )
                        logP += logPBetaMaskgPriorK( externalBeta(VS_IN_k,singleIdx_k) , w_ ,
                                                    arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) ) , ( 1./ sigmaRho(k,k) + xtxMultiplier(k) ) );
                    else
                        logP += logPBetaMaskgPriorK( externalBeta(VS_IN_k,singleIdx_k) , w_ ,
                                                    arma::inv_sympd( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) ,
                                                    ( 1./ sigmaRho(k,k) + xtxMultiplier(k) ) );
                }
                break;
            }
                
            case Beta_Type::independent :
            {
                for(unsigned int k=0; k<nOutcomes ; ++k)
                {
                    singleIdx_k(0) = k;
                    VS_IN_k = mask_( arma::find( mask_.col(1) == k ) , arma::zeros<arma::uvec>(1) );
                    
                    logP += Distributions::logPDFNormal( externalBeta(VS_IN_k,singleIdx_k) ,
                                                        w_ * arma::eye(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                }
                break;
            }
                
            case Beta_Type::reGroup :
            {
              for(unsigned int k=0; k<nOutcomes ; ++k)
                {
                  singleIdx_k(0) = k;
                  VS_IN_k = mask_( arma::find( mask_.col(1) == k ) , arma::zeros<arma::uvec>(1) );

                  //logP += Distributions::logPDFNormal( externalBeta(VS_IN_k,singleIdx_k) ,
                  //            w_ * arma::eye(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                  logP += Distributions::logPDFNormal( externalBeta(VS_IN_k,singleIdx_k) ,
                                                       arma::diagmat( arma::join_cols(w0_*arma::ones(nFixedPredictors),w_*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) ) );
                }
              break;
            }
                
            default:
                throw Bad_Beta_Type ( beta_type );
                
        }
    }
    return logP;
}

double SUR_Chain::logPBeta( const arma::mat&  externalBeta , const arma::umat& externalGamma , double w_ , double w0_ )
{
    arma::umat mask = createGammaMask( externalGamma );
    return logPBetaMask( externalBeta , mask , w_ , w0_  );
}

double SUR_Chain::logPBeta( )
{
    logP_beta = logPBetaMask( beta , gammaMask , w , w0 );
    return logP_beta;
}

double SUR_Chain::logPBeta( const arma::mat&  externalBeta )
{
    return logPBetaMask(  externalBeta , gammaMask , w , w0 );
}

// PREDICTIVE LIKELIHOODS
arma::mat SUR_Chain::predLikelihood( )
{
    predLik.set_size(nObservations, nOutcomes);
    arma::mat dataOutcome = data->cols( *outcomesIdx );
    
    for( unsigned int k=0; k<nOutcomes; ++k)
        for( unsigned int j=0; j<nObservations; ++j )
        {
            predLik(j,k) = std::exp( Distributions::logPDFNormal( dataOutcome(j,k) , (XB(j,k)+rhoU(j,k)) , sigmaRho(k,k)) );
        }
    
    return predLik;
}

// LOG LIKELIHOODS
double SUR_Chain::logLikelihood( )
{
    double logP = 0.;

    #ifdef _OPENMP
    #pragma omp parallel for default(shared) reduction(+:logP)
    #endif
    
    for( unsigned int k=0; k<nOutcomes; ++k)
    {
        logP += Distributions::logPDFNormal( data->col( (*outcomesIdx)(k) ) , (XB.col(k)+rhoU.col(k)) , sigmaRho(k,k));
    }
    
    logP /= temperature;
    
    log_likelihood = logP; // update internal state
    
    return logP;
}

double SUR_Chain::logLikelihood( const arma::umat&  externalGammaMask , const arma::mat& externalXB ,
                                const arma::mat& externalU , const arma::mat& externalRhoU , const arma::mat&  externalSigmaRho )
{
    double logP = 0.;

    #ifdef _OPENMP
    #pragma omp parallel for default(shared) reduction(+:logP)
    #endif
    
    for( unsigned int k=0; k<nOutcomes; ++k)
    {
        logP += Distributions::logPDFNormal( data->col( (*outcomesIdx)(k) ) , (externalXB.col(k) + externalRhoU.col(k)) ,  externalSigmaRho(k,k));
    }
    
    return logP/temperature;
}

double SUR_Chain::logLikelihood( arma::umat&  externalGammaMask , arma::mat& externalXB , arma::mat& externalU , arma::mat& externalRhoU , //gammaMask,XB,U,rhoU
                                const arma::mat&  externalBeta , const arma::umat& externalGamma , // beta , gamma
                                const arma::mat&  externalSigmaRho , const JunctionTree& externalJT ) // sigmaRho, jt
{
    externalGammaMask = createGammaMask(externalGamma);
    arma::uvec singleIdx_k(1), VS_IN_k;
    
    createQuantities(  externalGammaMask , externalXB , externalU , externalRhoU ,
                     externalGamma ,  externalBeta ,  externalSigmaRho , externalJT );
    
    double logP = 0.;

    #ifdef _OPENMP
    #pragma omp parallel for default(shared) reduction(+:logP)
    #endif
    
    for( unsigned int k=0; k<nOutcomes; ++k)
    {
        logP += Distributions::logPDFNormal( data->col( (*outcomesIdx)(k) ) , ( externalXB.col(k) + externalRhoU.col(k)) ,  externalSigmaRho(k,k));
    }
    
    return logP/temperature;
}


// *********************
// STEP FUNCTION - PERFORM ONE ITERATION FOR THE CHAIN
// *********************


// This function sample sigmas and rhos from their full conditionals and updates the relevant matrix rhoU to reflect thats
double SUR_Chain::sampleSigmaRhoGivenBeta( const arma::mat&  externalBeta , arma::mat& mutantSigmaRho , const JunctionTree& externalJT ,
                                          const arma::umat&  externalGammaMask , const arma::mat& externalXB , const arma::mat& externalU , arma::mat& mutantRhoU )
{
    double logP = 0.;
    
    mutantSigmaRho.zeros(nOutcomes,nOutcomes); // RESET THE WHOLE MATRIX !!!
    
    // hyperparameter of the posterior sampler
    arma::mat Sigma = ( externalU.t() * externalU ) / temperature; Sigma.diag() += tau;
    
    double thisSigmaTT;
    arma::uvec connectedNodes;
    arma::uvec conditioninIndexes;
    unsigned int nConditioninIndexes;
    
    double a,b;
    arma::mat rhoVar; // inverse matrix of the residual elements of Sigma in the component
    arma::rowvec rhoMean; // this is the partial Schur complement, needed for the sampler
    
    switch ( covariance_type )
    {
        case Covariance_Type::HIW :
        {
            arma::uvec singleIdx_l(1); // needed for convention with arma::submat
            std::vector<unsigned int> Prime_q,Res_q, Sep_q;
            unsigned int l;
            for( unsigned q=0; q < externalJT.perfectCliqueSequence.size(); ++q )
            {
                Sep_q = externalJT.perfectCliqueSequence[q]->getSeparator();
                Prime_q = externalJT.perfectCliqueSequence[q]->getNodes();
                Res_q.clear();
                std::set_difference(Prime_q.begin(), Prime_q.end(),
                                    Sep_q.begin(), Sep_q.end(),
                                    std::inserter(Res_q, Res_q.begin()));;
                
                for( unsigned int t=0; t<Res_q.size(); ++t )
                {
                    l = Res_q[t];
                    singleIdx_l(0) = l;
                    
                    // start computing interesting things
                    thisSigmaTT = Sigma(l,l);
                    
                    nConditioninIndexes = Sep_q.size() + t;
                    conditioninIndexes.zeros( nConditioninIndexes );
                    
                    if( nConditioninIndexes > 0 )
                    {
                        
                        if( Sep_q.size() > 0 )
                            conditioninIndexes(arma::span(0,Sep_q.size()-1)) = arma::conv_to<arma::uvec>::from( Sep_q );
                        
                        for( unsigned int inner_l=0; inner_l<t; ++inner_l )
                        {
                            conditioninIndexes( Sep_q.size() + inner_l ) = Res_q[inner_l];
                        }
                        /*test = */arma::inv_sympd( rhoVar , Sigma(conditioninIndexes,conditioninIndexes) ) ;
                        rhoMean = Sigma( singleIdx_l , conditioninIndexes ) * rhoVar ;
                        thisSigmaTT -= arma::as_scalar( rhoMean * Sigma( conditioninIndexes , singleIdx_l ) );
                    }
                    
                    // *** Diagonal Element
                    
                    // Compute parameters
                    a = 0.5 * ( nObservations/temperature + nu - nOutcomes + nConditioninIndexes + 1. ) ;
                    b = 0.5 * thisSigmaTT ;
                    
                    mutantSigmaRho(l,l) = randIGamma( a , b );
                    
                    logP += Distributions::logPDFIGamma( mutantSigmaRho(l,l), a , b );
                    
                    
                    // *** Off-Diagonal Element(s)
                    if( nConditioninIndexes > 0 )
                    {
                        
                        mutantSigmaRho( conditioninIndexes , singleIdx_l ) =
                        Distributions::randMvNormal( rhoMean.t() , mutantSigmaRho(l,l) * rhoVar );
                        
                        mutantSigmaRho( singleIdx_l , conditioninIndexes ) =
                        mutantSigmaRho( conditioninIndexes , singleIdx_l ).t();
                        
                        logP += Distributions::logPDFNormal( mutantSigmaRho( conditioninIndexes , singleIdx_l ) ,
                                                            rhoMean.t() , mutantSigmaRho(l,l) * rhoVar );
                    }
                    
                    // add zeros were set at the beginning with the SigmaRho reset, so no need to act now
                    
                } // end loop over elements inside components
            } // end loop over components
            break;
        }
            
        case Covariance_Type::IW :
        {
            arma::uvec singleIdx_k(1); // needed for convention with arma::submat
            for( unsigned int k=0; k<nOutcomes; ++k )
            {
                singleIdx_k(0) = k;
                
                // start computing interesting things
                thisSigmaTT = Sigma(k,k);
                
                nConditioninIndexes = k;
                conditioninIndexes.zeros( nConditioninIndexes );
                
                if( nConditioninIndexes > 0 )
                {
                    conditioninIndexes = arma::regspace<arma::uvec>(0, k-1);
                    
                    /*test = */arma::inv_sympd( rhoVar , Sigma(conditioninIndexes,conditioninIndexes) ) ;
                    rhoMean = Sigma( singleIdx_k , conditioninIndexes ) * rhoVar ;
                    thisSigmaTT -= arma::as_scalar( rhoMean * Sigma( conditioninIndexes , singleIdx_k ) );
                }
                
                // *** Diagonal Element
                
                // Compute parameters
                a = 0.5 * ( nObservations/temperature + nu - nOutcomes + nConditioninIndexes + 1. ) ;
                b = 0.5 * thisSigmaTT ;
                
                mutantSigmaRho(k,k) = randIGamma( a , b );
                
                logP += Distributions::logPDFIGamma( mutantSigmaRho(k,k), a , b );
                
                
                // *** Off-Diagonal Element(s)
                if( nConditioninIndexes > 0 )
                {
                    
                    mutantSigmaRho( conditioninIndexes , singleIdx_k ) =
                    Distributions::randMvNormal( rhoMean.t() , mutantSigmaRho(k,k) * rhoVar );
                    
                    mutantSigmaRho( singleIdx_k , conditioninIndexes ) =
                    mutantSigmaRho( conditioninIndexes , singleIdx_k ).t();
                    
                    logP += Distributions::logPDFNormal( mutantSigmaRho( conditioninIndexes , singleIdx_k ) ,
                                                        rhoMean.t() , mutantSigmaRho(k,k) * rhoVar );
                }
                
                // add zeros were set at the beginning with the SigmaRho reset, so no need to act now
                
            } // end loop over all elements
            break;
        }
            
        default:
            throw Bad_Covariance_Type ( covariance_type );
    }
    
    // modify useful quantities, only rhoU impacted
    //recompute rhoU as the rhos have changed
    mutantRhoU = createRhoU( externalU , mutantSigmaRho , externalJT );
    
    return logP;
}

double SUR_Chain::sampleBetaGivenSigmaRho( arma::mat& mutantBeta , const arma::mat&  externalSigmaRho , const JunctionTree& externalJT ,
                                          const arma::umat&  externalGammaMask , arma::mat& mutantXB , arma::mat& mutantU , arma::mat& mutantRhoU )
{
    double logP{0.}; // this is the log probability of the proposal
    // double logPrior{0.}; // this is the log prior of the new sample
    // logPrior is wasteful here because we end up not using it as a return value.
    // the prior is updated outside as this function is needed also in the global updates and we
    // don't want to update erroneously the state of a different chain
    
    mutantBeta.set_size(nVSPredictors + nFixedPredictors,nOutcomes); // resize to be sure
    mutantBeta.fill(0.);
    
    if(externalGammaMask.n_rows>0)
    {
        
        arma::vec mu_k; // beta samplers
        arma::mat W_k , iXtX; // indep prior uses W_k, gPrior uses iXtX
        double varianceFactor;
        
        arma::uvec singleIdx_k(1); // needed for convention with arma::submat
        
        arma::vec tmpVec;
        
        arma::uvec VS_IN_k;
        //bool test;
        
        arma::uvec xi = arma::conv_to<arma::uvec>::from(externalJT.perfectEliminationOrder);
        arma::vec xtxMultiplier(nOutcomes);
        arma::mat y_tilde = data->cols( (*outcomesIdx) ) - mutantRhoU ;
        
        // prepare posterior full conditional's hyperparameters
        y_tilde.each_row() /= ( externalSigmaRho.diag().t()) ; // divide each col by the corresponding element of sigma
        xtxMultiplier(xi(nOutcomes-1)) = 0;
        // y_tilde.col(xi(nOutcomes-1)) is already ok;
        
        for( unsigned int k=0; k < (nOutcomes-1); ++k)
        {
            xtxMultiplier(xi(k)) = 0;
            for(unsigned int l=k+1 ; l<nOutcomes ; ++l)
            {
                xtxMultiplier(xi(k)) += pow( externalSigmaRho(xi(l),xi(k)),2) /  externalSigmaRho(xi(l),xi(l));
                y_tilde.col(xi(k)) -= (  externalSigmaRho(xi(l),xi(k)) /  externalSigmaRho(xi(l),xi(l)) ) *
                ( mutantU.col(xi(l)) - mutantRhoU.col(xi(l)) +  externalSigmaRho(xi(l),xi(k)) * ( mutantU.col(xi(k)) - data->col( (*outcomesIdx)(xi(k)) ) ) );
            }
            
        }
        
        // actual sampling
        tmpVec.clear();
        
        // for( unsigned int j : externalJT.perfectEliminationOrder ) //shouldn't make a difference..
        for(unsigned int k=0; k<nOutcomes ; ++k)
        {
            
            VS_IN_k =  externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
            if(VS_IN_k.n_elem>0)
            {
                
                singleIdx_k(0) = k;
                
                switch ( beta_type )
                {
                    case Beta_Type::gprior :
                    {
                        varianceFactor = ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) );
                        
                        if( preComputedXtX )
                        {
                            arma::inv_sympd( iXtX , XtX(VS_IN_k,VS_IN_k) );
                            // W_k = iXtX * ( (w*temperature)/(w + temperature) ) / varianceFactor;
                        }else{
                            arma::inv_sympd( iXtX , data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) );
                            // W_k = iXtX * ( (w*temperature)/(w + temperature) ) / varianceFactor;
                        }
                        
                        mu_k = ( (w*temperature)/(w + temperature) / varianceFactor ) * iXtX *
                        ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * y_tilde.col(k) / temperature );
                        
                        tmpVec = Distributions::randMvNormal( mu_k , ( (w*temperature)/(w + temperature) / varianceFactor ) * iXtX );
                        logP += Distributions::logPDFNormal( tmpVec , mu_k , ( (w*temperature)/(w + temperature) / varianceFactor ) * iXtX );
                        
                        // logPrior += logPBetaMaskgPriorK( tmpVec , w , iXtX , varianceFactor );
                        
                        break;
                    }
                        
                    case Beta_Type::independent :
                    {
                        
                        if( preComputedXtX )
                            arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                        else
                            arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                        
                        mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * y_tilde.col(k) / temperature ) ;
                        
                        tmpVec = Distributions::randMvNormal( mu_k , W_k );
                        logP += Distributions::logPDFNormal( tmpVec , mu_k , W_k );
                        
                        break;
                    }
                        
                    case Beta_Type::reGroup :
                    {
                        
                        if( preComputedXtX )
                            //arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                            arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + arma::diagmat( arma::join_cols((1./w0)*arma::ones(nFixedPredictors),(1./w)*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) ) );
                        else
                          //arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                          arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) +  arma::diagmat( arma::join_cols((1./w0)*arma::ones(nFixedPredictors),(1./w)*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) ) );

                        mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * y_tilde.col(k) / temperature ) ;
                        /*
                        if( k==1 )
                        {
                            std::cout << "...Debug the update of w0 & w: " << w0 << "; " << w << std::endl;
                          std::cout << "...Debug the updated length of mu_1: " << arma::size(mu_k) << std::endl;
                          std::cout << "...Debug the update of b0: " << arma::conv_to<arma::rowvec>::from(beta.submat(0,0,2,0)) << std::endl;
                          std::cout << "...Debug the update of b: " << arma::conv_to<arma::rowvec>::from(beta.submat(nFixedPredictors-1+0,0,nFixedPredictors-1+2,0)) << std::endl;

                        }
                        */
                        tmpVec = Distributions::randMvNormal( mu_k , W_k );
                        logP += Distributions::logPDFNormal( tmpVec , mu_k , W_k );
                        
                        break;
                    }
                        
                    default:
                        throw Bad_Beta_Type ( beta_type );
                }
                
                mutantBeta(VS_IN_k,singleIdx_k) = tmpVec;
                
            }// end if VS_IN_k is non-empty
        }// end foreach outcome
    } // end if mask is non-empty
    
    // if not gPrior, update the prior value here at the end
    // if( beta_type != Beta_Type::gprior )
    //     logPrior = logPBetaMask( mutantBeta , externalGammaMask , w );
    
    // Now the beta have changed so X*B is changed as well as U, compute it to update it for the logLikelihood
    // finally as U changed, rhoU changes as well
    mutantXB = createXB(  externalGammaMask , mutantBeta );
    mutantU = createU( mutantXB );
    mutantRhoU = createRhoU( mutantU ,  externalSigmaRho , externalJT );
    
    return logP;
    
}

double SUR_Chain::sampleBetaKGivenSigmaRho( const unsigned int k , arma::mat& mutantBeta , const arma::mat&  externalSigmaRho , const JunctionTree& externalJT ,
                                           const arma::umat&  externalGammaMask , arma::mat& mutantXB , arma::mat& mutantU , arma::mat& mutantRhoU )
{
    double logP{0.};
    
    mutantBeta.col(k).fill( 0. );
    
    if(externalGammaMask.n_rows>0)
    {
        arma::uvec VS_IN_k;
        VS_IN_k =  externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
        
        if(VS_IN_k.n_elem>0)
        {
            arma::uvec singleIdx_k(1); // needed for convention with arma::submat
            singleIdx_k(0) = k;
            
            arma::vec mu_k; arma::mat W_k; // beta samplers
            arma::vec tmpVec;
            //bool test;
            
            arma::uvec xi = arma::conv_to<arma::uvec>::from(externalJT.perfectEliminationOrder);
            double xtxMultiplier;
            arma::vec y_tilde = data->col( (*outcomesIdx)(k) ) - mutantRhoU.col(k) ;
            
            // prepare posterior full conditional's hyperparameters
            xtxMultiplier = 0;
            y_tilde /=  externalSigmaRho(k,k);
            
            unsigned int k_idx = arma::as_scalar( arma::find( xi == k , 1 ) );
            
            for(unsigned int l=k_idx+1 ; l<nOutcomes ; ++l)
            {
                xtxMultiplier += pow( externalSigmaRho(xi(l),k),2) /  externalSigmaRho(xi(l),xi(l));
                y_tilde -= (  externalSigmaRho(xi(l),k) /  externalSigmaRho(xi(l),xi(l)) ) *
                ( mutantU.col(xi(l)) - mutantRhoU.col(xi(l)) +  externalSigmaRho(xi(l),k) * ( mutantU.col(k) - data->col( (*outcomesIdx)(k) ) ) );
            }
            
            // actual sampling
            tmpVec.clear();
            
            if( preComputedXtX )
            {
                switch ( beta_type )
                {
                    case Beta_Type::gprior :
                    {
                        W_k = (w*temperature)/(w + temperature) * arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) ) / ( 1./ externalSigmaRho(k,k) + xtxMultiplier );
                        break;
                    }
                        
                    case Beta_Type::independent :
                    {
                        arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                        break;
                    }
                        
                    case Beta_Type::reGroup :
                    {
                      //arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                      arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + arma::diagmat( arma::join_cols((1./w0)*arma::ones(nFixedPredictors),(1./w)*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) ) );
                      break;
                    }
                        
                    default:
                        throw Bad_Beta_Type ( beta_type );
                }
                
            }else{
                
                switch ( beta_type )
                {
                    case Beta_Type::gprior :
                    {
                        W_k = (w*temperature)/(w + temperature) * arma::inv_sympd( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) / ( 1./ externalSigmaRho(k,k) + xtxMultiplier );
                        break;
                    }
                        
                    case Beta_Type::independent :
                    {
                        arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                        break;
                    }

                    case Beta_Type::reGroup :
                    {
                        //arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                        arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + arma::diagmat( arma::join_cols((1./w0)*arma::ones(nFixedPredictors),(1./w)*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) ) );
                        break;
                    }
                        
                    default:
                        throw Bad_Beta_Type ( beta_type );
                }
            }
            
            mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * y_tilde / temperature ) ;
            
            tmpVec = Distributions::randMvNormal( mu_k , W_k );
            logP = Distributions::logPDFNormal( tmpVec , mu_k , W_k );
            mutantBeta(VS_IN_k,singleIdx_k) = tmpVec;
            
        } // end if VS_IN_k is non-empty
    } // end if gammaMask is non-empty
    
    // Now the beta have changed so X*B is changed as well as U, compute it to update it for the logLikelihood
    // finally as U changed, rhoU changes as well
    mutantXB = createXB(  externalGammaMask , mutantBeta );
    mutantU = createU( mutantXB );
    mutantRhoU = createRhoU( mutantU ,  externalSigmaRho , externalJT );
    
    return logP;
    
}

//logProbabilities of the above samplers (for the reverse moves)
// this function "simulate" a gibbs move and compute its proposal probability
double SUR_Chain::logPSigmaRhoGivenBeta( const arma::mat&  externalBeta , const arma::mat& mutantSigmaRho , const JunctionTree& externalJT ,
                                        const arma::umat&  externalGammaMask , const arma::mat& externalXB , const arma::mat& externalU , const arma::mat& mutantRhoU )
{
    double logP = 0.;
    
    // hyperparameter of the posterior sampler
    arma::mat Sigma = ( externalU.t() * externalU ) / temperature; Sigma.diag() += tau;
    
    double thisSigmaTT;
    arma::uvec connectedNodes;
    arma::uvec conditioninIndexes;
    unsigned int nConditioninIndexes;
    
    double a,b;
    arma::mat rhoVar; // inverse matrix of the residual elements of Sigma in the component
    arma::rowvec rhoMean; // this is the partial Schur complement, needed for the sampler
    
    switch ( covariance_type )
    {
        case Covariance_Type::HIW :
        {
            arma::uvec singleIdx_l(1); // needed for convention with arma::submat
            std::vector<unsigned int> Prime_q,Res_q, Sep_q;
            unsigned int l;
            
            for( unsigned q=0; q < externalJT.perfectCliqueSequence.size(); ++q )
            {
                Sep_q = externalJT.perfectCliqueSequence[q]->getSeparator();
                Prime_q = externalJT.perfectCliqueSequence[q]->getNodes();
                Res_q.clear();
                std::set_difference(Prime_q.begin(), Prime_q.end(),
                                    Sep_q.begin(), Sep_q.end(),
                                    std::inserter(Res_q, Res_q.begin()));;
                
                for( unsigned int t=0; t<Res_q.size(); ++t )
                {
                    l = Res_q[t];
                    singleIdx_l(0) = l;
                    
                    // start computing interesting things
                    thisSigmaTT = Sigma(l,l);
                    
                    nConditioninIndexes = Sep_q.size() + t;
                    conditioninIndexes.zeros( nConditioninIndexes );
                    
                    if( nConditioninIndexes > 0 )
                    {
                        if( Sep_q.size() > 0 )
                            conditioninIndexes(arma::span(0,Sep_q.size()-1)) = arma::conv_to<arma::uvec>::from( Sep_q );
                        
                        for( unsigned int inner_l=0; inner_l<t; ++inner_l )
                        {
                            conditioninIndexes( Sep_q.size() + inner_l ) = Res_q[inner_l];
                        }
                        /*test = */arma::inv_sympd( rhoVar , Sigma(conditioninIndexes,conditioninIndexes) ) ;
                        rhoMean = Sigma( singleIdx_l , conditioninIndexes ) * rhoVar ;
                        thisSigmaTT -= arma::as_scalar( rhoMean * Sigma( conditioninIndexes , singleIdx_l ) );
                    }
                    
                    // *** Diagonal Element
                    
                    // Compute parameters
                    a = 0.5 * ( nObservations/temperature + nu - nOutcomes + nConditioninIndexes + 1. ) ;
                    b = 0.5 * thisSigmaTT ;
                    
                    logP += Distributions::logPDFIGamma( mutantSigmaRho(l,l), a , b );
                    
                    
                    // *** Off-Diagonal Element(s)
                    if( nConditioninIndexes > 0 )
                    {
                        logP += Distributions::logPDFNormal( mutantSigmaRho( conditioninIndexes , singleIdx_l ) ,
                                                            rhoMean.t() , mutantSigmaRho(l,l) * rhoVar );
                    }
                    
                    // add zeros were set at the beginning with the SigmaRho reset, so no need to act now
                    
                } // end loop over elements inside components
            } // end loop over components
            break;
        }
            
        case Covariance_Type::IW :
        {
            arma::uvec singleIdx_k(1); // needed for convention with arma::submat
            for( unsigned int k=0; k<nOutcomes; ++k )
            {
                singleIdx_k(0) = k;
                
                // start computing interesting things
                thisSigmaTT = Sigma(k,k);
                
                nConditioninIndexes = k;
                conditioninIndexes.zeros( nConditioninIndexes );
                
                if( nConditioninIndexes > 0 )
                {
                    
                    conditioninIndexes = arma::regspace<arma::uvec>(0,k-1);
                    
                    /*test = */arma::inv_sympd( rhoVar , Sigma(conditioninIndexes,conditioninIndexes) ) ;
                    rhoMean = Sigma( singleIdx_k , conditioninIndexes ) * rhoVar ;
                    thisSigmaTT -= arma::as_scalar( rhoMean * Sigma( conditioninIndexes , singleIdx_k ) );
                }
                
                // *** Diagonal Element
                
                // Compute parameters
                a = 0.5 * ( nObservations/temperature + nu - nOutcomes + nConditioninIndexes + 1. ) ;
                b = 0.5 * thisSigmaTT ;
                
                logP += Distributions::logPDFIGamma( mutantSigmaRho(k,k), a , b );
                
                
                // *** Off-Diagonal Element(s)
                if( nConditioninIndexes > 0 )
                {
                    logP += Distributions::logPDFNormal( mutantSigmaRho( conditioninIndexes , singleIdx_k ) ,
                                                        rhoMean.t() , mutantSigmaRho(k,k) * rhoVar );
                }
                
                // add zeros were set at the beginning with the SigmaRho reset, so no need to act now
                
            } // end loop over elements inside components
            break;
        }
            
        default:
            throw Bad_Covariance_Type ( covariance_type );
    }
    return logP;
}

double SUR_Chain::logPBetaGivenSigmaRho( const arma::mat& mutantBeta , const arma::mat&  externalSigmaRho , const JunctionTree& externalJT ,
                                        const arma::umat& externalGammaMask , const arma::mat& mutantXB , const arma::mat& mutantU , const arma::mat& mutantRhoU )
{
    double logP{0.};
    
    if(externalGammaMask.n_rows>0)
    {
        arma::vec mu_k; arma::mat W_k; // beta samplers
        
        arma::uvec singleIdx_k(1); // needed for convention with arma::submat
        
        arma::uvec VS_IN_k;
        //bool test;
        
        arma::uvec xi = arma::conv_to<arma::uvec>::from( externalJT.perfectEliminationOrder );
        arma::vec xtxMultiplier(nOutcomes);
        arma::mat y_tilde = data->cols( *outcomesIdx ) - mutantRhoU ;
        
        // prepare posterior full conditional's hyperparameters
        y_tilde.each_row() /= ( externalSigmaRho.diag().t()) ; // divide each col by the corresponding element of sigma
        xtxMultiplier(xi(nOutcomes-1)) = 0;
        // y_tilde.col(xi(nOutcomes-1)) is already ok;
        
        for( unsigned int k=0; k < (nOutcomes-1); ++k)
        {
            xtxMultiplier(xi(k)) = 0;
            for(unsigned int l=k+1 ; l<nOutcomes ; ++l)
            {
                xtxMultiplier(xi(k)) += pow( externalSigmaRho(xi(l),xi(k)),2) /  externalSigmaRho(xi(l),xi(l));
                y_tilde.col(xi(k)) -= (  externalSigmaRho(xi(l),xi(k)) /  externalSigmaRho(xi(l),xi(l)) ) *
                ( U.col(xi(l)) - mutantRhoU.col(xi(l)) +  externalSigmaRho(xi(l),xi(k)) * ( mutantU.col(xi(k)) - data->col( (*outcomesIdx)(xi(k)) ) ) );
            }
            
        }
        
        // for( unsigned int j : externalJT.perfectEliminationOrder ) //shouldn't make a difference..
        for(unsigned int k=0; k<nOutcomes ; ++k)
        {
            
            VS_IN_k = externalGammaMask( arma::find( externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
            
            if(VS_IN_k.n_elem>0)
            {
                
                singleIdx_k(0) = k;
                
                if( preComputedXtX )
                {
                    switch ( beta_type )
                    {
                        case Beta_Type::gprior :
                        {
                            W_k = (w*temperature)/(w + temperature) * arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) ) / ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) );
                            break;
                        }
                            
                        case Beta_Type::independent :
                        {
                            arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                            break;
                        }
                            
                        case Beta_Type::reGroup :
                        {
                          //arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                          arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + arma::diagmat( arma::join_cols((1./w0)*arma::ones(nFixedPredictors),(1./w)*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) ) );
                              
                          break;
                        }
                            
                        default:
                            throw Bad_Beta_Type ( beta_type );
                    }
                    
                }else{
                    
                    switch ( beta_type )
                    {
                        case Beta_Type::gprior :
                        {
                            W_k = (w*temperature)/(w + temperature) * arma::inv_sympd( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) / ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) );
                            break;
                        }
                            
                        case Beta_Type::independent :
                        {
                            arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                            break;
                        }
                            
                        case Beta_Type::reGroup :
                        {
                          //arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                          arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + arma::diagmat( arma::join_cols((1./w0)*arma::ones(nFixedPredictors),(1./w)*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) ) );
                          break;
                        }
                            
                        default:
                            throw Bad_Beta_Type ( beta_type );
                    }
                }
                
                mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * y_tilde.col(k) / temperature ) ;
                
                logP += Distributions::logPDFNormal( mutantBeta(VS_IN_k,singleIdx_k) , mu_k , W_k );
            } // end if VS_IN_k is non-empty
        } // end for each outcome
    } // end ifmask is non-empty
    
    return logP;
}

double SUR_Chain::logPBetaKGivenSigmaRho( const unsigned int k , const arma::mat& mutantBeta , const arma::mat&  externalSigmaRho , const JunctionTree& externalJT ,
                                         const arma::umat&  externalGammaMask , const arma::mat& mutantXB , const arma::mat& mutantU , const arma::mat& mutantRhoU )
{
    double logP{0.};
    
    if(externalGammaMask.n_rows>0)
    {
        arma::uvec VS_IN_k;
        VS_IN_k =  externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
        
        if(VS_IN_k.n_elem>0)
        {
            arma::vec mu_k; arma::mat W_k; // beta samplers
            arma::vec tmpVec;
            
            arma::uvec singleIdx_k(1); // needed for convention with arma::submat
            singleIdx_k(0) = k;
            
            
            arma::uvec xi = arma::conv_to<arma::uvec>::from(externalJT.perfectEliminationOrder);
            double xtxMultiplier;
            arma::vec y_tilde = data->col( (*outcomesIdx)(k) ) - mutantRhoU.col(k) ;
            
            // prepare posterior full conditional's hyperparameters
            xtxMultiplier = 0;
            y_tilde /=  externalSigmaRho(k,k);
            
            unsigned int k_idx = arma::as_scalar( arma::find( xi == k , 1 ) );
            
            for(unsigned int l=k_idx+1 ; l<nOutcomes ; ++l)
            {
                xtxMultiplier += pow( externalSigmaRho(xi(l),k),2) /  externalSigmaRho(xi(l),xi(l));
                y_tilde -= (  externalSigmaRho(xi(l),k) /  externalSigmaRho(xi(l),xi(l)) ) *
                ( mutantU.col(xi(l)) - mutantRhoU.col(xi(l)) +  externalSigmaRho(xi(l),k) * ( mutantU.col(k) - data->col( (*outcomesIdx)(k) ) ) );
            }
            
            
            if( preComputedXtX )
            {
                switch ( beta_type )
                {
                    case Beta_Type::gprior :
                    {
                        W_k = (w*temperature)/(w + temperature) * arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) ) / ( 1./ externalSigmaRho(k,k) + xtxMultiplier );
                        break;
                    }
                        
                    case Beta_Type::independent :
                    {
                        arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                        break;
                    }
                    
                    case Beta_Type::reGroup :
                    {
                      //arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                      arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + arma::diagmat( arma::join_cols((1./w0)*arma::ones(nFixedPredictors),(1./w)*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) ) );

                      break;
                    }
                        
                    default:
                        throw Bad_Beta_Type ( beta_type );
                }
                
            }else{
                
                switch ( beta_type )
                {
                    case Beta_Type::gprior :
                    {
                        W_k = (w*temperature)/(w + temperature) * arma::inv_sympd( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) / ( 1./ externalSigmaRho(k,k) + xtxMultiplier );
                        break;
                    }
                        
                    case Beta_Type::independent :
                    {
                        arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                        break;
                    }
                        
                    case Beta_Type::reGroup :
                    {
                      //arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                      arma::inv_sympd( W_k ,  ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) +  arma::diagmat( arma::join_cols((1./w0)*arma::ones(nFixedPredictors),(1./w)*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) ) );
                      break;
                    }
                        
                    default:
                        throw Bad_Beta_Type ( beta_type );
                }
            }
            
            mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * y_tilde / temperature ) ;
            
            logP = Distributions::logPDFNormal( mutantBeta(VS_IN_k,singleIdx_k) , mu_k , W_k );
            
        }// end if VS_IN_k is non-empty
    } //end if mask is non-empty
    
    return logP;
}


// sample sigmaRho given Beta or Beta given sigmaRho
// simple interface to gibbs sampling that updates the internal states (given the internal states) from the full conditional
void SUR_Chain::sampleSigmaRhoGivenBeta()
{
    sampleSigmaRhoGivenBeta( beta , sigmaRho , jt , gammaMask , XB , U , rhoU );
}

void SUR_Chain::sampleBetaGivenSigmaRho()
{
    sampleBetaGivenSigmaRho( beta , sigmaRho , jt , gammaMask , XB , U , rhoU );
}

// note that since the above function sample using full-conditional (one outcome conditioned on all the others)
// I'm to update coefficients for one outcome at a time

// sampler for proposed updates on gamma
double SUR_Chain::gammaBanditProposal( arma::umat& mutantGamma , arma::uvec& updateIdx , unsigned int& outcomeUpdateIdx )
{
    
    double logProposalRatio;
    
    // decide on one outcome
    outcomeUpdateIdx = randIntUniform(0,nOutcomes-1);
    
    // Sample Zs (only for relevant outocome)
    for(unsigned int i=0; i<nVSPredictors; ++i)
    {
        banditZeta(i) = randBeta(banditAlpha(i,outcomeUpdateIdx),banditAlpha(i,outcomeUpdateIdx));
    }
    
    // Create mismatch (only for relevant outcome)
    for(unsigned int i=0; i<nVSPredictors; ++i)
    {
        mismatch(i) = (mutantGamma(i,outcomeUpdateIdx)==0)?(banditZeta(i)):(1.-banditZeta(i));   //mismatch
    }
    
    // Normalise
    // mismatch = arma::log(mismatch); //logscale ??? TODO
    // normalised_mismatch = mismatch - Utils::logspace_add(mismatch);
    
    normalised_mismatch = mismatch / arma::as_scalar(arma::sum(mismatch));
    
    if( randU01() < 0.5 )   // one deterministic update
    {
        // Decide which to update
        updateIdx = arma::zeros<arma::uvec>(1);
        updateIdx(0) = Distributions::randWeightedIndexSampleWithoutReplacement(nVSPredictors,normalised_mismatch); // sample the one
        
        // Update
        mutantGamma(updateIdx(0),outcomeUpdateIdx) = 1 - gamma(updateIdx(0),outcomeUpdateIdx); // deterministic, just switch

        // Compute logProposalRatio probabilities
        normalised_mismatch_backwards = mismatch;
        normalised_mismatch_backwards(updateIdx(0)) = 1. - normalised_mismatch_backwards(updateIdx(0)) ;
        
        // normalised_mismatch_backwards = normalised_mismatch_backwards - Utils::logspace_add(normalised_mismatch_backwards);
        normalised_mismatch_backwards = normalised_mismatch_backwards / arma::as_scalar(arma::sum(normalised_mismatch_backwards));
        
        logProposalRatio = ( std::log( normalised_mismatch_backwards(updateIdx(0)) ) ) -
        ( std::log( normalised_mismatch(updateIdx(0)) ) );
        
    }else{
        
        /*
         n_updates_bandit random (bern) updates
         Note that we make use of column indexing here for armadillo matrices
         */

        logProposalRatio = 0.;
        // Decide which to update
        updateIdx = arma::zeros<arma::uvec>(n_updates_bandit);
        updateIdx = Distributions::randWeightedIndexSampleWithoutReplacement(nVSPredictors,normalised_mismatch,n_updates_bandit); // sample n_updates_bandit indexes
        
        normalised_mismatch_backwards = mismatch; // copy for backward proposal
        
        // Update
        for(unsigned int i=0; i<n_updates_bandit; ++i)
        {
            mutantGamma(updateIdx(i),outcomeUpdateIdx) = randBernoulli(banditZeta(updateIdx(i))); // random update
            
            normalised_mismatch_backwards(updateIdx(i)) = 1.- normalised_mismatch_backwards(updateIdx(i));
            
            logProposalRatio += Distributions::logPDFBernoulli(gamma(updateIdx(i),outcomeUpdateIdx),banditZeta(updateIdx(i))) -
            Distributions::logPDFBernoulli(mutantGamma(updateIdx(i),outcomeUpdateIdx),banditZeta(updateIdx(i)));
        }

        // note that above I might be resampling a value equal to the current one, thus not updating da facto ... TODO
        
        // Compute logProposalRatio probabilities
        // normalised_mismatch_backwards = normalised_mismatch_backwards - Utils::logspace_add(normalised_mismatch_backwards);
        normalised_mismatch_backwards = normalised_mismatch_backwards / arma::as_scalar(arma::sum(normalised_mismatch_backwards));
        
        logProposalRatio += Distributions::logPDFWeightedIndexSampleWithoutReplacement(normalised_mismatch_backwards,updateIdx) -
        Distributions::logPDFWeightedIndexSampleWithoutReplacement(normalised_mismatch,updateIdx);
    }
    
    return logProposalRatio; // pass this to the outside
}

double SUR_Chain::gammaMC3Proposal( arma::umat& mutantGamma , arma::uvec& updateIdx  , unsigned int& outcomeUpdateIdx )
{
    updateIdx = arma::uvec(n_updates_MC3);
    
    // decide on one outcome
    outcomeUpdateIdx = randIntUniform(0,nOutcomes-1);
    
    for( unsigned int i=0; i<n_updates_MC3; ++i)
        updateIdx(i) = randIntUniform(0,nVSPredictors-1);    // note that I might be updating multiple times the same coeff
    
    for( auto i : updateIdx)
        mutantGamma(i,outcomeUpdateIdx) = ( randU01() < 0.5)? gamma(i,outcomeUpdateIdx) : 1-gamma(i,outcomeUpdateIdx); // could simply be ( 0.5 ? 1 : 0) ;
    
    return 0. ; // pass this to the outside, it's the (symmetric) logProposalRatio
}


// **************
// **** Methods that update the internal state of their parameter
// **************


// MH update, Normal in the log-space as tau is positive (with gamma prior)
void SUR_Chain::stepTau()
{
    double proposedTau = std::exp( std::log(tau) + randNormal(0.0, var_tau_proposal) );
    
    double proposedTauPrior = logPTau( proposedTau );
    double proposedSigmaRhoPrior = logPSigmaRho( sigmaRho, nu, proposedTau, jt);
    
    double logAccProb = (proposedTauPrior + proposedSigmaRhoPrior) - (logP_tau + logP_sigmaRho);
    
    if( randLogU01() < logAccProb )
    {
        tau = proposedTau;
        logP_tau = proposedTauPrior;
        logP_sigmaRho = proposedSigmaRhoPrior;
        
        ++tau_acc_count;
    }
}

// Here we have a full conditional available
void SUR_Chain::stepEta()
{
    double a = a_eta + 0.5*(arma::accu(jt.getAdjMat())/2. ) ; // divide by temperature if the prior on G is tempered
    double b = b_eta + ( (double)(nOutcomes * (nOutcomes-1.) * 0.5) - 0.5*(arma::accu(jt.getAdjMat())/2. ) ); // /temperature
    
    eta = randBeta( a , b );
    
    logPEta(); // update its prior value
    logPJT(); // update JT's pror value as it's impacted by the new eta
}

void SUR_Chain::stepJT()
{
    JunctionTree proposedJT;
    std::pair<bool,double> updated;
    
    double logProposalRatio, logAccProb = 0.;
    
    bool updateBeta = false; // should I update beta? Mostly debugging code
    
    arma::mat proposedSigmaRho, proposedRhoU;
    arma::mat proposedBeta, proposedXB, proposedU;
    
    double proposedSigmaRhoPrior, proposedJTPrior, proposedLikelihood;
    double proposedBetaPrior;
    
    // we will try to update the JT a few times in order to guarantee some mixing
    // in particular we will do n_updates_jt proposed moves
    // and each time we will select different nodes to join/split untill we find an actual possible move
    // i.e. not counting the trials where we select unsplittable/unjoinable x and y
    
    unsigned int count;
    // arma::uvec updateIdx(2);
    double val = randU01();
    
    for( unsigned int iter=0; iter < n_updates_jt; ++iter )
    {
        count = 0;
        
        jt.copyJT( proposedJT );
        
        do
        {
            if( val < 0.1 )   // 10% of the time, just shuffle the JT
            {
                proposedJT.randomJTPermutation();
                updated = { true , 0.0 }; // logProposal is symmetric
                
            }else if( val < 0.55 )  // 90% of the time, propose an update based on the JT sampler () half from sigle half from multiple e.u.
            {
                updated = proposedJT.propose_multiple_edge_update( );
            }else
            {
                updated = proposedJT.propose_single_edge_update( );
            }
            
        }while( !std::get<0>(updated) && count++ < 100 );
        
        logProposalRatio = std::get<1>(updated);
        
        // note for quantities below. The firt call to sampleXXX has the proposedQuantities set to the current value,
        // for them to be updated; the second call to logPXXX has them updated, needed for the backward probability
        // the main parameter of interest instead "changes to the current value" in the backward equation
        
        if( updateBeta )
        {
            // *********** heavy -- has a more vague graph though, as tau grows more ...
            
            proposedSigmaRho = sigmaRho;
            proposedRhoU = rhoU;
            
            logProposalRatio -= sampleSigmaRhoGivenBeta( beta , proposedSigmaRho , proposedJT ,
                                                        gammaMask , XB , U , proposedRhoU );
            logProposalRatio += logPSigmaRhoGivenBeta( beta , sigmaRho , jt ,
                                                      gammaMask , XB , U , proposedRhoU );
            
            proposedBeta = beta;
            proposedXB = XB;
            proposedU = U;
            
            logProposalRatio -= sampleBetaGivenSigmaRho( proposedBeta , proposedSigmaRho , proposedJT ,
                                                        gammaMask , proposedXB , proposedU , proposedRhoU );
            logProposalRatio += logPBetaGivenSigmaRho( beta , proposedSigmaRho , proposedJT ,
                                                      gammaMask , proposedXB , proposedU , proposedRhoU );
            
            proposedJTPrior = logPJT( proposedJT );
            proposedSigmaRhoPrior = logPSigmaRho( proposedSigmaRho , nu , tau , proposedJT );
            proposedBetaPrior = logPBeta( proposedBeta );
            proposedLikelihood = logLikelihood( gammaMask , proposedXB , proposedU , proposedRhoU , proposedSigmaRho );
            
            logAccProb = logProposalRatio +
            ( proposedJTPrior + proposedSigmaRhoPrior + proposedBetaPrior + proposedLikelihood ) -
            ( logP_jt + logP_sigmaRho + logP_beta + log_likelihood );
            
            // *********** end heavy
            
        }else{
            
            // *************** medium -- this seems good, a bit of a narrow graph wrt heavy (seems to be better) tau a bit shrunk
            
            proposedSigmaRho = sigmaRho;
            proposedRhoU = rhoU;
            
            logProposalRatio -= sampleSigmaRhoGivenBeta( beta , proposedSigmaRho , proposedJT ,
                                                        gammaMask , XB , U , proposedRhoU );
            logProposalRatio += logPSigmaRhoGivenBeta( beta , sigmaRho , jt ,
                                                      gammaMask , XB , U , proposedRhoU );
            
            proposedJTPrior = logPJT( proposedJT );
            proposedSigmaRhoPrior = logPSigmaRho( proposedSigmaRho , nu , tau , proposedJT );
            
            proposedLikelihood = logLikelihood( gammaMask , XB , U , proposedRhoU , proposedSigmaRho );
            
            logAccProb = logProposalRatio +
            ( proposedJTPrior + proposedSigmaRhoPrior + proposedLikelihood ) -
            ( logP_jt + logP_sigmaRho + log_likelihood );
            
            // *********** end medium
            
        }
        
        if( randLogU01() < logAccProb )
        {
            
            jt = proposedJT;
            sigmaRho = proposedSigmaRho;
            
            rhoU = proposedRhoU;
            
            logP_jt = proposedJTPrior;
            logP_sigmaRho = proposedSigmaRhoPrior;
            
            log_likelihood = proposedLikelihood;
            
            if(updateBeta){
                
                beta = proposedBeta;
                
                XB = proposedXB;
                U = proposedU;
                
                logP_beta = proposedBetaPrior;
            }
            
            jt_acc_count += 1./(double)n_updates_jt;
        }
    } // end for n_updates_jt
}

// MH update (log-normal)
// it might seems wastefull to recompute every logP_jk, but at the end you'd need to compute column J for the new value,
// recopute it for the old value to compare and then recompute the whole thing to replace the whole prior if accepted!
// The only solution would be to keep track of the whole matrix of values, funny idea but is it worth it memory-wise ?
void SUR_Chain::stepOneO()
{
    
    unsigned int k = randIntUniform(0,nOutcomes-1);
    arma::vec proposedO = o;
    
    double proposedOPrior, proposedGammaPrior, logAccProb;
    
    proposedO(k) = std::exp( std::log( o(k) ) + Distributions::randTruncNorm(0.0, var_o_proposal , -std::numeric_limits<double>::infinity() , -std::log( o(k) ) ) );
    
    if( arma::all( ( pi * proposedO(k) ) <= 1 ) )
    {
        proposedOPrior = logPO( proposedO );
        proposedGammaPrior = logPGamma( gamma, proposedO, pi);
        
        // A/R
        logAccProb = Distributions::logPDFTruncNorm( std::log( o(k) ) , std::log( proposedO(k) ) , var_o_proposal , -std::numeric_limits<double>::infinity() , -std::log( proposedO(k) ) ) -
        Distributions::logPDFTruncNorm( std::log( proposedO(k) ) , std::log( o(k) ) , var_o_proposal , -std::numeric_limits<double>::infinity() , -std::log( o(k) ) );
        logAccProb += (proposedOPrior + proposedGammaPrior) - (logP_o + logP_gamma);
        
        if( randLogU01() < logAccProb )
        {
            o(k) = proposedO(k);
            logP_o = proposedOPrior;
            logP_gamma = proposedGammaPrior;
            
            ++o_acc_count;
        }
    }
    
}

void SUR_Chain::stepO()
{
    
    arma::vec proposedO = o;
    
    double proposedOPrior, proposedGammaPrior, logAccProb;
    
    for( unsigned int k=0; k < nOutcomes ; ++k )
    {
        proposedO(k) = std::exp( std::log( o(k) ) + Distributions::randTruncNorm(0.0, var_o_proposal , -std::numeric_limits<double>::infinity() , -std::log( o(k) ) ) );
        
        if( arma::all( ( pi * proposedO(k) ) <= 1 ) )
        {
            proposedOPrior = logPO( proposedO );
            proposedGammaPrior = logPGamma( gamma, proposedO, pi);
            
            // A/R
            logAccProb = Distributions::logPDFTruncNorm( std::log( o(k) ) , std::log( proposedO(k) ) , var_o_proposal , -std::numeric_limits<double>::infinity() , -std::log( proposedO(k) ) ) -
            Distributions::logPDFTruncNorm( std::log( proposedO(k) ) , std::log( o(k) ) , var_o_proposal , -std::numeric_limits<double>::infinity() , -std::log( o(k) ) );
            logAccProb += (proposedOPrior + proposedGammaPrior) - (logP_o + logP_gamma);
            
            if( randLogU01() < logAccProb )
            {
                o(k) = proposedO(k);
                logP_o = proposedOPrior;
                logP_gamma = proposedGammaPrior;
                
                o_acc_count += o_acc_count / (double)nOutcomes;
            }else
                proposedO(k) = o(k);
        }else
            proposedO(k) = o(k);
    }
    
}

// MH update (log-normal) -- update one value at each iteration TODO worth doing more?
void SUR_Chain::stepOnePi()
{
    unsigned int j = randIntUniform(0,nVSPredictors-1);
    
    switch ( gamma_type )
    {
        case Gamma_Type::hotspot :
        {
            arma::vec proposedPi = pi;
            double proposedPiPrior, proposedGammaPrior, logAccProb;
            
            proposedPi(j) = std::exp( std::log( pi(j) ) + randNormal(0.0, var_pi_proposal) );
            
            if( arma::all( ( o * proposedPi(j) ) <= 1 ) )
            {
                proposedPiPrior = logPPi( proposedPi );
                proposedGammaPrior = logPGamma( gamma, o, proposedPi);
                
                // A/R
                logAccProb = (proposedPiPrior + proposedGammaPrior) - (logP_pi + logP_gamma);
                
                if( randLogU01() < logAccProb )
                {
                    pi(j) = proposedPi(j);
                    logP_pi = proposedPiPrior;
                    logP_gamma = proposedGammaPrior;
                    
                    ++pi_acc_count;
                }
            }
            break;
        }
            
        case Gamma_Type::hierarchical : // in this case it's conjugate
        {
            unsigned int k = arma::sum( gamma.row(j) );
            pi(j) = randBeta( a_pi + k , b_pi + nOutcomes - k );
            break;
        }
            
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
    
}

void SUR_Chain::stepPi()
{
    switch ( gamma_type )
    {
        case Gamma_Type::hotspot :
        {
            arma::vec proposedPi = pi;
            double proposedPiPrior, proposedGammaPrior, logAccProb;
            for( unsigned int j=0; j < nVSPredictors ; ++j )
            {
                proposedPi(j) = std::exp( std::log( pi(j) ) + randNormal(0.0, var_pi_proposal) );
                
                if( arma::all( ( o * proposedPi(j) ) <= 1 ) )
                {
                    proposedPiPrior = logPPi( proposedPi );
                    proposedGammaPrior = logPGamma( gamma, o, proposedPi);
                    
                    // A/R
                    logAccProb = (proposedPiPrior + proposedGammaPrior) - (logP_pi + logP_gamma);
                    
                    if( randLogU01() < logAccProb )
                    {
                        pi(j) = proposedPi(j);
                        logP_pi = proposedPiPrior;
                        logP_gamma = proposedGammaPrior;
                        
                        pi_acc_count += pi_acc_count / (double)nVSPredictors;
                    }else
                        proposedPi(j) = pi(j);
                }else
                    proposedPi(j) = pi(j);
            }
            break;
        }
            
        case Gamma_Type::hierarchical : // in this case it's conjugate
        {
            for( unsigned int j=0; j < nVSPredictors ; ++j )
            {
                unsigned int k = arma::sum( gamma.row(j) );
                pi(j) = randBeta( a_pi + k , b_pi + nOutcomes - k );
            }
            break;
        }
            
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
    
}

void SUR_Chain::stepW()
{
    switch ( beta_type )
    {
        case Beta_Type::gprior :
        {
            stepWMH();
            break;
        }
            
        case Beta_Type::independent :
        {
            stepWGibbs();
            break;
        }
            
        case Beta_Type::reGroup :
        {
          stepW0Gibbs();
          break;
        }
            
        default:
            throw Bad_Beta_Type ( beta_type );
    }
}


// Gibbs sampler available here again for w given all the current betas and the gammas -- TODO keep an eye on this
void SUR_Chain::stepWGibbs()
{
    double a = a_w + 0.5*( /*arma::accu(gamma) + intercept */ /*or*/ gammaMask.n_rows ); // divide by temperature if the prior on gamma is tempered
    double b = b_w + 0.5*( arma::accu( arma::square(arma::nonzeros(beta)) ) );   // all the beta_jk w/ gamma_jk=0 are 0 already // /temperature
    
    w = randIGamma( a , b );
    
    logPW(); // update its prior value
    logPBeta(); // update beta's log prior as it's impacted by the change in w
}

void SUR_Chain::stepW0Gibbs()
{
    double a = a_w + 0.5*( gammaMask.n_rows - nFixedPredictors*nOutcomes ); // divide by temperature if the prior on gamma is tempered
    double b = b_w + 0.5*( arma::accu( arma::square(beta.rows(nFixedPredictors,nFixedPredictors+nVSPredictors-1)) ) );   // all the beta_jk w/ gamma_jk=0 are 0 already // /temperature

    // std::cout << a_w << " -> " << a << "   ---   "<< b_w << " -> " << b << std::endl;
    // std::cout << arma::nonzeros(beta).t() << std::endl; std::cin >> w;

    w = randIGamma( a , b );
    logPW(); // update its prior value
    
    double a0 = a_w0 + 0.5*( nFixedPredictors*nOutcomes ); // divide by temperature if the prior on gamma is tempered
    double b0 = b_w0 + 0.5*( arma::accu( arma::square(beta.rows(0,nFixedPredictors-1)) ) );   // all the beta_jk w/ gamma_jk=0 are 0 already // /temperature
    w0 = randIGamma( a0 , b0 );
    logPW0(); // update its prior value

    logPBeta(); // update beta's log prior as it's impacted by the change in w
}

void SUR_Chain::stepWMH()
{
    double proposedW = std::exp( std::log(w) + randNormal(0.0, var_w_proposal) );
    double proposedW0 = std::exp( std::log(w0) + randNormal(0.0, var_w0_proposal) );
    
    double proposedWPrior = logPW( proposedW );
    double proposedW0Prior = logPW0( proposedW0 );
    double proposedBetaPrior = logPBetaMask( beta , gammaMask , proposedW , proposedW0 );
    
    double logAccProb = (proposedWPrior + proposedW0Prior + proposedBetaPrior) - (logP_w + logP_w0 + logP_beta);
    
    if( randLogU01() < logAccProb )
    {
        w = proposedW;
        w0 = proposedW0;
        logP_w = proposedWPrior;
        logP_w0 = proposedW0Prior;
        logP_beta = proposedBetaPrior;
        
        ++w_acc_count;
    }
}


void SUR_Chain::stepGamma()
{
    arma::umat proposedGamma = gamma;
    arma::uvec updateIdx;
    unsigned int outcomeUpdateIdx;
    
    double logProposalRatio = 0;
    
    // Update the proposed Gamma
    switch ( gamma_sampler_type )
    {
        case Gamma_Sampler_Type::bandit :
            logProposalRatio += gammaBanditProposal( proposedGamma , updateIdx , outcomeUpdateIdx );
            break;
            
        case Gamma_Sampler_Type::mc3 :
            logProposalRatio += gammaMC3Proposal( proposedGamma , updateIdx , outcomeUpdateIdx );
            break;
            
        default:
            break;
    }
    // given proposedGamma now, sample a new proposedBeta matrix and corresponging quantities
    arma::umat proposedGammaMask = createGammaMask( proposedGamma );
     
    // note only one outcome is updated
    // note for quantities below. The firt call to sampleXXX has the proposedQuantities set to the current value,
    // for them to be updated; the second call to logPXXX has them updated, needed for the backward probability
    // the main parameter of interest instead "changes to the current value" in the backward equation
    arma::mat proposedBeta = beta;
    
    arma::mat proposedXB = XB;
    arma::mat proposedU = U;
    arma::mat proposedRhoU = rhoU;
    
    logProposalRatio -= sampleBetaKGivenSigmaRho( outcomeUpdateIdx , proposedBeta , sigmaRho , jt ,
                                                 proposedGammaMask , proposedXB , proposedU , proposedRhoU );
    logProposalRatio += logPBetaKGivenSigmaRho( outcomeUpdateIdx , beta , sigmaRho , jt ,
                                               gammaMask , proposedXB , proposedU , proposedRhoU );
     
    // update log probabilities
    double proposedGammaPrior = logPGamma( proposedGamma );
    double proposedBetaPrior = logPBetaMask( proposedBeta , proposedGammaMask , w , w0 );
    double proposedLikelihood = logLikelihood( proposedGammaMask , proposedXB , proposedU , proposedRhoU , sigmaRho );
    
    double logAccProb = logProposalRatio +
    ( proposedGammaPrior + proposedBetaPrior + proposedLikelihood ) -
    ( logP_gamma + logP_beta + log_likelihood );
    
    if( randLogU01() < logAccProb )
    {
        gamma = proposedGamma;
        beta = proposedBeta;
        
        gammaMask = proposedGammaMask;
        XB = proposedXB;
        U = proposedU;
        rhoU = proposedRhoU;
        
        logP_gamma = proposedGammaPrior;
        logP_beta = proposedBetaPrior;
        log_likelihood = proposedLikelihood;
        
        // ++gamma_acc_count;
        gamma_acc_count += 1. ; // / updatedOutcomesIdx.n_elem
    }
     
    // after A/R, update bandit Related variables
    if( gamma_sampler_type == Gamma_Sampler_Type::bandit )
    {
        for(arma::uvec::iterator iter = updateIdx.begin(); iter != updateIdx.end(); ++iter)
        {
            // FINITE UPDATE
            if( banditAlpha(*iter,outcomeUpdateIdx) + banditBeta(*iter,outcomeUpdateIdx) <= banditLimit )
            {
                banditAlpha(*iter,outcomeUpdateIdx) += banditIncrement * gamma(*iter,outcomeUpdateIdx);
                banditBeta(*iter,outcomeUpdateIdx) += banditIncrement * (1-gamma(*iter,outcomeUpdateIdx));
            }
            
            // // CONTINUOUS UPDATE, alternative to the above, at most one has to be uncommented
            
            // banditAlpha(*iter,outcomeUpdateIdx) += banditIncrement * gamma(*iter,outcomeUpdateIdx);
            // banditBeta(*iter,outcomeUpdateIdx) += banditIncrement * (1-gamma(*iter,outcomeUpdateIdx));
            
            // // renormalise
            // if( banditAlpha(*iter,outcomeUpdateIdx) + banditBeta(*iter) > banditLimit )
            // {
            //     banditAlpha(*iter,outcomeUpdateIdx) = banditLimit * ( banditAlpha(*iter,outcomeUpdateIdx) / ( banditAlpha(*iter,outcomeUpdateIdx) + banditBeta(*iter,outcomeUpdateIdx) ));
            //     banditBeta(*iter,outcomeUpdateIdx) = banditLimit * (1. - ( banditAlpha(*iter,outcomeUpdateIdx) / ( banditAlpha(*iter,outcomeUpdateIdx) + banditBeta(*iter,outcomeUpdateIdx) )) );
            // }
            
        }
    }
}

void SUR_Chain::stepSigmaRhoAndBeta()
{
    sampleSigmaRhoGivenBeta();
    sampleBetaGivenSigmaRho();
    
    logPSigmaRho();
    logPBeta();
    logLikelihood();
}


// this updates all the internal states
void SUR_Chain::step()
{
    updateGammaMask();
    // update logP_gamma
    //logPGamma(); 
    
    // Update HyperParameters
    stepTau();
    stepW();
    
    switch ( gamma_type )
    {
        case Gamma_Type::hotspot :
            for( auto i=0; i<5; ++i)
            {
                stepOneO();
                stepOnePi();
            }
            break;
            
        case Gamma_Type::hierarchical :
            for( auto i=0; i<5; ++i)
                stepOnePi();
            break;
            
        case Gamma_Type::mrf :
            break; // nothing to do for this one yet
            
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
    
    // update log_likelihood
    //logLikelihood();
    
    if ( covariance_type == Covariance_Type::HIW )
    {
        stepEta();
        // Update JT
        if( internalIterationCounter >= jtStartIteration )
            stepJT();
    }
     
    // Update Sigmas, Rhos and Betas given all rest
    stepSigmaRhoAndBeta();
    
    // update logP_gamma; this might be redundant?
    logPGamma();
    // update gamma
    stepGamma();
     
    // increase iteration counter
    ++ internalIterationCounter;
    
    // update the MH proposal variance
    updateProposalVariances();
}

// update all the internal proposal RW variances based on their acceptance rate
// note that for o and pi, even though they're vectors, the proposal is independent and common for all elements
// see Roberts  & Rosenthal 2008
void SUR_Chain::updateProposalVariances()
{
    double delta, delta2;
    arma::vec deltaVec, delta2Vec;
    
    double adaptationFactor = 0.05;
    
    if( internalIterationCounter == 1 ) // init the mean and second moment
    {
        tauEmpiricalMean = std::log(tau);
        tauEmpiricalM2 = 0.;
        var_tau_proposal_init = var_tau_proposal;
        
        if ( gamma_type == Gamma_Type::hotspot )
        {
            oEmpiricalMean = arma::log(o);
            oEmpiricalM2 = arma::zeros<arma::vec>(nOutcomes);
            var_o_proposal_init = var_o_proposal;
            
            piEmpiricalMean = arma::log(pi);
            piEmpiricalM2 = arma::zeros<arma::vec>(nVSPredictors);
            var_pi_proposal_init = var_pi_proposal;
        }
        
        if( beta_type == Beta_Type::gprior )
        {
            wEmpiricalMean = w;
            wEmpiricalM2 = 0.;
            var_w_proposal_init = var_w_proposal;
        }
        
    }else if( internalIterationCounter > 1 )
    {
        // update running averages
        
        // tau
        delta = std::log(tau) - tauEmpiricalMean;
        tauEmpiricalMean = tauEmpiricalMean + ( delta / internalIterationCounter );
        delta2 = std::log(tau) - tauEmpiricalMean;
        tauEmpiricalM2 = tauEmpiricalM2 + delta * delta2;
        
        if ( gamma_type == Gamma_Type::hotspot )
        {
            // o
            deltaVec = arma::log(o) - oEmpiricalMean;
            oEmpiricalMean = oEmpiricalMean + ( deltaVec / internalIterationCounter );
            delta2Vec = arma::log(o) - oEmpiricalMean;
            oEmpiricalM2 = oEmpiricalM2 + deltaVec % delta2Vec ;
            
            // pi
            deltaVec = arma::log(pi) - piEmpiricalMean;
            piEmpiricalMean = piEmpiricalMean + ( deltaVec  / internalIterationCounter );
            delta2Vec = arma::log(pi) - piEmpiricalMean;
            piEmpiricalM2 = piEmpiricalM2 + deltaVec % delta2Vec ;
        }
        
        // w
        if( beta_type == Beta_Type::gprior )
        {
            delta = w - wEmpiricalMean;
            wEmpiricalMean = wEmpiricalMean + ( delta / internalIterationCounter );
            delta2 = w - wEmpiricalMean;
            wEmpiricalM2 = wEmpiricalM2 + delta * delta2;
        }
    }
    
    // Then if it's actually > n update the proposal variances
    
    if( internalIterationCounter > nObservations )  // update quantities and variance
    {
        
        // update proposal variances
        var_tau_proposal = adaptationFactor * var_tau_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * tauEmpiricalM2/(internalIterationCounter-1);
        if ( gamma_type == Gamma_Type::hotspot )
        {
            var_o_proposal = adaptationFactor * var_o_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * arma::mean( oEmpiricalM2/(internalIterationCounter-1) );
            var_pi_proposal = adaptationFactor * var_pi_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * arma::mean( piEmpiricalM2/(internalIterationCounter-1) );
        }
        if( beta_type == Beta_Type::gprior )
            var_w_proposal = adaptationFactor * var_w_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * wEmpiricalM2/(internalIterationCounter-1);
    }
}

// more complex functions could be defined outide through public methods
// but this as a baseline is good to have.


// *************************************
// Global operators between two chains
// *************************************
// asuming nu and other fixed hyperparameters are the same across chains, woudn;t make sense otherwise I think

void SUR_Chain::swapTau( std::shared_ptr<SUR_Chain>& that )
{
    double par = this->getTau();
    
    this->setTau( that->getTau() );
    that->setTau( par );
}

void SUR_Chain::swapEta( std::shared_ptr<SUR_Chain>& that )
{
    double par = this->getEta();
    
    this->setEta( that->getEta() );
    that->setEta( par );
}

void SUR_Chain::swapJT( std::shared_ptr<SUR_Chain>& that )
{
    JunctionTree par = this->getJT();
    
    this->setJT( that->getJT() );
    that->setJT( par );
}

void SUR_Chain::swapSigmaRho( std::shared_ptr<SUR_Chain>& that )
{
    arma::mat par = this->getSigmaRho();
    
    this->setSigmaRho( that->getSigmaRho() );
    that->setSigmaRho( par );
}

void SUR_Chain::swapO( std::shared_ptr<SUR_Chain>& that )
{
    arma::vec par = this->getO();
    
    this->setO( that->getO() );
    that->setO( par );
}

void SUR_Chain::swapPi( std::shared_ptr<SUR_Chain>& that )
{
    arma::vec par = this->getPi();
    
    this->setPi( that->getPi() );
    that->setPi( par );
}

void SUR_Chain::swapGamma( std::shared_ptr<SUR_Chain>& that )
{
    arma::umat par = this->getGamma();
    
    this->setGamma( that->getGamma() );
    that->setGamma( par );
}

void SUR_Chain::swapW( std::shared_ptr<SUR_Chain>& that )
{
    double par = this->getW();
    
    this->setW( that->getW() );
    that->setW( par );
}

void SUR_Chain::swapW0( std::shared_ptr<SUR_Chain>& that )
{
    double par = this->getW0();
    
    this->setW0( that->getW0() );
    that->setW0( par );
}

void SUR_Chain::swapBeta( std::shared_ptr<SUR_Chain>& that )
{
    arma::mat par = this->getBeta();
    
    this->setBeta( that->getBeta() );
    that->setBeta( par );
}

int SUR_Chain::exchangeGamma_step( std::shared_ptr<SUR_Chain>& that )
{
    // I'm exchanging the gammas AND the betas. So gammaMask, XB and U will follow and we will have to re-compute rhoU for both chains
    arma::umat swapGammaMask;
    arma::mat swapXB , swapU;
    
    arma::mat rhoU_1 = this->createRhoU( that->getU() , this->getSigmaRho() , this->getJT() );
    arma::mat rhoU_2 = that->createRhoU( this->getU() , that->getSigmaRho() , that->getJT() );
    
    double logLik_1 = this->logLikelihood( that->getGammaMask() , that->getXB() , that->getU() , rhoU_1 , this->getSigmaRho() );   // note that this and that lik are
    double logLik_2 = that->logLikelihood( this->getGammaMask() , this->getXB() , this->getU() , rhoU_2 , that->getSigmaRho() ); // important because of temperature
    
    double logPExchange = ( logLik_1 + logLik_2 ) -
    ( this->getLogLikelihood() + that->getLogLikelihood() );
    
    if( randLogU01() < logPExchange )
    {
        // parameters and priors
        this->swapGamma( that );
        this->swapBeta( that );
        
        // loglikelihood and related quantities
        swapGammaMask = this->getGammaMask() ;
        swapXB = this->getXB();
        swapU = this->getU();
        
        this->setGammaMask( that->getGammaMask() );
        this->setXB( that->getXB() );
        this->setU( that->getU() );
        
        that->setGammaMask( swapGammaMask );
        that->setXB( swapXB );
        that->setU( swapU );
        
        this->setRhoU( rhoU_1 );
        that->setRhoU( rhoU_2 );
        
        this->setLogLikelihood( logLik_1 );
        that->setLogLikelihood( logLik_2 );
        
        return 1;
    }else
        return 0;
}

int SUR_Chain::exchangeJT_step( std::shared_ptr<SUR_Chain>& that )
{
    // I'm exchanging the JT, sigmas and rhos.
    // So gammaMask, XB and U will stay the same and we will have to to re-compute rhoU for both chains
    arma::mat rhoU_1 = this->createRhoU( this->getU() , that->getSigmaRho() , that->getJT() );
    arma::mat rhoU_2 = that->createRhoU( that->getU() , this->getSigmaRho() , this->getJT() );
    
    double logLik_1 = this->logLikelihood( this->getGammaMask() , this->getXB() , this->getU() , rhoU_1 , that->getSigmaRho() ); // note that this and that lik are
    double logLik_2 = that->logLikelihood( that->getGammaMask() , that->getXB() , that->getU() , rhoU_2 , this->getSigmaRho() );   // important because of temperature
    
    double logPExchange = ( logLik_1 + logLik_2 ) -
    ( this->getLogLikelihood() + that->getLogLikelihood() );
    
    if( randLogU01() < logPExchange )
    {
        // parameters and priors
        this->swapJT( that );
        this->swapSigmaRho( that );
        
        // loglikelihood and related quantities
        this->setRhoU( rhoU_1 );
        that->setRhoU( rhoU_2 );
        
        this->setLogLikelihood( logLik_1 );
        that->setLogLikelihood( logLik_2 );
        
        return 1;
    }else
        return 0;
}

int SUR_Chain::adapt_crossOver_step( std::shared_ptr<SUR_Chain>& that )
{
    double pCrossOver;
    
    // Crossover operator hyper pars (see http://www3.stat.sinica.edu.tw/statistica/oldpdf/A10n21.pdf
    double pXO_0 = 0.1, pXO_1 = 0.2 , pXO_2 = 0.2;
    double p11 = pXO_0*pXO_0 + (1.-pXO_0)*(1.-pXO_0) ,p12 = 2.*pXO_0*(1.-pXO_0) ,p21= pXO_1*(1.-pXO_2) + pXO_2*(1.-pXO_1) ,p22 = pXO_1*pXO_2 + (1.-pXO_1)*(1.-pXO_2);
    
    unsigned int n11,n12,n21,n22;
    
    std::vector<arma::umat> gammaXO(2); gammaXO[0] = arma::umat(nVSPredictors,nOutcomes);  gammaXO[1] = arma::umat(nVSPredictors,nOutcomes);
    
    // Propose Crossover
    n11=0;n12=0;n21=0;n22=0;
    
    for(unsigned int j=0; j<nVSPredictors; ++j)
    {
        for(unsigned int k=0; k<nOutcomes; ++k)
        {
            if ( this->getGamma()(j,k) == that->getGamma()(j,k) )
            {
                gammaXO[0](j,k) = this->getGamma()(j,k);
                gammaXO[1](j,k) = this->getGamma()(j,k);
                
                gammaXO[0](j,k) = ( randU01() < pXO_0 )? 1-gammaXO[0](j,k) : gammaXO[0](j,k);
                gammaXO[1](j,k) = ( randU01() < pXO_0 )? 1-gammaXO[1](j,k) : gammaXO[1](j,k);
                
                if( gammaXO[0](j,k) == gammaXO[1](j,k) )
                    ++n11;
                else
                    ++n12;
            }
            else
            {
                gammaXO[0](j,k) = this->getGamma()(j,k);
                gammaXO[1](j,k) = that->getGamma()(j,k);
                
                gammaXO[0](j,k) = ( randU01() < pXO_1 )? 1-gammaXO[0](j,k) : gammaXO[0](j,k);
                gammaXO[1](j,k) = ( randU01() < pXO_2 )? 1-gammaXO[1](j,k) : gammaXO[1](j,k);
                
                if( gammaXO[0](j,k) == gammaXO[1](j,k) )
                    ++n21;
                else
                    ++n22;
            }
        }
    }
    
    pCrossOver = (n11 * std::log( p11 ) + n12 * std::log( p12 ) + n21 * std::log( p21 ) + n22 * std::log( p22 ) )-  // CrossOver proposal probability FORWARD
    (n11 * std::log( p11 ) + n12 * std::log( p21 ) + n21 * std::log( p12 ) + n22 * std::log( p22 ) );  // XO prop probability backward (note that ns stays the same but changes associated prob)
    
    // Propose betas that go with the new crossed-over states
    std::vector<arma::mat> betaXO(2);
    betaXO[0] = this->getBeta();
    betaXO[1] = that->getBeta();
    
    std::vector<arma::umat> gammaMask_XO(2);
    gammaMask_XO[0] = createGammaMask(gammaXO[0]);
    gammaMask_XO[1] = createGammaMask(gammaXO[1]);
    
    std::vector<arma::mat> XB_XO(2);
    XB_XO[0] = arma::mat( this->getXB() ); // copy into the new one
    XB_XO[1] = arma::mat( that->getXB() );
    std::vector<arma::mat> U_XO(2);
    U_XO[0] = arma::mat( this->getU() ); // copy into the new one
    U_XO[1] = arma::mat( that->getU() );
    std::vector<arma::mat> rhoU_XO(2);
    rhoU_XO[0] = arma::mat( this->getRhoU() ); // copy into the new one
    rhoU_XO[1] = arma::mat( that->getRhoU() );
    
    // propose for the first new chain
    pCrossOver -= this->sampleBetaGivenSigmaRho( betaXO[0] , this->getSigmaRho() , this->getJT() ,
                                                gammaMask_XO[0] , XB_XO[0] , U_XO[0] , rhoU_XO[0] );
    
    pCrossOver += this->logPBetaGivenSigmaRho( this->getBeta() , this->getSigmaRho() , this->getJT() ,
                                              this->getGammaMask() , XB_XO[0] , U_XO[0] , rhoU_XO[0] );
    
    // propose for the second new chain
    pCrossOver -= that->sampleBetaGivenSigmaRho( betaXO[1] , that->getSigmaRho() , that->getJT() ,
                                                gammaMask_XO[1] , XB_XO[1] , U_XO[1] , rhoU_XO[1] );
    
    pCrossOver += that->logPBetaGivenSigmaRho( that->getBeta() , that->getSigmaRho() , that->getJT() ,
                                              that->getGammaMask() , XB_XO[1] , U_XO[1] , rhoU_XO[1] );
    
    // log Posterior of the new chains
    double logLikFirst = this->logLikelihood( gammaMask_XO[0] , XB_XO[0] , U_XO[0] , rhoU_XO[0] , this->getSigmaRho() );
    double logLikSecond = that->logLikelihood( gammaMask_XO[1] , XB_XO[1] , U_XO[1] , rhoU_XO[1] , that->getSigmaRho() );
    
    double logPBetaFirst = this->logPBetaMask( betaXO[0] , gammaMask_XO[0] , this->getW() , this->getW0() );  // this and that are important for temperature and hyperparameters
    double logPBetaSecond = that->logPBetaMask( betaXO[1] , gammaMask_XO[1] , that->getW() , that->getW0() );
    
    double logPGammaFirst = this->logPGamma( gammaXO[0] );
    double logPGammaSecond = that->logPGamma( gammaXO[1] );
    
    pCrossOver +=   ( logLikFirst + logPBetaFirst + logPGammaFirst +
                     logLikSecond + logPBetaSecond + logPGammaSecond ) -
    ( this->getLogLikelihood() + this->getLogPBeta() + this->getLogPGamma() +
     that->getLogLikelihood() + that->getLogPBeta() + that->getLogPGamma() );
    
    if( randLogU01() < pCrossOver )
    {
        // -- first chain
        
        this->setGamma( gammaXO[0] , logPGammaFirst );
        this->setBeta( betaXO[0] , logPBetaFirst );
        
        this->setGammaMask( gammaMask_XO[0] );
        this->setXB( XB_XO[0] );
        this->setU( U_XO[0] );
        this->setRhoU( rhoU_XO[0] );
        
        this->setLogLikelihood( logLikFirst );
        
        // -- second chain
        
        that->setGamma( gammaXO[1] , logPGammaSecond );
        that->setBeta( betaXO[1] , logPBetaSecond );
        
        that->setGammaMask( gammaMask_XO[1] );
        that->setXB( XB_XO[1] );
        that->setU( U_XO[1] );
        that->setRhoU( rhoU_XO[1] );
        
        that->setLogLikelihood( logLikSecond );
        
        return 1;
    }else
        return 0;
    
    
}

int SUR_Chain::uniform_crossOver_step( std::shared_ptr<SUR_Chain>& that )
{
    double pCrossOver;
    
    std::vector<arma::umat> gammaXO(2); gammaXO[0] = arma::umat(nVSPredictors,nOutcomes);  gammaXO[1] = arma::umat(nVSPredictors,nOutcomes);
    
    // Propose Crossover
    for(unsigned int j=0; j<nVSPredictors; ++j)
    {
        for(unsigned int k=0; k<nOutcomes; ++k)
        {
            if( randU01() < 0.5 )
            {
                gammaXO[0](j,k) = this->getGamma()(j,k);
                gammaXO[1](j,k) = that->getGamma()(j,k);
                
            }else{
                gammaXO[0](j,k) = that->getGamma()(j,k);
                gammaXO[1](j,k) = this->getGamma()(j,k);
                
            }
        }
    }
    
    pCrossOver = 0; // XO prop probability symmetric now
    
    // Propose betas that go with the new crossed-over states
    std::vector<arma::mat> betaXO(2);
    betaXO[0] = this->getBeta();
    betaXO[1] = that->getBeta();
    
    std::vector<arma::umat> gammaMask_XO(2);
    gammaMask_XO[0] = createGammaMask(gammaXO[0]);
    gammaMask_XO[1] = createGammaMask(gammaXO[1]);
    
    std::vector<arma::mat> XB_XO(2);
    XB_XO[0] = arma::mat( this->getXB() ); // copy into the new one
    XB_XO[1] = arma::mat( that->getXB() );
    std::vector<arma::mat> U_XO(2);
    U_XO[0] = arma::mat( this->getU() ); // copy into the new one
    U_XO[1] = arma::mat( that->getU() );
    std::vector<arma::mat> rhoU_XO(2);
    rhoU_XO[0] = arma::mat( this->getRhoU() ); // copy into the new one
    rhoU_XO[1] = arma::mat( that->getRhoU() );
    
    // propose for the first new chain
    pCrossOver -= this->sampleBetaGivenSigmaRho( betaXO[0] , this->getSigmaRho() , this->getJT() ,
                                                gammaMask_XO[0] , XB_XO[0] , U_XO[0] , rhoU_XO[0] );
    
    pCrossOver += this->logPBetaGivenSigmaRho( this->getBeta() , this->getSigmaRho() , this->getJT() ,
                                              this->getGammaMask() , XB_XO[0] , U_XO[0] , rhoU_XO[0] );
    
    // propose for the second new chain
    pCrossOver -= that->sampleBetaGivenSigmaRho( betaXO[1] , that->getSigmaRho() , that->getJT() ,
                                                gammaMask_XO[1] , XB_XO[1] , U_XO[1] , rhoU_XO[1] );
    
    pCrossOver += that->logPBetaGivenSigmaRho( that->getBeta() , that->getSigmaRho() , that->getJT() ,
                                              that->getGammaMask() , XB_XO[1] , U_XO[1] , rhoU_XO[1] );
    
    // log Posterior of the new chains
    double logLikFirst = this->logLikelihood( gammaMask_XO[0] , XB_XO[0] , U_XO[0] , rhoU_XO[0] , this->getSigmaRho() );
    double logLikSecond = that->logLikelihood( gammaMask_XO[1] , XB_XO[1] , U_XO[1] , rhoU_XO[1] , that->getSigmaRho() );
    
    double logPBetaFirst = this->logPBetaMask( betaXO[0] , gammaMask_XO[0] , this->getW() , this->getW0() );  // this and that are important for temperature and hyperparameters
    double logPBetaSecond = that->logPBetaMask( betaXO[1] , gammaMask_XO[1] , that->getW() , that->getW0() );
    
    double logPGammaFirst = this->logPGamma( gammaXO[0] );
    double logPGammaSecond = that->logPGamma( gammaXO[1] );
    
    pCrossOver +=   ( logLikFirst + logPBetaFirst + logPGammaFirst +
                     logLikSecond + logPBetaSecond + logPGammaSecond ) -
    ( this->getLogLikelihood() + this->getLogPBeta() + this->getLogPGamma() +
     that->getLogLikelihood() + that->getLogPBeta() + that->getLogPGamma() );
    
    if( randLogU01() < pCrossOver )
    {
        // -- first chain
        
        this->setGamma( gammaXO[0] , logPGammaFirst );
        this->setBeta( betaXO[0] , logPBetaFirst );
        
        this->setGammaMask( gammaMask_XO[0] );
        this->setXB( XB_XO[0] );
        this->setU( U_XO[0] );
        this->setRhoU( rhoU_XO[0] );
        
        this->setLogLikelihood( logLikFirst );
        
        // -- second chain
        
        that->setGamma( gammaXO[1] , logPGammaSecond );
        that->setBeta( betaXO[1] , logPBetaSecond );
        
        that->setGammaMask( gammaMask_XO[1] );
        that->setXB( XB_XO[1] );
        that->setU( U_XO[1] );
        that->setRhoU( rhoU_XO[1] );
        
        that->setLogLikelihood( logLikSecond );
        
        return 1;
    }else
        return 0;
    
    
}

int SUR_Chain::block_crossOver_step( std::shared_ptr<SUR_Chain>& that , arma::mat& corrMatX , double threshold )
{
    double pCrossOver;
    
    std::vector<arma::umat> gammaXO(2); gammaXO[0] = arma::umat(nVSPredictors,nOutcomes);  gammaXO[1] = arma::umat(nVSPredictors,nOutcomes);
    
    // Propose Crossover
    
    // Select the ONE index to foor the block
    unsigned int predIdx = randIntUniform(0, nVSPredictors-1 ); // pred
    unsigned int outcIdx = randIntUniform(0, nOutcomes-1 ); // outcome
    
    arma::uvec covIdx;
    
    if( preComputedXtX ){
        covIdx = arma::find( arma::abs( corrMatX.row(predIdx) ) > threshold );  // this will include the original predIdx
    }else{
        // if not precomputed, compute the correlation only for the current row
        
        // #ifdef _OPENMP
        // #pragma omp parallel for
        // #endif
        // for( unsigned int j=0; j<nVSPredictors-1; ++j){
        //     tmpVec(j) = arma::as_scalar( arma::cor( data->col( (*VSPredictorsIdx)(predIdx) ) , data->col( (*VSPredictorsIdx(j)) ) );
        // }
        // covIdx = arma::find( arma::abs( tmpVec ) > threshold );
        // I'd rather leave parallelisation to armadillo and avoid the temp, but this version would work as well
        
        covIdx = arma::find( arma::abs( arma::cor( data->col( (*VSPredictorsIdx)(predIdx) ) , data->cols( (*VSPredictorsIdx) ) ) ) > threshold );
    }
    
    gammaXO[0] = this->getGamma();
    gammaXO[1] = that->getGamma();
    
    for(unsigned int j=0; j<covIdx.n_elem; ++j)
    {
        gammaXO[0](covIdx(j),outcIdx) = that->getGamma()(covIdx(j),outcIdx);
        gammaXO[1](covIdx(j),outcIdx) = this->getGamma()(covIdx(j),outcIdx);
    }
    
    pCrossOver = 0.;  // XO prop probability is weird, how do I compute it? Let's say is symmetric as is determnistic and both comes from the same corrMatX
    
    // Propose betas that go with the new crossed-over states
    std::vector<arma::mat> betaXO(2);
    betaXO[0] = this->getBeta();
    betaXO[1] = that->getBeta();
    
    std::vector<arma::umat> gammaMask_XO(2);
    gammaMask_XO[0] = createGammaMask(gammaXO[0]);
    gammaMask_XO[1] = createGammaMask(gammaXO[1]);
    
    std::vector<arma::mat> XB_XO(2);
    XB_XO[0] = arma::mat( this->getXB() ); // copy into the new one
    XB_XO[1] = arma::mat( that->getXB() );
    std::vector<arma::mat> U_XO(2);
    U_XO[0] = arma::mat( this->getU() ); // copy into the new one
    U_XO[1] = arma::mat( that->getU() );
    std::vector<arma::mat> rhoU_XO(2);
    rhoU_XO[0] = arma::mat( this->getRhoU() ); // copy into the new one
    rhoU_XO[1] = arma::mat( that->getRhoU() );
    
    // propose for the first new chain
    pCrossOver -= this->sampleBetaGivenSigmaRho( betaXO[0] , this->getSigmaRho() , this->getJT() ,
                                                gammaMask_XO[0] , XB_XO[0] , U_XO[0] , rhoU_XO[0] );
    
    pCrossOver += this->logPBetaGivenSigmaRho( this->getBeta() , this->getSigmaRho() , this->getJT() ,
                                              this->getGammaMask() , XB_XO[0] , U_XO[0] , rhoU_XO[0] );
    
    // propose for the second new chain
    pCrossOver -= that->sampleBetaGivenSigmaRho( betaXO[1] , that->getSigmaRho() , that->getJT() ,
                                                gammaMask_XO[1] , XB_XO[1] , U_XO[1] , rhoU_XO[1] );
    
    pCrossOver += that->logPBetaGivenSigmaRho( that->getBeta() , that->getSigmaRho() , that->getJT() ,
                                              that->getGammaMask() , XB_XO[1] , U_XO[1] , rhoU_XO[1] );
    
    // log Posterior of the new chains
    double logLikFirst = this->logLikelihood( gammaMask_XO[0] , XB_XO[0] , U_XO[0] , rhoU_XO[0] , this->getSigmaRho() );
    double logLikSecond = that->logLikelihood( gammaMask_XO[1] , XB_XO[1] , U_XO[1] , rhoU_XO[1] , that->getSigmaRho() );
    
    double logPBetaFirst = this->logPBetaMask( betaXO[0] , gammaMask_XO[0] , this->getW() , this->getW0() );  // this and that are important for temperature and hyperparameters
    double logPBetaSecond = that->logPBetaMask( betaXO[1] , gammaMask_XO[1] , that->getW() , that->getW0() );
    
    double logPGammaFirst = this->logPGamma( gammaXO[0] );
    double logPGammaSecond = that->logPGamma( gammaXO[1] );
    
    pCrossOver +=   ( logLikFirst + logPBetaFirst + logPGammaFirst +
                     logLikSecond + logPBetaSecond + logPGammaSecond ) -
    ( this->getLogLikelihood() + this->getLogPBeta() + this->getLogPGamma() +
     that->getLogLikelihood() + that->getLogPBeta() + that->getLogPGamma() );
    
    if( randLogU01() < pCrossOver )
    {
        // -- first chain
        
        this->setGamma( gammaXO[0] , logPGammaFirst );
        this->setBeta( betaXO[0] , logPBetaFirst );
        
        this->setGammaMask( gammaMask_XO[0] );
        this->setXB( XB_XO[0] );
        this->setU( U_XO[0] );
        this->setRhoU( rhoU_XO[0] );
        
        this->setLogLikelihood( logLikFirst );
        
        // -- second chain
        
        that->setGamma( gammaXO[1] , logPGammaSecond );
        that->setBeta( betaXO[1] , logPBetaSecond );
        
        that->setGammaMask( gammaMask_XO[1] );
        that->setXB( XB_XO[1] );
        that->setU( U_XO[1] );
        that->setRhoU( rhoU_XO[1] );
        
        that->setLogLikelihood( logLikSecond );
        
        return 1;
    }else
        return 0;
    
}

void SUR_Chain::swapAll( std::shared_ptr<SUR_Chain>& thatChain )
{
    
    // HARD SWAP cause swapping "this" is not an option
    // swap quantities
    arma::umat swapGammaMask;
    arma::mat swapMat;
    
    swapGammaMask = this->getGammaMask() ;
    this->setGammaMask( thatChain->getGammaMask() );
    thatChain->setGammaMask( swapGammaMask );
    
    swapMat = this->getXB();
    this->setXB( thatChain->getXB() );
    thatChain->setXB( swapMat );
    
    swapMat = this->getU();
    this->setU( thatChain->getU() );
    thatChain->setU( swapMat );
    
    swapMat = this->getRhoU();
    this->setRhoU( thatChain->getRhoU() );
    thatChain->setRhoU( swapMat );
    
    // parameters and priors
    this->swapTau( thatChain );
    
    if ( covariance_type == Covariance_Type::HIW )
    {
        this->swapEta( thatChain );
        this->swapJT( thatChain );
    }
    
    this->swapSigmaRho( thatChain );
    
    if ( gamma_type == Gamma_Type::hotspot )
    {
        this->swapO( thatChain );
        this->swapPi( thatChain );
    }
    else if ( gamma_type == Gamma_Type::hierarchical )
    {
        this->swapPi( thatChain );
    }
    
    this->swapGamma( thatChain );
    
    this->swapW( thatChain );
    this->swapW0( thatChain );
    this->swapBeta( thatChain );
    
    // recompute likelihood
    this->logLikelihood();
    thatChain->logLikelihood();
    
}

int SUR_Chain::globalStep( std::shared_ptr<SUR_Chain>& that )
{
    
    unsigned int globalType {6};  // the default skips the global step
    switch ( covariance_type )
    {
        case Covariance_Type::HIW :
            globalType = randIntUniform(0,5);
            break;
            
        case Covariance_Type::IW :
            globalType = randIntUniform(0,4);
            break;
            
        default:
            break;
    }
    
    switch(globalType){
            
            // -- Exchange and CrossOver
        case 0:
            return this -> exchangeGamma_step( that );
            break;
            
        case 1:
            return this -> adapt_crossOver_step( that );
            break;
            
        case 2:
            return this -> uniform_crossOver_step( that );
            break;
            
        case 3:
            return this -> block_crossOver_step( that , corrMatX , 0.25 );
            break;
            
        case 4:
            return this -> exchangeAll_step( that );
            break;
            
        case 5:
            return this -> exchangeJT_step( that );
            break;
            
        default:
            break;
    }
    
    return 0;
}

int SUR_Chain::exchangeAll_step( std::shared_ptr<SUR_Chain>& thatChain )
{
    
    double logPExchange = ( this->getLogLikelihood() * this->getTemperature() -
                           thatChain->getLogLikelihood() * thatChain->getTemperature() ) *
    ( 1. / thatChain->getTemperature() - 1. / this->getTemperature() );
    //  no priors because that is not tempered so it cancels out
    
    if( randLogU01() < logPExchange )
    {
        // Swap all the states
        this -> swapAll( thatChain );
        
        return 1;
    }else
        return 0;
}

// end NB

// *******************************
// Other Methods
// *******************************

// update relavant quantities
arma::umat SUR_Chain::createGammaMask( const arma::umat& gamma )
{
    
    // CREATE HERE THE GAMMA "MASK"
    // INITIALISE THE INDEXES FOR THE GAMMA MASK
    arma::umat mask = arma::zeros<arma::umat>(nFixedPredictors*nOutcomes,2); //this is just an initialisation
    if( nFixedPredictors )
    {
      for( unsigned int j=0; j<nFixedPredictors; ++j)
      {
        for(unsigned int k=0 ; k<nOutcomes ; ++k)  //add gammas for the fixed variables
        {
          mask(j*nOutcomes+k,0) = j; mask(j*nOutcomes+k,1) = k;
        }
      }
    }
    
    for(unsigned int k=0 ; k<nOutcomes ; ++k)  //add the other gammas
    {
        arma::uvec tmpUVec = arma::find(gamma.col(k) != 0);
        unsigned int tmpIdx = mask.n_rows;
        
        if( tmpUVec.n_elem > 0 )
        {
            mask.insert_rows( tmpIdx , arma::zeros<arma::umat>( tmpUVec.n_elem , 2 ));
            mask.submat( tmpIdx, 0, mask.n_rows-1 , 0 ) = tmpUVec + nFixedPredictors ;
            mask.submat( tmpIdx, 1, mask.n_rows-1 , 1 ).fill(k);
        }
    }
    // Gamma mask done
    
    return mask;
}


void SUR_Chain::updateGammaMask()
{
    // CREATE HERE THE GAMMA "MASK"
    // INITIALISE THE INDEXES FOR THE GAMMA MASK
    gammaMask.zeros(nFixedPredictors*nOutcomes,2); //this is just an initialisation
    if( nFixedPredictors )
    {
      for( unsigned int j=0; j<nFixedPredictors; ++j)
      {
        for(unsigned int k=0 ; k<nOutcomes ; ++k)  //add gammas for the fixed variables
        {
          gammaMask(j*nOutcomes+k,0) = j; gammaMask(j*nOutcomes+k,1) = k;
        }
      }
    }
    
    for(unsigned int k=0 ; k<nOutcomes ; ++k)   //add the other gammas
    {
        arma::uvec tmpUVec = arma::find(gamma.col(k) != 0);
        unsigned int tmpIdx = gammaMask.n_rows;
        
        if( tmpUVec.n_elem > 0 )
        {
            gammaMask.insert_rows( tmpIdx , arma::zeros<arma::umat>( tmpUVec.n_elem , 2 ));
            gammaMask.submat( tmpIdx, 0, gammaMask.n_rows-1 , 0 ) = tmpUVec + nFixedPredictors ;
            gammaMask.submat( tmpIdx, 1, gammaMask.n_rows-1 , 1 ).fill(k);
        }
    }
    // Gamma mask updated
}

arma::mat SUR_Chain::createXB( const arma::umat&  externalGammaMask , const arma::mat&  externalBeta )
{
    arma::uvec singleIdx_k(1), VS_IN_k;
    arma::mat externalXB = arma::zeros<arma::mat>(nObservations,nOutcomes);
    
    if(externalGammaMask.n_rows > 0)
    {
        for(unsigned int k=0; k<nOutcomes; ++k)
        {
            singleIdx_k(0) = k;
            VS_IN_k =  externalGammaMask( arma::find(  externalGammaMask.col(1) == k ) , arma::zeros<arma::uvec>(1) );
            externalXB.col(k) = (data->cols( (*predictorsIdx)(VS_IN_k) ) *  externalBeta.submat(VS_IN_k,singleIdx_k) );
        }
    }
    return externalXB;
}

void SUR_Chain::updateXB()
{
    arma::uvec singleIdx_k(1), VS_IN_k;
    XB.set_size(nObservations,nOutcomes); // reset without initialising nor preserving data
    XB.fill(0.);
    
    if(gammaMask.n_rows > 0)
    {
        for(unsigned int k=0; k<nOutcomes; ++k)
        {
            singleIdx_k(0) = k;
            VS_IN_k = gammaMask( arma::find( gammaMask.col(1) == k ) , arma::zeros<arma::uvec>(1) );
            XB.col(k) = (data->cols( (*predictorsIdx)(VS_IN_k) ) * beta.submat(VS_IN_k,singleIdx_k) );
        }
    }
}

arma::mat SUR_Chain::createU( const arma::mat& externalXB )
{
    arma::mat externalU = data->cols( *outcomesIdx ) - externalXB;
    return externalU;
}

void SUR_Chain::updateU()
{
    U = data->cols( *outcomesIdx ) - XB;
}

arma::mat SUR_Chain::createRhoU( const arma::mat& externalU , const arma::mat&  externalSigmaRho , const JunctionTree& externalJT )
{
    
    arma::mat externalRhoU = arma::zeros<arma::mat>(nObservations,nOutcomes);
    
    switch ( covariance_type )
    {
        case Covariance_Type::HIW :
        {
            arma::uvec xi = arma::conv_to<arma::uvec>::from(externalJT.perfectEliminationOrder);
            
            for( unsigned int k=1; k < nOutcomes; ++k)
            {
                for(unsigned int l=0 ; l<k ; ++l)
                {
                    if(  externalSigmaRho(xi(k),xi(l)) != 0 )
                        externalRhoU.col(xi(k)) += externalU.col(xi(l)) *  externalSigmaRho(xi(k),xi(l));
                }
            }
            break;
        }
            
        case Covariance_Type::IW :
        {
            for( unsigned int k=1; k < nOutcomes; ++k)
            {
                for(unsigned int l=0 ; l<k ; ++l)
                {
                    if(  externalSigmaRho(k,l) != 0 )
                        externalRhoU.col(k) += externalU.col(l) *  externalSigmaRho(k,l);
                }
            }
            break;
        }
            
        default:
            throw Bad_Covariance_Type ( covariance_type );
    }
    
    return externalRhoU;
}


void SUR_Chain::updateRhoU()
{
    rhoU.zeros(nObservations,nOutcomes);
    
    switch ( covariance_type )
    {
        case Covariance_Type::HIW :
        {
            arma::uvec xi = arma::conv_to<arma::uvec>::from(jt.perfectEliminationOrder);
            for( unsigned int k=1; k < nOutcomes; ++k)
            {
                for(unsigned int l=0 ; l<k ; ++l)
                {
                    if(  sigmaRho(xi(k),xi(l)) != 0 )
                        rhoU.col(xi(k)) += U.col(xi(l)) * sigmaRho(xi(k),xi(l));
                }
            }
            break;
        }
            
        case Covariance_Type::IW :
        {
            for( unsigned int k=1; k < nOutcomes; ++k)
            {
                for(unsigned int l=0 ; l<k ; ++l)
                {
                    if(  sigmaRho(k,l) != 0 )
                        rhoU.col(k) += U.col(l) * sigmaRho(k,l);
                }
            }
            break;
        }
            
        default:
            break;
    }
    
}

void SUR_Chain::createQuantities( arma::umat&  externalGammaMask , arma::mat& externalXB , arma::mat& externalU , arma::mat& externalRhoU ,
                                 const arma::umat& externalGamma , const arma::mat&  externalBeta ,
                                 const arma::mat&  externalSigmaRho , const JunctionTree& externalJT )
{
    externalGammaMask = createGammaMask( externalGamma );
    externalXB = createXB(  externalGammaMask ,  externalBeta );
    externalU = createU( externalXB );
    externalRhoU = createRhoU( externalU ,  externalSigmaRho , externalJT );
}

void SUR_Chain::updateQuantities()
{
    updateGammaMask();
    updateXB();
    updateU();
    updateRhoU();
}



// Bandit-sampling related methods
void SUR_Chain::banditInit()// initialise all the private memebers
{
    banditZeta = arma::vec(nVSPredictors);
    
    banditAlpha = arma::mat(nVSPredictors,nOutcomes);
    banditAlpha.fill( 0.5 );
    
    banditBeta = arma::mat(nVSPredictors,nOutcomes);
    banditBeta.fill( 0.5 );
    
    mismatch = arma::vec(nVSPredictors);
    normalised_mismatch = arma::vec(nVSPredictors);
    normalised_mismatch_backwards = arma::vec(nVSPredictors);
    
    n_updates_bandit = 4; // this needs to be low as its O(n_updates!)
    
    banditLimit = (double)nObservations;
    banditIncrement = 1.;
}

// MC3 init
void SUR_Chain::MC3Init()
{
    n_updates_MC3 = std::ceil( nVSPredictors/40 ); //arbitrary number, should I use something different?
}
