#include "HRR_Chain.h"

/*******************************
 * HRR HAS NO PARALLEL OMP FOR NOW
 * even if the likelihood would seem embarassingly parallel
 * the are problems in adding any omp pragma...
 * Note that contrarily from SUR here the model still have sigma in the prior for beta for computaional convenience
 *******************************/


// *******************************
// Constructors
// *******************************

HRR_Chain::HRR_Chain( std::shared_ptr<arma::mat> data_, std::shared_ptr<arma::mat> mrfG_, unsigned int nObservations_,
                     unsigned int nOutcomes_, unsigned int nVSPredictors_, unsigned int nFixedPredictors_,
                     std::shared_ptr<arma::uvec> outcomesIdx_, std::shared_ptr<arma::uvec> VSPredictorsIdx_,
                     std::shared_ptr<arma::uvec> fixedPredictorsIdx_, std::shared_ptr<arma::umat> missingDataArrayIdx_, std::shared_ptr<arma::uvec> completeCases_,
                     Gamma_Sampler_Type gamma_sampler_type_ , Gamma_Type gamma_type_ ,
                     Beta_Type beta_type_ , Covariance_Type covariance_type_ , bool output_CPO , int maxThreads , int tick ,
                     double externalTemperature ):
data(data_), mrfG(mrfG_), outcomesIdx(outcomesIdx_), VSPredictorsIdx(VSPredictorsIdx_), fixedPredictorsIdx(fixedPredictorsIdx_),
missingDataArrayIdx(missingDataArrayIdx_), completeCases(completeCases_),
nObservations(nObservations_), nOutcomes(nOutcomes_), nVSPredictors(nVSPredictors_), nFixedPredictors(nFixedPredictors_),
temperature(externalTemperature),internalIterationCounter(0),
covariance_type(covariance_type_),gamma_type(gamma_type_),beta_type(beta_type_),gamma_sampler_type(gamma_sampler_type_)
{
    
    if( covariance_type != Covariance_Type::IG )
        throw Bad_Covariance_Type ( covariance_type );
    
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
    
    wInit();
    w0Init();
    
    sigmaABInit();
    
    logLikelihood();
    if( output_CPO && (temperature = 1.) ){
        predLikelihood();
    }
    
}


HRR_Chain::HRR_Chain( Utils::SUR_Data& surData,
                     Gamma_Sampler_Type gamma_sampler_type_ , Gamma_Type gamma_type_ ,
                     Beta_Type beta_type_ , Covariance_Type covariance_type_ , bool output_CPO , int maxThreads , int tick ,
                     double externalTemperature ):
HRR_Chain(surData.data,surData.mrfG,surData.nObservations,surData.nOutcomes,surData.nVSPredictors,surData.nFixedPredictors,
surData.outcomesIdx,surData.VSPredictorsIdx,surData.fixedPredictorsIdx,surData.missingDataArrayIdx,surData.completeCases,
          gamma_sampler_type_,gamma_type_,beta_type_,covariance_type_,output_CPO,maxThreads,tick,externalTemperature){ }

HRR_Chain::HRR_Chain( Utils::SUR_Data& surData, double externalTemperature ):
HRR_Chain(surData.data,surData.mrfG,surData.nObservations,surData.nOutcomes,surData.nVSPredictors,surData.nFixedPredictors,
surData.outcomesIdx,surData.VSPredictorsIdx,surData.fixedPredictorsIdx,surData.missingDataArrayIdx,surData.completeCases,
          Gamma_Sampler_Type::bandit , Gamma_Type::hotspot , Beta_Type::independent , Covariance_Type::IG , false ,
          externalTemperature){ }

// *******************************
// Getters and Setters
// *******************************

// data
void HRR_Chain::setXtX()
{
    
    // Compute XtX
    if( (nFixedPredictors+nVSPredictors) < 100000 )  // kinda arbitrary value, how can we assess a more sensible one?
    {
        preComputedXtX = true;
        XtX = data->cols( *predictorsIdx ).t() * data->cols( *predictorsIdx );
        corrMatX = arma::cor( data->submat(arma::regspace<arma::uvec>(0,nObservations-1), *VSPredictorsIdx ) );  // this is only for values to be selected
    }else{
        
        preComputedXtX = false;
        XtX.clear();          // if not precomputed, these two are just reset
        corrMatX.clear();
    }
}

// gPrior
void HRR_Chain::gPriorInit() // g Prior can only be init at the start, so no proper set method
{
    if( internalIterationCounter > 0 )
        throw std::runtime_error(std::string("gPrior can only be initialised at the start of the MCMC"));
    
    // set the boot to true
    beta_type = Beta_Type::gprior;
    
    // re-initialise the w parameter in line with the new prior, as w now has a different meaning
    wInit( (double)nObservations , 0.5*nOutcomes + nOutcomes -1. , 0.5*nObservations*nOutcomes ); // these values are taken from Lewin 2016
    
    // update internals
    logPW();
    log_likelihood = logLikelihood();
    
}

// usefull quantities to keep track of
arma::umat& HRR_Chain::getGammaMask(){ return gammaMask; }
void HRR_Chain::setGammaMask( arma::umat  externalGammaMask )
{
    gammaMask =  externalGammaMask ;
}

arma::urowvec& HRR_Chain::getModelSize() const
{
    static arma::urowvec modelSize;
    modelSize = nFixedPredictors + arma::sum( gamma , 0 ); // 0 is to get the sum of the elements in each column
    return modelSize;
}

// MCMC related tuning parameters
double HRR_Chain::getTemperature() const{ return temperature; }
void HRR_Chain::setTemperature( double temp_ )
{
    log_likelihood = log_likelihood * temperature / temp_ ; // re-temper the log_likelihood first
    temperature = temp_ ;
}

unsigned int HRR_Chain::getinternalIterationCounter() const{ return internalIterationCounter ; }

Gamma_Sampler_Type HRR_Chain::getGammaSamplerType(){ return gamma_sampler_type ; }
void HRR_Chain::setGammaSamplerType( Gamma_Sampler_Type gamma_sampler_type_ )
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
unsigned int HRR_Chain::getNUpdatesBandit() const{ return n_updates_bandit ; }
void HRR_Chain::setNUpdatesBandit( unsigned int n_updates_bandit_ ){ n_updates_bandit = n_updates_bandit_ ; }

arma::mat& HRR_Chain::getBanditZeta(){ return banditZeta; }
void HRR_Chain::setBanditZeta( arma::mat banditZeta_ ){ banditZeta = banditZeta_ ; }

arma::mat& HRR_Chain::getBanditAlpha(){ return banditAlpha ; }
void HRR_Chain::setBanditAlpha( arma::mat banditAlpha_ ){ banditAlpha = banditAlpha_ ; }

arma::mat& HRR_Chain::getBanditBeta(){ return banditBeta ; }
void HRR_Chain::setBanditBeta( arma::mat banditBeta_ ){ banditBeta = banditBeta_ ; }

arma::vec& HRR_Chain::getBanditMismatch(){ return mismatch; }
void HRR_Chain::setBanditMismatch( arma::vec mismatch_ ){ mismatch = mismatch_ ; }

arma::vec& HRR_Chain::getBanditNormalisedMismatch(){ return normalised_mismatch ; }
void HRR_Chain::setBanditNormalisedMismatch( arma::vec normalised_mismatch_ ){ normalised_mismatch = normalised_mismatch_ ; }

arma::vec& HRR_Chain::getBanditNormalisedMismatchBackwards(){ return normalised_mismatch_backwards ; }
void HRR_Chain::setBanditNormalisedMismatchBackwards( arma::vec normalised_mismatch_backwards_ ){ normalised_mismatch_backwards = normalised_mismatch_backwards_ ; }

// Parameter states etc

void HRR_Chain::sigmaABInit()
{
    a_sigma = 1.;     // their mean is gonna be b/(a-1)
    b_sigma = 1.;     // with variance b^2/(a-1)^2/(a-2)
}

double HRR_Chain::getSigmaA() const{ return a_sigma ; }
void HRR_Chain::setSigmaA( double a_sigma_ )
{
    a_sigma = a_sigma_ ;
    logLikelihood();
}

double HRR_Chain::getSigmaB() const{ return b_sigma ; }
void HRR_Chain::setSigmaB( double b_sigma_ )
{
    b_sigma = b_sigma_ ;
    logLikelihood();
}

void HRR_Chain::setSigmaAB( double a_sigma_ , double b_sigma_ )
{
    if ( covariance_type != Covariance_Type::IG )
        throw Bad_Covariance_Type( covariance_type );
    
    a_sigma = a_sigma_ ;
    b_sigma = b_sigma_ ;
    logLikelihood();
}


// o_k
arma::vec& HRR_Chain::getO(){ return o ; }
void HRR_Chain::setO( arma::vec& o_ )
{
    o = o_ ;
    logPO();
}

void HRR_Chain::setO( arma::vec& o_ , double logP_o_ )
{
    o = o_ ;
    logP_o = logP_o_ ;
}

double HRR_Chain::getOA() const{ return a_o ; }
void HRR_Chain::setOA( double a_o_ )
{
    a_o = a_o_ ;
    logPO();
}

double HRR_Chain::getOB() const{ return b_o ; }
void HRR_Chain::setOB( double b_o_ )
{
    b_o = b_o_ ;
    logPO();
}

void HRR_Chain::setOAB( double a_o_ , double b_o_ )
{
    if ( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type( gamma_type );
    
    a_o = a_o_ ;
    b_o = b_o_ ;
    logPO();
}

double HRR_Chain::getVarOProposal() const{ return var_o_proposal ; }
void HRR_Chain::setVarOProposal( double var_o_proposal_ ){ var_o_proposal = var_o_proposal_ ; }

double HRR_Chain::getOAccRate() const{ return o_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double HRR_Chain::getLogPO() const{ return logP_o ; }
// no setter for this, dedicated setter below

// pi_j
arma::vec& HRR_Chain::getPi(){ return pi ; }
void HRR_Chain::setPi( arma::vec& pi_ )
{
    pi = pi_ ;
    logPPi();
}

void HRR_Chain::setPi( arma::vec& pi_ , double logP_pi_ )
{
    pi = pi_ ;
    logP_pi = logP_pi_ ;
}

double HRR_Chain::getPiA() const{ return a_pi ; }
void HRR_Chain::setPiA( double a_pi_ )
{
    a_pi = a_pi_ ;
    logPPi();
}

double HRR_Chain::getPiB() const{ return b_pi ; }
void HRR_Chain::setPiB( double b_pi_ )
{
    b_pi = b_pi_ ;
    logPPi();
}

void HRR_Chain::setPiAB( double a_pi_ , double b_pi_ )
{
    if ( gamma_type != Gamma_Type::hotspot && gamma_type != Gamma_Type::hierarchical )
        throw Bad_Gamma_Type( gamma_type );
    
    a_pi = a_pi_ ;
    b_pi = b_pi_ ;
    logPPi();
}

double HRR_Chain::getVarPiProposal() const{ return var_pi_proposal ; }
void HRR_Chain::setVarPiProposal( double var_pi_proposal_ ){ var_pi_proposal = var_pi_proposal_ ; }

double HRR_Chain::getPiAccRate() const{ return pi_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double HRR_Chain::getLogPPi() const{ return logP_pi ; }
// no setter for this, dedicated setter below

// GAMMA
arma::umat& HRR_Chain::getGamma(){ return gamma ; }
void HRR_Chain::setGamma( arma::umat& externalGamma )
{
    gamma = externalGamma ;
    logPGamma();
    log_likelihood = logLikelihood( gammaMask , gamma ); // update internal state
}

void HRR_Chain::setGamma( arma::umat& externalGamma , double logP_gamma_ )
{
    gamma = externalGamma ;
    logP_gamma = logP_gamma_ ;
    log_likelihood = logLikelihood( gammaMask , gamma ); // update internal state
}

double HRR_Chain::getGammaD() const{ return mrf_d ; }
void HRR_Chain::setGammaD( double mrf_d_ )
{
    if( gamma_type != Gamma_Type::mrf )
        throw Bad_Gamma_Type( gamma_type );
    
    mrf_d = mrf_d_ ;
    logPGamma();
}

double HRR_Chain::getGammaE() const{ return mrf_e ; }
void HRR_Chain::setGammaE( double mrf_e_ )
{
    if( gamma_type != Gamma_Type::mrf )
        throw Bad_Gamma_Type( gamma_type );
    
    mrf_e = mrf_e_ ;
    logPGamma();
}

void HRR_Chain::setGammaDE( double mrf_d_ , double mrf_e_ )
{
    if( gamma_type != Gamma_Type::mrf )
        throw Bad_Gamma_Type( gamma_type );
    
    mrf_d = mrf_d_ ;
    mrf_e = mrf_e_ ;
    logPGamma();
}

unsigned int HRR_Chain::getNUpdatesMC3() const{ return n_updates_MC3 ; }
void HRR_Chain::setNUpdatesMC3( unsigned int n_updates_MC3_ ){ n_updates_MC3 = n_updates_MC3_ ; }

double HRR_Chain::getGammaAccRate() const{ return gamma_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double HRR_Chain::getLogPGamma() const{ return logP_gamma ; }
// no setter for this, dedicated setter below

// W
double HRR_Chain::getW() const{ return w ; }
void HRR_Chain::setW( double w_ )
{
    w = w_ ;
    logPW();
    logLikelihood();
}

void HRR_Chain::setW( double w_ , double logP_w_ )
{
    w = w_ ;
    logP_w = logP_w_ ;
    logLikelihood();
}

double HRR_Chain::getWA() const{ return a_w ; }
void HRR_Chain::setWA( double a_w_ )
{
    a_w = a_w_ ;
    logPW();
}

double HRR_Chain::getWB() const{ return b_w ; }
void HRR_Chain::setWB( double b_w_ )
{
    b_w = b_w_ ;
    logPW();
}

void HRR_Chain::setWAB( double a_w_ , double b_w_ )
{
    a_w = a_w_ ;
    b_w = b_w_ ;
    logPW();
}

double HRR_Chain::getVarWProposal() const{ return var_w_proposal ; }
void HRR_Chain::setVarWProposal( double var_w_proposal_ ){ var_w_proposal = var_w_proposal_ ; }

double HRR_Chain::getWAccRate() const{ return w_acc_count/(double)internalIterationCounter ; }

double HRR_Chain::getLogPW() const{ return logP_w; }
// no setter for this, dedicated setter below

// W0
double HRR_Chain::getW0() const{ return w0 ; }
void HRR_Chain::setW0( double w0_ )
{
  w0 = w0_ ;
  logPW0();
  logLikelihood();
}

void HRR_Chain::setW0( double w0_ , double logP_w0_ )
{
  w0 = w0_ ;
  logP_w0 = logP_w0_ ;
  logLikelihood();
}

double HRR_Chain::getW0A() const{ return a_w0 ; }
void HRR_Chain::setW0A( double a_w0_ )
{
  a_w0 = a_w0_ ;
  logPW0();
}

double HRR_Chain::getW0B() const{ return b_w0 ; }
void HRR_Chain::setW0B( double b_w0_ )
{
  b_w0 = b_w0_ ;
  logPW0();
}

void HRR_Chain::setW0AB( double a_w0_ , double b_w0_ )
{
  a_w0 = a_w0_ ;
  b_w0 = b_w0_ ;
  logPW0();
}

double HRR_Chain::getVarW0Proposal() const{ return var_w0_proposal ; }
void HRR_Chain::setVarW0Proposal( double var_w0_proposal_ ){ var_w0_proposal = var_w0_proposal_ ; }

double HRR_Chain::getW0AccRate() const{ return w0_acc_count/(double)internalIterationCounter ; }

double HRR_Chain::getLogPW0() const{ return logP_w0; }
// no setter for this, dedicated setter below

arma::mat& HRR_Chain::getBeta() //const
{
    // Init to zero
    static arma::mat beta = arma::zeros<arma::mat>(nFixedPredictors+nVSPredictors,nOutcomes); // note the static here
    beta.fill(0.);
    updateGammaMask();
    
    for( unsigned int k=0; k<nOutcomes; ++k)
    {
        arma::uvec VS_IN = gammaMask( arma::find( gammaMask.col(1) == k ) , arma::zeros<arma::uvec>(1) );
        
        if(VS_IN.n_elem>0)
        {
            arma::uvec singleIdx_k = { k };
            arma::mat W_k;
            
            if( preComputedXtX )
            {
                switch ( beta_type )
                {
                    case Beta_Type::gprior :
                    {
                        //W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( XtX(VS_IN,VS_IN) );
                        W_k = (w+temperature)/(w*temperature) * XtX(VS_IN,VS_IN);
                        break;
                    }
                        
                    case Beta_Type::independent :
                    {
                        //W_k = arma::inv_sympd( XtX(VS_IN,VS_IN)/temperature + 1./w * arma::eye<arma::mat>(VS_IN.n_elem,VS_IN.n_elem) );
                        W_k = XtX(VS_IN,VS_IN) + 1./w * arma::eye<arma::mat>(VS_IN.n_elem,VS_IN.n_elem);
                        break;
                    }
                        
                    case Beta_Type::reGroup :
                    {
                        //W_k = arma::inv_sympd( XtX(VS_IN,VS_IN)/temperature + 1./w * arma::eye<arma::mat>(VS_IN.n_elem,VS_IN.n_elem) );
                        W_k = XtX(VS_IN,VS_IN) + arma::diagmat( arma::join_cols(1./w0*arma::ones(nFixedPredictors),1./w*arma::ones(VS_IN.n_elem-nFixedPredictors)) );
                        break;
                    }
                        
                    default:
                        throw Bad_Beta_Type ( beta_type );
                }
            }
            else
            {
                switch ( beta_type )
                {
                    case Beta_Type::gprior :
                    {
                        //W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN) ).t() * data->cols( (*predictorsIdx)(VS_IN) ) ) );
                        W_k = (w+temperature)/(w*temperature) * data->cols( (*predictorsIdx)(VS_IN) ).t() * data->cols( (*predictorsIdx)(VS_IN) );
                        break;
                    }
                        
                    case Beta_Type::independent :
                    {
                        //W_k = arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN) ).t() * data->cols( (*predictorsIdx)(VS_IN) ) )/temperature + 1./w * arma::eye<arma::mat>(VS_IN.n_elem,VS_IN.n_elem) );
                        W_k = data->cols( (*predictorsIdx)(VS_IN) ).t() * data->cols( (*predictorsIdx)(VS_IN) ) + 1./w * arma::eye<arma::mat>(VS_IN.n_elem,VS_IN.n_elem);
                        break;
                    }
                        
                    case Beta_Type::reGroup :
                    {
                      // W_k = arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN) ).t() * data->cols( (*predictorsIdx)(VS_IN) ) )/temperature + 1./w * arma::eye<arma::mat>(VS_IN.n_elem,VS_IN.n_elem) );
                      W_k = data->cols( (*predictorsIdx)(VS_IN) ).t() * data->cols( (*predictorsIdx)(VS_IN) ) + arma::diagmat( arma::join_cols(1./w0*arma::ones(nFixedPredictors),1./w*arma::ones(VS_IN.n_elem-nFixedPredictors)) );
                      break;
                    }
                        
                    default:
                        throw Bad_Beta_Type ( beta_type );
                }
            }
            
            W_k = arma::inv_sympd( W_k ); // might add settings argument arma::inv_sympd( W_k, inv_opts::allow_approx )
            arma::vec mu_k = W_k * ( data->cols( VS_IN ).t() * data->col( (*outcomesIdx)(k) ) / temperature ); // we divide by temp later
            
            beta.submat(VS_IN,singleIdx_k) = Distributions::randMvNormal( mu_k , W_k );
        }
        
    }
    
    return beta;
}

// PREDICTIVE-LIKELIHOOD OF INDIVIDUALS FOR THE SSUR MODEL
arma::mat HRR_Chain::getPredLikelihood() { return predLik ; }
void HRR_Chain::setPredLikelihood( arma::mat predLik_ ){ predLik = predLik_ ; }

// LOG-LIKELIHOOD FOR THE HRR MODEL
double HRR_Chain::getLogLikelihood() const{ return log_likelihood ; }
void HRR_Chain::setLogLikelihood( double log_likelihood_ ){ log_likelihood = log_likelihood_ ; }

double HRR_Chain::getJointLogPrior() const
{
    return logP_o + logP_pi + logP_gamma + logP_w  + logP_w0;
}

double HRR_Chain::getJointLogPosterior() const
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
void HRR_Chain::oInit( arma::vec& o_init , double a_o_ , double b_o_ , double var_o_proposal_ )
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

void HRR_Chain::oInit()
{
    if( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    arma::vec init = arma::ones<arma::vec>(nOutcomes) / std::max( 500. , (double)nVSPredictors ) ;
    oInit( init , 2. , std::max(500., (double)nVSPredictors)-2. , 0.005 );
}

void HRR_Chain::oInit( arma::vec& o_init )
{
    if( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    oInit( o_init , 2. , std::max(500., (double)nVSPredictors)-2. , 0.005 );
}


// this is only for the hotspot
void HRR_Chain::piInit( arma::vec& pi_init , double a_pi_ , double b_pi_ , double var_pi_proposal_ )
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
void HRR_Chain::piInit( arma::vec& pi_init , double a_pi_ , double b_pi_ )
{
    if ( gamma_type != Gamma_Type::hierarchical )
        throw Bad_Gamma_Type ( gamma_type );
    
    pi = pi_init;
    a_pi = a_pi_;
    b_pi = b_pi_;
    logPPi();
}

// this is either hotspot or hierarchical
void HRR_Chain::piInit()
{
    arma::vec init = arma::ones<arma::vec>(nVSPredictors);
    switch ( gamma_type )
    {
        case Gamma_Type::hotspot :
            piInit( init , 2. , 1. , 0.02 );
            break;
            
        case Gamma_Type::hierarchical :
            for(unsigned i=0; i<nOutcomes; i++)
            {
                init(i) = randBeta(1., (double)nOutcomes-1.);
            }
            piInit( init , 1. , (double)nOutcomes-1. );
            break;
            
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
}

void HRR_Chain::piInit( arma::vec& pi_init )
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


void HRR_Chain::mrfGInit()
{
    if( gamma_type != Gamma_Type::mrf )
        throw Bad_Gamma_Type ( gamma_type );
    
    mrf_G = arma::zeros<arma::mat>(2,2);
    mrf_d = -3. ;
    mrf_e = 0.2 ;
}
/*
void HRR_Chain::mrfGInit( arma::mat& mrf_G_ )
{
    if( gamma_type != Gamma_Type::mrf )
        throw Bad_Gamma_Type ( gamma_type );
    
    mrf_G = mrf_G_ ;
    mrf_d = -3. ;
    mrf_e = 0.2 ;
}
*/
void HRR_Chain::gammaInit( arma::umat& gamma_init )
{
    gamma = gamma_init;
    gamma_acc_count = 0.;
    logPGamma();
    updateGammaMask();
}

void HRR_Chain::gammaInit()
{
    arma::umat init = arma::zeros<arma::umat>(nVSPredictors,nOutcomes);
    gammaInit( init );
}

void HRR_Chain::wInit( double w_init , double a_w_ , double b_w_ , double var_w_proposal_ )
{
    w = w_init;
    a_w = a_w_;
    b_w = b_w_;
    
    var_w_proposal = var_w_proposal_ ;
    
    w_acc_count = 0.;
    
    logPW();
}

void HRR_Chain::wInit( double w_init , double a_w_ , double b_w_ )
{
    wInit( w_init , a_w_ , b_w_ , 0.02 );
}


void HRR_Chain::wInit( double w_init )
{
    wInit( w_init , 2. , 5. , 0.02 );
}

void HRR_Chain::wInit()
{
    wInit( 1. , 2. , 5. , 0.02 );
}

void HRR_Chain::w0Init( double w0_init , double a_w0_ , double b_w0_ , double var_w0_proposal_ )
{
  w0 = w0_init;
  a_w0 = a_w0_;
  b_w0 = b_w0_;

  var_w0_proposal = var_w0_proposal_ ;

  w0_acc_count = 0.;
    
  logPW0();
}

void HRR_Chain::w0Init( double w0_init , double a_w0_ , double b_w0_ )
{
  w0Init( w0_init , a_w0_ , b_w0_ , 0.02 );
}


void HRR_Chain::w0Init( double w0_init )
{
  w0Init( w0_init , 2. , 5. , 0.02 );
}

void HRR_Chain::w0Init()
{
  w0Init( 1. , 2. , 5. , 0.02 );
}

// *****************************
// Methods for Log Probabilities
// *****************************

// LOG PRIORS
// logP for all parameters in 3 versions
// empty for re-computing (and updating) the logP given all current values and hyperparameter values
// with one argument for computing with a different value given the current hyperapameter values
// with full arguments for computing the same logP but with different values and different hyperameter values

// o_k
double HRR_Chain::logPO( const arma::vec& o_ , double a_o_ , double b_o_ )
{
    if ( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    double logP = 0.;
    for(unsigned int k=0; k<nOutcomes; ++k)
        logP += Distributions::logPDFBeta( o_(k) , a_o_, b_o_ );
    
    return logP;
}

double HRR_Chain::logPO( )
{
    if ( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    logP_o = logPO( o , a_o , b_o );
    return logP_o;
}

double HRR_Chain::logPO( const arma::vec& o_ )
{
    if ( gamma_type != Gamma_Type::hotspot )
        throw Bad_Gamma_Type ( gamma_type );
    
    return logPO( o_ , a_o , b_o );
}


// pi_j
double HRR_Chain::logPPi( arma::vec& pi_ , double a_pi_ , double b_pi_ )
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

double HRR_Chain::logPPi( )
{
    if ( gamma_type != Gamma_Type::hotspot && gamma_type != Gamma_Type::hierarchical )
        throw Bad_Gamma_Type ( gamma_type );
    
    logP_pi = logPPi( pi , a_pi , b_pi );
    return logP_pi;
}
double HRR_Chain::logPPi( arma::vec& pi_ )
{
    if ( gamma_type != Gamma_Type::hotspot && gamma_type != Gamma_Type::hierarchical )
        throw Bad_Gamma_Type ( gamma_type );
    
    return logPPi( pi_ , a_pi , b_pi );
}

// GAMMA
// this is the hotspot prior
double HRR_Chain::logPGamma( const arma::umat& externalGamma , const arma::vec& o_ , const arma::vec& pi_ )
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
double HRR_Chain::logPGamma( const arma::umat& externalGamma , const arma::vec& pi_ )
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
double HRR_Chain::logPGamma( const arma::umat& externalGamma , double d , double e )
{
    if( gamma_type != Gamma_Type::mrf )
        throw Bad_Gamma_Type ( gamma_type );
    
    arma::mat externalMRFG = mrfG->cols( arma::linspace<arma::uvec>(0,2,3) );
    
    double logP = 0.;
    // calculate the linear and quadratic parts in MRF by using all edges of G
    arma::vec gammaVec = arma::conv_to< arma::vec >::from(arma::vectorise(externalGamma));
    double quad_mrf = 0.;
    double linear_mrf = 0.;
    int count_linear_mrf = 0;
    for( unsigned i=0; i < (externalMRFG).n_rows; ++i )
    {
        if( (externalMRFG)(i,0) != (externalMRFG)(i,1) ){
            quad_mrf += e * 2.0 * gammaVec( (externalMRFG)(i,0) ) * gammaVec( (externalMRFG)(i,1) ) * (externalMRFG)(i,2);
        }else{
                if( gammaVec( (externalMRFG)(i,0) ) == 1 ){
                    linear_mrf += d * gammaVec( (externalMRFG)(i,0) ) * (externalMRFG)(i,2);
                    count_linear_mrf ++;
                }
        }
    }
    logP = arma::as_scalar( linear_mrf + d * (arma::accu( externalGamma ) - count_linear_mrf) + e * 2.0 * quad_mrf );
    
    return logP;
}

// these below are general interfaces
double HRR_Chain::logPGamma( )
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
            logP_gamma = logPGamma( gamma , mrf_d , mrf_e );
            break;
        }
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
    return logP_gamma;
}

double HRR_Chain::logPGamma( const arma::umat& externalGamma )
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
            logP = logPGamma( externalGamma , mrf_d , mrf_e );
            break;
        }
        default:
            throw Bad_Gamma_Type ( gamma_type );
    }
    
    return logP;
}

// W
double HRR_Chain::logPW( double w_ , double a_w_ , double b_w_ )
{
    return Distributions::logPDFIGamma( w_ , a_w_, b_w_ );
}

double HRR_Chain::logPW( )
{
    logP_w = logPW( w , a_w , b_w );
    return logP_w;
}

double HRR_Chain::logPW( double w_ )
{
    return logPW( w_ , a_w , b_w );
}

// W0
double HRR_Chain::logPW0( double w0_ , double a_w0_ , double b_w0_ )
{
  return Distributions::logPDFIGamma( w0_ , a_w0_, b_w0_ );
}

double HRR_Chain::logPW0( )
{
  logP_w0 = logPW0( w0 , a_w0 , b_w0 );
  return logP_w0;
}

double HRR_Chain::logPW0( double w0_ )
{
  return logPW0( w0_ , a_w0 , b_w0 );
}

// PREDICTIVE LIKELIHOOD of the HRR model is not implemented
arma::mat HRR_Chain::predLikelihood( )
{
    return predLik;
}

// LOG LIKELIHOODS
double HRR_Chain::logLikelihood( )
{
    
    double logP {0};
    predLik.set_size(nObservations, nOutcomes);
    
    // yMean is needed if y is not standardized. This might cause MCMC stopped
    arma::mat yMean = data->cols( *outcomesIdx );
    for( unsigned k=0; k<nOutcomes; k++)
    {
        yMean.col(k) = arma::mean(yMean.col(k)) * arma::ones(nObservations);
    }

    #ifdef _OPENMP
    #pragma omp parallel for default(shared) reduction(+:logP)
    #endif
    
    for( unsigned int k=0; k<nOutcomes; ++k)
    {
        arma::uvec VS_IN_k = {}; // be sure it's empty by default
        if(gammaMask.n_rows>0)
            VS_IN_k = gammaMask( arma::find(  gammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
        
        arma::mat W_k;
        if( preComputedXtX )
        {
            switch ( beta_type )
            {
                case Beta_Type::gprior :
                {
                    //W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) );
                    W_k = (w+temperature)/(w*temperature) * XtX(VS_IN_k,VS_IN_k);
                    break;
                }
                    
                case Beta_Type::independent :
                {
                    W_k = XtX(VS_IN_k,VS_IN_k) + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem);
                    break;
                }
                    
                case Beta_Type::reGroup :
                {
                  //W_k = arma::inv_sympd( XtX(VS_IN_k,VS_IN_k)/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                  W_k = XtX(VS_IN_k,VS_IN_k) + arma::diagmat( arma::join_cols(1./w0*arma::ones(nFixedPredictors),1./w*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) );
                  break;
                }
                    
                default:
                    throw Bad_Beta_Type ( beta_type );
            }
        }
        else
        {
            switch ( beta_type )
            {
                case Beta_Type::gprior :
                {
                    //W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) );
                    W_k = (w+temperature)/(w*temperature) * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) );
                    break;
                }
                    
                case Beta_Type::independent :
                {
                    //W_k = arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) )/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                    W_k = data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem);
                    break;
                }
                    
                case Beta_Type::reGroup :
                {
                    //W_k = arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) )/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                    W_k = data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) + arma::diagmat( arma::join_cols(1./w0*arma::ones(nFixedPredictors),1./w*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) );
                    break;
                }
                    
                default:
                    throw Bad_Beta_Type ( beta_type );
            }
        }
        
        W_k = arma::inv_sympd( W_k );
        arma::vec mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * (data->col( (*outcomesIdx)(k) ) - yMean.col(k)) / temperature ); // we divide by temp later
        //arma::vec mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * (data->col( (*outcomesIdx)(k) )) ); // we divide by temp later
        
        double a_sigma_k = a_sigma + 0.5*(double)nObservations/temperature;
        double b_sigma_k = b_sigma + 0.5* arma::as_scalar( (( data->col((*outcomesIdx)(k)) - yMean.col(k) ).t() * (data->col((*outcomesIdx)(k)) - yMean.col(k) )) - ( mu_k.t() * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * ( data->col((*outcomesIdx)(k)) - yMean.col(k) ) ) )/temperature;
        //double b_sigma_k = b_sigma + 0.5* arma::as_scalar( (( data->col((*outcomesIdx)(k)) ).t() * (data->col((*outcomesIdx)(k)) )) - ( mu_k.t() * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * ( data->col((*outcomesIdx)(k)) ) ) )/temperature;
        
        double sign, tmp;
        arma::log_det(tmp, sign, W_k );
        logP += 0.5*tmp;
        
        // arma::log_det(tmp, sign, w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        logP -= 0.5 * (double)VS_IN_k.n_elem * log(w);
        
        logP += a_sigma*log(b_sigma) - a_sigma_k*log(b_sigma_k);
        
        logP += std::lgamma(a_sigma_k) - std::lgamma(a_sigma);
        
        // posterior predictive - t distribution after shifting and scaling by some quantities; from the multivariate t distribution p(y_tilde |y)
        if( output_CPO && (temperature = 1.) ){
            for( unsigned int j=0; j<nObservations; ++j )
            {
                double mu_scale, W_scale, t1, t2;
                
                mu_scale = (data->col(k))(j) - arma::as_scalar(arma::conv_to<arma::vec>::from(data->row(j))((*predictorsIdx)(VS_IN_k)).t() *mu_k);
                W_scale = b_sigma_k/a_sigma_k * ( 1. + arma::as_scalar(arma::conv_to<arma::vec>::from(data->row(j))((*predictorsIdx)(VS_IN_k)).t() * W_k * arma::conv_to<arma::vec>::from(data->row(j))((*predictorsIdx)(VS_IN_k))) );
                
                t1 = std::lgamma(a_sigma_k+0.5) - 0.5*std::log(2.*a_sigma_k*M_PI) - 0.5*std::log( W_scale ) - std::lgamma(a_sigma_k);
                t2 = -(a_sigma_k+0.5) * std::log( 1. + mu_scale*mu_scale/W_scale/2./a_sigma_k );
                
                predLik(j,k) = std::exp( t1+t2 );
            }
        }
        
    }
    
    logP += -log(M_PI)*((double)nObservations*(double)nOutcomes*0.5); // normalising constant remaining from the likelhood
    log_likelihood = logP; // update internal state
    
    return logP;
}

double HRR_Chain::logLikelihood( const arma::umat&  externalGammaMask )
{
    double logP {0};
    
    // yMean is needed if y is not standardized
    arma::mat yMean = data->cols( *outcomesIdx );
    for( unsigned k=0; k<nOutcomes; k++)
    {
        yMean.col(k) = arma::mean(yMean.col(k)) * arma::ones(nObservations);
    }
    
    #ifdef _OPENMP
    #pragma omp parallel for default(shared) reduction(+:logP)
    #endif
    
    for( unsigned int k=0; k<nOutcomes; ++k)
    {
        arma::uvec VS_IN_k = {};
        if(externalGammaMask.n_rows>0)
            VS_IN_k = externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
        
        arma::mat W_k;
        if( preComputedXtX )
        {
            switch ( beta_type )
            {
                case Beta_Type::gprior :
                {
                    //W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) );
                    W_k = (w+temperature)/(w*temperature) * XtX(VS_IN_k,VS_IN_k);
                    break;
                }
                    
                case Beta_Type::independent :
                {
                    //W_k = arma::inv_sympd( XtX(VS_IN_k,VS_IN_k)/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                    W_k = XtX(VS_IN_k,VS_IN_k) + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem);
                    break;
                }
                    
                case Beta_Type::reGroup :
            {
                //W_k = arma::inv_sympd( XtX(VS_IN_k,VS_IN_k)/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                W_k = XtX(VS_IN_k,VS_IN_k) + arma::diagmat( arma::join_cols(1./w0*arma::ones(nFixedPredictors),1./w*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) );
                break;
            }
                    
                default:
                    throw Bad_Beta_Type ( beta_type );
            }
        }
        else
        {
            switch ( beta_type )
            {
                case Beta_Type::gprior :
                {
                    //W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) );
                    W_k = (w+temperature)/(w*temperature) * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) );
                    break;
                }
                    
                case Beta_Type::independent :
                {
                    //W_k = arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) )/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                    W_k = data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem);
                    break;
                }
                    
                case Beta_Type::reGroup :
                {
                    //W_k = arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) )/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                    W_k = data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) + arma::diagmat( arma::join_cols(1./w0*arma::ones(nFixedPredictors),1./w*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) );
                    break;
                }
                    
                default:
                    throw Bad_Beta_Type ( beta_type );
            }
        }
        
        W_k = arma::inv_sympd( W_k );
        arma::vec mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * (data->col( (*outcomesIdx)(k) ) - yMean.col(k)) / temperature ); // we divide by temp later
        //arma::vec mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * (data->col( (*outcomesIdx)(k) )) ); // we divide by temp later
        
        double a_sigma_k = a_sigma + 0.5*(double)nObservations/temperature;
        double b_sigma_k = b_sigma + 0.5* arma::as_scalar( (( data->col((*outcomesIdx)(k)) - yMean.col(k) ).t() * (data->col((*outcomesIdx)(k)) - yMean.col(k) )) - ( mu_k.t() * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * ( data->col((*outcomesIdx)(k)) - yMean.col(k) ) ) )/temperature;
        //double b_sigma_k = b_sigma + 0.5* arma::as_scalar( (( data->col((*outcomesIdx)(k)) ).t() * (data->col((*outcomesIdx)(k)) )) - ( mu_k.t() * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * ( data->col((*outcomesIdx)(k)) ) ) )/temperature;
        
        double sign, tmp; //sign is needed for the implementation, but we 'assume' that all the matrices are (semi-)positive-definite (-> det>=0)
        arma::log_det(tmp, sign, W_k );
        logP += 0.5*tmp;
        
        // arma::log_det(tmp, sign, w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        logP -= 0.5 * (double)VS_IN_k.n_elem * log(w);
        
        logP += a_sigma*log(b_sigma) - a_sigma_k*log(b_sigma_k);
        
        logP += std::lgamma(a_sigma_k) - std::lgamma(a_sigma);
        /*
        if(logP - log(M_PI)*((double)nObservations*(double)nOutcomes*0.5) > 0.){
            std::cout << "...Debug positive logP=" << logP <<
            "; k=" << k <<
            "; beta_type=" << (beta_type == Beta_Type::independent) <<
            "; tmp=" << tmp <<
            "; sign=" << sign <<
            "; w=" << w <<
            "; a_sigma=" << a_sigma <<
            "; b_sigma=" << b_sigma <<
            "; a_sigma_k=" << a_sigma_k <<
            "; b_sigma_k=" << b_sigma_k << std::endl;
        } */
    }
    
    logP += -log(M_PI)*((double)nObservations*(double)nOutcomes*0.5); // normalising constant remaining from the likelhood
    return logP;
    
}

double HRR_Chain::logLikelihood( arma::umat& externalGammaMask , const arma::umat& externalGamma ) // gammaMask , gamma
{
    
    double logP{0};
    externalGammaMask = createGammaMask(externalGamma);
    
    // yMean is needed if y is not standardized
    arma::mat yMean = data->cols( *outcomesIdx );
    for( unsigned k=0; k<nOutcomes; k++)
    {
        yMean.col(k) = arma::mean(yMean.col(k)) * arma::ones(nObservations);
    }
    
    #ifdef _OPENMP
    #pragma omp parallel for default(shared) reduction(+:logP)
    #endif
    
    for( unsigned int k=0; k<nOutcomes; ++k)
    {
        arma::uvec VS_IN_k = {};
        if(externalGammaMask.n_rows>0)
            VS_IN_k = externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
        
        arma::mat W_k;
        if( preComputedXtX )
        {
            switch ( beta_type )
            {
                case Beta_Type::gprior :
                {
                    //W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) );
                    W_k = (w+temperature)/(w*temperature) * XtX(VS_IN_k,VS_IN_k);
                    break;
                }
                    
                case Beta_Type::independent :
                {
                    //W_k = arma::inv_sympd( XtX(VS_IN_k,VS_IN_k)/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                    W_k = XtX(VS_IN_k,VS_IN_k) + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem);
                    break;
                }
                    
                case Beta_Type::reGroup :
                {
                  //W_k = arma::inv_sympd( XtX(VS_IN_k,VS_IN_k)/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                  W_k = XtX(VS_IN_k,VS_IN_k) + arma::diagmat( arma::join_cols(1./w0*arma::ones(nFixedPredictors),1./w*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) );
                  break;
                }
                    
                default:
                    throw Bad_Beta_Type ( beta_type );
            }
        }
        else
        {
            switch ( beta_type )
            {
                case Beta_Type::gprior :
                {
                    //W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) );
                    W_k = (w+temperature)/(w*temperature) * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) );
                    break;
                }
                    
                case Beta_Type::independent :
                {
                    //W_k = arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) )/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                    W_k = data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem);
                    break;
                }
                    
                case Beta_Type::reGroup :
                {
                  //W_k = arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) )/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                  W_k = data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) +  arma::diagmat( arma::join_cols(1./w0*arma::ones(nFixedPredictors),1./w*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) );
                  break;
                }
                    
                default:
                    throw Bad_Beta_Type ( beta_type );
            }
        }
        
        W_k = arma::inv_sympd( W_k );
        arma::vec mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * (data->col( (*outcomesIdx)(k) ) - yMean.col(k))/temperature ); // we divide by temp later
        //arma::vec mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * (data->col( (*outcomesIdx)(k) )) ); // we divide by temp later
        
        double a_sigma_k = a_sigma + 0.5*(double)nObservations/temperature;
        double b_sigma_k = b_sigma + 0.5* arma::as_scalar( (( data->col((*outcomesIdx)(k)) - yMean.col(k) ).t() * (data->col((*outcomesIdx)(k)) - yMean.col(k) )) - ( mu_k.t() * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * ( data->col((*outcomesIdx)(k)) - yMean.col(k) ) ) )/temperature;
        //double b_sigma_k = b_sigma + 0.5* arma::as_scalar( (( data->col((*outcomesIdx)(k)) ).t() * (data->col((*outcomesIdx)(k)) )) - ( mu_k.t() * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * ( data->col((*outcomesIdx)(k)) ) ) )/temperature;
        
        double sign, tmp; //sign is needed for the implementation, but we 'assume' that all the matrices are (semi-)positive-definite (-> det>=0)
        arma::log_det(tmp, sign, W_k );
        logP += 0.5*tmp;
        
        // arma::log_det(tmp, sign, w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        logP -= 0.5 * (double)VS_IN_k.n_elem * log(w);
        
        logP += a_sigma*log(b_sigma) - a_sigma_k*log(b_sigma_k);
        
        logP += std::lgamma(a_sigma_k) - std::lgamma(a_sigma);
    }
    
    logP += -log(M_PI)*((double)nObservations*(double)nOutcomes*0.5); // normalising constant remaining from the likelhood
    return logP;
    
}

double HRR_Chain::logLikelihood( const arma::umat& externalGammaMask , const double externalW, const double externalW0 , const double externalA_sigma, const double externalB_sigma)
{
    
    double logP{0};
    
    // yMean is needed if y is not standardized
    arma::mat yMean = data->cols( *outcomesIdx );
    for( unsigned k=0; k<nOutcomes; k++)
    {
        yMean.col(k) = arma::mean(yMean.col(k)) * arma::ones(nObservations);
    }
   
    #ifdef _OPENMP
    #pragma omp parallel for default(shared) reduction(+:logP)
    #endif
    
    for( unsigned int k=0; k<nOutcomes; ++k)
    {
        arma::uvec VS_IN_k = {};
        if(externalGammaMask.n_rows>0)
            VS_IN_k = externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
        
        arma::mat W_k;
        if( preComputedXtX )
        {
            switch ( beta_type )
            {
                case Beta_Type::gprior :
                {
                    //W_k = (externalW*temperature)/(externalW+temperature) * arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) );
                    W_k = (externalW+temperature)/(externalW*temperature) * XtX(VS_IN_k,VS_IN_k);
                    break;
                }
                    
                case Beta_Type::independent :
                {
                    W_k = XtX(VS_IN_k,VS_IN_k) + 1./externalW * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem);
                    break;
                }
                    
                case Beta_Type::reGroup :
                {
                    // W_k = arma::inv_sympd( XtX(VS_IN_k,VS_IN_k)/temperature + 1./externalW * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                    W_k = XtX(VS_IN_k,VS_IN_k) + arma::diagmat( arma::join_cols(1./externalW0*arma::ones(nFixedPredictors),1./externalW*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) );
                    break;
                }
                    
                default:
                    throw Bad_Beta_Type ( beta_type );
            }
        }
        else
        {
            switch ( beta_type )
            {
                case Beta_Type::gprior :
                {
                    //W_k = (externalW*temperature)/(externalW+temperature) * arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) ) );
                    W_k = (externalW+temperature)/(externalW*temperature) * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) );
                    break;
                }
                    
                case Beta_Type::independent :
                {
                    W_k = data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) + 1./externalW * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem);
                    break;
                }
                    
                case Beta_Type::reGroup :
                {
                  //W_k = arma::inv_sympd( ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) )/temperature + 1./externalW * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
                  W_k = data->cols( (*predictorsIdx)(VS_IN_k) ).t() * data->cols( (*predictorsIdx)(VS_IN_k) ) + arma::diagmat( arma::join_cols(1./externalW0*arma::ones(nFixedPredictors),1./externalW*arma::ones(VS_IN_k.n_elem-nFixedPredictors)) );
                  break;
                }
                    
                default:
                    throw Bad_Beta_Type ( beta_type );
            }
        }
        
        W_k = arma::inv_sympd( W_k );
        arma::vec mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * (data->col( (*outcomesIdx)(k) ) - yMean.col(k)) / temperature ); // we divide by temp later
        //arma::vec mu_k = W_k * ( data->cols( (*predictorsIdx)(VS_IN_k) ).t() * (data->col( (*outcomesIdx)(k) )) ); // we divide by temp later
        
        double a_sigma_k = externalA_sigma + 0.5*(double)nObservations/temperature;
        double b_sigma_k = externalB_sigma + 0.5* arma::as_scalar( (( data->col((*outcomesIdx)(k)) - yMean.col(k) ).t() * (data->col((*outcomesIdx)(k)) - yMean.col(k) )) - ( mu_k.t() * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * ( data->col((*outcomesIdx)(k)) - yMean.col(k) ) ) )/temperature;
        //double b_sigma_k = externalB_sigma + 0.5* arma::as_scalar( (( data->col((*outcomesIdx)(k)) ).t() * (data->col((*outcomesIdx)(k)) )) - ( mu_k.t() * data->cols( (*predictorsIdx)(VS_IN_k) ).t() * ( data->col((*outcomesIdx)(k)) ) ) )/temperature;
        
        double sign, tmp; //sign is needed for the implementation, but we 'assume' that all the matrices are (semi-)positive-definite (-> det>=0)
        arma::log_det(tmp, sign, W_k );
        logP += 0.5*tmp;
        
        // arma::log_det(tmp, sign, w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        logP -= 0.5 * (double)VS_IN_k.n_elem * log(externalW);
        
        logP += externalA_sigma*log(externalB_sigma) - a_sigma_k*log(b_sigma_k);
        
        logP += std::lgamma(a_sigma_k) - std::lgamma(externalA_sigma);
    }
    
    logP += -log(M_PI)*((double)nObservations*(double)nOutcomes*0.5); // initialise with the normalising constant remaining from the likelhood
    return logP;
}
// *********************
// STEP FUNCTION - PERFORM ONE ITERATION FOR THE CHAIN
// *********************

// sampler for proposed updates on gamma
double HRR_Chain::gammaBanditProposal( arma::umat& mutantGamma , arma::uvec& updateIdx , unsigned int& outcomeUpdateIdx )
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

double HRR_Chain::gammaMC3Proposal( arma::umat& mutantGamma , arma::uvec& updateIdx , unsigned int& outcomeUpdateIdx )
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

// MH update (log-normal) -- update one value at each iteration TODO worth doing more?
void HRR_Chain::stepOneO()
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

void HRR_Chain::stepO()
{
    
    arma::vec proposedO = o;
    
    double proposedOPrior, proposedGammaPrior, logAccProb;
    
    for( unsigned int k=0; k<nOutcomes ; ++k )
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
void HRR_Chain::stepOnePi()
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

void HRR_Chain::stepPi()
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

// MH sampling for slab's variance w, also for random effects' w0 if needed
void HRR_Chain::stepW()
{
    double proposedW = std::exp( std::log(w) + randNormal(0.0, var_w_proposal) );
    
    double proposedWPrior = logPW( proposedW );
    double proposedLikelihood = logLikelihood( gammaMask , proposedW , proposedW , a_sigma , b_sigma );
    
    double logAccProb = (proposedWPrior + proposedLikelihood) - (logP_w + log_likelihood);
    
    if( randLogU01() < logAccProb )
    {
        w = proposedW;
        logP_w = proposedWPrior;
        log_likelihood = proposedLikelihood;
        
        ++w_acc_count;
    }
    
    // update w0 if there are random effects
    if( nFixedPredictors > 0 )
    {
        double proposedW0 = std::exp( std::log(w0) + randNormal(0.0, var_w0_proposal) );
        double proposedW0Prior = logPW0( proposedW0 );
        proposedLikelihood = logLikelihood( gammaMask , w , proposedW0 , a_sigma , b_sigma );

        double logAccProb = (proposedW0Prior + proposedLikelihood) - (logP_w0 + log_likelihood);

        if( randLogU01() < logAccProb )
          {
            w0 = proposedW0;
            logP_w0 = proposedW0Prior;
            log_likelihood = proposedLikelihood;
          }
    }
}
/*
void HRR_Chain::stepW0()
{
  //double proposedW = std::exp( std::log(w) + randNormal(0.0, var_w_proposal) ); // no need to sample a new w
  double proposedW0 = std::exp( std::log(w0) + randNormal(0.0, var_w0_proposal) );

  //double proposedWPrior = logPW( proposedW );
  double proposedW0Prior = logPW0( proposedW0 );
  double proposedLikelihood = logLikelihood( gammaMask , w , proposedW0 , a_sigma , b_sigma );

  //double logAccProb = (proposedWPrior + proposedW0Prior + proposedLikelihood) - (logP_w + logP_w0 + log_likelihood);
  double logAccProb = (proposedW0Prior + proposedLikelihood) - (logP_w0 + log_likelihood);

  if( randLogU01() < logAccProb )
    {
      //w = proposedW;
      w0 = proposedW0;
      //logP_w = proposedWPrior;
      logP_w0 = proposedW0Prior;
      log_likelihood = proposedLikelihood;

      ++w_acc_count;
      //++w0_acc_count;
    }
} */

void HRR_Chain::stepGamma()
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
    // update log probabilities
    double proposedGammaPrior = logPGamma( proposedGamma );
    double proposedLikelihood = logLikelihood( proposedGammaMask );
    
    double logAccProb = logProposalRatio +
    ( proposedGammaPrior + proposedLikelihood ) -
    ( logP_gamma + log_likelihood );
    
    if( randLogU01() < logAccProb )
    {
        gamma = proposedGamma;
        gammaMask = proposedGammaMask;
        
        logP_gamma = proposedGammaPrior;
        log_likelihood = proposedLikelihood;
        
        // ++gamma_acc_count;
        gamma_acc_count += 1.; // / updatedOutcomesIdx.n_elem;
    }
    
    //updateGammaMask();
    
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


// this updates all the internal states
void HRR_Chain::step()
{
    updateGammaMask();
    
    // update logP_gamma
    //logPGamma();
    
    // Update HyperParameters
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
    
    // update log_likelihood; this might be redundant?
    //logLikelihood();
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
void HRR_Chain::updateProposalVariances()
{
    double delta, delta2;
    arma::vec deltaVec, delta2Vec;
    
    double adaptationFactor = 0.05;
    
    if( internalIterationCounter == 1 ) // init the mean and second moment
    {
        if ( gamma_type == Gamma_Type::hotspot )
        {
            oEmpiricalMean = arma::log(o);
            oEmpiricalM2 = arma::zeros<arma::vec>(nOutcomes);
            var_o_proposal_init = var_o_proposal;
            
            piEmpiricalMean = arma::log(pi);
            piEmpiricalM2 = arma::zeros<arma::vec>(nVSPredictors);
            var_pi_proposal_init = var_pi_proposal;
        }
        
        wEmpiricalMean = std::log(w);
        wEmpiricalM2 = 0.;
        var_w_proposal_init = var_w_proposal;
        
        w0EmpiricalMean = std::log(w0);
        w0EmpiricalM2 = 0.;
        var_w0_proposal_init = var_w0_proposal;
        
    }else if( internalIterationCounter > 1 )
    {
        // update running averages
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
        delta = std::log(w) - wEmpiricalMean;
        wEmpiricalMean = wEmpiricalMean + ( delta / internalIterationCounter );
        delta2 = std::log(w) - wEmpiricalMean;
        wEmpiricalM2 = wEmpiricalM2 + delta * delta2;
        
        // w0
               delta = std::log(w0) - w0EmpiricalMean;
               w0EmpiricalMean = w0EmpiricalMean + ( delta / internalIterationCounter );
               delta2 = std::log(w0) - w0EmpiricalMean;
               w0EmpiricalM2 = w0EmpiricalM2 + delta * delta2;
    }
    
    // Then if it's actually > n update the proposal variances
    
    if( internalIterationCounter > nObservations  )  // update quantities and variance
    {
        
        // update proposal variances
        if ( gamma_type == Gamma_Type::hotspot )
        {
            var_o_proposal = adaptationFactor * var_o_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * arma::mean( oEmpiricalM2/(internalIterationCounter-1) );
            var_pi_proposal = adaptationFactor * var_pi_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * arma::mean( piEmpiricalM2/(internalIterationCounter-1) );
        }
        var_w_proposal = adaptationFactor * var_w_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * wEmpiricalM2/(internalIterationCounter-1);
        var_w0_proposal = adaptationFactor * var_w0_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * w0EmpiricalM2/(internalIterationCounter-1);
    }
}

// more complex functions could be defined outide through public methods
// but this as a baseline is good to have.


// *******************************
// Global operators between two chains
// *******************************
// asuming nu and other fixed hyperparameters are the same across chains, woudn;t make sense otherwise I think

void HRR_Chain::swapO( std::shared_ptr<HRR_Chain>& that )
{
    arma::vec par = this->getO();
    
    this->setO( that->getO() );
    that->setO( par );
}

void HRR_Chain::swapPi( std::shared_ptr<HRR_Chain>& that )
{
    arma::vec par = this->getPi();
    
    this->setPi( that->getPi() );
    that->setPi( par );
}

void HRR_Chain::swapGamma( std::shared_ptr<HRR_Chain>& that )
{
    arma::umat par = this->getGamma();
    
    this->setGamma( that->getGamma() );
    that->setGamma( par );
}

void HRR_Chain::swapW( std::shared_ptr<HRR_Chain>& that )
{
    double par = this->getW();
    
    this->setW( that->getW() );
    that->setW( par );
}

void HRR_Chain::swapW0( std::shared_ptr<HRR_Chain>& that )
{
  double par = this->getW0();

  this->setW0( that->getW0() );
  that->setW0( par );
}

int HRR_Chain::exchangeGamma_step( std::shared_ptr<HRR_Chain>& that )
{
    // I'm exchanging the gammas AND the betas. So gammaMask, XB and U will follow and we will have to re-compute rhoU for both chains
    //arma::umat swapGammaMask;
    arma::mat swapXB , swapU;
    
    double logLik_1 = this->logLikelihood( that->getGammaMask() );   // note that this and that lik are
    double logLik_2 = that->logLikelihood( this->getGammaMask() ); // important because of temperature
    
    double logPExchange = ( logLik_1 + logLik_2 ) -
    ( this->getLogLikelihood() + that->getLogLikelihood() );
    
    if( randLogU01() < logPExchange )
    {
        // parameters and priors
        this->swapGamma( that );
        
        // loglikelihood and related quantities
        /*swapGammaMask = this->getGammaMask() ;
        this->setGammaMask( that->getGammaMask() );
        that->setGammaMask( swapGammaMask );*/
        this->gammaMask.swap(that->gammaMask);
        
        this->setLogLikelihood( logLik_1 );
        that->setLogLikelihood( logLik_2 );
        
        return 1;
    }else
        return 0;
}


int HRR_Chain::adapt_crossOver_step( std::shared_ptr<HRR_Chain>& that )
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
    std::vector<arma::umat> gammaMask_XO(2);
    gammaMask_XO[0] = createGammaMask(gammaXO[0]);
    gammaMask_XO[1] = createGammaMask(gammaXO[1]);
    
    // log Posterior of the new chains
    double logLikFirst = this->logLikelihood( gammaMask_XO[0] );
    double logLikSecond = that->logLikelihood( gammaMask_XO[1] );
    
    double logPGammaFirst = this->logPGamma( gammaXO[0] );
    double logPGammaSecond = that->logPGamma( gammaXO[1] );
    
    pCrossOver +=   ( logLikFirst + logPGammaFirst +
                     logLikSecond + logPGammaSecond ) -
    ( this->getLogLikelihood() + this->getLogPGamma() +
     that->getLogLikelihood() + that->getLogPGamma() );
    
    if( randLogU01() < pCrossOver )
    {
        // -- first chain
        
        this->setGamma( gammaXO[0] , logPGammaFirst );
        this->setGammaMask( gammaMask_XO[0] );
        
        this->setLogLikelihood( logLikFirst );
        
        // -- second chain
        
        that->setGamma( gammaXO[1] , logPGammaSecond );
        that->setGammaMask( gammaMask_XO[1] );
        
        that->setLogLikelihood( logLikSecond );
        
        return 1;
    }else
        return 0;
    
    
}

int HRR_Chain::uniform_crossOver_step( std::shared_ptr<HRR_Chain>& that )
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
    
    std::vector<arma::umat> gammaMask_XO(2);
    gammaMask_XO[0] = createGammaMask(gammaXO[0]);
    gammaMask_XO[1] = createGammaMask(gammaXO[1]);
    
    // log Posterior of the new chains
    double logLikFirst = this->logLikelihood( gammaMask_XO[0] );
    double logLikSecond = that->logLikelihood( gammaMask_XO[1] );
    
    double logPGammaFirst = this->logPGamma( gammaXO[0] );
    double logPGammaSecond = that->logPGamma( gammaXO[1] );
    
    pCrossOver +=   ( logLikFirst + logPGammaFirst +
                     logLikSecond + logPGammaSecond ) -
    ( this->getLogLikelihood()+ this->getLogPGamma() +
     that->getLogLikelihood() + that->getLogPGamma() );
    
    if( randLogU01() < pCrossOver )
    {
        // -- first chain
        
        this->setGamma( gammaXO[0] , logPGammaFirst );
        this->setGammaMask( gammaMask_XO[0] );
        
        this->setLogLikelihood( logLikFirst );
        
        // -- second chain
        
        that->setGamma( gammaXO[1] , logPGammaSecond );
        that->setGammaMask( gammaMask_XO[1] );
        
        that->setLogLikelihood( logLikSecond );
        
        return 1;
    }else
        return 0;
    
    
}

int HRR_Chain::block_crossOver_step( std::shared_ptr<HRR_Chain>& that , arma::mat& corrMatX , double threshold )
{
    double pCrossOver;
    
    std::vector<arma::umat> gammaXO(2); gammaXO[0] = arma::umat(nVSPredictors,nOutcomes);  gammaXO[1] = arma::umat(nVSPredictors,nOutcomes);
    
    // Propose Crossover
    
    // Select the ONE index to foor the block
    unsigned int predIdx = randIntUniform(0, nVSPredictors-1 ); // pred
    unsigned int outcIdx = randIntUniform(0, nOutcomes-1 ); // outcome
    
    arma::uvec covIdx = arma::find( arma::abs( corrMatX.row(predIdx) ) > threshold );  // this will include the original predIdx
    
    gammaXO[0] = this->getGamma();
    gammaXO[1] = that->getGamma();
    
    for(unsigned int j=0; j<covIdx.n_elem; ++j)
    {
        gammaXO[0](covIdx(j),outcIdx) = that->getGamma()(covIdx(j),outcIdx);
        gammaXO[1](covIdx(j),outcIdx) = this->getGamma()(covIdx(j),outcIdx);
    }
    
    pCrossOver = 0.;  // XO prop probability is weird, how do I compute it? Let's say is symmetric as is determnistic and both comes from the same corrMatX
    
    std::vector<arma::umat> gammaMask_XO(2);
    gammaMask_XO[0] = createGammaMask(gammaXO[0]);
    gammaMask_XO[1] = createGammaMask(gammaXO[1]);
    
    // log Posterior of the new chains
    double logLikFirst = this->logLikelihood( gammaMask_XO[0] );
    double logLikSecond = that->logLikelihood( gammaMask_XO[1] );
    
    double logPGammaFirst = this->logPGamma( gammaXO[0] );
    double logPGammaSecond = that->logPGamma( gammaXO[1] );
    
    pCrossOver +=   ( logLikFirst + logPGammaFirst +
                     logLikSecond + logPGammaSecond ) -
    ( this->getLogLikelihood() + this->getLogPGamma() +
     that->getLogLikelihood() + that->getLogPGamma() );
    
    if( randLogU01() < pCrossOver )
    {
        // -- first chain
        
        this->setGamma( gammaXO[0] , logPGammaFirst );
        this->setGammaMask( gammaMask_XO[0] );
        
        this->setLogLikelihood( logLikFirst );
        
        // -- second chain
        
        that->setGamma( gammaXO[1] , logPGammaSecond );
        that->setGammaMask( gammaMask_XO[1] );
        
        that->setLogLikelihood( logLikSecond );
        
        return 1;
    }else
        return 0;
    
}

void HRR_Chain::swapAll( std::shared_ptr<HRR_Chain>& thatChain )
{
    
    // HARD SWAP cause swapping "this" is not an option
    // swap quantities
    arma::umat swapGammaMask;
    arma::mat swapMat;
    
    swapGammaMask = this->getGammaMask() ;
    this->setGammaMask( thatChain->getGammaMask() );
    thatChain->setGammaMask( swapGammaMask );
    //this->gammaMask.swap(thatChain->gammaMask); // Not sure if function swap() sloves an swapping issue?
    
    // parameters and priors
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
    
    // recompute likelihood
    this->logLikelihood();
    thatChain->logLikelihood();
    
}

int HRR_Chain::globalStep( std::shared_ptr<HRR_Chain>& that )
{
    
    unsigned int globalType = randIntUniform(0,5);
    
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
            
        default:
            break;
    }
    
    return 0;
}

int HRR_Chain::exchangeAll_step( std::shared_ptr<HRR_Chain>& thatChain )
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


// *******************************
// Other Methods
// *******************************

// update relavant quantities
arma::umat HRR_Chain::createGammaMask( const arma::umat& gamma )
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



void HRR_Chain::updateGammaMask()
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

// Bandit-sampling related methods
void HRR_Chain::banditInit()// initialise all the private memebers
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
void HRR_Chain::MC3Init()
{
    n_updates_MC3 = std::ceil( nVSPredictors/40 ); //arbitrary number, should I use something different?
}
