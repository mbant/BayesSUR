#include "HESS_Chain.h"

/*******************************
 * HESS HAS NO PARALLEL OMP FOR NOW
 * even if the likelihood would seem embarassingly parallel
 * the are problems in adding any omp pragma...
 * Note that contrarily from SUR here the model still have sigma in the prior for beta for computaional convenience
*******************************/


// *******************************
// Constructors
// *******************************

// empty, everything is default, except the data
HESS_Chain::HESS_Chain( arma::mat& externalY , arma::mat& externalX , double externalTemperature = 1. ) // Y and X  (and temperature that'll have a default value of 1.)
{

    setData( std::make_shared<arma::mat>(externalY), std::make_shared<arma::mat>(externalX) );

    temperature = externalTemperature;

    internalIterationCounter = 0;

    gammaSamplerType = "bandit";
    gPrior = false;

    banditInit();
    // MC3Init();

    oInit();
    piInit();
    gammaInit();
    wInit();

    sigmaABInit();

    updateGammaMask();

    logLikelihood();

}

// full, every parameter object is initialised from existing objects
HESS_Chain::HESS_Chain( arma::mat& externalY , arma::mat& externalX , // Y and X
        arma::vec& o_init , arma::vec& pi_init , arma::umat& gamma_init , double w_init , // o, pi, gamma, w
        double externalTemperature = 1. ) // temperature
{
    std::string banditSampler = "bandit";
    bool gPrior_= false;
    HESS_Chain( externalY , externalX , o_init , pi_init , 
        gamma_init , w_init , banditSampler , gPrior_ , externalTemperature );
}


// full, every parameter object is initialised from existing objects, plus the type of gamma sampler
HESS_Chain::HESS_Chain( arma::mat& externalY , arma::mat& externalX , // Y and X
        arma::vec& o_init , arma::vec& pi_init , arma::umat& gamma_init , double w_init , // o, pi, gamma, w
        std::string gammaSamplerType_ , bool gPrior_ , double externalTemperature = 1. ) // gamma sampler type , gPrior , temperature
{

    setData( std::make_shared<arma::mat>(externalY), std::make_shared<arma::mat>(externalX) );

    temperature = externalTemperature;

    internalIterationCounter = 0;

    gPrior = gPrior_;

    gammaSamplerType = gammaSamplerType_;
    if( gammaSamplerType == "B" || gammaSamplerType == "bandit" || gammaSamplerType == "Bandit" || gammaSamplerType == "b" )
        banditInit();
    else
       MC3Init();    

    oInit( o_init );
    piInit( pi_init );
    gammaInit( gamma_init );
    wInit( w_init );

    sigmaABInit();

    updateGammaMask();

    logLikelihood();
}

// *******************************
// Getters and Setters
// *******************************

// data
std::shared_ptr<arma::mat> HESS_Chain::getY() const{ return Y ; }

std::shared_ptr<arma::mat> HESS_Chain::getX() const{ return X ; }

void HESS_Chain::setData( std::shared_ptr<arma::mat> externalY , std::shared_ptr<arma::mat> externalX )
{ 
    Y = externalY ;
    
    n = (*Y).n_rows;
    s = (*Y).n_cols;

    X = externalX ;

    if( arma::all( (*X).col(0) == 1. ) )
    {
        // if the interceot is there, remove it temporarilty
        (*X).shed_col(0);
    }

    p = (*X).n_cols;
    corrMatX = arma::cor( *X );  // this need to be on pxp values, not with the intercept

    // standardise data
    (*Y).each_col( [](arma::vec& a){ a = ( a - arma::mean(a) ) / arma::stddev(a); } );
    (*X).each_col( [](arma::vec& a){ a = ( a - arma::mean(a) ) / arma::stddev(a); } );
    
    // add intercept (again)
    (*X).insert_cols(0, arma::ones<arma::vec>(n) );

    // Compute XtX
    if( p < 3000 )  // kinda arbitrary value, how can we assess a more sensible one?
    {
        preComputedXtX = true;
        XtX = (*X).t() * (*X);
    }else{
        preComputedXtX = false;
        XtX = arma::mat();
    }


    // this should never be called once the model is set. 
    // In case it is, all the liklihood-impacted states are reset
    gamma = arma::zeros<arma::umat>(p,s);
    updateGammaMask();
    
    log_likelihood = logLikelihood();

}

arma::mat& HESS_Chain::getXtX(){ return XtX ; }

unsigned int HESS_Chain::getN() const{ return n ; }
unsigned int HESS_Chain::getP() const{ return p ; }
unsigned int HESS_Chain::getS() const{ return s ; }
// no setters as they are linked to the data

// gPrior
void HESS_Chain::gPriorInit() // g Prior can only be init at the start, so no proper set method
{
    if( internalIterationCounter > 0 )
        throw std::runtime_error(std::string("gPrior can only be initialised at the start of the MCMC"));

    // set the boot to true
    gPrior = true;

    // re-initialise the w parameter in line with the new prior, as w now has a different meaning
    wInit( (double)n , 0.5*s + s -1. , 0.5*n*s ); // these values are taken from Lewin 2016
    
    // update internals
    logPW();
    log_likelihood = logLikelihood();

}

bool HESS_Chain::getGPrior() const { return gPrior ; }

// usefull quantities to keep track of
arma::umat& HESS_Chain::getGammaMask(){ return gammaMask; }
void HESS_Chain::setGammaMask( arma::umat  externalGammaMask )
{ 
    gammaMask =  externalGammaMask ; 
}

// MCMC related tuning parameters
double HESS_Chain::getTemperature() const{ return temperature; }
void HESS_Chain::setTemperature( double temp_ )
{ 
    log_likelihood = log_likelihood * temperature / temp_ ; // re-temper the log_likelihood first
    temperature = temp_ ;
}

unsigned int HESS_Chain::getinternalIterationCounter() const{ return internalIterationCounter ; }

std::string HESS_Chain::getGammaSamplerType(){ return gammaSamplerType ; }

void HESS_Chain::setGammaSamplerType( std::string& gammaSamplerType_ )
{
    if( gammaSamplerType != gammaSamplerType_ )
    {
        gammaSamplerType = gammaSamplerType_ ; 
        if( gammaSamplerType == "B" || gammaSamplerType == "bandit" || gammaSamplerType == "Bandit" || gammaSamplerType == "b" )
        {
            banditInit();

        }else if( gammaSamplerType == "MC3" || gammaSamplerType == "mc3" )
        {
            MC3Init();

        }else{
            throw std::runtime_error(std::string("Sampler type can only be Bandit or MC3"));
        }
    }
}

// Bandit-sampling related quantities
unsigned int HESS_Chain::getNUpdatesBandit() const{ return n_updates_bandit ; }
void HESS_Chain::setNUpdatesBandit( unsigned int n_updates_bandit_ ){ n_updates_bandit = n_updates_bandit_ ; }

arma::mat& HESS_Chain::getBanditZeta(){ return banditZeta; }
void HESS_Chain::setBanditZeta( arma::mat banditZeta_ ){ banditZeta = banditZeta_ ; }

arma::mat& HESS_Chain::getBanditAlpha(){ return banditAlpha ; }
void HESS_Chain::setBanditAlpha( arma::mat banditAlpha_ ){ banditAlpha = banditAlpha_ ; }

arma::mat& HESS_Chain::getBanditBeta(){ return banditBeta ; }
void HESS_Chain::setBanditBeta( arma::mat banditBeta_ ){ banditBeta = banditBeta_ ; }

arma::vec& HESS_Chain::getBanditMismatch(){ return mismatch; }
void HESS_Chain::setBanditMismatch( arma::vec mismatch_ ){ mismatch = mismatch_ ; }

arma::vec& HESS_Chain::getBanditNormalisedMismatch(){ return normalised_mismatch ; }
void HESS_Chain::setBanditNormalisedMismatch( arma::vec normalised_mismatch_ ){ normalised_mismatch = normalised_mismatch_ ; }

arma::vec& HESS_Chain::getBanditNormalisedMismatchBackwards(){ return normalised_mismatch_backwards ; }
void HESS_Chain::setBanditNormalisedMismatchBackwards( arma::vec normalised_mismatch_backwards_ ){ normalised_mismatch_backwards = normalised_mismatch_backwards_ ; }

// Parameter states etc

void HESS_Chain::sigmaABInit()
{
    a_sigma = 1.;     // their mean is gonna be b/(a-1)
    b_sigma = 1.;     // with variance b^2/(a-1)^2/(a-2)
}

double HESS_Chain::getSigmaA() const{ return a_sigma ; }
void HESS_Chain::setSigmaA( double a_sigma_ )
{ 
    a_sigma = a_sigma_ ; 
    logLikelihood();
}

double HESS_Chain::getSigmaB() const{ return b_sigma ; }
void HESS_Chain::setSigmaB( double b_sigma_ )
{
    b_sigma = b_sigma_ ; 
    logLikelihood();
}

// o_k
arma::vec& HESS_Chain::getO(){ return o ; }
void HESS_Chain::setO( arma::vec& o_ )
{
    o = o_ ;
    logPO();
}

void HESS_Chain::setO( arma::vec& o_ , double logP_o_ )
{
    o = o_ ;
    logP_o = logP_o_ ;
}

double HESS_Chain::getOA() const{ return a_o ; }
void HESS_Chain::setOA( double a_o_ )
{ 
    a_o = a_o_ ;
    logPO();
}

double HESS_Chain::getOB() const{ return b_o ; }
void HESS_Chain::setOB( double b_o_ )
{
    b_o = b_o_ ;
    logPO();
}

double HESS_Chain::getVarOProposal() const{ return var_o_proposal ; }
void HESS_Chain::setVarOProposal( double var_o_proposal_ ){ var_o_proposal = var_o_proposal_ ; }

double HESS_Chain::getOAccRate() const{ return o_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double HESS_Chain::getLogPO() const{ return logP_o ; }
// no setter for this, dedicated setter below

// pi_j
arma::vec& HESS_Chain::getPi(){ return pi ; }
void HESS_Chain::setPi( arma::vec& pi_ )
{
    pi = pi_ ;
    logPPi();
}

void HESS_Chain::setPi( arma::vec& pi_ , double logP_pi_ )
{
    pi = pi_ ;
    logP_pi = logP_pi_ ;
}

double HESS_Chain::getPiA() const{ return a_pi ; }
void HESS_Chain::setPiA( double a_pi_ )
{
    a_pi = a_pi_ ; 
    logPPi();
}

double HESS_Chain::getPiB() const{ return b_pi ; }
void HESS_Chain::setPiB( double b_pi_ )
{ 
    b_pi = b_pi_ ;
    logPPi();
}

double HESS_Chain::getVarPiProposal() const{ return var_pi_proposal ; }
void HESS_Chain::setVarPiProposal( double var_pi_proposal_ ){ var_pi_proposal = var_pi_proposal_ ; }

double HESS_Chain::getPiAccRate() const{ return pi_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double HESS_Chain::getLogPPi() const{ return logP_pi ; }
// no setter for this, dedicated setter below

// GAMMA
arma::umat& HESS_Chain::getGamma(){ return gamma ; }
void HESS_Chain::setGamma( arma::umat& externalGamma )
{
    gamma = externalGamma ; 
    logPGamma();
    log_likelihood = logLikelihood( gammaMask , gamma ); // update internal state
}

void HESS_Chain::setGamma( arma::umat& externalGamma , double logP_gamma_ )
{
    gamma = externalGamma ; 
    logP_gamma = logP_gamma_ ;
    log_likelihood = logLikelihood( gammaMask , gamma ); // update internal state
}

unsigned int HESS_Chain::getNUpdatesMC3() const{ return n_updates_MC3 ; }
void HESS_Chain::setNUpdatesMC3( unsigned int n_updates_MC3_ ){ n_updates_MC3 = n_updates_MC3_ ; }

double HESS_Chain::getGammaAccRate() const{ return gamma_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double HESS_Chain::getLogPGamma() const{ return logP_gamma ; }
// no setter for this, dedicated setter below

// W
double HESS_Chain::getW() const{ return w ; }
void HESS_Chain::setW( double w_ )
{
    w = w_ ;
    logPW();
    logLikelihood();
}

void HESS_Chain::setW( double w_ , double logP_w_ )
{
    w = w_ ;
    logP_w = logP_w_ ;
    logLikelihood();
}

double HESS_Chain::getWA() const{ return a_w ; }
void HESS_Chain::setWA( double a_w_ )
{ 
    a_w = a_w_ ;
    logPW();
}

double HESS_Chain::getWB() const{ return b_w ; }
void HESS_Chain::setWB( double b_w_ )
{ 
    b_w = b_w_ ;
    logPW();
}

double HESS_Chain::getVarWProposal() const{ return var_w_proposal ; }
void HESS_Chain::setVarWProposal( double var_w_proposal_ ){ var_w_proposal = var_w_proposal_ ; }

double HESS_Chain::getWAccRate() const{ return w_acc_count/(double)internalIterationCounter ; }

double HESS_Chain::getLogPW() const{ return logP_w; }
// no setter for this, dedicated setter below

// LOG-LIKELIHOOD FOR THE HESS MODEL
double HESS_Chain::getLogLikelihood() const{ return log_likelihood ; }
void HESS_Chain::setLogLikelihood( double log_likelihood_ ){ log_likelihood = log_likelihood_ ; }

double HESS_Chain::getJointLogPrior() const
{
    return logP_o + logP_pi + logP_gamma + logP_w;
}

double HESS_Chain::getJointLogPosterior() const
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
void HESS_Chain::oInit( arma::vec& o_init , double a_o_ , double b_o_ , double var_o_proposal_ )
{
    o = o_init;
    a_o = a_o_;
    b_o = b_o_;
    var_o_proposal = var_o_proposal_;
    o_acc_count = 0.;
    
    logPO();
}

void HESS_Chain::oInit()
{
    arma::vec init = arma::ones<arma::vec>(s) / std::max( 500. , (double)p ) ;
    oInit( init , 2. , std::max(500.,(double)p)-2. , 0.005 );
}

void HESS_Chain::oInit( arma::vec& o_init )
{
    oInit( o_init , 2. , std::max(500.,(double)p)-2. , 0.005 );
}

void HESS_Chain::piInit( arma::vec& pi_init , double a_pi_ , double b_pi_ , double var_pi_proposal_ )
{
    pi = pi_init;
    a_pi = a_pi_;
    b_pi = b_pi_;
    var_pi_proposal = var_pi_proposal_;
    pi_acc_count = 0.;
    
    logPPi();
}

void HESS_Chain::piInit()
{
    arma::vec init = arma::ones<arma::vec>(p) ;
    piInit( init , 2. , 1. , 0.02 );
}

void HESS_Chain::piInit( arma::vec& pi_init )
{
    piInit( pi_init , 2. , 1. , 0.02 );
}

void HESS_Chain::gammaInit( arma::umat& gamma_init )
{
    gamma = gamma_init; 
    gamma_acc_count = 0.;
    logPGamma();
    updateGammaMask();
}

void HESS_Chain::gammaInit()
{
    arma::umat init = arma::zeros<arma::umat>(p,s);
    gammaInit( init );
}

void HESS_Chain::wInit( double w_init , double a_w_ , double b_w_ , double var_w_proposal_ )
{
    w = w_init;
    a_w = a_w_;
    b_w = b_w_;

    var_w_proposal = var_w_proposal_ ;

    w_acc_count = 0.;
    
    logPW();
}

void HESS_Chain::wInit( double w_init , double a_w_ , double b_w_ )
{
    wInit( w_init , a_w_ , b_w_ , 0.02 );
}


void HESS_Chain::wInit( double w_init )
{
    wInit( w_init , 2. , 5. , 0.02 );
}

void HESS_Chain::wInit()
{
    wInit( 1. , 2. , 5. , 0.02 );
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
double HESS_Chain::logPO( const arma::vec& o_ , double a_o_ , double b_o_ )
{
    double logP = 0.;
    for(unsigned int k=0; k<s; ++k)
        logP += Distributions::logPDFBeta( o_(k) , a_o_, b_o_ );

    return logP;
}

double HESS_Chain::logPO( )
{
    logP_o = logPO( o , a_o , b_o );
    return logP_o;
}

double HESS_Chain::logPO( const arma::vec& o_ )
{
    return logPO( o_ , a_o , b_o );
}


// pi_j
double HESS_Chain::logPPi( arma::vec& pi_ , double a_pi_ , double b_pi_ )
{
    double logP = 0.;
    for(unsigned int j=0; j<p; ++j)
        logP += Distributions::logPDFGamma( pi_(j) , a_pi_, b_pi_ );

    return logP;    
}

double HESS_Chain::logPPi( )
{
    logP_pi = logPPi( pi , a_pi , b_pi );
    return logP_pi;
}
double HESS_Chain::logPPi( arma::vec& pi_ )
{
    return logPPi( pi_ , a_pi , b_pi );
}

// GAMMA
double HESS_Chain::logPGamma( const arma::umat& externalGamma , const arma::vec& o_ , const arma::vec& pi_ )
{
    double logP = 0.;
    for(unsigned int j=0; j<p; ++j)
    {
        for(unsigned int k=0; k<s; ++k)
        {
            if( ( o_(k) * pi_(j) ) > 1 )
                return -std::numeric_limits<double>::infinity();

            logP += Distributions::logPDFBernoulli( externalGamma(j,k), o_(k) * pi_(j) );
        }
    }
    return logP;
}

double HESS_Chain::logPGamma( )
{
    logP_gamma = logPGamma( gamma , o , pi );
    return logP_gamma;
}

double HESS_Chain::logPGamma( const arma::umat& externalGamma )
{
    return logPGamma( externalGamma , o , pi );
}

// W
double HESS_Chain::logPW( double w_ , double a_w_ , double b_w_ )
{
    return Distributions::logPDFIGamma( w_ , a_w_, b_w_ );
}

double HESS_Chain::logPW( )
{
    logP_w = logPW( w , a_w , b_w );
    return logP_w;
}

double HESS_Chain::logPW( double w_ )
{
    return logPW( w_ , a_w , b_w );
}

// LOG LIKELIHOODS
double HESS_Chain::logLikelihood( )
{
    double logP, sign, tmp; //sign is needed for the implementation, but we 'assume' that all the matrices are (semi-)positive-definite (-> det>=0)

    logP = -log(M_PI)*((double)n*(double)s*0.5); // initialise with the normalising constant remaining from the likelhood

    arma::uvec VS_IN_k;
    arma::mat W_k;
    arma::vec mu_k;
    double a_sigma_k, b_sigma_k;

    // HERE WOULD BE A CANDIDATE FOR private(VS_IN_k, W_k, mu_k, a_sigma_k, b_sigma_k, tmp, sign) reduction(+:logP)
    // but it introduces random arma::inv_sympd errors and problems with the logP =/
    for( unsigned int k=0; k<s; ++k)
    {
        VS_IN_k = gammaMask( arma::find(  gammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
        
        if( preComputedXtX )
        {
            if( gPrior )
                W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) );
            else
                W_k = arma::inv_sympd( XtX(VS_IN_k,VS_IN_k)/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        }
        else
        {
            if( gPrior )
                W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) ) );
            else
                W_k = arma::inv_sympd( ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) )/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        }

        mu_k = W_k * ( (*X).cols(VS_IN_k).t() * (*Y).col(k) ); // we divide by temp later

        a_sigma_k = a_sigma + 0.5*(double)n/temperature;
        b_sigma_k = b_sigma + 0.5* arma::as_scalar( ((*Y).col(k).t() * (*Y).col(k)) - ( mu_k.t() * (*X).cols(VS_IN_k).t() * (*Y).col(k) ) )/temperature;

        arma::log_det(tmp, sign, W_k );
        logP += 0.5*tmp; 

        // arma::log_det(tmp, sign, w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        logP -= 0.5 * (double)VS_IN_k.n_elem * log(w); 

        logP += a_sigma*log(b_sigma) - a_sigma_k*log(b_sigma_k);

        logP += std::lgamma(a_sigma_k) - std::lgamma(a_sigma);
    }

    log_likelihood = logP; // update internal state

    return logP;
}

double HESS_Chain::logLikelihood( const arma::umat&  externalGammaMask )
{
    double logP, sign, tmp; //sign is needed for the implementation, but we 'assume' that all the matrices are (semi-)positive-definite (-> det>=0)

    logP = -log(M_PI)*((double)n*(double)s*0.5); // initialise with the normalising constant remaining from the likelhood

    arma::uvec VS_IN_k;
    arma::mat W_k;
    arma::vec mu_k;
    double a_sigma_k, b_sigma_k;

    // HERE WOULD BE A CANDIDATE FOR private(VS_IN_k, W_k, mu_k, a_sigma_k, b_sigma_k, tmp, sign) reduction(+:logP)
    // but it introduces random arma::inv_sympd errors and problems with the logP =/
    for( unsigned int k=0; k<s; ++k)
    {
        VS_IN_k = externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );

        if( preComputedXtX )
        {
            if( gPrior )
                W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) );
            else
                W_k = arma::inv_sympd( XtX(VS_IN_k,VS_IN_k)/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        }
        else
        {
            if( gPrior )
                W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) ) );
            else
                W_k = arma::inv_sympd( ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) )/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        }

        mu_k = W_k * ( (*X).cols(VS_IN_k).t() * (*Y).col(k) ); // we divide by temp later

        a_sigma_k = a_sigma + 0.5*(double)n/temperature;
        b_sigma_k = b_sigma + 0.5* arma::as_scalar( ((*Y).col(k).t() * (*Y).col(k)) - ( mu_k.t() * (*X).cols(VS_IN_k).t() * (*Y).col(k) ) )/temperature;

        arma::log_det(tmp, sign, W_k );
        logP += 0.5*tmp; 

        // arma::log_det(tmp, sign, w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        logP -= 0.5 * (double)VS_IN_k.n_elem * log(w); 

        logP += a_sigma*log(b_sigma) - a_sigma_k*log(b_sigma_k);

        logP += std::lgamma(a_sigma_k) - std::lgamma(a_sigma);
    }

    return logP;

}                                 

double HESS_Chain::logLikelihood( arma::umat& externalGammaMask , const arma::umat& externalGamma ) // gammaMask , gamma 
{

    double logP, sign, tmp; //sign is needed for the implementation, but we 'assume' that all the matrices are (semi-)positive-definite (-> det>=0)

    externalGammaMask = createGammaMask(externalGamma);

    logP = -log(M_PI)*((double)n*(double)s*0.5); // initialise with the normalising constant remaining from the likelhood

    arma::uvec VS_IN_k;
    arma::mat W_k;
    arma::vec mu_k;
    double a_sigma_k, b_sigma_k;

    // HERE WOULD BE A CANDIDATE FOR private(VS_IN_k, W_k, mu_k, a_sigma_k, b_sigma_k, tmp, sign) reduction(+:logP)
    // but it introduces random arma::inv_sympd errors and problems with the logP =/
    for( unsigned int k=0; k<s; ++k)
    {
        VS_IN_k = externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );

        if( preComputedXtX )
        {
            if( gPrior )
                W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) );
            else
                W_k = arma::inv_sympd( XtX(VS_IN_k,VS_IN_k)/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        }
        else
        {
            if( gPrior )
                W_k = (w*temperature)/(w+temperature) * arma::inv_sympd( ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) ) );
            else
                W_k = arma::inv_sympd( ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) )/temperature + 1./w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        }

        mu_k = W_k * ( (*X).cols(VS_IN_k).t() * (*Y).col(k) ); // we divide by temp later

        a_sigma_k = a_sigma + 0.5*(double)n/temperature;
        b_sigma_k = b_sigma + 0.5* arma::as_scalar( ((*Y).col(k).t() * (*Y).col(k)) - ( mu_k.t() * (*X).cols(VS_IN_k).t() * (*Y).col(k) ) )/temperature;

        arma::log_det(tmp, sign, W_k );
        logP += 0.5*tmp; 

        // arma::log_det(tmp, sign, w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        logP -= 0.5 * (double)VS_IN_k.n_elem * log(w); 

        logP += a_sigma*log(b_sigma) - a_sigma_k*log(b_sigma_k);

        logP += std::lgamma(a_sigma_k) - std::lgamma(a_sigma);
    }

    return logP;

}   

double HESS_Chain::logLikelihood( const arma::umat& externalGammaMask , const double externalW , const double externalA_sigma, const double externalB_sigma)
{

    double logP, sign, tmp; //sign is needed for the implementation, but we 'assume' that all the matrices are (semi-)positive-definite (-> det>=0)
    logP = -log(M_PI)*((double)n*(double)s*0.5); // initialise with the normalising constant remaining from the likelhood

    arma::uvec VS_IN_k;
    arma::mat W_k;
    arma::vec mu_k;
    double a_sigma_k, b_sigma_k;

    // HERE WOULD BE A CANDIDATE FOR private(VS_IN_k, W_k, mu_k, a_sigma_k, b_sigma_k, tmp, sign) reduction(+:logP)
    // but it introduces random arma::inv_sympd errors and problems with the logP =/
    for( unsigned int k=0; k<s; ++k)
    {
        VS_IN_k = externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );

        if( preComputedXtX )
        {
            if( gPrior )
                W_k = (externalW*temperature)/(externalW+temperature) * arma::inv_sympd( XtX(VS_IN_k,VS_IN_k) );
            else
                W_k = arma::inv_sympd( XtX(VS_IN_k,VS_IN_k)/temperature + 1./externalW * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        }
        else
        {
            if( gPrior )
                W_k = (externalW*temperature)/(externalW+temperature) * arma::inv_sympd( ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) ) );
            else
                W_k = arma::inv_sympd( ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) )/temperature + 1./externalW * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        }


        mu_k = W_k * ( (*X).cols(VS_IN_k).t() * (*Y).col(k) );

        a_sigma_k = externalA_sigma + 0.5*(double)n/temperature;
        b_sigma_k = externalB_sigma + 0.5* arma::as_scalar( (*Y).col(k).t() * (*Y).col(k) - ( mu_k.t() * ( (*X).cols(VS_IN_k).t() * (*Y).col(k) ) ) )/temperature;

        arma::log_det(tmp, sign, W_k );
        logP += 0.5*tmp; 

        // arma::log_det(tmp, sign, w * arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        logP -= 0.5 * (double)VS_IN_k.n_elem * log(externalW); 

        logP += externalA_sigma*log(externalB_sigma) - a_sigma_k*log(b_sigma_k);

        logP += std::lgamma(a_sigma_k) - std::lgamma(externalA_sigma);
    }

    return logP;
}  
// *********************
// STEP FUNCTION - PERFORM ONE ITERATION FOR THE CHAIN
// *********************

// sampler for proposed updates on gamma
double HESS_Chain::gammaBanditProposal( arma::umat& mutantGamma , arma::uvec& updateIdx , unsigned int& outcomeIdx )
{

    double logProposalRatio;

    // decide on one outcome
    outcomeIdx = Distributions::randIntUniform(0,s-1);

    // Sample Zs (only for relevant outocome)
    for(unsigned int i=0; i<p; ++i)
    {
            banditZeta(i) = Distributions::randBeta(banditAlpha(i,outcomeIdx),banditAlpha(i,outcomeIdx));
    }

    // Create mismatch (only for relevant outcome)
    for(unsigned int i=0; i<p; ++i)
    {
        mismatch(i) = (mutantGamma(i,outcomeIdx)==0)?(banditZeta(i)):(1.-banditZeta(i));   //mismatch
    }

    // Normalise
    // mismatch = arma::log(mismatch); //logscale ??? TODO
    // normalised_mismatch = mismatch - Utils::logspace_add(mismatch);

    normalised_mismatch = mismatch / arma::as_scalar(arma::sum(mismatch));

    if( Distributions::randU01() < 0.5 )   // one deterministic update
    {

        // Decide which to update
        updateIdx = arma::zeros<arma::uvec>(1);
        updateIdx(0) = Distributions::randWeightedIndexSampleWithoutReplacement(p,normalised_mismatch); // sample the one

        // Update
        mutantGamma(updateIdx(0),outcomeIdx) = 1 - gamma(updateIdx(0),outcomeIdx); // deterministic, just switch

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
        updateIdx = Distributions::randWeightedIndexSampleWithoutReplacement(p,normalised_mismatch,n_updates_bandit); // sample n_updates_bandit indexes

        normalised_mismatch_backwards = mismatch; // copy for backward proposal

        // Update
        for(unsigned int i=0; i<n_updates_bandit; ++i)
        {
            mutantGamma(updateIdx(i),outcomeIdx) = Distributions::randBernoulli(banditZeta(updateIdx(i))); // random update

            normalised_mismatch_backwards(updateIdx(i)) = 1.- normalised_mismatch_backwards(updateIdx(i));

            logProposalRatio += Distributions::logPDFBernoulli(gamma(updateIdx(i),outcomeIdx),banditZeta(updateIdx(i))) -
                Distributions::logPDFBernoulli(mutantGamma(updateIdx(i),outcomeIdx),banditZeta(updateIdx(i)));
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
                
double HESS_Chain::gammaMC3Proposal( arma::umat& mutantGamma , arma::uvec& updateIdx , unsigned int& outcomeIdx )
{
    updateIdx = arma::uvec(n_updates_MC3);

    // decide on one outcome
    outcomeIdx = Distributions::randIntUniform(0,s-1);

    for( unsigned int i=0; i<n_updates_MC3; ++i)
        updateIdx(i) = Distributions::randIntUniform(0,p-1);    // note that I might be updating multiple times the same coeff

    for( auto i : updateIdx)
    mutantGamma(i,outcomeIdx) = ( Distributions::randU01() < 0.5)? gamma(i,outcomeIdx) : 1-gamma(i,outcomeIdx); // could simply be ( 0.5 ? 1 : 0) ;

    return 0. ; // pass this to the outside, it's the (symmetric) logProposalRatio
}



// **************
// **** Methods that update the internal state of their parameter
// **************

// MH update (log-normal) -- update one value at each iteration TODO worth doing more?
void HESS_Chain::stepOneO()
{
    
    unsigned int k = Distributions::randIntUniform(0,s-1);
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

        if( Distributions::randLogU01() < logAccProb )
        {
            o(k) = proposedO(k);
            logP_o = proposedOPrior;
            logP_gamma = proposedGammaPrior;

            ++o_acc_count;
        }
    }

}

void HESS_Chain::stepO()
{
    
    arma::vec proposedO = o;

    double proposedOPrior, proposedGammaPrior, logAccProb;

    for( unsigned int k=0; k < s ; ++k )
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

            if( Distributions::randLogU01() < logAccProb )
            {
                o(k) = proposedO(k);
                logP_o = proposedOPrior;
                logP_gamma = proposedGammaPrior;

                o_acc_count += o_acc_count / (double)s;
            }else
                proposedO(k) = o(k);
        }else
            proposedO(k) = o(k);
    }    

}

// MH update (log-normal) -- update one value at each iteration TODO worth doing more?
void HESS_Chain::stepOnePi()
{
    unsigned int j = Distributions::randIntUniform(0,p-1);
    arma::vec proposedPi = pi;

    double proposedPiPrior, proposedGammaPrior, logAccProb;

    proposedPi(j) = std::exp( std::log( pi(j) ) + Distributions::randNormal(0.0, var_pi_proposal) );

    if( arma::all( ( o * proposedPi(j) ) <= 1 ) )
    {
        proposedPiPrior = logPPi( proposedPi );
        proposedGammaPrior = logPGamma( gamma, o, proposedPi);

        // A/R
        logAccProb = (proposedPiPrior + proposedGammaPrior) - (logP_pi + logP_gamma);

        if( Distributions::randLogU01() < logAccProb )
        {
            pi(j) = proposedPi(j);
            logP_pi = proposedPiPrior;
            logP_gamma = proposedGammaPrior;

            ++pi_acc_count;
        }
    }

}

void HESS_Chain::stepPi()
{
    arma::vec proposedPi = pi;

    double proposedPiPrior, proposedGammaPrior, logAccProb;

    for( unsigned int j=0; j < p ; ++j )
    {
        proposedPi(j) = std::exp( std::log( pi(j) ) + Distributions::randNormal(0.0, var_pi_proposal) );

        if( arma::all( ( o * proposedPi(j) ) <= 1 ) )
        {
            proposedPiPrior = logPPi( proposedPi );
            proposedGammaPrior = logPGamma( gamma, o, proposedPi);

            // A/R
            logAccProb = (proposedPiPrior + proposedGammaPrior) - (logP_pi + logP_gamma);

            if( Distributions::randLogU01() < logAccProb )
            {
                pi(j) = proposedPi(j);
                logP_pi = proposedPiPrior;
                logP_gamma = proposedGammaPrior;

                pi_acc_count += pi_acc_count / (double)p;
            }else
                proposedPi(j) = pi(j);
        }else
            proposedPi(j) = pi(j);
    }
}

// Gibbs sampler NOT available here
// so MH sampler
void HESS_Chain::stepW()
{
    double proposedW = std::exp( std::log(w) + Distributions::randNormal(0.0, var_w_proposal) );

    double proposedWPrior = logPW( proposedW );
    double proposedLikelihood = logLikelihood( gammaMask , proposedW , a_sigma , b_sigma );

    double logAccProb = (proposedWPrior + proposedLikelihood) - (logP_w + log_likelihood);

    if( Distributions::randLogU01() < logAccProb )
    {
        w = proposedW;
        logP_w = proposedWPrior;
        log_likelihood = proposedLikelihood;

        ++w_acc_count;
    }
}

void HESS_Chain::stepGamma()
{
    arma::umat proposedGamma = gamma;
    arma::uvec updateIdx;
    unsigned int outcomeIdx;

    double logProposalRatio = 0;

    // Update the proposed Gamma
    if( gammaSamplerType == "B" || gammaSamplerType == "bandit" || gammaSamplerType == "Bandit" || gammaSamplerType == "b" )
    {
        logProposalRatio += gammaBanditProposal( proposedGamma , updateIdx , outcomeIdx );

    }else if( gammaSamplerType == "MC3" || gammaSamplerType == "mc3" )
    {
        logProposalRatio += gammaMC3Proposal( proposedGamma , updateIdx , outcomeIdx );

    }else{
        logProposalRatio += gammaBanditProposal( proposedGamma , updateIdx , outcomeIdx ); // default
    }


    // given proposedGamma now, sample a new proposedBeta matrix and corresponging quantities
    arma::umat proposedGammaMask = createGammaMask( proposedGamma );

    // now only one outcome is updated
    // arma::uvec updatedOutcomesIdx = arma::unique( arma::floor( updateIdx / p )); // every p I get to the new column, and as columns correspond to outcomes ... 
    
    // update log probabilities
    double proposedGammaPrior = logPGamma( proposedGamma );
    double proposedLikelihood = logLikelihood( proposedGammaMask );
    
    double logAccProb = logProposalRatio +
                ( proposedGammaPrior + proposedLikelihood ) - 
                ( logP_gamma + log_likelihood );

    if( Distributions::randLogU01() < logAccProb )
    {
        gamma = proposedGamma;
        gammaMask = proposedGammaMask;
        
        logP_gamma = proposedGammaPrior;
        log_likelihood = proposedLikelihood;

        // ++gamma_acc_count;
        gamma_acc_count += 1.; // / updatedOutcomesIdx.n_elem;
    }

    // after A/R, update bandit Related variables
    if(!( gammaSamplerType == "MC3" || gammaSamplerType == "mc3" ))
    {
        for(arma::uvec::iterator iter = updateIdx.begin(); iter != updateIdx.end(); ++iter)
        {
            // FINITE UPDATE
            if( banditAlpha(*iter,outcomeIdx) + banditBeta(*iter,outcomeIdx) <= banditLimit )
            {
                banditAlpha(*iter,outcomeIdx) += banditIncrement * gamma(*iter,outcomeIdx);
                banditBeta(*iter,outcomeIdx) += banditIncrement * (1-gamma(*iter,outcomeIdx));
            }

            // // CONTINUOUS UPDATE
            // banditAlpha(*iter,outcomeIdx) += banditIncrement * gamma(*iter,outcomeIdx);
            // banditBeta(*iter,outcomeIdx) += banditIncrement * (1-gamma(*iter,outcomeIdx));

            // // then renormalise them
            // if( banditAlpha(*iter,outcomeIdx) + banditBeta(*iter) > banditLimit )
            // {
            //     banditAlpha(*iter,outcomeIdx) = banditLimit * ( banditAlpha(*iter,outcomeIdx) / ( banditAlpha(*iter,outcomeIdx) + banditBeta(*iter,outcomeIdx) ));
            //     banditBeta(*iter,outcomeIdx) = banditLimit * (1. - ( banditAlpha(*iter,outcomeIdx) / ( banditAlpha(*iter,outcomeIdx) + banditBeta(*iter,outcomeIdx) )) );
            // }
        }
    }
}


// this updates all the internal states
void HESS_Chain::step()
{
    // Update HyperParameters
    stepW();
    
    for( auto i=0; i<5; ++i){
        stepOneO();
        stepOnePi();
    }
    
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
void HESS_Chain::updateProposalVariances()
{
    double delta, delta2;
    arma::vec deltaVec, delta2Vec;

    double adaptationFactor = 0.05;

    if( internalIterationCounter == 1 ) // init the mean and second moment
    {
        wEmpiricalMean = std::log(w);
        oEmpiricalMean = arma::log(o);
        piEmpiricalMean = arma::log(pi);

        wEmpiricalM2 = 0.;
        oEmpiricalM2 = arma::zeros<arma::vec>(s);
        piEmpiricalM2 = arma::zeros<arma::vec>(p);

        var_w_proposal_init = var_w_proposal;
        var_o_proposal_init = var_o_proposal;
        var_pi_proposal_init = var_pi_proposal;

    }else if( internalIterationCounter > 1 ) 
    {
        // update running averages

        // tau
        delta = std::log(w) - wEmpiricalMean;
        wEmpiricalMean = wEmpiricalMean + ( delta / internalIterationCounter );
        delta2 = std::log(w) - wEmpiricalMean;
        wEmpiricalM2 = wEmpiricalM2 + delta * delta2;

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

    // Then if it's actually > n update the proposal variances

    if( internalIterationCounter > n  )  // update quantities and variance
    {

        // update proposal variances
        var_w_proposal = adaptationFactor * var_w_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * wEmpiricalM2/(internalIterationCounter-1);
        var_o_proposal = adaptationFactor * var_o_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * arma::mean( oEmpiricalM2/(internalIterationCounter-1) );
        var_pi_proposal = adaptationFactor * var_pi_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * arma::mean( piEmpiricalM2/(internalIterationCounter-1) );

    }
}

// more complex functions could be defined outide through public methods
// but this as a baseline is good to have.


// *******************************
// Global operators between two chains
// *******************************
// asuming nu and other fixed hyperparameters are the same across chains, woudn;t make sense otherwise I think 

void HESS_Chain::swapO( std::shared_ptr<HESS_Chain>& that )
{
    arma::vec par = this->getO();

    this->setO( that->getO() );
    that->setO( par );
}

void HESS_Chain::swapPi( std::shared_ptr<HESS_Chain>& that )
{
    arma::vec par = this->getPi();

    this->setPi( that->getPi() );
    that->setPi( par );
}

void HESS_Chain::swapGamma( std::shared_ptr<HESS_Chain>& that )
{
    arma::umat par = this->getGamma();

    this->setGamma( that->getGamma() );
    that->setGamma( par );
}

void HESS_Chain::swapW( std::shared_ptr<HESS_Chain>& that )
{
    double par = this->getW();

    this->setW( that->getW() );
    that->setW( par );
}

int HESS_Chain::exchangeGamma_step( std::shared_ptr<HESS_Chain>& that )
{
    // I'm exchanging the gammas AND the betas. So gammaMask, XB and U will follow and we will have to re-compute rhoU for both chains
    arma::umat swapGammaMask;
    arma::mat swapXB , swapU;

    double logLik_1 = this->logLikelihood( that->getGammaMask() );   // note that this and that lik are
    double logLik_2 = that->logLikelihood( this->getGammaMask() ); // important because of temperature

    double logPExchange = ( logLik_1 + logLik_2 ) -
                        ( this->getLogLikelihood() + that->getLogLikelihood() );
 
    if( Distributions::randLogU01() < logPExchange )
    {
        // parameters and priors
        this->swapGamma( that );

        // loglikelihood and related quantities
        swapGammaMask = this->getGammaMask() ;
        this->setGammaMask( that->getGammaMask() );
        that->setGammaMask( swapGammaMask );

        this->setLogLikelihood( logLik_1 );
        that->setLogLikelihood( logLik_2 );

        return 1;
    }else
        return 0;
}


int HESS_Chain::adapt_crossOver_step( std::shared_ptr<HESS_Chain>& that )
{
    double pCrossOver;

    // Crossover operator hyper pars (see http://www3.stat.sinica.edu.tw/statistica/oldpdf/A10n21.pdf
	double pXO_0 = 0.1, pXO_1 = 0.2 , pXO_2 = 0.2;
	double p11 = pXO_0*pXO_0 + (1.-pXO_0)*(1.-pXO_0) ,p12 = 2.*pXO_0*(1.-pXO_0) ,p21= pXO_1*(1.-pXO_2) + pXO_2*(1.-pXO_1) ,p22 = pXO_1*pXO_2 + (1.-pXO_1)*(1.-pXO_2);

    unsigned int n11,n12,n21,n22;

    std::vector<arma::umat> gammaXO(2); gammaXO[0] = arma::umat(p,s);  gammaXO[1] = arma::umat(p,s); 

    // Propose Crossover
    n11=0;n12=0;n21=0;n22=0;

    for(unsigned int j=0; j<p; ++j)
    {
        for(unsigned int k=0; k<s; ++k)
        {
            if ( this->getGamma()(j,k) == that->getGamma()(j,k) )
            {
                gammaXO[0](j,k) = this->getGamma()(j,k);
                gammaXO[1](j,k) = this->getGamma()(j,k);

                gammaXO[0](j,k) = ( Distributions::randU01() < pXO_0 )? 1-gammaXO[0](j,k) : gammaXO[0](j,k);
                gammaXO[1](j,k) = ( Distributions::randU01() < pXO_0 )? 1-gammaXO[1](j,k) : gammaXO[1](j,k);

                if( gammaXO[0](j,k) == gammaXO[1](j,k) )
                    ++n11;
                else
                    ++n12;
            }
            else
            {
                gammaXO[0](j,k) = this->getGamma()(j,k);
                gammaXO[1](j,k) = that->getGamma()(j,k);

                gammaXO[0](j,k) = ( Distributions::randU01() < pXO_1 )? 1-gammaXO[0](j,k) : gammaXO[0](j,k);
                gammaXO[1](j,k) = ( Distributions::randU01() < pXO_2 )? 1-gammaXO[1](j,k) : gammaXO[1](j,k);

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
    
    if( Distributions::randLogU01() < pCrossOver )
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

int HESS_Chain::uniform_crossOver_step( std::shared_ptr<HESS_Chain>& that )
{
    double pCrossOver;

    std::vector<arma::umat> gammaXO(2); gammaXO[0] = arma::umat(p,s);  gammaXO[1] = arma::umat(p,s); 

    // Propose Crossover
    for(unsigned int j=0; j<p; ++j)
    {
        for(unsigned int k=0; k<s; ++k)
        {
            if( Distributions::randU01() < 0.5 )
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
    
    if( Distributions::randLogU01() < pCrossOver )
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

int HESS_Chain::block_crossOver_step( std::shared_ptr<HESS_Chain>& that , arma::mat& corrMatX , double threshold )
{
    double pCrossOver;

    std::vector<arma::umat> gammaXO(2); gammaXO[0] = arma::umat(p,s);  gammaXO[1] = arma::umat(p,s); 

    // Propose Crossover

    // Select the ONE index to foor the block
    unsigned int predIdx = Distributions::randIntUniform(0, p-1 ); // pred
    unsigned int outcIdx = Distributions::randIntUniform(0, s-1 ); // outcome

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
    
    if( Distributions::randLogU01() < pCrossOver )
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

void HESS_Chain::swapAll( std::shared_ptr<HESS_Chain>& thatChain )
{
    
    // HARD SWAP cause swapping "this" is not an option
    // swap quantities
    arma::umat swapGammaMask;
    arma::mat swapMat;

    swapGammaMask = this->getGammaMask() ;
    this->setGammaMask( thatChain->getGammaMask() );
    thatChain->setGammaMask( swapGammaMask );
    
    // parameters and priors
    this->swapO( thatChain );
    this->swapPi( thatChain );
    this->swapGamma( thatChain );

    this->swapW( thatChain );

    // recompute likelihood
    this->logLikelihood();
    thatChain->logLikelihood();

}

int HESS_Chain::globalStep( std::shared_ptr<HESS_Chain>& that )
{

    unsigned int globalType = Distributions::randIntUniform(0,3);

    // std::cout << globalType << std::flush;

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
        
        default: 
            break;
    }

    return 0;
}

// *******************************
// Other Methods
// *******************************

// update relavant quantities
arma::umat HESS_Chain::createGammaMask( const arma::umat& gamma )
{

    // CREATE HERE THE GAMMA "MASK"
    // INITIALISE THE INDEXES FOR THE GAMMA MASK
    arma::umat mask = arma::zeros<arma::umat>(s,2); //this is just an initialisation
    arma::uvec tmpUVec;
    unsigned int tmpIdx;
    for(unsigned int k=0 ; k<s ; ++k)  //add gammas for the intercepts
    {
        mask(k,0) = 0; mask(k,1) = k;
    }

    for(unsigned int k=0 ; k<s ; ++k)  //add the other gammas
    {
        tmpUVec = arma::find(gamma.col(k) != 0);
        tmpIdx = mask.n_rows;

        if( tmpUVec.n_elem > 0 )
        {
            mask.insert_rows( tmpIdx , arma::zeros<arma::umat>( tmpUVec.n_elem , 2 ));
            mask.submat( tmpIdx, 0, mask.n_rows-1 , 0 ) = tmpUVec + 1 ; // +1 cause gamma doesn't have the intercept
            mask.submat( tmpIdx, 1, mask.n_rows-1 , 1 ).fill(k);
        }
    }
    // Gamma mask done

    return mask;
}


void HESS_Chain::updateGammaMask()
{
    // CREATE HERE THE GAMMA "MASK"
    // INITIALISE THE INDEXES FOR THE GAMMA MASK
    gammaMask.zeros(s,2); //this is just an initialisation  -- 2*s means at least all the intercept plus one other covariate for each equation
    arma::uvec tmpUVec;
    unsigned int tmpIdx;
    for(unsigned int k=0 ; k<s ; ++k)  //add gammas for the intercepts
    {
        gammaMask(k,0) = 0; gammaMask(k,1) = k;
    }

    for(unsigned int k=0 ; k<s ; ++k)  //add gammas for the intercepts
    {
        tmpUVec = arma::find(gamma.col(k) != 0);
        tmpIdx = gammaMask.n_rows;

        if( tmpUVec.n_elem > 0 )
        {
            gammaMask.insert_rows( tmpIdx , arma::zeros<arma::umat>( tmpUVec.n_elem , 2 ));
            gammaMask.submat( tmpIdx, 0, gammaMask.n_rows-1 , 0 ) = tmpUVec + 1 ;
            gammaMask.submat( tmpIdx, 1, gammaMask.n_rows-1 , 1 ).fill(k);
        }
    }
    // Gamma mask updated
}

// Bandit-sampling related methods
void HESS_Chain::banditInit()// initialise all the private memebers
{
    banditZeta = arma::vec(p);

    banditAlpha = arma::mat(p,s);
    banditAlpha.fill( 0.5 );
    
    banditBeta = arma::mat(p,s);
    banditBeta.fill( 0.5 );

    mismatch = arma::vec(p);
    normalised_mismatch = arma::vec(p);
    normalised_mismatch_backwards = arma::vec(p);

    n_updates_bandit = 4; // this needs to be low as its O(n_updates!)

    banditLimit = (double)n;
    banditIncrement = 1.;
}

// MC3 init
void HESS_Chain::MC3Init()
{
    n_updates_MC3 = std::ceil( p/40 ); //arbitrary number, should I use something different?
}