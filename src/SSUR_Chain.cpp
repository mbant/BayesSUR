#include "SSUR_Chain.h"

// *******************************
// Constructors
// *******************************

// empty, everything is default, except the data
SSUR_Chain::SSUR_Chain( arma::mat& externalY , arma::mat& externalX , double externalTemperature = 1. )
{


    setData( std::make_shared<arma::mat>(externalY), std::make_shared<arma::mat>(externalX) );

    temperature = externalTemperature;

    internalIterationCounter = 0;
    jtStartIteration = 0;

    gammaSamplerType = "Bandit";

    banditInit();
    MC3Init();    

    tauInit();
    etaInit();
    jtInit();
    oInit();
    piInit();
    gammaInit();
    wInit();

    betaInit();
    sigmaRhoInit();

    updateQuantities();

    logLikelihood();

    // init for sigma rho and beta to reasonable values -- one step of gibbs
    stepSigmaRhoAndBeta();

}

SSUR_Chain::SSUR_Chain( arma::mat& externalY , arma::mat& externalX , std::string& gammaSamplerType_ , double externalTemperature = 1. ):
    SSUR_Chain(externalY, externalX, externalTemperature){ gammaSamplerType = gammaSamplerType_ ; }

// full, every parameter object is initialised from existing objects
SSUR_Chain::SSUR_Chain( arma::mat& externalY , arma::mat& externalX , // Y and X
        double tau_init , double eta_init , JunctionTree& jt_init , arma::mat& sigmaRho_init , // tau, eta, jt, sigmaRho 
        arma::vec& o_init , arma::vec& pi_init , arma::umat& gamma_init , double w_init , arma::mat& beta_init , // o, pi, gamma, w, beta
        double externalTemperature = 1. ) // temperature
{
    std::string banditSampler = "bandit";
    SSUR_Chain( externalY , externalX , tau_init , eta_init , jt_init , sigmaRho_init , o_init , pi_init , 
        gamma_init , w_init , beta_init , banditSampler , externalTemperature );
}


// full, every parameter object is initialised from existing objects, plus the type of gamma sampler
SSUR_Chain::SSUR_Chain( arma::mat& externalY , arma::mat& externalX , // Y and X
        double tau_init , double eta_init , JunctionTree& jt_init , arma::mat& sigmaRho_init , // tau, eta, jt, sigmaRho 
        arma::vec& o_init , arma::vec& pi_init , arma::umat& gamma_init , double w_init , arma::mat& beta_init , // o, pi, gamma, w, beta
        std::string gammaSamplerType_ , double externalTemperature = 1. ) // gamma sampler type , temperature
{

    setData( std::make_shared<arma::mat>(externalY), std::make_shared<arma::mat>(externalX) );

    temperature = externalTemperature;

    internalIterationCounter = 0;
    jtStartIteration = 0;

    gammaSamplerType = gammaSamplerType_;
    if( gammaSamplerType == "B" || gammaSamplerType == "bandit" || gammaSamplerType == "Bandit" || gammaSamplerType == "b" )
        banditInit();
    MC3Init();    

    tauInit( tau_init );
    etaInit( eta_init );
    jtInit( jt_init );
    oInit( o_init );
    piInit( pi_init );
    gammaInit( gamma_init );
    wInit( w_init );

    betaInit( beta_init );
    sigmaRhoInit( sigmaRho_init );

    updateQuantities();

    logLikelihood();
}

// *******************************
// Getters and Setters
// *******************************

// data
std::shared_ptr<arma::mat> SSUR_Chain::getY() const{ return Y ; }

std::shared_ptr<arma::mat> SSUR_Chain::getX() const{ return X ; }

void SSUR_Chain::setData( std::shared_ptr<arma::mat> externalY , std::shared_ptr<arma::mat> externalX )
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
    XtX = (*X).t() * (*X);


    // this should never be called once the model is set. 
    // In case it is, all the liklihood-impacted states are reset
    gamma = arma::zeros<arma::umat>(p,s);
    beta = arma::zeros<arma::mat>(p+1,s);
    jt = JunctionTree();
    sigmaRho = arma::eye(s,s);
    
    log_likelihood = -std::numeric_limits<double>::infinity();

}

arma::mat& SSUR_Chain::getXtX(){ return XtX ; }

unsigned int SSUR_Chain::getN() const{ return n ; }
unsigned int SSUR_Chain::getP() const{ return p ; }
unsigned int SSUR_Chain::getS() const{ return s ; }
// no setters as they are linked to the data

// usefull quantities to keep track of
arma::umat& SSUR_Chain::getGammaMask(){ return gammaMask; }
void SSUR_Chain::setGammaMask( arma::umat  externalGammaMask ){ gammaMask =  externalGammaMask ; }

arma::mat& SSUR_Chain::getXB(){ return XB; }
void SSUR_Chain::setXB( arma::mat externalXB ){ XB = externalXB ; }

arma::mat& SSUR_Chain::getU(){ return U ; }
void SSUR_Chain::setU( arma::mat externalU ){ U = externalU ; }

arma::mat& SSUR_Chain::getRhoU(){ return rhoU ; }
void SSUR_Chain::setRhoU( arma::mat externalRhoU ){ rhoU = externalRhoU ; }

// MCMC related tuning parameters
double SSUR_Chain::getTemperature() const{ return temperature; }
void SSUR_Chain::setTemperature( double temp_ )
{ 
    log_likelihood = log_likelihood * temperature / temp_ ; // re-temper the log_likelihood first
    temperature = temp_ ;
}

unsigned int SSUR_Chain::getinternalIterationCounter() const{ return internalIterationCounter ; }

unsigned int SSUR_Chain::getJTStartIteration() const{ return jtStartIteration; }
void SSUR_Chain::setJTStartIteration( const unsigned int jts ){ jtStartIteration = jts; }

std::string SSUR_Chain::getGammaSamplerType(){ return gammaSamplerType ; }
void SSUR_Chain::setGammaSamplerType( std::string& gammaSamplerType_ ){ gammaSamplerType = gammaSamplerType_ ; }

// Bandit-sampling related quantities
unsigned int SSUR_Chain::getNUpdatesBandit() const{ return n_updates_bandit ; }
void SSUR_Chain::setNUpdatesBandit( unsigned int n_updates_bandit_ ){ n_updates_bandit = n_updates_bandit_ ; }

arma::mat& SSUR_Chain::getBanditZeta(){ return banditZeta; }
void SSUR_Chain::setBanditZeta( arma::mat banditZeta_ ){ banditZeta = banditZeta_ ; }

arma::mat& SSUR_Chain::getBanditAlpha(){ return banditAlpha ; }
void SSUR_Chain::setBanditAlpha( arma::mat banditAlpha_ ){ banditAlpha = banditAlpha_ ; }

arma::mat& SSUR_Chain::getBanditBeta(){ return banditBeta ; }
void SSUR_Chain::setBanditBeta( arma::mat banditBeta_ ){ banditBeta = banditBeta_ ; }

arma::vec& SSUR_Chain::getBanditMismatch(){ return mismatch; }
void SSUR_Chain::setBanditMismatch( arma::vec mismatch_ ){ mismatch = mismatch_ ; }

arma::vec& SSUR_Chain::getBanditNormalisedMismatch(){ return normalised_mismatch ; }
void SSUR_Chain::setBanditNormalisedMismatch( arma::vec normalised_mismatch_ ){ normalised_mismatch = normalised_mismatch_ ; }

arma::vec& SSUR_Chain::getBanditNormalisedMismatchBackwards(){ return normalised_mismatch_backwards ; }
void SSUR_Chain::setBanditNormalisedMismatchBackwards( arma::vec normalised_mismatch_backwards_ ){ normalised_mismatch_backwards = normalised_mismatch_backwards_ ; }

// Parameter states etc

// TAU
double SSUR_Chain::getTau() const{ return tau; }
void SSUR_Chain::setTau( double tau_ )
{
    tau = tau_;
    logPTau(); // this updates the internal logP variable
}

void SSUR_Chain::setTau( double tau_ , double logP_tau_)
{
    tau = tau_;
    logP_tau = logP_tau_ ;
}

double SSUR_Chain::getTauA() const{ return a_tau; }
void SSUR_Chain::setTauA( double a_tau_ )
{ 
    a_tau = a_tau_ ;
    logPTau(); 
}

double SSUR_Chain::getTauB() const{ return b_tau; }
void SSUR_Chain::setTauB( double b_tau_ )
{ 
    b_tau = b_tau_ ; 
    logPTau();
}

double SSUR_Chain::getVarTauProposal() const{ return var_tau_proposal ; }
void SSUR_Chain::setVarTauProposal( double var_tau_proposal_ ){ var_tau_proposal = var_tau_proposal_ ; }

double SSUR_Chain::getTauAccRate() const{ return tau_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double SSUR_Chain::getLogPTau() const{ return logP_tau ; }
// no setter for this, dedicated setter below

// ETA
double SSUR_Chain::getEta() const{ return eta ; }
void SSUR_Chain::setEta( double eta_ )
{
    eta = eta_ ;
    logPEta(); // this updates the internal logP variable
}

void SSUR_Chain::setEta( double eta_ , double logP_eta_ )
{
    eta = eta_ ;
    logP_eta = logP_eta_ ;
}

double SSUR_Chain::getEtaA() const{ return a_eta ; }
void SSUR_Chain::setEtaA( double a_eta_ )
{ 
    a_eta = a_eta_ ;
    logPEta();
}

double SSUR_Chain::getEtaB() const{ return b_eta ; }
void SSUR_Chain::setEtaB( double b_eta_ )
{
    b_eta = b_eta_ ; 
    logPEta();
}

double SSUR_Chain::getLogPEta() const{ return logP_eta ; }
// no setter for this, dedicated setter below

// JT
JunctionTree& SSUR_Chain::getJT(){ return jt; }  // need to check this works correctly TODO (i.e. correct behaviour is returning a copy of jt)
arma::sp_umat SSUR_Chain::getGAdjMat() const{ return jt.getAdjMat(); }  // need to check this works correctly TODO (i.e. correct behaviour is returning a copy of jt)
void SSUR_Chain::setJT( JunctionTree& externalJT )
{ 
    jt = externalJT ; // old jt gets destroyed as there's no reference to him anymore, new jt "points" to the foreign object
    logPJT();
} 

void SSUR_Chain::setJT( JunctionTree& externalJT , double logP_jt_ )
{ 
    jt = externalJT ; // old jt gets destroyed as there's no reference to him anymore, new jt "points" to the foreign object
    logP_jt = logP_jt_ ;
} 

unsigned int SSUR_Chain::getNUpdatesJT() const{ return n_updates_jt; }
void SSUR_Chain::setNUpdatesJT( unsigned int n_updates_jt_ ){ n_updates_jt = n_updates_jt_ ; }

double SSUR_Chain::getJTAccRate() const{ return jt_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double SSUR_Chain::getLogPJT() const{ return logP_jt ; }
// no setter for this, dedicated setter below

// sigmas and rhos
arma::mat& SSUR_Chain::getSigmaRho(){ return sigmaRho ; }
void SSUR_Chain::setSigmaRho( arma::mat&  externalSigmaRho )
{ 
    sigmaRho =  externalSigmaRho ; 
    logPSigmaRho();
    // THIS DO NOT ACT ON LOGLIK, SO IT'S IMPORTANT TO REMEMBER TO CHANGE THAT
}

void SSUR_Chain::setSigmaRho( arma::mat&  externalSigmaRho , double logP_sigmaRho_ )
{ 
    sigmaRho =  externalSigmaRho ; 
    logP_sigmaRho = logP_sigmaRho_ ;
    // THIS DO NOT ACT ON LOGLIK, SO IT'S IMPORTANT TO REMEMBER TO CHANGE THAT
}

double SSUR_Chain::getNu() const{ return nu; } 
void SSUR_Chain::setNu( double nu_ )
{ 
    nu = nu_ ; 
    logPSigmaRho();
} 

double SSUR_Chain::getLogPSigmaRho(){ return logP_sigmaRho ; }
// no setter for this, dedicated setter below

// o_k
arma::vec& SSUR_Chain::getO(){ return o ; }
void SSUR_Chain::setO( arma::vec& o_ )
{
    o = o_ ;
    logPO();
}

void SSUR_Chain::setO( arma::vec& o_ , double logP_o_ )
{
    o = o_ ;
    logP_o = logP_o_ ;
}

double SSUR_Chain::getOA() const{ return a_o ; }
void SSUR_Chain::setOA( double a_o_ )
{ 
    a_o = a_o_ ;
    logPO();
}

double SSUR_Chain::getOB() const{ return b_o ; }
void SSUR_Chain::setOB( double b_o_ )
{
    b_o = b_o_ ;
    logPO();
}

double SSUR_Chain::getVarOProposal() const{ return var_o_proposal ; }
void SSUR_Chain::setVarOProposal( double var_o_proposal_ ){ var_o_proposal = var_o_proposal_ ; }

double SSUR_Chain::getOAccRate() const{ return o_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double SSUR_Chain::getLogPO() const{ return logP_o ; }
// no setter for this, dedicated setter below

// pi_j
arma::vec& SSUR_Chain::getPi(){ return pi ; }
void SSUR_Chain::setPi( arma::vec& pi_ )
{
    pi = pi_ ;
    logPPi();
}

void SSUR_Chain::setPi( arma::vec& pi_ , double logP_pi_ )
{
    pi = pi_ ;
    logP_pi = logP_pi_ ;
}

double SSUR_Chain::getPiA() const{ return a_pi ; }
void SSUR_Chain::setPiA( double a_pi_ )
{
    a_pi = a_pi_ ; 
    logPPi();
}

double SSUR_Chain::getPiB() const{ return b_pi ; }
void SSUR_Chain::setPiB( double b_pi_ )
{ 
    b_pi = b_pi_ ;
    logPPi();
}

double SSUR_Chain::getVarPiProposal() const{ return var_pi_proposal ; }
void SSUR_Chain::setVarPiProposal( double var_pi_proposal_ ){ var_pi_proposal = var_pi_proposal_ ; }

double SSUR_Chain::getPiAccRate() const{ return pi_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double SSUR_Chain::getLogPPi() const{ return logP_pi ; }
// no setter for this, dedicated setter below

// GAMMA
arma::umat& SSUR_Chain::getGamma(){ return gamma ; }
void SSUR_Chain::setGamma( arma::umat& externalGamma )
{
    gamma = externalGamma ; 
    logPGamma();
}

void SSUR_Chain::setGamma( arma::umat& externalGamma , double logP_gamma_ )
{
    gamma = externalGamma ; 
    logP_gamma = logP_gamma_ ;
}

unsigned int SSUR_Chain::getNUpdatesMC3() const{ return n_updates_MC3 ; }
void SSUR_Chain::setNUpdatesMC3( unsigned int n_updates_MC3_ ){ n_updates_MC3 = n_updates_MC3_ ; }

double SSUR_Chain::getGammaAccRate() const{ return gamma_acc_count/(double)internalIterationCounter ; }
// no setter for this, is updated internally

double SSUR_Chain::getLogPGamma() const{ return logP_gamma ; }
// no setter for this, dedicated setter below

// W
double SSUR_Chain::getW() const{ return w ; }
void SSUR_Chain::setW( double w_ )
{
    w = w_ ;
    logPW();
}

void SSUR_Chain::setW( double w_ , double logP_w_ )
{
    w = w_ ;
    logP_w = logP_w_ ;
}

double SSUR_Chain::getWA() const{ return a_w ; }
void SSUR_Chain::setWA( double a_w_ )
{ 
    a_w = a_w_ ;
    logPW();
}

double SSUR_Chain::getWB() const{ return b_w ; }
void SSUR_Chain::setWB( double b_w_ )
{ 
    b_w = b_w_ ;
    logPW();
}

double SSUR_Chain::getLogPW() const{ return logP_w; }
// no setter for this, dedicated setter below

// BETA
arma::mat& SSUR_Chain::getBeta(){ return beta ; }
void SSUR_Chain::setBeta( arma::mat&  externalBeta )
{
    beta =  externalBeta ;
    logPBeta();
    // THIS DO NOT ACT ON LOGLIK, SO IT'S IMPORTANT TO REMEMBER TO CHANGE THAT
}

void SSUR_Chain::setBeta( arma::mat&  externalBeta , double logP_beta_ )
{
    beta =  externalBeta ;
    logP_beta = logP_beta_ ;
    // THIS DO NOT ACT ON LOGLIK, SO IT'S IMPORTANT TO REMEMBER TO CHANGE THAT
}

double SSUR_Chain::getLogPBeta() const{ return logP_beta ; }
// no setter for this, dedicated setter below

// LOG-LIKELIHOOD FOR THE SSUR MODEL
double SSUR_Chain::getLogLikelihood() const{ return log_likelihood ; }
void SSUR_Chain::setLogLikelihood( double log_likelihood_ ){ log_likelihood = log_likelihood_ ; }

double SSUR_Chain::getJointLogPrior() const
{
    return logP_tau + logP_eta + logP_jt + logP_sigmaRho + logP_o + logP_pi + logP_gamma + logP_w + logP_beta;
}

double SSUR_Chain::getJointLogPosterior() const
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
void SSUR_Chain::tauInit( double tau_init , double a_tau_ , double b_tau_ , double var_tau_proposal_ )
{
    tau = tau_init ;

    a_tau = a_tau_; 
    b_tau = b_tau_;

    var_tau_proposal = var_tau_proposal_;

    tau_acc_count = 0.;
    
    logPTau();
}

void SSUR_Chain::tauInit()
{
    tauInit( 1. , 0.1 , 10. , 0.02 );
}

void SSUR_Chain::tauInit( double tau_init )
{
    tauInit( tau_init , 0.1 , 10. , 0.02 );
}


void SSUR_Chain::etaInit( double eta_init , double a_eta_ , double b_eta_ )
{
    eta = eta_init;

    a_eta = a_eta_;
    b_eta = b_eta_;

    logPEta();
}

void SSUR_Chain::etaInit()
{
    etaInit( 0.1 , /*0.5*s*(s-1)*0.1*/ 1. , /*0.5*s*(s-1)*(1.-0.1)*/ 1. ); // 0.1  is essentially the expected active proportion
}

void SSUR_Chain::etaInit( double eta_init )
{
    etaInit( eta_init , /*0.5*s*(s-1)*0.1*/ 1. , /*0.5*s*(s-1)*(1.-0.1)*/ 1. ); // should I use eta_init in the prior parameters?
}


// **** care, due to the need to pass a whole object, for constuctor purposes
// this one function differently than the rest (I need to repeat the code in both function essentally..)
void SSUR_Chain::jtInit( JunctionTree& jt_init )
{
    jt = jt_init ;

    jt_acc_count = 0.;
    n_updates_jt = 5; // default value, should I pick something different?

    logPJT();
}

void SSUR_Chain::jtInit()
{
    jt = JunctionTree( s , "empty" ); //empty constructor, diagonal adj matrix of dimension s
    // jt = JunctionTree( s , "full" ); // full constructor, dense adj matrix of dimension s
 
    jt_acc_count = 0.;
    n_updates_jt = 5; // default value, should I pick something different?

    logPJT();
}
// end different from the rest


void SSUR_Chain::sigmaRhoInit( arma::mat& sigmaRho_init , double nu_ )
{
    sigmaRho = sigmaRho_init;
    nu = nu_;
    logPSigmaRho();
}

void SSUR_Chain::sigmaRhoInit()
{
    arma::mat init = arma::eye<arma::mat>(s,s);
    sigmaRhoInit( init , s+2. );
}

void SSUR_Chain::sigmaRhoInit( arma::mat& sigmaRho_init )
{
    sigmaRhoInit( sigmaRho_init , s+2. );
}

void SSUR_Chain::oInit( arma::vec& o_init , double a_o_ , double b_o_ , double var_o_proposal_ )
{
    o = o_init;
    a_o = a_o_;
    b_o = b_o_;
    var_o_proposal = var_o_proposal_;
    o_acc_count = 0.;
    
    logPO();
}

void SSUR_Chain::oInit()
{
    arma::vec init = arma::ones<arma::vec>(s) / std::max( 500. , (double)p ) ;
    oInit( init , 2. , std::max(500.,(double)p)-2. , 0.005 );
}

void SSUR_Chain::oInit( arma::vec& o_init )
{
    oInit( o_init , 2. , std::max(500.,(double)p)-2. , 0.005 );
}

void SSUR_Chain::piInit( arma::vec& pi_init , double a_pi_ , double b_pi_ , double var_pi_proposal_ )
{
    pi = pi_init;
    a_pi = a_pi_;
    b_pi = b_pi_;
    var_pi_proposal = var_pi_proposal_;
    pi_acc_count = 0.;
    
    logPPi();
}

void SSUR_Chain::piInit()
{
    arma::vec init = arma::ones<arma::vec>(p) ;
    piInit( init , 2. , 1. , 0.02 );
}

void SSUR_Chain::piInit( arma::vec& pi_init )
{
    piInit( pi_init , 2. , 1. , 0.02 );
}

void SSUR_Chain::gammaInit( arma::umat& gamma_init )
{
    gamma = gamma_init; 
    gamma_acc_count = 0.;
    logPGamma();
    updateGammaMask();
}

void SSUR_Chain::gammaInit()
{
    arma::umat init = arma::zeros<arma::umat>(p,s);
    gammaInit( init );
}

void SSUR_Chain::wInit( double w_init , double a_w_ , double b_w_ )
{
    w = w_init;
    a_w = a_w_;
    b_w = b_w_;
    
    logPW();
}

void SSUR_Chain::wInit( double w_init )
{
    wInit( w_init , 2. , 5. );
}

void SSUR_Chain::wInit()
{
    wInit( 1. , 2. , 5. );
}

void SSUR_Chain::betaInit( arma::mat& beta_init )
{
    beta = beta_init;
    logPBeta();
}

void SSUR_Chain::betaInit()
{
    arma::mat init = arma::zeros<arma::mat>(p+1,s);
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
double SSUR_Chain::logPTau( double tau_ , double a_tau_ , double b_tau_ )
{
    return Distributions::logPDFGamma( tau_ , a_tau_, b_tau_ );
}

double SSUR_Chain::logPTau( )
{
    logP_tau = logPTau( tau , a_tau, b_tau );
    return logP_tau;
}

double SSUR_Chain::logPTau( double tau_ )
{
    return logPTau( tau_ , a_tau, b_tau );
}

// ETA
double SSUR_Chain::logPEta( double eta_ , double a_eta_ , double b_eta_ )
{
    return Distributions::logPDFBeta( eta_ , a_eta_ , b_eta_ );
}

double SSUR_Chain::logPEta( )
{
    logP_eta = logPEta( eta , a_eta, b_eta );
    return logP_eta;
}

double SSUR_Chain::logPEta( double eta_ )
{
     return logPEta( eta_ , a_eta, b_eta );
}

// JT
double SSUR_Chain::logPJT( const JunctionTree& externalJT , double eta_ )
{
    double logP = 0.;
    for(unsigned int k=0; k<(s-1); ++k)
    {
        for(unsigned int l=k+1; l<s; ++l)
        {
            logP += Distributions::logPDFBernoulli( externalJT.adjacencyMatrix(k,l) , eta_ );
        }
    }

    return logP;
}

double SSUR_Chain::logPJT( )
{
    logP_jt = logPJT( jt , eta );
    return logP_jt;
}

double SSUR_Chain::logPJT( const JunctionTree& externalJT )
{
    return logPJT( externalJT , eta );
}


// sigma + rhos
double SSUR_Chain::logPSigmaRho( const arma::mat&  externalSigmaRho , double nu_ , double tau_ , const JunctionTree& externalJT )
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
            a = 0.5 * ( nu_ - s + nConditioninIndexes + 1. ) ;
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

double SSUR_Chain::logPSigmaRho( )
{
    logP_sigmaRho = logPSigmaRho( sigmaRho , nu , tau , jt );
    return logP_sigmaRho;
}

double SSUR_Chain::logPSigmaRho( const arma::mat&  externalSigmaRho )
{
    return logPSigmaRho(  externalSigmaRho , nu , tau , jt );
}

// o_k
double SSUR_Chain::logPO( const arma::vec& o_ , double a_o_ , double b_o_ )
{
    double logP = 0.;
    for(unsigned int k=0; k<s; ++k)
        logP += Distributions::logPDFBeta( o_(k) , a_o_, b_o_ );

    return logP;
}

double SSUR_Chain::logPO( )
{
    logP_o = logPO( o , a_o , b_o );
    return logP_o;
}

double SSUR_Chain::logPO( const arma::vec& o_ )
{
    return logPO( o_ , a_o , b_o );
}


// pi_j
double SSUR_Chain::logPPi( arma::vec& pi_ , double a_pi_ , double b_pi_ )
{
    double logP = 0.;
    for(unsigned int j=0; j<p; ++j)
        logP += Distributions::logPDFGamma( pi_(j) , a_pi_, b_pi_ );

    return logP;    
}

double SSUR_Chain::logPPi( )
{
    logP_pi = logPPi( pi , a_pi , b_pi );
    return logP_pi;
}
double SSUR_Chain::logPPi( arma::vec& pi_ )
{
    return logPPi( pi_ , a_pi , b_pi );
}

// GAMMA
double SSUR_Chain::logPGamma( const arma::umat& externalGamma , const arma::vec& o_ , const arma::vec& pi_ )
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

double SSUR_Chain::logPGamma( )
{
    logP_gamma = logPGamma( gamma , o , pi );
    return logP_gamma;
}

double SSUR_Chain::logPGamma( const arma::umat& externalGamma )
{
    return logPGamma( externalGamma , o , pi );
}

// W
double SSUR_Chain::logPW( double w_ , double a_w_ , double b_w_ )
{
    return Distributions::logPDFIGamma( w_ , a_w_, b_w_ );
}

double SSUR_Chain::logPW( )
{
    logP_w = logPW( w , a_w , b_w );
    return logP_w;
}

double SSUR_Chain::logPW( double w_ )
{
    return logPW( w_ , a_w , b_w );
}

// BETA
double SSUR_Chain::logPBeta( const arma::mat&  externalBeta , const arma::umat& externalGamma , double w_ )
{
    arma::uvec VS_IN_j;
    arma::uvec singleIdx_j(1);

    arma::umat mask = createGammaMask( externalGamma );

    double logP = 0.;

    for(unsigned int j=0; j<s ; ++j)
    {
        singleIdx_j(0) = j;
        VS_IN_j = mask( arma::find( mask.col(1) == j) , arma::zeros<arma::uvec>(1) );

        logP += Distributions::logPDFNormal(  externalBeta(VS_IN_j,singleIdx_j) ,
                    w_ * arma::eye(VS_IN_j.n_elem,VS_IN_j.n_elem) );

    }

    return logP;
}

double SSUR_Chain::logPBetaMask( const arma::mat&  externalBeta , const arma::umat& mask_ , double w_ )
{
    arma::uvec VS_IN_j;
    arma::uvec singleIdx_j(1);

    double logP = 0.;

    for(unsigned int j=0; j<s ; ++j)
    {
        singleIdx_j(0) = j;
        VS_IN_j = mask_( arma::find( mask_.col(1) == j) , arma::zeros<arma::uvec>(1) );

        logP += Distributions::logPDFNormal(  externalBeta(VS_IN_j,singleIdx_j) ,
                    w_ * arma::eye(VS_IN_j.n_elem,VS_IN_j.n_elem) );

    }

    return logP;
}

double SSUR_Chain::logPBeta( )
{
    logP_beta = logPBetaMask( beta , gammaMask , w );
    return logP_beta;
}

double SSUR_Chain::logPBeta( const arma::mat&  externalBeta )
{
    return logPBetaMask(  externalBeta , gammaMask , w );
}

// LOG LIKELIHOODS
double SSUR_Chain::logLikelihood( )
{
    double logP = 0.;
    #pragma omp parallel for default(shared) reduction(+:logP)
    for( unsigned int k=0; k<s; ++k)
    {
        logP += Distributions::logPDFNormal( (*Y).col(k) , (XB.col(k)+rhoU.col(k)) , sigmaRho(k,k));
    }

    logP /= temperature;

    log_likelihood = logP; // update internal state

    return logP;
}

double SSUR_Chain::logLikelihood( const arma::umat&  externalGammaMask , const arma::mat& externalXB ,
                                  const arma::mat& externalU , const arma::mat& externalRhoU , const arma::mat&  externalSigmaRho )
{
        double logP = 0.;
		#pragma omp parallel for default(shared) reduction(+:logP)
		for( unsigned int k=0; k<s; ++k)
		{
			logP += Distributions::logPDFNormal( (*Y).col(k) , (externalXB.col(k) + externalRhoU.col(k)) ,  externalSigmaRho(k,k));
		}

		return logP/temperature;
}                                  

double SSUR_Chain::logLikelihood( arma::umat&  externalGammaMask , arma::mat& externalXB , arma::mat& externalU , arma::mat& externalRhoU , //gammaMask,XB,U,rhoU
                        const arma::mat&  externalBeta , const arma::umat& externalGamma , // beta , gamma 
                        const arma::mat&  externalSigmaRho , const JunctionTree& externalJT ) // sigmaRho, jt
{
    externalGammaMask = createGammaMask(externalGamma);
    arma::uvec singleIdx_k(1), VS_IN_k;

    createQuantities(  externalGammaMask , externalXB , externalU , externalRhoU ,
                      externalGamma ,  externalBeta ,  externalSigmaRho , externalJT );

    double logP = 0.;
    #pragma omp parallel for default(shared) reduction(+:logP)
    for( unsigned int k=0; k<s; ++k)
    {
        logP += Distributions::logPDFNormal( (*Y).col(k) , ( externalXB.col(k) + externalRhoU.col(k)) ,  externalSigmaRho(k,k));
    }

    return logP/temperature;
}


// *********************
// STEP FUNCTION - PERFORM ONE ITERATION FOR THE CHAIN
// *********************


// This function sample sigmas and rhos from their full conditionals and updates the relevant matrix rhoU to reflect thats
double SSUR_Chain::sampleSigmaRhoGivenBeta( const arma::mat&  externalBeta , arma::mat& mutantSigmaRho , const JunctionTree& externalJT ,
                const arma::umat&  externalGammaMask , const arma::mat& externalXB , const arma::mat& externalU , arma::mat& mutantRhoU )
{
    double logP = 0.;

    mutantSigmaRho.zeros(s,s); // RESET THE WHOLE MATRIX !!!

    // hyperparameter of the posterior sampler
    arma::mat Sigma = ( externalU.t() * externalU ) / temperature; Sigma.diag() += tau;

    double thisSigmaTT;
    arma::uvec connectedNodes;
    arma::uvec conditioninIndexes;
    unsigned int nConditioninIndexes;
    arma::uvec singleIdx_l(1); // needed for convention with arma::submat

    double a,b;
    arma::mat rhoVar; // inverse matrix of the residual elements of Sigma in the component 
    arma::rowvec rhoMean; // this is the partial Schur complement, needed for the sampler 

    //bool test;

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
            a = 0.5 * ( n/temperature + nu - s + nConditioninIndexes + 1. ) ;
            b = 0.5 * thisSigmaTT ;

            mutantSigmaRho(l,l) = Distributions::randIGamma( a , b );

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

    // modify useful quantities, only rhoU impacted
    //recompute crossU as the rhos have changed
    mutantRhoU = createRhoU( externalU , mutantSigmaRho , externalJT );

    return logP;
}

double SSUR_Chain::sampleBetaGivenSigmaRho( arma::mat& mutantBeta , const arma::mat&  externalSigmaRho , const JunctionTree& externalJT ,
                const arma::umat&  externalGammaMask , arma::mat& mutantXB , arma::mat& mutantU , arma::mat& mutantRhoU )
{
    double logP = 0.;

    arma::vec mu_k; arma::mat W_k; // beta samplers

    arma::uvec singleIdx_k(1); // needed for convention with arma::submat

    arma::vec tmpVec;

    arma::uvec VS_IN_k;
    //bool test;

    arma::uvec xi = arma::conv_to<arma::uvec>::from(externalJT.perfectEliminationOrder);
    arma::vec xtxMultiplier(s);
    arma::mat y_tilde = (*Y) - mutantRhoU ;

    // prepare posterior full conditional's hyperparameters
    y_tilde.each_row() /= ( externalSigmaRho.diag().t()) ; // divide each col by the corresponding element of sigma
    xtxMultiplier(xi(s-1)) = 0;
    // y_tilde.col(xi(s-1)) is already ok;

    for( unsigned int k=0; k < (s-1); ++k)
    {
        xtxMultiplier(xi(k)) = 0;
        for(unsigned int l=k+1 ; l<s ; ++l)
        {
            xtxMultiplier(xi(k)) += pow( externalSigmaRho(xi(l),xi(k)),2) /  externalSigmaRho(xi(l),xi(l));
            y_tilde.col(xi(k)) -= (  externalSigmaRho(xi(l),xi(k)) /  externalSigmaRho(xi(l),xi(l)) ) * 
                ( mutantU.col(xi(l)) - mutantRhoU.col(xi(l)) +  externalSigmaRho(xi(l),xi(k)) * ( mutantU.col(xi(k)) - (*Y).col(xi(k)) ) );
        }

    }

    // actual sampling
    tmpVec.clear();

    // for( unsigned int j : externalJT.perfectEliminationOrder ) //shouldn't make a difference..
    for(unsigned int k=0; k<s ; ++k)
    {

        VS_IN_k =  externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
        singleIdx_k(0) = k;

        // /*test = */arma::inv_sympd( W_k ,  ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        /*test = */arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        mu_k = W_k * ( (*X).cols(VS_IN_k).t() * y_tilde.col(k) / temperature ) ;

        tmpVec = Distributions::randMvNormal( mu_k , W_k );
        logP += Distributions::logPDFNormal( tmpVec , mu_k , W_k );

        // to be sure
        mutantBeta.col(k).fill( 0. );

        for(unsigned int j=0 ; j<VS_IN_k.n_elem ; ++j)
            mutantBeta(VS_IN_k(j),k) = tmpVec(j);

    }

    // Now the beta have changed so X*B is changed as well as U, compute it to update it for the logLikelihood
    // finally as U changed, rhoU changes as well
    mutantXB = createXB(  externalGammaMask , mutantBeta );
    mutantU = createU( mutantXB );
    mutantRhoU = createRhoU( mutantU ,  externalSigmaRho , externalJT );

    return logP;

}   

double SSUR_Chain::sampleBetaKGivenSigmaRho( const unsigned int k , arma::mat& mutantBeta , const arma::mat&  externalSigmaRho , const JunctionTree& externalJT ,
                const arma::umat&  externalGammaMask , arma::mat& mutantXB , arma::mat& mutantU , arma::mat& mutantRhoU )
{
    double logP;

    arma::vec mu_k; arma::mat W_k; // beta samplers

    arma::uvec singleIdx_k(1); // needed for convention with arma::submat

    arma::vec tmpVec;

    arma::uvec VS_IN_k;
    //bool test;

    arma::uvec xi = arma::conv_to<arma::uvec>::from(externalJT.perfectEliminationOrder);
    double xtxMultiplier;
    arma::vec y_tilde = (*Y).col(k) - mutantRhoU.col(k) ;

    // prepare posterior full conditional's hyperparameters
    xtxMultiplier = 0;
    y_tilde /=  externalSigmaRho(k,k);
    // y_tilde.col(xi(s-1)) is already ok;
    // xtxMultiplier(xi(s-1)) = 0; is already ok

    unsigned int k_idx = arma::as_scalar( arma::find( xi == k , 1 ) );

    for(unsigned int l=k_idx+1 ; l<s ; ++l)
    {
        xtxMultiplier += pow( externalSigmaRho(xi(l),k),2) /  externalSigmaRho(xi(l),xi(l));
        y_tilde -= (  externalSigmaRho(xi(l),k) /  externalSigmaRho(xi(l),xi(l)) ) * 
            ( mutantU.col(xi(l)) - mutantRhoU.col(xi(l)) +  externalSigmaRho(xi(l),k) * ( mutantU.col(k) - (*Y).col(k) ) );
    }

    // actual sampling
    tmpVec.clear();

    VS_IN_k =  externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
    singleIdx_k(0) = k;

    // /*test = */arma::inv_sympd( W_k ,  ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
    /*test = */arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
    mu_k = W_k * ( (*X).cols(VS_IN_k).t() * y_tilde / temperature ) ;

    tmpVec = Distributions::randMvNormal( mu_k , W_k );
    logP = Distributions::logPDFNormal( tmpVec , mu_k , W_k );

    // to be sure
    mutantBeta.col(k).fill( 0. );

    for(unsigned int j=0 ; j<VS_IN_k.n_elem ; ++j)
        mutantBeta(VS_IN_k(j),k) = tmpVec(j);

    // Now the beta have changed so X*B is changed as well as U, compute it to update it for the logLikelihood
    // finally as U changed, rhoU changes as well
    mutantXB = createXB(  externalGammaMask , mutantBeta );
    mutantU = createU( mutantXB );
    mutantRhoU = createRhoU( mutantU ,  externalSigmaRho , externalJT );

    return logP;

}  

//logProbabilities of the above samplers (for the reverse moves)
// this function "simulate" a gibbs move and compute its proposal probability
double SSUR_Chain::logPSigmaRhoGivenBeta( const arma::mat&  externalBeta , const arma::mat& mutantSigmaRho , const JunctionTree& externalJT ,
                                            const arma::umat&  externalGammaMask , const arma::mat& externalXB , const arma::mat& externalU , const arma::mat& mutantRhoU )
{
    double logP = 0.;

    // hyperparameter of the posterior sampler
    arma::mat Sigma = ( externalU.t() * externalU ) / temperature; Sigma.diag() += tau;

    double thisSigmaTT;
    arma::uvec connectedNodes;
    arma::uvec conditioninIndexes;
    unsigned int nConditioninIndexes;
    arma::uvec singleIdx_l(1); // needed for convention with arma::submat

    double a,b;
    arma::mat rhoVar; // inverse matrix of the residual elements of Sigma in the component 
    arma::rowvec rhoMean; // this is the partial Schur complement, needed for the sampler 

    //bool test;

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
            a = 0.5 * ( n/temperature + nu - s + nConditioninIndexes + 1. ) ;
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

    return logP;
}

double SSUR_Chain::logPBetaGivenSigmaRho( const arma::mat& mutantBeta , const arma::mat&  externalSigmaRho , const JunctionTree& externalJT ,
                const arma::umat& externalGammaMask , const arma::mat& mutantXB , const arma::mat& mutantU , const arma::mat& mutantRhoU )
{
    double logP = 0.;

    arma::vec mu_k; arma::mat W_k; // beta samplers

    arma::uvec singleIdx_k(1); // needed for convention with arma::submat

    arma::uvec VS_IN_k;
    //bool test;

    arma::uvec xi = arma::conv_to<arma::uvec>::from( externalJT.perfectEliminationOrder );
    arma::vec xtxMultiplier(s);
    arma::mat y_tilde = (*Y) - mutantRhoU ;

    // prepare posterior full conditional's hyperparameters
    y_tilde.each_row() /= ( externalSigmaRho.diag().t()) ; // divide each col by the corresponding element of sigma
    xtxMultiplier(xi(s-1)) = 0;
    // y_tilde.col(xi(s-1)) is already ok;

    for( unsigned int k=0; k < (s-1); ++k)
    {
        xtxMultiplier(xi(k)) = 0;
        for(unsigned int l=k+1 ; l<s ; ++l)
        {
            xtxMultiplier(xi(k)) += pow( externalSigmaRho(xi(l),xi(k)),2) /  externalSigmaRho(xi(l),xi(l));
            y_tilde.col(xi(k)) -= (  externalSigmaRho(xi(l),xi(k)) /  externalSigmaRho(xi(l),xi(l)) ) * 
                ( U.col(xi(l)) - mutantRhoU.col(xi(l)) +  externalSigmaRho(xi(l),xi(k)) * ( mutantU.col(xi(k)) - (*Y).col(xi(k)) ) );
        }

    }

    // for( unsigned int j : externalJT.perfectEliminationOrder ) //shouldn't make a difference..
    for(unsigned int k=0; k<s ; ++k)
    {

        VS_IN_k = externalGammaMask( arma::find( externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
        singleIdx_k(0) = k;

        // /*test = */arma::inv_sympd( W_k ,  ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        /*test = */arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
        mu_k = W_k * ( (*X).cols(VS_IN_k).t() * y_tilde.col(k) / temperature ) ;

        logP += Distributions::logPDFNormal( mutantBeta(VS_IN_k,singleIdx_k) , mu_k , W_k );
    }

    return logP;
}   

double SSUR_Chain::logPBetaKGivenSigmaRho( const unsigned int k , const arma::mat& mutantBeta , const arma::mat&  externalSigmaRho , const JunctionTree& externalJT ,
                const arma::umat&  externalGammaMask , const arma::mat& mutantXB , const arma::mat& mutantU , const arma::mat& mutantRhoU )
{
    double logP;

    arma::vec mu_k; arma::mat W_k; // beta samplers

    arma::uvec singleIdx_k(1); // needed for convention with arma::submat

    arma::vec tmpVec;

    arma::uvec VS_IN_k;
    //bool test;

    arma::uvec xi = arma::conv_to<arma::uvec>::from(externalJT.perfectEliminationOrder);
    double xtxMultiplier;
    arma::vec y_tilde = (*Y).col(k) - mutantRhoU.col(k) ;

    // prepare posterior full conditional's hyperparameters
    xtxMultiplier = 0;
    y_tilde /=  externalSigmaRho(k,k);
    // y_tilde.col(xi(s-1)) is already ok;
    // xtxMultiplier(xi(s-1)) = 0; is already ok

    unsigned int k_idx = arma::as_scalar( arma::find( xi == k , 1 ) );

    for(unsigned int l=k_idx+1 ; l<s ; ++l)
    {
        xtxMultiplier += pow( externalSigmaRho(xi(l),k),2) /  externalSigmaRho(xi(l),xi(l));
        y_tilde -= (  externalSigmaRho(xi(l),k) /  externalSigmaRho(xi(l),xi(l)) ) * 
            ( mutantU.col(xi(l)) - mutantRhoU.col(xi(l)) +  externalSigmaRho(xi(l),k) * ( mutantU.col(k) - (*Y).col(k) ) );
    }

    VS_IN_k =  externalGammaMask( arma::find(  externalGammaMask.col(1) == k) , arma::zeros<arma::uvec>(1) );
    singleIdx_k(0) = k;

    // /*test = */arma::inv_sympd( W_k ,  ( (*X).cols(VS_IN_k).t() * (*X).cols(VS_IN_k) ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier(k) ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
    /*test = */arma::inv_sympd( W_k ,  ( XtX(VS_IN_k,VS_IN_k) / temperature ) * ( 1./ externalSigmaRho(k,k) + xtxMultiplier ) + (1./w)*arma::eye<arma::mat>(VS_IN_k.n_elem,VS_IN_k.n_elem) );
    mu_k = W_k * ( (*X).cols(VS_IN_k).t() * y_tilde / temperature ) ;

    logP = Distributions::logPDFNormal( mutantBeta(VS_IN_k,singleIdx_k) , mu_k , W_k );

    return logP;

}  


// sample sigmaRho given Beta or Beta given sigmaRho
// simple interface to gibbs sampling that updates the internal states (given the internal states) from the full conditional
void SSUR_Chain::sampleSigmaRhoGivenBeta()
{
    sampleSigmaRhoGivenBeta( beta , sigmaRho , jt , gammaMask , XB , U , rhoU );
}

void SSUR_Chain::sampleBetaGivenSigmaRho()
{
    sampleBetaGivenSigmaRho( beta , sigmaRho , jt , gammaMask , XB , U , rhoU );
}


// sampler for proposed updates on gamma
double SSUR_Chain::gammaBanditProposal( arma::umat& mutantGamma , arma::uvec& updateIdx )
{

    double logProposalRatio;

    // Sample Zs
    for(unsigned int i=0; i<p; ++i)
    {
        for(unsigned int j=0; j<s; ++j)
        {
            banditZeta(i,j) = Distributions::randBeta(banditAlpha(i,j),banditAlpha(i,j));
        }
    }

    // Create mismatch
    for(unsigned int i=0; i<(p*s); ++i)
    {
        mismatch(i) = (mutantGamma(i)==0)?(banditZeta(i)):(1.-banditZeta(i));   //mismatch
    }

    // Normalise
    // mismatch = arma::log(mismatch); //logscale ??? TODO
    // normalised_mismatch = mismatch - Utils::logspace_add(mismatch);

    normalised_mismatch = mismatch / arma::as_scalar(arma::sum(mismatch));

    if( Distributions::randU01() < 0.5 )   // one deterministic update
    {

        // Decide which to update
        updateIdx = arma::zeros<arma::uvec>(1);
        updateIdx(0) = Distributions::randWeightedIndexSampleWithoutReplacement(p*s,normalised_mismatch); // sample the one

        // Update
        mutantGamma(updateIdx(0)) = 1 - gamma(updateIdx(0)); // deterministic, just switch

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
        updateIdx = Distributions::randWeightedIndexSampleWithoutReplacement(p*s,normalised_mismatch,n_updates_bandit); // sample n_updates_bandit indexes

        normalised_mismatch_backwards = mismatch; // copy for backward proposal

        // Update
        for(unsigned int i=0; i<n_updates_bandit; ++i)
        {
            mutantGamma(updateIdx(i)) = Distributions::randBernoulli(banditZeta(updateIdx(i))); // random update

            normalised_mismatch_backwards(updateIdx(i)) = 1.- normalised_mismatch_backwards(updateIdx(i));

            logProposalRatio += Distributions::logPDFBernoulli(gamma(updateIdx(i)),banditZeta(updateIdx(i))) -
                Distributions::logPDFBernoulli(mutantGamma(updateIdx(i)),banditZeta(updateIdx(i)));
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
                
double SSUR_Chain::gammaMC3Proposal( arma::umat& mutantGamma , arma::uvec& updateIdx )
{
    updateIdx = arma::uvec(n_updates_MC3);

    for( unsigned int i=0; i<n_updates_MC3; ++i)
        updateIdx(i) = Distributions::randIntUniform(0,(p*s)-1);    // note that I might be updating multiple times the same coeff

    for( auto i : updateIdx)
    mutantGamma(i) = ( Distributions::randU01() < 0.5)? gamma(i) : 1-gamma(i); // could simply be ( 0.5 ? 1 : 0) ;

    return 0. ; // pass this to the outside, it's the (symmetric) logProposalRatio
}


// **************
// **** Methods that update the internal state of their parameter
// **************


// MH update, Normal in the log-space as tau is positive (with gamma prior)
void SSUR_Chain::stepTau()
{
    double proposedTau = std::exp( std::log(tau) + Distributions::randNormal(0.0, var_tau_proposal) );

    double proposedTauPrior = logPTau( proposedTau );
    double proposedSigmaRhoPrior = logPSigmaRho( sigmaRho, nu, proposedTau, jt);

    double logAccProb = (proposedTauPrior + proposedSigmaRhoPrior) - (logP_tau + logP_sigmaRho);

    if( Distributions::randLogU01() < logAccProb )
    {
        tau = proposedTau;
        logP_tau = proposedTauPrior;
        logP_sigmaRho = proposedSigmaRhoPrior;

        ++tau_acc_count;
    }
}

// Here we have a full conditional available
void SSUR_Chain::stepEta()
{
    double a = a_eta + 0.5*(arma::accu(jt.getAdjMat())/2. ) ; // divide by temperature if the prior on G is tempered
    double b = b_eta + ( (double)(s * (s-1.) * 0.5) - 0.5*(arma::accu(jt.getAdjMat())/2. ) ); // /temperature

    eta = Distributions::randBeta( a , b );

    logPEta(); // update its prior value
    logPJT(); // update JT's pror value as it's impacted by the new eta
}

void SSUR_Chain::stepJT()
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
    double val = Distributions::randU01();
    
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

        if( Distributions::randLogU01() < logAccProb )
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

// MH update (log-normal) -- update one value at each iteration TODO worth doing more?
void SSUR_Chain::stepOneO()
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

void SSUR_Chain::stepO()
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
void SSUR_Chain::stepOnePi()
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

void SSUR_Chain::stepPi()
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

// Gibbs sampler available here again for w given all the current betas and the gammas -- TODO keep an eye on this
void SSUR_Chain::stepW()
{
    double a = a_w + 0.5*( /*arma::accu(gamma) + intercept */ /*or*/ gammaMask.n_rows ); // divide by temperature if the prior on gamma is tempered
    double b = b_w + 0.5*( arma::accu( arma::square(arma::nonzeros(beta)) ) );   // all the beta_jk w/ gamma_jk=0 are 0 already // /temperature

    // std::cout << a_w << " -> " << a << "   ---   "<< b_w << " -> " << b << std::endl; 
    // std::cout << arma::nonzeros(beta).t() << std::endl; std::cin >> w;

    w = Distributions::randIGamma( a , b );

    logPW(); // update its prior value
    logPBeta(); // update beta's log prior as it's impacted by the change in w
}

void SSUR_Chain::stepGamma()
{
    arma::umat proposedGamma = gamma;
    arma::uvec updateIdx;

    double logProposalRatio = 0;

    // Update the proposed Gamma
    if( gammaSamplerType == "B" || gammaSamplerType == "bandit" || gammaSamplerType == "Bandit" || gammaSamplerType == "b" )
    {
        logProposalRatio += gammaBanditProposal( proposedGamma , updateIdx );

    }else if( gammaSamplerType == "MC3" || gammaSamplerType == "mc3" )
    {
        logProposalRatio += gammaMC3Proposal( proposedGamma , updateIdx );

    }else{
        logProposalRatio += gammaBanditProposal( proposedGamma , updateIdx ); // default
    }

    // given proposedGamma now, sample a new proposedBeta matrix and corresponging quantities
    arma::umat proposedGammaMask = createGammaMask( proposedGamma );
    arma::uvec updatedOutcomesIdx = arma::unique( arma::floor( updateIdx / p )); // every p I get to the new column, and as columns correspond to outcomes ... 
    
    // note for quantities below. The firt call to sampleXXX has the proposedQuantities set to the current value,
        // for them to be updated; the second call to logPXXX has them updated, needed for the backward probability
        // the main parameter of interest instead "changes to the current value" in the backward equation
    arma::mat proposedBeta = beta;

    arma::mat proposedXB = XB;
    arma::mat proposedU = U;
    arma::mat proposedRhoU = rhoU;

    for(unsigned int k : updatedOutcomesIdx)
    {
        logProposalRatio -= sampleBetaKGivenSigmaRho( k , proposedBeta , sigmaRho , jt ,
                                proposedGammaMask , proposedXB , proposedU , proposedRhoU );
        logProposalRatio += logPBetaKGivenSigmaRho( k , beta , sigmaRho , jt ,
                                gammaMask , proposedXB , proposedU , proposedRhoU );
    }

    // update log probabilities
    double proposedGammaPrior = logPGamma( proposedGamma );
    double proposedBetaPrior = logPBetaMask( proposedBeta , proposedGammaMask , w );
    double proposedLikelihood = logLikelihood( proposedGammaMask , proposedXB , proposedU , proposedRhoU , sigmaRho );
    
    double logAccProb = logProposalRatio +
                ( proposedGammaPrior + proposedBetaPrior + proposedLikelihood ) - 
                ( logP_gamma + logP_beta + log_likelihood );

    if( Distributions::randLogU01() < logAccProb )
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
        gamma_acc_count += 1. / updatedOutcomesIdx.n_elem;
    }

    // after A/R, update bandit Related variables
    if(!( gammaSamplerType == "MC3" || gammaSamplerType == "mc3" ))
    {
        for(arma::uvec::iterator iter = updateIdx.begin(); iter != updateIdx.end(); ++iter)
        {
            // FINITE UPDATE
            if( banditAlpha(*iter) + banditBeta(*iter) <= banditLimit )
            {
                banditAlpha(*iter) += banditIncrement * gamma(*iter);
                banditBeta(*iter) += banditIncrement * (1-gamma(*iter));
            }

            // // CONTINUOUS UPDATE
            // banditAlpha(*iter) += banditIncrement * gamma(*iter);
            // banditBeta(*iter) += banditIncrement * (1-gamma(*iter));

            // // then renormalise them
            // if( banditAlpha(*iter) + banditBeta(*iter) > banditLimit )
            // {
            //     banditAlpha(*iter) = banditLimit * ( banditAlpha(*iter) / ( banditAlpha(*iter) + banditBeta(*iter) ));
            //     banditBeta(*iter) = banditLimit * (1. - ( banditAlpha(*iter) / ( banditAlpha(*iter) + banditBeta(*iter) )) );
            // }
        }
    }
}

void SSUR_Chain::stepSigmaRhoAndBeta()
{
    sampleSigmaRhoGivenBeta();
    sampleBetaGivenSigmaRho();
    
    logPSigmaRho();
    logPBeta();
    logLikelihood();
}


// this updates all the internal states
void SSUR_Chain::step()
{

    // Update HyperParameters
    stepTau();
    stepEta();
    stepW();
    
    for( auto i=0; i<5; ++i){
        stepOneO();
        stepOnePi();
    }
    
    // Update JT
    if( internalIterationCounter >= jtStartIteration )
        stepJT();

    // update gamma
    stepGamma();

    // Update Sigmas, Rhos and Betas given all rest
    stepSigmaRhoAndBeta();

    // increase iteration counter
    ++ internalIterationCounter;

    // update the MH proposal variance
    updateProposalVariances();
}

// update all the internal proposal RW variances based on their acceptance rate
// note that for o and pi, even though they're vectors, the proposal is independent and common for all elements
// see Roberts  & Rosenthal 2008
void SSUR_Chain::updateProposalVariances()
{
    double delta, delta2;
    arma::vec deltaVec, delta2Vec;

    double adaptationFactor = 0.05;

    if( internalIterationCounter == 1 ) // init the mean and second moment
    {
        tauEmpiricalMean = std::log(tau);
        oEmpiricalMean = arma::log(o);
        piEmpiricalMean = arma::log(pi);

        tauEmpiricalM2 = 0.;
        oEmpiricalM2 = arma::zeros<arma::vec>(s);
        piEmpiricalM2 = arma::zeros<arma::vec>(p);

        var_tau_proposal_init = var_tau_proposal;
        var_o_proposal_init = var_o_proposal;
        var_pi_proposal_init = var_pi_proposal;

    }else if( internalIterationCounter > 1 ) 
    {
        // update running averages

        // tau
        delta = std::log(tau) - tauEmpiricalMean;
        tauEmpiricalMean = tauEmpiricalMean + ( delta / internalIterationCounter );
        delta2 = std::log(tau) - tauEmpiricalMean;
        tauEmpiricalM2 = tauEmpiricalM2 + delta * delta2;

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
        var_tau_proposal = adaptationFactor * var_tau_proposal_init + (1. - adaptationFactor) * (2.38*2.38) * tauEmpiricalM2/(internalIterationCounter-1);
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

void SSUR_Chain::swapTau( std::shared_ptr<SSUR_Chain>& that )
{
    double par = this->getTau();

    this->setTau( that->getTau() );
    that->setTau( par );
}

void SSUR_Chain::swapEta( std::shared_ptr<SSUR_Chain>& that )
{
    double par = this->getEta();

    this->setEta( that->getEta() );
    that->setEta( par );
}

void SSUR_Chain::swapJT( std::shared_ptr<SSUR_Chain>& that )
{
    JunctionTree par = this->getJT();

    this->setJT( that->getJT() );
}

void SSUR_Chain::swapSigmaRho( std::shared_ptr<SSUR_Chain>& that )
{
    arma::mat par = this->getSigmaRho();

    this->setSigmaRho( that->getSigmaRho() );
    that->setSigmaRho( par );
}

void SSUR_Chain::swapO( std::shared_ptr<SSUR_Chain>& that )
{
    arma::vec par = this->getO();

    this->setO( that->getO() );
    that->setO( par );
}

void SSUR_Chain::swapPi( std::shared_ptr<SSUR_Chain>& that )
{
    arma::vec par = this->getPi();

    this->setPi( that->getPi() );
    that->setPi( par );
}

void SSUR_Chain::swapGamma( std::shared_ptr<SSUR_Chain>& that )
{
    arma::umat par = this->getGamma();

    this->setGamma( that->getGamma() );
    that->setGamma( par );
}

void SSUR_Chain::swapW( std::shared_ptr<SSUR_Chain>& that )
{
    double par = this->getW();

    this->setW( that->getW() );
    that->setW( par );
}

void SSUR_Chain::swapBeta( std::shared_ptr<SSUR_Chain>& that )
{
    arma::mat par = this->getBeta();

    this->setBeta( that->getBeta() );
    that->setBeta( par );
}

int SSUR_Chain::exchangeGamma_step( std::shared_ptr<SSUR_Chain>& that )
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
 
    if( Distributions::randLogU01() < logPExchange )
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

int SSUR_Chain::exchangeJT_step( std::shared_ptr<SSUR_Chain>& that )
{
    // I'm exchanging the JT, sigmas and rhos.
    // So gammaMask, XB and U will stay the same and we will have to to re-compute rhoU for both chains
    arma::mat rhoU_1 = this->createRhoU( this->getU() , that->getSigmaRho() , that->getJT() );
    arma::mat rhoU_2 = that->createRhoU( that->getU() , this->getSigmaRho() , this->getJT() );

    double logLik_1 = this->logLikelihood( this->getGammaMask() , this->getXB() , this->getU() , rhoU_1 , that->getSigmaRho() ); // note that this and that lik are
    double logLik_2 = that->logLikelihood( that->getGammaMask() , that->getXB() , that->getU() , rhoU_2 , this->getSigmaRho() );   // important because of temperature

    double logPExchange = ( logLik_1 + logLik_2 ) -
                        ( this->getLogLikelihood() + that->getLogLikelihood() );
 
    if( Distributions::randLogU01() < logPExchange )
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

int SSUR_Chain::adapt_crossOver_step( std::shared_ptr<SSUR_Chain>& that )
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

    double logPBetaFirst = this->logPBetaMask( betaXO[0] , gammaMask_XO[0] , this->getW() );  // this and that are important for temperature and hyperparameters
    double logPBetaSecond = that->logPBetaMask( betaXO[1] , gammaMask_XO[1] , that->getW() );
    
    double logPGammaFirst = this->logPGamma( gammaXO[0] );
    double logPGammaSecond = that->logPGamma( gammaXO[1] );

    pCrossOver +=   ( logLikFirst + logPBetaFirst + logPGammaFirst +
                        logLikSecond + logPBetaSecond + logPGammaSecond ) -
                    ( this->getLogLikelihood() + this->getLogPBeta() + this->getLogPGamma() +
                        that->getLogLikelihood() + that->getLogPBeta() + that->getLogPGamma() );
    
    if( Distributions::randLogU01() < pCrossOver )
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

int SSUR_Chain::uniform_crossOver_step( std::shared_ptr<SSUR_Chain>& that )
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

    double logPBetaFirst = this->logPBetaMask( betaXO[0] , gammaMask_XO[0] , this->getW() );  // this and that are important for temperature and hyperparameters
    double logPBetaSecond = that->logPBetaMask( betaXO[1] , gammaMask_XO[1] , that->getW() );
    
    double logPGammaFirst = this->logPGamma( gammaXO[0] );
    double logPGammaSecond = that->logPGamma( gammaXO[1] );

    pCrossOver +=   ( logLikFirst + logPBetaFirst + logPGammaFirst +
                        logLikSecond + logPBetaSecond + logPGammaSecond ) -
                    ( this->getLogLikelihood() + this->getLogPBeta() + this->getLogPGamma() +
                        that->getLogLikelihood() + that->getLogPBeta() + that->getLogPGamma() );
    
    if( Distributions::randLogU01() < pCrossOver )
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

int SSUR_Chain::block_crossOver_step( std::shared_ptr<SSUR_Chain>& that , arma::mat& corrMatX , double threshold )
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

    double logPBetaFirst = this->logPBetaMask( betaXO[0] , gammaMask_XO[0] , this->getW() );  // this and that are important for temperature and hyperparameters
    double logPBetaSecond = that->logPBetaMask( betaXO[1] , gammaMask_XO[1] , that->getW() );
    
    double logPGammaFirst = this->logPGamma( gammaXO[0] );
    double logPGammaSecond = that->logPGamma( gammaXO[1] );

    pCrossOver +=   ( logLikFirst + logPBetaFirst + logPGammaFirst +
                        logLikSecond + logPBetaSecond + logPGammaSecond ) -
                    ( this->getLogLikelihood() + this->getLogPBeta() + this->getLogPGamma() +
                        that->getLogLikelihood() + that->getLogPBeta() + that->getLogPGamma() );
    
    if( Distributions::randLogU01() < pCrossOver )
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

void SSUR_Chain::swapAll( std::shared_ptr<SSUR_Chain>& thatChain )
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
    this->swapEta( thatChain );
    
    this->swapJT( thatChain );
    this->swapSigmaRho( thatChain );

    this->swapO( thatChain );
    this->swapPi( thatChain );
    this->swapGamma( thatChain );

    this->swapW( thatChain );
    this->swapBeta( thatChain );

    // recompute likelihood
    this->logLikelihood();
    thatChain->logLikelihood();

}

int SSUR_Chain::globalStep( std::shared_ptr<SSUR_Chain>& that )
{

    unsigned int globalType = Distributions::randIntUniform(0,5);

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
        
        case 4: 
            return this -> exchangeJT_step( that );
            break;
        
        case 5: 
            return this -> exchangeAll_step( that );
            break;
        
        default: 
            break;
    }

    return 0;
}

int SSUR_Chain::exchangeAll_step( std::shared_ptr<SSUR_Chain>& thatChain )
{

    double logPExchange = ( this->getLogLikelihood() * this->getTemperature() -
                            thatChain->getLogLikelihood() * thatChain->getTemperature() ) * 
					( 1. / thatChain->getTemperature() - 1. / this->getTemperature() );
                    //  no priors because that is not tempered so it cancels out

    if( Distributions::randLogU01() < logPExchange )
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
arma::umat SSUR_Chain::createGammaMask( const arma::umat& gamma )
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


void SSUR_Chain::updateGammaMask()
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

arma::mat SSUR_Chain::createXB( const arma::umat&  externalGammaMask , const arma::mat&  externalBeta )
{
    arma::uvec singleIdx_k(1), VS_IN_k;
    arma::mat externalXB(n,s);

    for(unsigned int k=0; k<s; ++k)
    {
        singleIdx_k(0) = k;
        VS_IN_k =  externalGammaMask( arma::find(  externalGammaMask.col(1) == k ) , arma::zeros<arma::uvec>(1) );
        externalXB.col(k) = ((*X).cols(VS_IN_k) *  externalBeta.submat(VS_IN_k,singleIdx_k) );
    }
    return externalXB;
}

void SSUR_Chain::updateXB()
{
    arma::uvec singleIdx_k(1), VS_IN_k;
    XB.set_size(n,s); // reset without initialising nor preserving data

    for(unsigned int k=0; k<s; ++k)
    {
        singleIdx_k(0) = k;
        VS_IN_k = gammaMask( arma::find( gammaMask.col(1) == k ) , arma::zeros<arma::uvec>(1) );
        XB.col(k) = ((*X).cols(VS_IN_k) * beta.submat(VS_IN_k,singleIdx_k) );
    }

}

arma::mat SSUR_Chain::createU( const arma::mat& externalXB )
{
    arma::mat externalU = (*Y) - externalXB;
    return externalU;
}

void SSUR_Chain::updateU()
{
    U = (*Y) - XB;
}

arma::mat SSUR_Chain::createRhoU( const arma::mat& externalU , const arma::mat&  externalSigmaRho , const JunctionTree& externalJT )
{
    arma::uvec xi = arma::conv_to<arma::uvec>::from(externalJT.perfectEliminationOrder);
    arma::mat externalRhoU = arma::zeros<arma::mat>(n,s);

    for( unsigned int k=1; k < s; ++k)
    {
        for(unsigned int l=0 ; l<k ; ++l)
        {
            if(  externalSigmaRho(xi(k),xi(l)) != 0 )
                externalRhoU.col(xi(k)) += externalU.col(xi(l)) *  externalSigmaRho(xi(k),xi(l));
        }
    }

    return externalRhoU;
}


void SSUR_Chain::updateRhoU()
{
    arma::uvec xi = arma::conv_to<arma::uvec>::from(jt.perfectEliminationOrder);
    rhoU.zeros(n,s);

    for( unsigned int k=1; k < s; ++k)
    {
        for(unsigned int l=0 ; l<k ; ++l)
        {
            if(  sigmaRho(xi(k),xi(l)) != 0 )
                rhoU.col(xi(k)) += U.col(xi(l)) * sigmaRho(xi(k),xi(l));
        }
    }
}

void SSUR_Chain::createQuantities( arma::umat&  externalGammaMask , arma::mat& externalXB , arma::mat& externalU , arma::mat& externalRhoU ,
                                    const arma::umat& externalGamma , const arma::mat&  externalBeta ,
                                    const arma::mat&  externalSigmaRho , const JunctionTree& externalJT )
{
     externalGammaMask = createGammaMask( externalGamma );
    externalXB = createXB(  externalGammaMask ,  externalBeta );
    externalU = createU( externalXB );
    externalRhoU = createRhoU( externalU ,  externalSigmaRho , externalJT );
}

void SSUR_Chain::updateQuantities()
{
    updateGammaMask();
    updateXB();
    updateU();
    updateRhoU();
}



// Bandit-sampling related methods
void SSUR_Chain::banditInit()// initialise all the private memebers
{

    banditZeta = arma::mat(p,s);

    banditAlpha = arma::mat(p,s);
    banditAlpha.fill( 0.5 );
    
    banditBeta = arma::mat(p,s);
    banditBeta.fill( 0.5 );

    mismatch = arma::vec(p*s);
    normalised_mismatch = arma::vec(p*s);
    normalised_mismatch_backwards = arma::vec(p*s);

    n_updates_bandit = 4; // this needs to be low as its O(n_updates!)

    banditLimit = (double)n;
    banditIncrement = 1.;

}

// MC3 init
void SSUR_Chain::MC3Init()
{
    n_updates_MC3 = std::ceil( p/10 ); //arbitrary nunmber, should I use something different?
}
