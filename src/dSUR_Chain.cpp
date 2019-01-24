#include "dSUR_Chain.h"

// *******************************
// Constructors
// *******************************

#include "dSUR_Chain.h"

// *******************************
// Constructors
// *******************************

dSUR_Chain::dSUR_Chain ( std::shared_ptr<arma::mat> data, unsigned int nObservations, 
            unsigned int nOutcomes, unsigned int nVSPredictors, unsigned int nFixedPredictors,
            std::shared_ptr<arma::uvec> outcomesIdx, std::shared_ptr<arma::uvec> VSPredictorsIdx,
            std::shared_ptr<arma::uvec> fixedPredictorIdx, std::shared_ptr<arma::umat> NAArrayIdx, std::shared_ptr<arma::uvec> completeCases, 
            std::string gammaSamplerType_ = "Bandit", bool usingGprior = false, double externalTemperature = 1. ):
    SSUR_Chain ( data, nObservations, nOutcomes, nVSPredictors, nFixedPredictors,
            outcomesIdx, VSPredictorsIdx, fixedPredictorIdx, NAArrayIdx, completeCases,
            gammaSamplerType_, usingGprior, externalTemperature )
    {

        // after the SSUR_Chain constructor has been executed run the following to "correct"
        jtInit(); 

        betaInit(beta);
        sigmaRhoInit(sigmaRho);

        updateQuantities();

        logLikelihood();
    }


dSUR_Chain::dSUR_Chain( Utils::SUR_Data& surData, std::string gammaSamplerType_ = "Bandit", bool usingGprior = false, double externalTemperature = 1. ):
    dSUR_Chain(surData.data,surData.nObservations,surData.nOutcomes,surData.nVSPredictors,surData.nFixedPredictors,
        surData.outcomesIdx,surData.VSPredictorsIdx,surData.fixedPredictorsIdx,surData.missingDataArrayIdx,surData.completeCases,
        gammaSamplerType_,usingGprior,externalTemperature){ }

dSUR_Chain::dSUR_Chain( Utils::SUR_Data& surData, double externalTemperature = 1. ):
    dSUR_Chain(surData.data,surData.nObservations,surData.nOutcomes,surData.nVSPredictors,surData.nFixedPredictors,
        surData.outcomesIdx,surData.VSPredictorsIdx,surData.fixedPredictorsIdx,surData.missingDataArrayIdx,surData.completeCases,
        "Bandit",false,externalTemperature){ }


// ******************************
// Init Methods
// ******************************

void dSUR_Chain::jtInit( JunctionTree& jt_init )
{
    jt = jt_init ;

    jt_acc_count = 0.;
    n_updates_jt = 5; // default value, should I pick something different?

    logPJT();
}

void dSUR_Chain::jtInit()
{
    jt = JunctionTree( nOutcomes , "full" ); // full constructor, dense adj matrix of dimension s
 
    jt_acc_count = 0.;
    n_updates_jt = 5; // default value, should I pick something different?

    logPJT();
}


// *********************
// STEP FUNCTION - PERFORM ONE ITERATION FOR THE CHAIN
// *********************


// This function sample sigmas and rhos from their full conditionals and updates the relevant matrix rhoU to reflect thats
double dSUR_Chain::sampleSigmaRhoGivenBeta( const arma::mat&  externalBeta , arma::mat& mutantSigmaRho ,
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
    arma::uvec singleIdx_k(1); // needed for convention with arma::submat

    double a,b;
    arma::mat rhoVar; // inverse matrix of the residual elements of Sigma in the component 
    arma::rowvec rhoMean; // this is the partial Schur complement, needed for the sampler 

    //bool test;

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

        mutantSigmaRho(k,k) = Distributions::randIGamma( a , b );

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

    } // end loop over elements inside components

    // modify useful quantities, only rhoU impacted
    //recompute crossU as the rhos have changed
    mutantRhoU = createRhoU( externalU , mutantSigmaRho );

    return logP;
}

//logProbabilities of the above samplers (for the reverse moves)
// this function "simulate" a gibbs move and compute its proposal probability
double dSUR_Chain::logPSigmaRhoGivenBeta( const arma::mat&  externalBeta , const arma::mat& mutantSigmaRho , 
                                            const arma::umat&  externalGammaMask , const arma::mat& externalXB , const arma::mat& externalU , const arma::mat& mutantRhoU )
{
    double logP = 0.;

    // hyperparameter of the posterior sampler
    arma::mat Sigma = ( externalU.t() * externalU ) / temperature; Sigma.diag() += tau;

    double thisSigmaTT;
    arma::uvec connectedNodes;
    arma::uvec conditioninIndexes;
    unsigned int nConditioninIndexes;
    arma::uvec singleIdx_k(1); // needed for convention with arma::submat

    double a,b;
    arma::mat rhoVar; // inverse matrix of the residual elements of Sigma in the component 
    arma::rowvec rhoMean; // this is the partial Schur complement, needed for the sampler 

    //bool test;


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

    return logP;
}

// this updates all the internal states
void dSUR_Chain::step()
{

    // Update HyperParameters
    stepTau();
    stepW();
    
    for( auto i=0; i<5; ++i){
        stepOneO();
        stepOnePi();
    }

    // update gamma
    stepGamma();

    // Update Sigmas, Rhos and Betas given all rest
    stepSigmaRhoAndBeta();

    // increase iteration counter
    ++ internalIterationCounter;

    // update the MH proposal variance
    updateProposalVariances();
}


// *******************************
// Global operators between two chains
// *******************************
// asuming nu and other fixed hyperparameters are the same across chains, woudn;t make sense otherwise I think 


void dSUR_Chain::swapTau( std::shared_ptr<dSUR_Chain>& that )
{
    double par = this->getTau();

    this->setTau( that->getTau() );
    that->setTau( par );
}

void dSUR_Chain::swapSigmaRho( std::shared_ptr<dSUR_Chain>& that )
{
    arma::mat par = this->getSigmaRho();

    this->setSigmaRho( that->getSigmaRho() );
    that->setSigmaRho( par );
}

void dSUR_Chain::swapO( std::shared_ptr<dSUR_Chain>& that )
{
    arma::vec par = this->getO();

    this->setO( that->getO() );
    that->setO( par );
}

void dSUR_Chain::swapPi( std::shared_ptr<dSUR_Chain>& that )
{
    arma::vec par = this->getPi();

    this->setPi( that->getPi() );
    that->setPi( par );
}

void dSUR_Chain::swapGamma( std::shared_ptr<dSUR_Chain>& that )
{
    arma::umat par = this->getGamma();

    this->setGamma( that->getGamma() );
    that->setGamma( par );
}

void dSUR_Chain::swapW( std::shared_ptr<dSUR_Chain>& that )
{
    double par = this->getW();

    this->setW( that->getW() );
    that->setW( par );
}

void dSUR_Chain::swapBeta( std::shared_ptr<dSUR_Chain>& that )
{
    arma::mat par = this->getBeta();

    this->setBeta( that->getBeta() );
    that->setBeta( par );
}


int dSUR_Chain::exchangeGamma_step( std::shared_ptr<dSUR_Chain>& that )
{
    // I'm exchanging the gammas AND the betas. So gammaMask, XB and U will follow and we will have to re-compute rhoU for both chains
    arma::umat swapGammaMask;
    arma::mat swapXB , swapU;

    arma::mat rhoU_1 = this->createRhoU( that->getU() , this->getSigmaRho() );
    arma::mat rhoU_2 = that->createRhoU( this->getU() , that->getSigmaRho() );

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

int dSUR_Chain::adapt_crossOver_step( std::shared_ptr<dSUR_Chain>& that )
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

    pCrossOver += this->logPBetaGivenSigmaRho( this->getBeta() , this->getSigmaRho() , this->getJT(),   
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

int dSUR_Chain::uniform_crossOver_step( std::shared_ptr<dSUR_Chain>& that )
{
    double pCrossOver;

    std::vector<arma::umat> gammaXO(2); gammaXO[0] = arma::umat(nVSPredictors,nOutcomes);  gammaXO[1] = arma::umat(nVSPredictors,nOutcomes); 

    // Propose Crossover
    for(unsigned int j=0; j<nVSPredictors; ++j)
    {
        for(unsigned int k=0; k<nOutcomes; ++k)
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
    pCrossOver -= this->sampleBetaGivenSigmaRho( betaXO[0] , this->getSigmaRho() , this->getJT(), 
                                gammaMask_XO[0] , XB_XO[0] , U_XO[0] , rhoU_XO[0] );

    pCrossOver += this->logPBetaGivenSigmaRho( this->getBeta() , this->getSigmaRho() , this->getJT(),   
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

int dSUR_Chain::block_crossOver_step( std::shared_ptr<dSUR_Chain>& that , arma::mat& corrMatX , double threshold )
{
    double pCrossOver;

    std::vector<arma::umat> gammaXO(2); gammaXO[0] = arma::umat(nVSPredictors,nOutcomes);  gammaXO[1] = arma::umat(nVSPredictors,nOutcomes); 

    // Propose Crossover

    // Select the ONE index to foor the block
    unsigned int predIdx = Distributions::randIntUniform(0, nVSPredictors-1 ); // pred
    unsigned int outcIdx = Distributions::randIntUniform(0, nOutcomes-1 ); // outcome

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
    pCrossOver -= this->sampleBetaGivenSigmaRho( betaXO[0] , this->getSigmaRho() , this->getJT(), 
                                gammaMask_XO[0] , XB_XO[0] , U_XO[0] , rhoU_XO[0] );

    pCrossOver += this->logPBetaGivenSigmaRho( this->getBeta() , this->getSigmaRho() , this->getJT(),   
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

void dSUR_Chain::swapAll( std::shared_ptr<dSUR_Chain>& thatChain )
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

int dSUR_Chain::globalStep( std::shared_ptr<dSUR_Chain>& that )
{

    unsigned int globalType = Distributions::randIntUniform(0,4);

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
            return this -> exchangeAll_step( that );
            break;
        
        default: 
            break;
    }

    return 0;
}

int dSUR_Chain::exchangeAll_step( std::shared_ptr<dSUR_Chain>& thatChain )
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

// *******************************
// Other Methods
// *******************************

arma::mat dSUR_Chain::createRhoU( const arma::mat& externalU , const arma::mat&  externalSigmaRho )
{
    arma::mat externalRhoU = arma::zeros<arma::mat>(nObservations,nOutcomes);

    for( unsigned int k=1; k < nOutcomes; ++k)
    {
        for(unsigned int l=0 ; l<k ; ++l)
        {
            if(  externalSigmaRho(k,l) != 0 )
                externalRhoU.col(k) += externalU.col(l) *  externalSigmaRho(k,l);
        }
    }

    return externalRhoU;
}


void dSUR_Chain::updateRhoU()
{
    rhoU.zeros(nObservations,nOutcomes);

    for( unsigned int k=1; k < nOutcomes; ++k)
    {
        for(unsigned int l=0 ; l<k ; ++l)
        {
            if(  sigmaRho(k,l) != 0 )
                rhoU.col(k) += U.col(l) * sigmaRho(k,l);
        }
    }
}

void dSUR_Chain::createQuantities( arma::umat&  externalGammaMask , arma::mat& externalXB , arma::mat& externalU , arma::mat& externalRhoU ,
                                    const arma::umat& externalGamma , const arma::mat&  externalBeta ,
                                    const arma::mat&  externalSigmaRho )
{
     externalGammaMask = createGammaMask( externalGamma );
    externalXB = createXB(  externalGammaMask ,  externalBeta );
    externalU = createU( externalXB );
    externalRhoU = createRhoU( externalU ,  externalSigmaRho );
}

void dSUR_Chain::updateQuantities()
{
    updateGammaMask();
    updateXB();
    updateU();
    updateRhoU();
}