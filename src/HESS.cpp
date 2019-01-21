#include "global.h"
#include "utils.h"
#include "distr.h"
#include "HESS.h"
#include <armadillo>
#include <cmath>
#include <omp.h>

extern omp_lock_t RNGlock; //defined in global.h

namespace Model
{

	double logPrior(const arma::mat &omega, const arma::umat &gamma, const arma::vec& a_0, const arma::vec& b_0)
	{
		unsigned int p = gamma.n_rows;
		unsigned int s = gamma.n_cols;
		double logP = 0;

		for(unsigned int j=0; j<p; ++j)
		{
			for(unsigned int l=0; l<s; ++l)
			{
				logP += Distributions::logPDFBeta( omega(j,l) , a_0(j), b_0(j) );
				logP += Distributions::logPDFBernoulli( gamma(j,l), omega(j,l) );
			}	
		}

		return logP;
	}

	double logLikelihood(const arma::mat&Y, const arma::mat &X, const arma::umat &gamma, 
						 const double a_r_0, const double b_r_0, const arma::mat &W_0, double temp)
	{
		unsigned int n = Y.n_rows;
		unsigned int p = gamma.n_rows;
		unsigned int s = Y.n_cols;

		arma::uvec VS_IN;  // Variable Selection -- IN variables

		// std::cout << gamma.t() << find(gamma != 0).t() << (find(gamma != 0)+1).t() << VS_IN.t() << std::endl<< std::endl;

		// define variables
		arma::mat X_gamma;
		arma::mat XtX;
		// arma::vec hat_B;
		arma::mat W_n;
		arma::vec tilde_B;
		double a_r_n, b_r_n;
		double logP, sign, tmp; //sign is needed for the implementation, but we 'assume' that all the matrices are (semi-)positive-definite (-> det>=0)
		// tmp is a temporary container for the log determinants needed

		// int_lik_num2 = 0.5*log(det(W_n) )                            + sum(alpha_0*log(beta_0)) + sum(lgamma(alpha_n))
		// int_lik_den2 = 0.5*log(det( (diag(s) %x% W_0)[vidx,vidx] ) ) + sum(alpha_n*log(beta_n)) + sum(lgamma(alpha_0)) + log(2*pi)*(n*s/2)

		logP = -log(M_PI)*((double)n*(double)s*0.5); // initialise with the normalising constant remaining from the likelhood

		for(unsigned int l=0; l<s; ++l)
		{

			VS_IN = find(gamma.col(l) != 0)+1; VS_IN.insert_rows(0,1); // note that insert_rows works by inserting 1 row in position 0, and set it to zero since I'm passing an int and not a type arma::uvec(1)!!

			X_gamma = X.cols(VS_IN);				// predictors for the l^th outcome
			XtX = arma::trans(X_gamma) * X_gamma;

			// hat_B = arma::inv_sympd(XtX) * X_gamma.t() * Y.col(l);
			W_n = arma::inv_sympd( XtX + arma::inv_sympd( W_0(VS_IN,VS_IN) ) );
			tilde_B = W_n * ( X_gamma.t()*Y.col(l)  /* + W_0.i() * ZERO  */ );

			a_r_n = a_r_0 + 0.5*(double)n;

			b_r_n = b_r_0 + 0.5* arma::as_scalar( Y.col(l).t() * Y.col(l) - ( tilde_B.t() * arma::inv_sympd(W_n) * tilde_B ) );

			arma::log_det(tmp, sign, W_n );
			logP += 0.5*tmp; 

			arma::log_det(tmp, sign, W_0(VS_IN,VS_IN) );
			logP += -0.5*tmp; 

			logP += a_r_0*log(b_r_0) - a_r_n*log(b_r_n);

			logP += std::lgamma(a_r_n) - std::lgamma(a_r_0);
		}

		return logP/temp;

	}

	void banditProposal
	(
	arma::umat &gamma, arma::mat &zeta, arma::umat &gamma_prop, arma::mat& alpha_z, arma::mat& beta_z,
	arma::vec& mismatch, arma::vec& normalised_mismatch, arma::vec& normalised_mismatch_backwards,
	const arma::mat&Y, const arma::mat &X, arma::mat& omega_curr, double& logPrior_curr, double &logLik_curr,
	const double a_r_0, const double b_r_0, const arma::mat& W_0, const arma::vec& a_0, const arma::vec& b_0,
	double& accCount, unsigned int nUpdates, double temp
	)
	{

		unsigned int p = gamma.n_rows;
		unsigned int s = gamma.n_cols;
		double log_proposal = 0.;

		double banditIncrement = 1.;
		double banditLimit = 500.;

		// Sample Zs
		for(unsigned int i=0; i<p; ++i)
		{
			for(unsigned int j=0; j<s; ++j)
			{
				zeta(i,j) = Distributions::randBeta(alpha_z(i,j),beta_z(i,j));
			}
		}

		// Create mismatch
		for(unsigned int i=0; i<(p*s); ++i)
		{
			mismatch(i) = (gamma(i)==0)?(zeta(i)):(1.-zeta(i));   //mismatch
		}

		// Normalise
		// mismatch = arma::log(mismatch); //logscale
		// normalised_mismatch = mismatch - Utils::logspace_add(mismatch);


		normalised_mismatch = mismatch / arma::as_scalar(arma::sum(mismatch));

		arma::uvec update_idx;


		if( Distributions::randU01() < 0.5 )   // one deterministic update
		{

			// Decide which to update
			update_idx = arma::zeros<arma::uvec>(1);

			update_idx = Distributions::randWeightedIndexSampleWithoutReplacement(p*s,normalised_mismatch); // sample the one

			// Update
			gamma_prop(update_idx) = 1 - gamma(update_idx); // deterministic, just switch

			// Compute log_proposal probabilities
			normalised_mismatch_backwards = mismatch;
			normalised_mismatch_backwards(update_idx) = 1. - normalised_mismatch_backwards(update_idx) ;
			// normalised_mismatch_backwards = normalised_mismatch_backwards - Utils::logspace_add(normalised_mismatch_backwards);
			normalised_mismatch_backwards = normalised_mismatch_backwards / arma::sum(normalised_mismatch_backwards);

			log_proposal = arma::as_scalar( arma::log( normalised_mismatch_backwards(update_idx) ) ) -
						 arma::as_scalar( arma::log( normalised_mismatch(update_idx) ) );

		}else{

			/*
			L(nUpdates) random (bern) updates
			Note that we make use of column indexing here for armadillo matrices
			*/

			log_proposal = 0.;
			// Decide which to update
			update_idx = arma::zeros<arma::uvec>(nUpdates);
			update_idx = Distributions::randWeightedIndexSampleWithoutReplacement(p*s,normalised_mismatch,nUpdates); // sample nUpdates indexes

			normalised_mismatch_backwards = mismatch; // copy for backward proposal


			// Update
			for(unsigned int i=0; i<nUpdates; ++i)
			{
				gamma_prop(update_idx(i)) = Distributions::randBernoulli(zeta(update_idx(i))); // random update

				if(gamma_prop(update_idx(i)) != gamma(update_idx(i)))
				{
					normalised_mismatch_backwards(update_idx(i)) = 1.- normalised_mismatch_backwards(update_idx(i));
					log_proposal += Distributions::logPDFBernoulli(gamma(update_idx(i)),zeta(update_idx(i))) -
						Distributions::logPDFBernoulli(gamma_prop(update_idx(i)),zeta(update_idx(i)));
				}
			}

			// Compute log_proposal probabilities
			// normalised_mismatch_backwards = normalised_mismatch_backwards - Utils::logspace_add(normalised_mismatch_backwards);
			normalised_mismatch_backwards = normalised_mismatch_backwards / arma::sum(normalised_mismatch_backwards);

			log_proposal += Distributions::logPDFWeightedIndexSampleWithoutReplacement(normalised_mismatch_backwards,update_idx) -
				Distributions::logPDFWeightedIndexSampleWithoutReplacement(normalised_mismatch,update_idx);
		}

		// A/R step

		double logLik_prop = logLikelihood(Y,X, gamma_prop, a_r_0, b_r_0, W_0, temp);
		double logPrior_prop = logPrior(omega_curr, gamma, a_0, b_0)/temp;

		double logAccProb = (logLik_prop - logLik_curr) + (logPrior_prop - logPrior_curr) + log_proposal;

		if( Distributions::randLogU01() < logAccProb )
		{
			accCount += update_idx.n_elem;/*/(double)nUpdates*/;
			gamma = gamma_prop;
			logLik_curr = logLik_prop;
			logPrior_curr = logPrior_prop;

		}

		for(arma::uvec::iterator iter = update_idx.begin(); iter != update_idx.end(); ++iter)
		{
			if( alpha_z(*iter) + beta_z(*iter) < banditLimit )
			{
				alpha_z(*iter) += banditIncrement * gamma(*iter);
				beta_z(*iter) += banditIncrement * (1-gamma(*iter));
			}
		}


	}


	void bandit_SUR_MCMC_step(const arma::mat&Y, const arma::mat &X, arma::mat& omega_curr, arma::umat& gamma_curr, double& logPrior_curr, double &logLik_curr, 
					const double a_r_0, const double b_r_0, const arma::mat& W_0, 
					const arma::vec& a_0, const arma::vec& b_0, double& accCount, unsigned int nUpdates, double temp,
					arma::mat& zeta, arma::mat& alpha_z, arma::mat& beta_z, arma::vec& mismatch,
					arma::vec& normalised_mismatch, arma::vec& normalised_mismatch_backwards)
	{
		
		// Initialise proposal
		// arma::mat omega_prop;
		arma::umat gamma_prop = gamma_curr;

		// unsigned int n = Y.n_rows;
		unsigned int p = gamma_curr.n_rows;
		unsigned int s = Y.n_cols;
		arma::imat predUpdate;       // which predictor are we updating?
		arma::uword /*arma::ivec*/ outUpdate; 		 // which outcome?
		
		double logPrior_prop, logLik_prop;
		double logAccProb;

		predUpdate = arma::join_rows( Distributions::randIntUniform(nUpdates, 0,p-1) , Distributions::randIntUniform(nUpdates,0,s-1) );    // note that I might be updating multiple times the same coeff
		
		double a,b;

		// ## Update Omega
		// FULL CONDITIONAL (no sense in making it parallel for now as the RNG is still sequential)
		// note that a_0 and b_0 are vectors of size p, meaning there's one prior for each beta, regardless of the associated Y
		for(unsigned int k=0; k<s; ++k)
		{
			for(unsigned int j=0; j<p; ++j)
			{
				a = a_0(j) + gamma_curr(j,k); b =  b_0(j) + 1. - gamma_curr(j,k);
				omega_curr(j,k) = Distributions::randBeta( (a+1.)/temp - 1. , b/temp );
			}
		}
		
		// Update the prior
		logPrior_curr = logPrior(omega_curr, gamma_curr, a_0, b_0)/temp; // update the prior probability


		// ## Update Gamma -- Local move
		banditProposal(gamma_curr, zeta, gamma_prop, alpha_z, beta_z,
				mismatch, normalised_mismatch, normalised_mismatch_backwards,
				Y, X, omega_curr, logPrior_curr, logLik_curr, a_r_0, b_r_0, W_0, a_0, b_0,
				accCount, nUpdates, temp);

	}

		void MC3_SUR_MCMC_step(const arma::mat&Y, const arma::mat &X, arma::mat& omega_curr, arma::umat& gamma_curr, double& logPrior_curr, double &logLik_curr, 
					const double a_r_0, const double b_r_0, const arma::mat& W_0, 
					const arma::vec& a_0, const arma::vec& b_0, double& accCount, unsigned int nUpdates, double temp)
	{
		
		// Initialise proposal
		// arma::mat omega_prop;
		arma::umat gamma_prop = gamma_curr;

		// unsigned int n = Y.n_rows;
		unsigned int p = gamma_curr.n_rows;
		unsigned int s = Y.n_cols;
		arma::imat predUpdate;       // which predictor are we updating?
		arma::uword /*arma::ivec*/ outUpdate; 		 // which outcome?
		
		double logPrior_prop, logLik_prop;
		double logAccProb;

		predUpdate = arma::join_rows( Distributions::randIntUniform(nUpdates, 0,p-1) , Distributions::randIntUniform(nUpdates,0,s-1) );    // note that I might be updating multiple times the same coeff
		
		double a,b;

		// ## Update Omega
		// FULL CONDITIONAL (no sense in making it parallel for now as the RNG is still sequential)
		// note that a_0 and b_0 are vectors of size p, meaning there's one prior for each beta, regardless of the associated Y
		for(unsigned int k=0; k<s; ++k)
		{
			for(unsigned int j=0; j<p; ++j)
			{
				a = a_0(j) + gamma_curr(j,k); b =  b_0(j) + 1. - gamma_curr(j,k);
				omega_curr(j,k) = Distributions::randBeta( (a+1.)/temp - 1. , b/temp );
			}
		}
		
		// Update the prior
		logPrior_curr = logPrior(omega_curr, gamma_curr, a_0, b_0)/temp; // update the prior probability


		// ## Update Gamma -- Local move
		predUpdate = Distributions::randIntUniform(nUpdates,0,p-1);    // note that I might be updating multiple times the same coeff

		for(unsigned int j=0; j<nUpdates ; ++j)
		{
			gamma_prop = gamma_curr;

			// Decide which l (l=0,...s-1) to sample
			outUpdate = Distributions::randIntUniform(0,s-1);    //maybe I need to update more often here that I have a matrix?

			/// Sample Gamma (see http://www.genetics.org/content/suppl/2011/09/16/genetics.111.131425.DC1/131425SI.pdf pag 5)
			gamma_prop(predUpdate(j),outUpdate) = Distributions::randBernoulli( pow(omega_curr(predUpdate(j),outUpdate),1./temp) / ( pow(omega_curr(predUpdate(j),outUpdate),1./temp) + pow( 1. - omega_curr(predUpdate(j),outUpdate) , 1./temp) ) ); 
			
			if( gamma_prop(predUpdate(j),outUpdate) != gamma_curr(predUpdate(j),outUpdate) )
			{
				// Compute log Probability for the current point
				logLik_prop = logLikelihood(Y,X, gamma_prop, a_r_0, b_r_0, W_0, temp);

				// A/R   
				logAccProb = (logLik_prop - logLik_curr);

				if( Distributions::randLogU01() < logAccProb )
				{
					accCount += 1./(double)nUpdates;
					gamma_curr(predUpdate(j),outUpdate) = gamma_prop(predUpdate(j),outUpdate);				
					logLik_curr = logLik_prop;
					logPrior_curr = logPrior(omega_curr, gamma_curr, a_0, b_0)/temp; // update the prior probability
				}
			}
		} // finish nUpdates on Gamma
	}


	void exchange_step(arma::ucube& gamma_state, arma::cube& omega_state, 
						arma::vec& logPrior_state, arma::vec& logLik_state,    				// states
						const arma::vec& temperature, const unsigned int nChains, const unsigned int nGlobalUpdates,// tuning pars
						double& accCountGlobalUpdates, unsigned int& countGlobalUpdates)
	{


		double pExchange;
		unsigned int chainIdx, firstChain, secondChain;

		arma::mat omegaSwap;
		arma::umat gammaSwap; // TODO, there's really no need to use a temp matrix when only a temp double is needed to swap in a loop...
		double logPriorSwap;
		double logLikSwap;


		// Select the chains to swap
		chainIdx = (nChains>2) ? Distributions::randIntUniform(1, (nChains)*(nChains-1)/2 ) : 1;   // (nChains-1)*(nChains-2)/2 is the number of possible chain combinations with nChains

		for(unsigned int c=1; c<nChains; ++c)
		{
			for(unsigned int r=0; r<c; ++r)
			{
				if( (--chainIdx) == 0 ){
					firstChain = r;
					secondChain = c;
					break;
				}
			}
		}

		countGlobalUpdates++;

		// Swap probability
		pExchange = logLik_state(secondChain) * ( temperature(secondChain)/temperature(firstChain) - 1.) +   // the plus is correct, don't doubt
					logLik_state(firstChain) * ( temperature(firstChain)/temperature(secondChain) - 1.) ;

		// The prior doesn't come into play as it's not tempered and hence there's no effect from it

		// DEBUG CHECK
		if(firstChain == 0 || secondChain == 0)
		{
			// std::cout << std::exp(pExchange) << std::endl <<std::flush;
			// std::cout << "Current -- " << logLik_state(firstChain)  << "  --  " << logLik_state(secondChain) << "  ===  " << logLik_state(firstChain) + logLik_state(secondChain)  << std::endl <<std::flush;
			// std::cout << "Proposed -- " << logLik_state(secondChain)*temperature(secondChain)/temperature(firstChain)  << "  --  " << logLik_state(firstChain)*temperature(firstChain)/temperature(secondChain)  <<
			//  "  ===  " << logLik_state(secondChain)*temperature(secondChain)/temperature(firstChain) + logLik_state(firstChain)*temperature(firstChain)/temperature(secondChain) << std::endl <<std::flush;
			// std::cout << pExchange - (
			// logLik_state(secondChain)*temperature(secondChain)/temperature(firstChain) + 
			// logLik_state(firstChain)*temperature(firstChain)/temperature(secondChain) - 
			// logLik_state(firstChain) - logLik_state(secondChain) ) << std::endl <<std::flush;
			// std::cout << temperature(firstChain)  << "  --  " << temperature(secondChain) << std::endl <<std::flush;
			// int jnk; std::cin >> jnk;
			// return ;
		}


		// A/R
		if( Distributions::randLogU01() < pExchange )
		{

			omegaSwap = omega_state.slice(secondChain);
			gammaSwap = gamma_state.slice(secondChain);
			logPriorSwap = logPrior_state(secondChain);
			logLikSwap = logLik_state(secondChain)*temperature(secondChain)/temperature(firstChain);

			omega_state.slice(secondChain) = omega_state.slice(firstChain);
			gamma_state.slice(secondChain) = gamma_state.slice(firstChain);

			logPrior_state(secondChain) = logPrior_state(firstChain);
			logLik_state(secondChain) =  logLik_state(firstChain)*temperature(firstChain)/temperature(secondChain);

			omega_state.slice(firstChain) = omegaSwap;
			gamma_state.slice(firstChain) = gammaSwap;
			logPrior_state(firstChain) = logPriorSwap;
			logLik_state(firstChain) = logLikSwap;

			accCountGlobalUpdates++;	

		}


	} //end generic Exchange step

	void nearExchange_step(arma::ucube& gamma_state, arma::cube& omega_state,
						arma::vec& logPrior_state, arma::vec& logLik_state,    				// states
						const arma::vec& temperature, const unsigned int nChains, const unsigned int nGlobalUpdates,// tuning pars
						double& accCountGlobalUpdates, unsigned int& countGlobalUpdates)
	{

		double pExchange;
		unsigned int chainIdx, firstChain, secondChain;

		arma::mat omegaSwap;
		arma::umat gammaSwap;
		double logPriorSwap;
		double logLikSwap;

		// Select the chains to swap

		if( nChains>2 )
		{
			firstChain = Distributions::randIntUniform(1, nChains-2 );  // so not the first (0) or last (nChains-1) indexes
			secondChain = ( Distributions::randU01() < 0.5 ) ? firstChain-1 : firstChain+1 ; // then select a neighbour
		}else{
			// with just 2 chains
			firstChain = 0;
			secondChain = 1;
		}


		// if( firstChain == 0 || secondChain == 0 )   // if the swap involves the main chain
		countGlobalUpdates++;

		// Swap probability
		pExchange = 	logLik_state(secondChain) * ( temperature(secondChain)/temperature(firstChain) - 1.) +   // the plus is correct, don't doubt
						logLik_state(firstChain) * ( temperature(firstChain)/temperature(secondChain) - 1.) ;

		// A/R
		if( Distributions::randLogU01() < pExchange )
		{

			omegaSwap = omega_state.slice(secondChain);
			gammaSwap = gamma_state.slice(secondChain);
						logPriorSwap = logPrior_state(secondChain);
			logLikSwap = logLik_state(secondChain)*temperature(secondChain)/temperature(firstChain);

			omega_state.slice(secondChain) = omega_state.slice(firstChain);
			gamma_state.slice(secondChain) = gamma_state.slice(firstChain);

			logPrior_state(secondChain) = logPrior_state(firstChain);
			logLik_state(secondChain) =  logLik_state(firstChain)*temperature(firstChain)/temperature(secondChain);

			omega_state.slice(firstChain) = omegaSwap;
			gamma_state.slice(firstChain) = gammaSwap;
			logPrior_state(firstChain) = logPriorSwap;
			logLik_state(firstChain) = logLikSwap;

			// if (firstChain == 0 || secondChain == 0 )
			accCountGlobalUpdates++;

		}
	} //end Near Exchange

	void allExchange_step(arma::ucube& gamma_state, arma::cube& omega_state,
						arma::vec& logPrior_state, arma::vec& logLik_state,    				// states
						const arma::vec& temperature, const unsigned int nChains, const unsigned int nGlobalUpdates,// tuning pars
						double& accCountGlobalUpdates, unsigned int& countGlobalUpdates)
	{

		arma::vec pExchange( (nChains)*(nChains-1)/2 +1 );
		unsigned int swapIdx, firstChain, secondChain;

		arma::umat indexTable( pExchange.n_elem, 2);
		unsigned int tabIndex = 0;
		indexTable(tabIndex,0) = 0; indexTable(tabIndex,1) = 0;
		tabIndex++;

		for(unsigned int c=1; c<nChains; ++c)
		{
			for(unsigned int r=0; r<c; ++r)
			{
				indexTable(tabIndex,0) = r; indexTable(tabIndex,1) = c;
				tabIndex++;
			}
		}

		arma::mat omegaSwap;
		arma::umat gammaSwap;
		double logPriorSwap;
		double logLikSwap;

		countGlobalUpdates++;

		// Select the chains to swap
		tabIndex = 0;
		pExchange(tabIndex) = 0.; // these are log probabilities, remember!
		tabIndex++;

		#pragma omp parallel for private(tabIndex, firstChain, secondChain)
		for(tabIndex = 1; tabIndex <= ((nChains)*(nChains-1)/2); ++tabIndex)
		{

			firstChain = indexTable(tabIndex,0);
			secondChain  = indexTable(tabIndex,1);

			// Swap probability
			// pExchange(tabIndex) = logLik_state(firstChain) * temperature(firstChain) / temperature(secondChain) +
			// 	logLik_state(secondChain) * temperature(secondChain) / temperature(firstChain) -
			// 	logLik_state(firstChain) - logLik_state(secondChain);

			pExchange(tabIndex) = 	logLik_state(secondChain) * ( temperature(secondChain)/temperature(firstChain) - 1.) +   // the plus is correct, don't doubt
				logLik_state(firstChain) * ( temperature(firstChain)/temperature(secondChain) - 1.) ;


		}

		// normalise and cumulate the weights
		double logSumWeights = Utils::logspace_add(pExchange); // normaliser
		arma::vec cumulPExchange = arma::cumsum( arma::exp( pExchange - logSumWeights ) ); // this should sum to one

		// Now select which swap happens
		double val = Distributions::randU01();

		swapIdx = 0;
		while( val > cumulPExchange(swapIdx) )
		{
			swapIdx++;
		}

		if( swapIdx != 0 )
		{

			firstChain = indexTable(swapIdx,0);
			secondChain  = indexTable(swapIdx,1);

			accCountGlobalUpdates++;

			omegaSwap = omega_state.slice(secondChain);
			gammaSwap = gamma_state.slice(secondChain);
			logPriorSwap = logPrior_state(secondChain);
			logLikSwap = logLik_state(secondChain)*temperature(secondChain)/temperature(firstChain);

			omega_state.slice(secondChain) = omega_state.slice(firstChain);
			gamma_state.slice(secondChain) = gamma_state.slice(firstChain);

			logPrior_state(secondChain) = logPrior_state(firstChain);
			logLik_state(secondChain) =  logLik_state(firstChain)*temperature(firstChain)/temperature(secondChain);

			omega_state.slice(firstChain) = omegaSwap;
			gamma_state.slice(firstChain) = gammaSwap;
			logPrior_state(firstChain) = logPriorSwap;
			logLik_state(firstChain) = logLikSwap;
		}
		//else swapIdx = 0 means no swap at all

	} //end ALL Exchange

	void uniformCrossOver_step(arma::ucube& gamma_state, arma::cube& omega_state,
						arma::vec& logPrior_state, arma::vec& logLik_state,    	// states
						const double a_r_0, const double b_r_0, const arma::mat& W_0,
						const arma::vec& a_0, const arma::vec& b_0,const arma::mat& Y,const arma::mat& X,	// Prior pars, data
						const arma::vec& temperature, const unsigned int nChains, const unsigned int nGlobalUpdates, // hyper tuning pars
						double& accCountGlobalUpdates, unsigned int& countGlobalUpdates)
	{


		unsigned int p = X.n_cols-1; // there's the intercept
		unsigned int s = Y.n_cols;

		double pCrossOver;

		double logPriorFirst;
		double logPriorSecond;
		double logLikFirst;
		double logLikSecond;

		arma::ucube gammaXO(p,s,2);
		arma::cube omegaXO(p,s,2);
		unsigned int chainIdx, firstChain, secondChain;

		// Select the chains to XO 
		chainIdx = (nChains>2) ? Distributions::randIntUniform(1, (nChains)*(nChains-1)/2 ) : 1;   // (nChains-1)*(nChains-2)/2 is the number of possible chain combinations with nChains

		for(unsigned int c=1; c<nChains; ++c)
		{
			for(unsigned int r=0; r<c; ++r)
			{
				if( (--chainIdx) == 0 ){
					firstChain = r;
					secondChain = c;
					break;
				}
			}
		}

		countGlobalUpdates++; // a global update is happening

		// Propose Crossover
		for(unsigned int j=0; j<p; ++j)
		{
			for(unsigned int l=0; l<s; ++l)
			{
				if( Distributions::randU01() < 0.5 )
				{
					gammaXO(j,l,0) = gamma_state(j,l,firstChain);
					gammaXO(j,l,1) = 1-gammaXO(j,l,0);

					omegaXO(j,l,0) = omega_state(j,l,firstChain);
					omegaXO(j,l,1) = omega_state(j,l,secondChain);

				}else{
					gammaXO(j,l,0) = gamma_state(j,l,secondChain);
					gammaXO(j,l,1) = 1-gammaXO(j,l,0);

					omegaXO(j,l,0) = omega_state(j,l,secondChain);
					omegaXO(j,l,1) = omega_state(j,l,firstChain);
				}
			}
		}

		// Probability of acceptance

		logPriorFirst = logPrior(omega_state.slice(firstChain), gammaXO.slice(0), a_0, b_0)/temperature(firstChain);
		logPriorSecond = logPrior(omega_state.slice(secondChain), gammaXO.slice(1), a_0, b_0)/temperature(secondChain);
		logLikFirst = logLikelihood(Y, X, gammaXO.slice(0), a_r_0, b_r_0, W_0, temperature(firstChain));
		logLikSecond = logLikelihood(Y, X, gammaXO.slice(1), a_r_0, b_r_0, W_0, temperature(secondChain));

		pCrossOver = (	logPriorFirst + logLikFirst - logPrior_state(firstChain) - logLik_state(firstChain) ) - 
						(	logPriorSecond + logLikSecond - logPrior_state(secondChain) - logLik_state(secondChain) );

		pCrossOver += 0.;  // XO prop probability simmetric now

		// A/R
		if( Distributions::randLogU01() < pCrossOver )
		{
			gamma_state.slice(firstChain) = gammaXO.slice(0);
			omega_state.slice(firstChain) = omegaXO.slice(0);
			logLik_state(firstChain) = logLikFirst;

			gamma_state.slice(secondChain) = gammaXO.slice(1);
			omega_state.slice(secondChain) = omegaXO.slice(1);
			logLik_state(secondChain) = logLikSecond;

			accCountGlobalUpdates++;

		} // end CrossOver
	}


	void blockCrossOver_step(arma::ucube& gamma_state, arma::cube& omega_state,
						arma::vec& logPrior_state, arma::vec& logLik_state,    	// states
						const double a_r_0, const double b_r_0, const arma::mat& W_0,
						const arma::vec& a_0, const arma::vec& b_0,const arma::mat& Y,const arma::mat& X,
						const arma::mat& covariatesCorrelation, const arma::vec& temperature,	// Prior pars, data
						const double threshold, const unsigned int nChains, const unsigned int nGlobalUpdates,		// hyper tuning pars
						double& accCountGlobalUpdates, unsigned int& countGlobalUpdates)
	{

		unsigned int p = X.n_cols-1; // there's the intercept
		unsigned int s = Y.n_cols;

		double pCrossOver;

		double logPriorFirst;
		double logPriorSecond;
		double logLikFirst;
		double logLikSecond;

		arma::ucube gammaXO(p,s,2);
		arma::cube omegaXO(p,s,2);

		unsigned int chainIdx, firstChain, secondChain;

		// Select the chains to XO
		chainIdx = (nChains>2) ? Distributions::randIntUniform(1, (nChains)*(nChains-1)/2 ) : 1;   // (nChains-1)*(nChains-2)/2 is the number of possible chain combinations with nChains

		for(unsigned int c=1; c<nChains; ++c)
		{
			for(unsigned int r=0; r<c; ++r)
			{
				if( (--chainIdx) == 0 ){
					firstChain = r;
					secondChain = c;
					break;
				}
			}
		}

		countGlobalUpdates++; // a global update is happening

		// Propose Crossover

		// Select the ONE index to foor the block
		unsigned int predIdx = Distributions::randIntUniform(0, p-1 ); // pred
		unsigned int outcIdx = Distributions::randIntUniform(0, s-1 ); // outcome


		arma::uvec covIdx = arma::find( arma::abs( covariatesCorrelation.row(predIdx) ) > threshold );  // this will include the original predIdx

		gammaXO.slice(0) = gamma_state.slice(firstChain);
		gammaXO.slice(1) = gamma_state.slice(secondChain);

		omegaXO.slice(0) = omega_state.slice(firstChain);
		omegaXO.slice(1) = omega_state.slice(secondChain);

		for(unsigned int j=0; j<covIdx.n_elem; ++j)
		{
			gammaXO(covIdx(j),outcIdx,0) = gamma_state(covIdx(j),outcIdx,secondChain);
			gammaXO(covIdx(j),outcIdx,1) = gamma_state(covIdx(j),outcIdx,firstChain);

			omegaXO(covIdx(j),outcIdx,0) = omega_state(covIdx(j),outcIdx,secondChain);
			omegaXO(covIdx(j),outcIdx,1) = omega_state(covIdx(j),outcIdx,firstChain);

		}
		// Probability of acceptance
		logPriorFirst = logPrior(omega_state.slice(firstChain), gammaXO.slice(0), a_0, b_0)/temperature(firstChain);
		logPriorSecond = logPrior(omega_state.slice(secondChain), gammaXO.slice(1), a_0, b_0)/temperature(secondChain);
		logLikFirst = logLikelihood(Y, X, gammaXO.slice(0), a_r_0, b_r_0, W_0, temperature(firstChain));
		logLikSecond = logLikelihood(Y, X, gammaXO.slice(1), a_r_0, b_r_0, W_0, temperature(secondChain));

		pCrossOver = (	logPriorFirst + logLikFirst - logPrior_state(firstChain) - logLik_state(firstChain) ) - 
						(	logPriorSecond + logLikSecond - logPrior_state(secondChain) - logLik_state(secondChain) );

		pCrossOver += 0.;  // XO prop probability is weird, how do I compute it? Let's say is symmetric as is determnistic and both comes from the same covariatesCorrelation

		// A/R
		if( Distributions::randLogU01() < pCrossOver )
		{
			gamma_state.slice(firstChain) = gammaXO.slice(0);
			omega_state.slice(firstChain) = omegaXO.slice(0);
			logLik_state(firstChain) = logLikFirst;

			gamma_state.slice(secondChain) = gammaXO.slice(1);
			omega_state.slice(secondChain) = omegaXO.slice(1);
			logLik_state(secondChain) = logLikSecond;

			accCountGlobalUpdates++;

		} // end CrossOver
	}

	void adapCrossOver_step(arma::ucube& gamma_state, arma::cube& omega_state, 
						arma::vec& logPrior_state, arma::vec& logLik_state,    	// states
						const double a_r_0, const double b_r_0, const arma::mat& W_0,
						const arma::vec& a_0, const arma::vec& b_0,const arma::mat& Y,const arma::mat& X,	// Prior pars, data
						const arma::vec& temperature, double pXO_0, double pXO_1, double pXO_2,
						double p11, double p12, double p21, double p22,									// tuning pars
						const unsigned int nChains, const unsigned int nGlobalUpdates,					// hyper tuning pars
						double& accCountGlobalUpdates, unsigned int& countGlobalUpdates) 
	{

		unsigned int p = gamma_state.n_rows;
		unsigned int s = Y.n_cols;

		double pCrossOver;
		unsigned int n11,n12,n21,n22;

		double logPriorFirst;
		double logPriorSecond;
		double logLikFirst;
		double logLikSecond;
		
		arma::ucube gammaXO(p,s,2);
		unsigned int chainIdx, firstChain, secondChain;


		// Select the chains to XO
		chainIdx = (nChains>2) ? Distributions::randIntUniform(1, (nChains)*(nChains-1)/2 ) : 1;   // (nChains-1)*(nChains-2)/2 is the number of possible chain combinations with nChains

		for(unsigned int c=1; c<nChains; ++c)
		{
			for(unsigned int r=0; r<c; ++r)
			{
				if( (--chainIdx) == 0 ){
					firstChain = r;
					secondChain = c;
					break;
				}
			}
		}
		// if( firstChain == 0 || secondChain == 0 )
		countGlobalUpdates++;

		// Propose Crossover
		n11=0;n12=0;n21=0;n22=0;

		for(unsigned int j=0; j<p; ++j)
		{ 
			for(unsigned int l=0; l<s; ++l)
			{
				if ( gamma_state(j,l,firstChain) == gamma_state(j,l,secondChain) )
				{
					gammaXO(j,l,0) = gamma_state(j,l,firstChain);
					gammaXO(j,l,1) = gamma_state(j,l,secondChain);

					gammaXO(j,l,0) = ( Distributions::randU01() < pXO_0 )? 1-gammaXO(j,l,0) : gammaXO(j,l,0);
					gammaXO(j,l,1) = ( Distributions::randU01() < pXO_0 )? 1-gammaXO(j,l,1) : gammaXO(j,l,1);
					
					if( gammaXO(j,l,0) == gammaXO(j,l,1) )
						++n11;
					else
						++n12;
				}
				else
				{
					gammaXO(j,l,0) = gamma_state(j,l,firstChain);
					gammaXO(j,l,1) = gamma_state(j,l,secondChain);

					gammaXO(j,l,0) = ( Distributions::randU01() < pXO_1 )? 1-gammaXO(j,l,0) : gammaXO(j,l,0);
					gammaXO(j,l,1) = ( Distributions::randU01() < pXO_2 )? 1-gammaXO(j,l,1) : gammaXO(j,l,1);

					if( gammaXO(j,l,0) == gammaXO(j,l,1) )
						++n21;
					else
						++n22;
				}
			}
		}
		// Probability of acceptance

		logPriorFirst = logPrior(omega_state.slice(firstChain), gammaXO.slice(0), a_0, b_0)/temperature(firstChain);
		logPriorSecond = logPrior(omega_state.slice(secondChain), gammaXO.slice(1), a_0, b_0)/temperature(secondChain);
		logLikFirst = logLikelihood(Y, X, gammaXO.slice(0), a_r_0, b_r_0, W_0, temperature(firstChain));
		logLikSecond = logLikelihood(Y, X, gammaXO.slice(1), a_r_0, b_r_0, W_0, temperature(secondChain));

		pCrossOver = (	logPriorFirst + logLikFirst - logPrior_state(firstChain) - logLik_state(firstChain) ) - 
						(	logPriorSecond + logLikSecond - logPrior_state(secondChain) - logLik_state(secondChain) );

		pCrossOver += (n11 * log( p11 ) + n12 * log( p12 ) + n21 * log( p21 ) + n22 * log( p22 ) )-  // CrossOver proposal probability FORWARD
						(n11 * log( p11 ) + n12 * log( p21 ) + n21 * log( p12 ) + n22 * log( p22 ) );  // XO prop probability backward (note that ns stays the same but changes associated prob)

		// A/R
		if( Distributions::randLogU01() < pCrossOver )
		{
			gamma_state.slice(firstChain) = gammaXO.slice(0);
			logPrior_state(firstChain) = logPriorFirst;
			logLik_state(firstChain) = logLikFirst;

			gamma_state.slice(secondChain) = gammaXO.slice(1);
			logPrior_state(secondChain) = logPriorSecond;
			logLik_state(secondChain) = logLikSecond;

			// If I'm crossOver'ing on the first chain update also the mcmc_* variables
			// if(firstChain == 0)
			// {
			accCountGlobalUpdates++;
			// }
		}

	} // end CrossOver


} // end namespace