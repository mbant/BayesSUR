#include "global.h"
#include "distr.h"
#include "utils.h"

#ifndef CCODE
	using Rcpp::Rcout;
	#define BOOST_DISABLE_ASSERTS
#else
	#define Rcout std::cout
#endif

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

using namespace Rcpp;

    // [[Rcpp::export]]
	double randU01()
	{
		return R::runif( 0., 1. );
	}

    // [[Rcpp::export]]
	double randLogU01()
	{
		return log( R::runif( 0., 1. ) );
	}

    // [[Rcpp::export]]
	int randIntUniform(const int a,const int b)
	{
		return ceil( R::runif( a-1, b ) );
	}

    // [[Rcpp::export]]
	double randExponential(const double lambda)
	{
		return R::rexp( lambda );
	}

    // [[Rcpp::export]]
	arma::vec randVecExponential(const unsigned int n, const double lambda)
	{
		arma::vec res(n);
		for(unsigned int i=0; i<n; ++i)
		{
			res(i) = R::rexp( lambda );
		}
		return res;
	}

    // [[Rcpp::export]]
	unsigned int randBinomial(const unsigned int n, const double p) // slow but safe (CARE, n here is the binomial parameters, return value is always ONE integer)
	{
		return R::rbinom( n, p );
        
	}

    // [[Rcpp::export]]
	arma::uvec randMultinomial(unsigned int n, const arma::vec prob)
	{

	  unsigned int K = prob.n_elem;
	  arma::uvec rN = arma::zeros<arma::uvec>(K);
	  double p_tot = sum(prob);
	  double pp;

	  for(unsigned int k = 0 ; k < (K-1) ; ++k)
	  {
	    if(prob(k)>0) {
	    	pp = prob(k) / p_tot;
	    	rN(k) = ((pp < 1.) ? R::rbinom(n,  pp) : n);
	    	n -= rN(k);
	    }else{
	    	rN(k) = 0;
	    }


	    if(n <= 0) /* we have all*/ return rN;
	    p_tot -= prob(k); /* i.e. = sum(prob[(k+1):K]) */
	  }
	  rN(K-1) = n - sum(rN);
	  return rN;

	}

    // [[Rcpp::export]]
	double randNormal(const double m=0., const double sigmaSquare=1.) // random normal interface, parameters mean and variance
	{
		if( sigmaSquare< 0 )
			throw Distributions::negativeParameters();

    	return R::rnorm( m, sigmaSquare );
	}

    // [[Rcpp::export]]
	arma::vec randVecNormal(const unsigned int n, const double m=0., const double sigmaSquare=1.) // n-sample normal, parameters mean and variance
	{

		if( sigmaSquare< 0 )
			throw Distributions::negativeParameters();
        
    	arma::vec res(n);
    	for(unsigned int i=0; i<n; ++i)
		{
			res(i) = R::rnorm( m, sigmaSquare );
		}
		return res;
	}

    // [[Rcpp::export]]
	double randT(const double nu)
	{
    	return R::rt( nu );
	}

    // [[Rcpp::export]]
	arma::vec randVecT(const unsigned int n, const double nu)
	{
    	arma::vec res(n);
    	for(unsigned int i=0; i<n; ++i)
		{
			res(i) = R::rt( nu );
		}
		return res;
	}

    // [[Rcpp::export]]
	arma::vec randMvT(const double &nu, const arma::vec &m, const arma::mat &Sigma)
	{
		unsigned int d = m.n_elem;
		arma::rowvec res = randVecT(d,nu).t() * arma::chol(Sigma);

		return res.t() + m;
	}

    // [[Rcpp::export]]
	double randGamma(double shape, double scale)   // shape scale parametrisation
	{
		//check
		if(shape <= 0 || scale <= 0 )
		{
			throw Distributions::negativeParameters(); // THROW EXCPTION
		}

		return R::rgamma( shape, scale );
	}

    // [[Rcpp::export]]
	double randIGamma(double shape, double scale)
	{
		//check
		if(shape <= 0 || scale <= 0 )
		{
			Rcout << " Negative parameter in the gamma sampler " << '\n';
			throw Distributions::negativeParameters(); // THROW EXCPTION
		}

		return 1./R::rgamma(shape, 1./scale);
        //return 1./Rcpp::rgamma(1, shape, 1./scale)[0];
	}

    // [[Rcpp::export]]
	arma::mat randWishart(double df, const arma::mat& S)   // unsigned int df is obsolete, I see no reason to keep it
	{
		// Dimension of returned wishart
		unsigned int m = S.n_rows;

		// Z composition:
		// sqrt chisqs on diagonal (with different parameters, so no need to create the distribution object here)
		// random normals below diagonal
        
		//std::normal_distribution<> normal01(0.,1.);
        
		// misc above diagonal
		arma::mat Z(m,m);

		// Fill the diagonal
		for(unsigned int i = 0; i < m; i++){
			Z(i,i) = sqrt( randGamma( (df-i)/2.,2. ) );    // (note it's df-1:m+1)
		}

		// Fill the lower matrix with random normals
		for(unsigned int j = 0; j < m; j++){
			for(unsigned int i = j+1; i < m; i++){
		  		Z(i,j) = R::rnorm(0.,1.);
			}
		}

		// Lower triangle * chol decomp
		arma::mat C = arma::trimatl(Z).t() * arma::chol(S);

		// Return random wishart
		return C.t()*C;
	}

    // [[Rcpp::export]]
	double randBeta(double a, double b)
	{
		return R::rbeta( a, b );
	}

    // [[Rcpp::export]]
	unsigned int randBernoulli(double pi)
	{
		return R::rbinom( 1, pi );
	}


namespace Distributions{

    arma::mat randIWishart(double df, const arma::mat& S)
    {

          return arma::inv_sympd( randWishart(df,S.i()) );   // return the inverse of the correspondent Wishart variate ... is this even fast enough?
    }

    arma::vec randMvNormal(const arma::vec &m, const arma::mat &Sigma) // random normal interface to arma::randn
    {
        unsigned int d = m.n_elem;
        //check
        if(Sigma.n_rows != d || Sigma.n_cols != d )
        {
            Rcout << " Dimension not matching in the multivariate normal sampler";
            throw dimensionsNotMatching();
        }

            arma::mat A;
            arma::vec eigval;
            arma::mat eigvec;
            arma::rowvec res;

            if( arma::chol(A,Sigma) )
            {
                res = randVecNormal(d).t() * A ;
            }
            else
            {
                if( eig_sym(eigval, eigvec, Sigma) )
                {
                    res = (eigvec * arma::diagmat(arma::sqrt(eigval)) * randVecNormal(d)).t();
                }else{
                    Rcout << "randMvNorm failing because of singular Sigma matrix" << '\n';
                    throw negativeDefiniteParameters();
                }
            }

        return res.t() + m;
    }

    arma::mat randMN(const arma::mat &M, const arma::mat &rowCov, const arma::mat &colCov)
    {
        arma::mat C = arma::chol( arma::kron(colCov,rowCov) );
        arma::mat z = randVecNormal( (unsigned int)(M.n_cols * M.n_rows) ).t() * C;
        z.reshape( arma::size(M) );
        return (z + M);
    }


	double randTruncNorm(double m, double sd,double lower, double upper) // Naive, but it'll do for now -- notice now parameters are mean and standard deviation!
	{
		double ret = randNormal(m,sd);

		while( ret < lower || ret > upper)
			ret = randNormal(m,sd);

		return ret;
	}

	arma::uvec randSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const arma::uvec& population, // population to draw from
	    unsigned int sampleSize        // size of each sample
	) // output, sample is a zero-offset indices to selected items, output is the subsampled populaiton.
	{
		arma::uvec samples(sampleSize);

	    int t = 0; // total input records dealt with
	    unsigned int m = 0; // number of items selected so far
	    double u;

	    while (m < sampleSize)
	    {
	        u = randU01(); // call a uniform(0,1) random number generator

	        if ( (populationSize - t)*u >= sampleSize - m )
	        {
	            t++;
	        }
	        else
	        {
	            samples(m) = t;
	            t++; m++;
	        }
	    }

	    return population(samples);
	}

	std::vector<unsigned int> randSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const std::vector<unsigned int>& population, // population to draw from
	    unsigned int sampleSize        // size of each sample
	) // output, sample is a zero-offset indices to selected items, output is the subsampled populaiton.
	{
		std::vector<unsigned int> samplesIndexes(sampleSize);

	    int t = 0; // total input records dealt with
	    unsigned int m = 0; // number of items selected so far
	    double u;

	    while (m < sampleSize)
	    {
	        u = randU01(); // call a uniform(0,1) random number generator

	        if ( (populationSize - t)*u >= sampleSize - m )
	        {
	            t++;
	        }
	        else
	        {
	            samplesIndexes[m] = t;
	            t++; m++;
	        }
	    }

		std::vector<unsigned int> res(sampleSize);
		m = 0;
		for( auto i : samplesIndexes )
		{
			res[m++] = population[i];
		}
		
	    return res;
	}

	// IMPLEMENTATION FROM Efraimidistr and Spirakis 2006
	// (probably efficient when n is close to N rather than in our case, but is there a better alternative?)
	//  even for n=1, we'd still need to compute the cumulative sum of all the weights if we want to use bisection, or what I use below which is O(N) anyway...
	arma::uvec randWeightedSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const arma::vec& weights,	   // probability for each element
	    unsigned int sampleSize,        // size of each sample
	    const arma::uvec& population // population to draw from
	) // sample is a zero-offset indices to selected items, output is the subsampled population.
	{
	    arma::vec score = randVecExponential(populationSize,1.)/weights;
	    arma::uvec result = population( (arma::sort_index(weights,"ascend")) );

	    return result.subvec(0,sampleSize-1);
	}

	// overload with sampleSize equal to one
	arma::uword randWeightedSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const arma::vec& weights,	   // probability for each element
	    const arma::uvec& population // population to draw from
	) // sample is a zero-offset indices to selected items, output is the subsampled population.
	{
	    double u = randU01();
	    double tmp = weights(0);
	    int t = 0;

	    while(u > tmp)
	    {
	    	tmp += weights(++t);
	    }

	    return population( t );
	}


	// Versions that return indexes only
	arma::uvec randWeightedIndexSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const arma::vec& weights,	   // (log) probability for each element
	    unsigned int sampleSize         // size of each sample
	) // sample is a zero-offset indices to selected items, output is the subsampled population.
	{
		// note I can do everything in the log scale as the ordering won't change!
	    arma::vec score = randVecExponential(populationSize,1.) - weights;
	    arma::uvec result = arma::sort_index(score,"ascend");

	    return result.subvec(0,sampleSize-1);
	}

	// Overload with equal weights
	arma::uvec randWeightedIndexSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    unsigned int sampleSize         // size of each sample
	) // sample is a zero-offset indices to selected items, output is the subsampled population.
	{
		// note I can do everything in the log scale as the ordering won't change!
	    arma::vec score = randVecExponential(populationSize,1.);
	    arma::uvec result = arma::sort_index(score,"ascend");
        
	    return result.subvec(0,sampleSize-1);
	}

	// overload with sampleSize equal to one
	arma::uword randWeightedIndexSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const arma::vec& weights     // probability for each element
	) // sample is a zero-offset indices to selected items, output is the subsampled population.
	{
		// note I can do everything in the log scale as the ordering won't change!

	    double u = randU01();
	    double tmp = weights(0);
	    unsigned int t = 0;

	    while(u > tmp)
	    {
	    	// tmp = Utils::logspace_add(tmp,logWeights(++t));
	    	tmp += weights(++t);
	    }

	    return t;
	}



	/// ################### NOW LOG PDFs


	// logPDF rand Weighted Indexes (need to implement the one for the original starting vector?)
	double logPDFWeightedIndexSampleWithoutReplacement(const arma::vec& weights, const arma::uvec& indexes)
	{
		// arma::vec logP_permutation = arma::zeros<arma::vec>((int)std::tgamma(indexes.n_elem+1));  //too big of a vector
		double logP_permutation = 0.; double tmp;

		std::vector<unsigned int> v = arma::conv_to<std::vector<unsigned int>>::from(arma::sort(indexes));
		// vector should be sorted at the beginning.

		arma::uvec current_permutation;
		arma::vec current_weights;

	    do {
	        current_permutation = arma::conv_to<arma::uvec>::from(v);
	        current_weights = weights;
			tmp = 0.;

			while( current_permutation.n_elem > 0 )
			{
			   tmp += log(current_weights(current_permutation(0)));
			   current_permutation.shed_row(0);
			   current_weights = current_weights/arma::sum(current_weights(current_permutation));   // this will gets array weights that do not sum to 1 in total, but will only use relevant elements
    		}

			logP_permutation = Utils::logspace_add(logP_permutation,tmp);

	    } while (std::next_permutation(v.begin(), v.end()));

		return logP_permutation;
	}


	double logPDFIGamma(double x, double a, double b)
	{
		if( x < 0 || b < 0 || a < 0 )
			return -std::numeric_limits<double>::infinity();
		else
			return a*log(b) -std::lgamma(a) + (-a-1.)*log(x) -b/x;
	}

	double logPDFGamma(double x, double a, double b)
	{
		if( x < 0 || b < 0 || a < 0 )
			return -std::numeric_limits<double>::infinity();
		else
			return -a*log(b) -std::lgamma(a) + (a-1.)*log(x) -x/b;
	}


	double logPDFIWishart(const arma::mat& X, double nu, const arma::mat& Sigma)
	{
		unsigned int p = X.n_rows;
		double ret = -0.5*(double)p*nu*log(2.) - lMvGamma(p,nu) - 0.5*arma::trace( Sigma * arma::inv_sympd(X) );
        double tmp = 0.;
        double sign;
		arma::log_det(tmp, sign, X );
		ret += -0.5*( (double)p + nu + 1. )*tmp;

		arma::log_det(tmp, sign, Sigma );
		ret += +0.5*nu*tmp;

		return ret;
	}


	double logPDFMN(const arma::mat& X, const arma::mat& rowCov, const arma::mat colCov)
	{
		unsigned int n = X.n_rows;
		unsigned int m = X.n_cols;

		double ret = -0.5*arma::trace( arma::inv_sympd(colCov) * X.t() * arma::inv_sympd(rowCov) * X ) -
					(double)n*(double)m*0.5*log(2.*M_PI);
        double tmp = 0.;
        double sign;
		arma::log_det(tmp, sign, colCov );
		ret += -0.5*(double)n*tmp;

		arma::log_det(tmp, sign, rowCov );
		ret += -0.5*(double)m*tmp;

		return ret;
	}


 	double logPDFNormal(const double& x, const double& m,const  double& sigmaSquare)
	{

		return -0.5*log(2.*M_PI) -0.5*log(sigmaSquare) -(0.5/sigmaSquare)*arma::as_scalar( pow(x-m,2) );

	}


	double logPDFNormal(const arma::vec& x, const  arma::mat& Sigma)  // zeroMean
	{
		unsigned int k = Sigma.n_cols;

		double sign, tmp;
		arma::log_det(tmp, sign, Sigma ); //sign is not importantas det SHOULD be > 0 as for positive definiteness!

		return -0.5*(double)k*log(2.*M_PI) -0.5*tmp -0.5* arma::as_scalar( (x).t() * arma::inv_sympd(Sigma) * (x) );

	}


	double logPDFNormal(const arma::vec& x, const arma::vec& m,const  arma::mat& Sigma)
	{
		unsigned int k = Sigma.n_cols;

		double sign, tmp;
		arma::log_det(tmp, sign, Sigma ); //sign is not importantas det SHOULD be > 0 as for positive definiteness!

		return -0.5*(double)k*log(2.*M_PI) -0.5*tmp -0.5* arma::as_scalar( (x-m).t() * arma::inv_sympd(Sigma) * (x-m) );

	}

	double logPDFNormal(const arma::vec& x, const arma::vec& m, const double& Sigma)
	{
		//this is more a log likelihood rather than logPDF here, since the input vector is indep realisations with same sigma and (possibly) different means
		// we rely on amradillo for parallelisation wrt to individuals
		unsigned int n = x.n_elem;

		return -0.5*(double)n*log(2.*M_PI) -0.5*n*log(Sigma) -0.5/Sigma * arma::as_scalar( (x-m).t() * (x-m) );

	}

	double logPDFNormal(arma::vec& x, arma::vec& m, const arma::mat& rowCov ,const arma::mat& colCov)   // vectorised version of a matrix normal
	{
		unsigned int k = rowCov.n_rows;
		unsigned int d = colCov.n_rows;

		double logP = -0.5*(double)k*log(2.*M_PI) - 0.5 * arma::as_scalar( ( (x-m).t() * arma::inv_sympd( arma::kron( colCov , rowCov ) ) * (x-m) ) );

		double sign, tmp;
		arma::log_det(tmp, sign, rowCov );
		logP += -0.5*(double)d*tmp;

		arma::log_det(tmp, sign, colCov );
		logP += -0.5*(double)k*tmp;

		return logP;
	}

	double lBeta(double a,double b){    //log beta function
		return std::lgamma(a) + std::lgamma(b) - std::lgamma(a+b);
	}

	double logPDFBeta(double x, double a, double b)
	{
		if( x <= 0. || x >= 1. )
			return -std::numeric_limits<double>::infinity();
		else
			return -lBeta(a,b) + (a-1)*log(x) + (b-1)*log(1-x);
	}

	double logPDFBernoulli(unsigned int x, double pi)
	{
		if( x > 1 ) // remember x is UNSIGNED int here
			return -std::numeric_limits<double>::infinity();
		else
			return x*log(pi) + (1-x)*log(1.-pi);
	}

	double logPDFBernoulli(const arma::uvec& x, double pi)
	{
		double n = (double)x.n_elem;
		double k = arma::sum(x);

		return k*log(pi) + (n-k)*log(1.-pi);
	}

	// could be int but we'd need to cast them to doubles anyway...
	double log_binomial_coefficient(double n, double k) { return -log(n+1)-lBeta(n-k+1,k+1); }

	double logPDFBinomial(unsigned int k, unsigned int n, double pi)
	{
		if ( n < k )
			return -std::numeric_limits<double>::infinity();
		else
			return log_binomial_coefficient(n, k) + k*log(pi) + (n-k)*log(1.-pi) ;
			// return log(boost::math::binomial_coefficient<double>(n, k) ) + k*log(pi) + (n-k)*log(1.-pi) ; // need double template https://www.boost.org/doc/libs/1_69_0/libs/math/doc/html/math_toolkit/factorials/sf_binomial.html
	}

	double CDFNormal(double x, double m, double sd)
	{
		return 0.5 * std::erfc(-((x-m)/sd) * M_SQRT1_2); // ... right?
	}

	// BELOW not needed and introduces extra dependency on
	//#include <boost/math/special_functions/erf.hpp>
	// hence excluded
	// double invCDFNormal(double x, double m, double sd)
	// {
	// 	return sqrt(2.) * boost::math::erf_inv(2.*((x-m)/sd)-1.);; // ... right?
	// }


	double logPDFTruncNorm(double x, double m, double sd, double lower, double upper)
	{
		double truncNormConst = log( CDFNormal(upper,m,sd) - CDFNormal(lower,m,sd) );
		return -log(sqrt(2.*M_PI)) -log(sd) -0.5*(x-m)*(x-m)/(sd*sd) - truncNormConst;
	}

	double lMvGamma(unsigned int n, double a)					// Multivariate GAMMA FUNCTION! NOT THE PDF/CDF
	{
		double lmvG = 0;
		for(unsigned int j=0; j<n ; ++j)
		{
			lmvG += std::lgamma( a + 0.5*(1.-(double)j+1.) );   //last +1 cause indexes start from 0
		}

		return (n*(n-1.)*0.25)*log(M_PI) + lmvG;
	}
    
    /*double logPDFt( double x, double d )
    {
        
        return std::lgamma( (d+1.)/2. ) - 0.5*log(d*M_PI) - std::lgamma( d/2. ) - (d+1.)/2.*log( 1+x*x/d );
        
    }*/



}
