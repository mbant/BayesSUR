#ifndef DISTR
#define DISTR

#include <armadillo>
#include <tgmath.h>
#include <boost/math/special_functions/erf.hpp> // can I do this?
#include <random>

namespace Distributions{

	double randExponential(const double lambda);
	arma::vec randExponential(const unsigned int n, const double lambda);

	unsigned int randBinomial(const unsigned int n, const double p);
	arma::uvec randMultinomial(unsigned int n, const arma::vec prob);

	double randNormal(const double m, const double sigmaSquare); // random normal interface to arma::randn
	arma::vec randNormal(const int n, const double m, const double sigmaSquare);
	arma::vec randMvNormal(const arma::vec &m, const arma::mat &Sigma);

	double randT(const double nu);
	arma::vec randT(const unsigned int n, const double nu);
	arma::vec randMvT(const double &nu, const arma::vec &m, const arma::mat &Sigma);


	arma::mat randWishart(double df, const arma::mat& S);
	arma::mat randIWishart(double df, const arma::mat& S);
	arma::mat randMN(const arma::mat &M, const arma::mat &rowCov, const arma::mat &colCov);
	double randBeta(double a, double b);
	unsigned int randBernoulli(double pi);
	double randTruncNorm(double m, double sd,double lower, double upper);
	double randU01();
	double randLogU01();
	int randIntUniform(const int a,const int b);
	arma::ivec randIntUniform(const unsigned int n, const int a,const int b);

	double randIGamma(double a, double b);

	double logPDFIWishart(const arma::mat& M, double nu, const arma::mat& Sigma);
	double logPDFMN(const arma::mat& X, const arma::mat& rowCov, const arma::mat colCov); //assuming zero mean

	double logPDFNormal(const double& x, const double& m, const double& sigmaSquare);
	double logPDFNormal(const arma::vec& x, const  arma::mat& Sigma); //zero mean
	double logPDFNormal(const arma::vec& x, const arma::vec& m, const arma::mat& Sigma);
	double logPDFNormal(const arma::vec& x, const arma::vec& m,const  double& Sigma);
	double logPDFNormal(arma::vec& x, arma::vec& m, const arma::mat& rowCov, const arma::mat& colCov);

	double lBeta(double a,double b);
	double logPDFBeta(double x, double a, double b);
	double logPDFBernoulli(unsigned int x, double pi);
	double logPDFTruncNorm(double x, double m, double sd, double lower, double upper);
	double logPDFGamma(double x, double a, double b);
	double logPDFIGamma(double x, double a, double b);

	double lMvGamma(unsigned int n, double a);

	double CDFNormal(double x, double m=0., double sd=1.);
	double invCDFNormal(double x, double m=0., double sd=1.);

	arma::uvec randSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const arma::uvec& population, // population to draw from
	    unsigned int sampleSize        // size of each sample
	); // output, sample is a zero-offset indices to selected items, output is the subsampled populaiton.

	std::vector<unsigned int> randSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const std::vector<unsigned int>& population, // population to draw from
	    unsigned int sampleSize        // size of each sample
	); // output, sample is a zero-offset indices to selected items, output is the subsampled populaiton.

	arma::uvec randWeightedSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const arma::vec& weights,	   // probability for each element
	    unsigned int sampleSize,        // size of each sample
	    const arma::uvec& population // population to draw from
	); // sample is a zero-offset indices to selected items, output is the subsampled population.


	// overload with sampleSize equal to one
	arma::uword randWeightedSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const arma::vec& weights,	   // probability for each element
	    const arma::uvec& population // population to draw from
	); // sample is a zero-offset indices to selected items, output is the subsampled population.


	// Versions that return indexes only
	arma::uvec randWeightedIndexSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const arma::vec& weights,	   // probability for each element
	    unsigned int sampleSize        // size of each sample
	); // sample is a zero-offset indices to selected items, output is the subsampled population.

	// overload with equal weights
	arma::uvec randWeightedIndexSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    unsigned int sampleSize         // size of each sample
	); // sample is a zero-offset indices to selected items, output is the subsampled population.

	// overload with sampleSize equal to one
	arma::uword randWeightedIndexSampleWithoutReplacement
	(
	    unsigned int populationSize,    // size of set sampling from
	    const arma::vec& weights   // probability for each element
	); // sample is a zero-offset indices to selected items, output is the subsampled population.


	// logPDF rand Weighted Indexes (need to implement the one for the original starting vector?)
	double logPDFWeightedIndexSampleWithoutReplacement(const arma::vec& weights, const arma::uvec& indexes);


}



#endif
