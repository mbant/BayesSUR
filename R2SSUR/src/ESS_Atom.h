#ifndef ESS_ATOM_H
#define ESS_ATOM_H

#include <iostream>  // for std::cout
#include <string>
#include <vector>
#include <memory>

#include "utils.h"
#include "distr.h"
#include "junction_tree.h"

/******************************************************
 * Parameter types, to be used in all derived classes 
 * ****************************************************/

class Bad_Covariance_Type : public std::exception{

  virtual const char* what() const throw()
  {
    return "The covariance type here is not valid";
  }

};

class Bad_Gamma_Type : public std::exception{

  virtual const char* what() const throw()
  {
    return "The Gamma type here is not valid";
  }

};

class Bad_Beta_Type : public std::exception{

  virtual const char* what() const throw()
  {
    return "The Beta type here is not valid";
  }

};

class Bad_Gamma_Sampler_Type : public std::exception{

  virtual const char* what() const throw()
  {
    return "The Gamma_Sampler type here is not valid";
  }

};

enum class Covariance_Type {
    sparse=1, dense
};

enum class Gamma_Type {
    hotspot=1, hierarchical, mrf
};

enum class Beta_Type {
    independent=1, gprior
};

enum class Gamma_Sampler_Type {
    bandit=1, mc3
};

/************************************
 * Interface class that works with
 * the ESS_Sampler class
 * CRTP used for global exchanges
 ***********************************/

class ESS_Base {}; // non-templated base-class needed for static assert checks

template<typename T>
class ESS_Atom : public ESS_Base
{

    public:

        // NEED TO IMPLEMENT A CONSTRUCTOR AS BELOW
        // ESS_Atom( Utils::SUR_Data& surData , double temperature_ );

        virtual void step() = 0;

        virtual int globalStep( std::shared_ptr<T>& ) = 0;

        virtual double logLikelihood() = 0;

        virtual double getLogLikelihood() const = 0;
        virtual double getJointLogPrior() const = 0;
        virtual double getJointLogPosterior() const = 0;

};

#endif