#ifndef PARAM_TYPES_H
#define PARAM_TYPES_H

#ifdef CCODE
	#include <iostream>
#endif

/******************************************************
 * Parameter types, to be used in all derived classes 
 * ****************************************************/

enum class Covariance_Type {
    HIW=1, IW, IG
};

enum class Gamma_Type {
    hotspot=1, hierarchical, mrf
};

enum class Beta_Type {
    independent=1, gprior, reGroup
};

enum class Gamma_Sampler_Type {
    bandit=1, mc3
};

class Bad_Covariance_Type : public std::exception{

  public:
    Bad_Covariance_Type( Covariance_Type ct ): ct_(ct) {}

    virtual const char* what() const throw()
    {
      switch ( ct_ )
      {
        case Covariance_Type::IW :
          return "The DENSE (IW) COVARIANCE type is not valid here";
      
        case Covariance_Type::HIW :
          return "The SPARSE (HIW) COVARIANCE type is not valid here";
      
        case Covariance_Type::IG :
          return "The INDEPENDENT (IG) COVARIANCES type is not valid here";

        default:
          return "The covariance type here is not valid -- unknown type";
      }
    }

  private:
    Covariance_Type ct_;

};

class Bad_Gamma_Type : public std::exception{

  public:
    Bad_Gamma_Type( Gamma_Type gt ): gt_(gt){}

    virtual const char* what() const throw()
    {
      switch ( gt_ )
      {
        case Gamma_Type::hotspot :
          return "The HOTSPOT GAMMA type is not valid here";
      
        case Gamma_Type::hierarchical :
          return "The simple HIERARCHICAL GAMMA type is not valid here";
      
        case Gamma_Type::mrf :
          return "The Markov Random Field GAMMA type is not valid here";
      
        default:
          return "The gamma type here is not valid -- unknown type";
      }
    }

  private:
    Gamma_Type gt_;

};

class Bad_Beta_Type : public std::exception{

  public:
    Bad_Beta_Type( Beta_Type bt ): bt_(bt){}

    virtual const char* what() const throw()
    {
      switch ( bt_ )
      {
        case Beta_Type::independent :
          return "The INDEPENDENT BETA type is not valid here";
      
        case Beta_Type::gprior :
          return "The GPRIOR BETA type is not valid here";

        case Beta_Type::reGroup :
        return "The REGROUP BETA type is not valid here";
      
        default:
          return "The beta type here is not valid -- unknown type";
      }
    }

  private:
    Beta_Type bt_;

};

class Bad_Gamma_Sampler_Type : public std::exception{

  public:
    Bad_Gamma_Sampler_Type( Gamma_Sampler_Type gst ): gst_(gst){}

    virtual const char* what() const throw()
    {
      switch ( gst_ )
      {
        case Gamma_Sampler_Type::bandit :
          return "The BANDIT GAMMA SAMPLER type is not valid here";
      
        case Gamma_Sampler_Type::mc3 :
          return "The MC3 GAMMA SAMPLER type is not valid here";
      
        default:
          return "The GAMMA SAMPLER type here is not valid -- unknown type";
      }
    }

  private:
    Gamma_Sampler_Type gst_;
    
};

#endif
