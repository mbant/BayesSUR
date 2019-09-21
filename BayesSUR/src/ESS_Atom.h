#ifndef ESS_ATOM_H
#define ESS_ATOM_H

#ifdef CCODE
	#include <iostream>
#endif

#include <string>
#include <vector>
#include <memory>

#include "utils.h"
#include "distr.h"
#include "junction_tree.h"

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
