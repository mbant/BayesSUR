#include "junction_tree.h"

#include <vector>
#include <iostream>
#include <string>
#include <armadillo>
#include <cmath>
#include <limits>
#include <omp.h>

// #include "global.h"
// #include "utils.h"
// #include "distr.h"

extern omp_lock_t RNGlock; //defined in global.h
extern std::vector<std::mt19937_64> rng;


int main()
{
    arma::sp_umat A(3,3);
    A.print();

    A(1,1) = 1; A(0,2) = 1; A(2,0) = 1;
    A.print();

    A = arma::trimatl(A);
    A.print();

    std::cout << "Create a JT .. " <<std::endl;
    JunctionTree jt( 10 );

    std::cout << ".. Done! Now Print it .." << std::endl;
    jt.print();

    std::cout << std::flush << std::endl << std::endl << 
        "Assign-Reference a new JT .. " <<std::endl;
    JunctionTree jt2 = jt;   // this should point to the same components and modifying one should modify the other

    std::cout << ".. Done! Now Print it .." << std::endl;
    jt2.print();

    std::cout << std::flush << std::endl << std::endl << 
        "Copy into a new JT .. " <<std::endl;
    JunctionTree jt3; //use with care as this is just a defult constructor
    jt.copyJT( jt3 ); //  this should copy all components and modifying one should NOT modify the other

    std::cout << ".. Done! Now Print it .." << std::endl;
    jt3.print();

    std::cout << "Create a FULL JT .. " <<std::endl;
    JunctionTree jt4( 10 , "full" );

    std::cout << ".. Done! Now Print it .." << std::endl;
    jt4.print();

    // Init Random Device  ####################################
    omp_init_lock(&RNGlock);  // RNG lock for the parallel part

    std::random_device r;
	unsigned int nThreads = omp_get_max_threads();
	rng.reserve(nThreads);  // reserve the correct space for the vector of rng engines
	std::seed_seq seedSeq;	// and declare the seedSequence
	std::vector<unsigned int> seedInit(8);

	// seed all the engines
	for(unsigned int i=0; i<nThreads; ++i)
	{
		rng[i].seed(r());
	}
    // #################################### End init

    std::cout << std::flush << std::endl << std::endl << 
        "Update it a few times into a new JT .. " <<std::endl;
    std::pair<bool,double> updated;
    unsigned int i=0;

    while( ++i < 20 )
    {
        jt2.copyJT( jt3 );
        updated = jt3.propose_single_edge_update( );

        if( jt3.perfectCliqueSequence[0]->getSeparator().size()>0 ||
            (jt3.perfectCliqueSequence[0]->getParent()) ||
            jt3.perfectEliminationOrder.size() < 10 )
            {
                std::cout << std::endl << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  -----  " << i << std::endl << std::endl;
                std::cout << std::endl << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  -----  " << i << std::endl << std::endl;
                std::cout << std::endl << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  -----  " << i << std::endl << std::endl;
                break;
            }

        if( Distributions::randU01() < std::get<1>(updated) )
        {
            jt2 = jt3;
        }
    }

    std::cout << ".. Done! Now Print it .." << std::endl;
    jt2.print();

    
    std::cout << std::flush << std::endl << std::endl << 
        "Now copy this newer more complex one into a new JT .. " <<std::endl;
    jt2.copyJT( jt3 ); //  this should copy all components and modifying one should NOT modify the other

    std::cout << ".. Done! Now Print it .." << std::endl;
    jt3.print();


    arma::uvec w = {0,0,3,3,1,3,2,1};
    std::cout << arma::index_max(w) <<std::endl;


    return 0;
}