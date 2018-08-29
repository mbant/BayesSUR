#include "global.h"

#include <omp.h>
#include <vector>
#include <random>

omp_lock_t RNGlock; // GLOBAL NAMESPACE CAUSE I NEED IT NO MATTER WHERE 
    //use with 
    // omp_set_lock(&RNGlock);
    // omp_unset_lock(&RNGlock);

std::vector<std::mt19937_64> rng;
