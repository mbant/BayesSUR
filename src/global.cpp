#include "global.h"

#ifdef _OPENMP
  omp_lock_t RNGlock; // GLOBAL NAMESPACE CAUSE I NEED IT NO MATTER WHERE
#endif
  //use with 
  // omp_set_lock(&RNGlock);
  // omp_unset_lock(&RNGlock);
//std::vector<std::mt19937_64> rng;
//Rcpp::RNGScope scope;
