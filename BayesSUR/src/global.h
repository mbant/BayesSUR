#ifndef GLOBAL
#define GLOBAL

  #include <armadillo>  // include this globally
  #include <random>
  #include <vector>

  #ifdef _OPENMP
    #include <omp.h>
  #endif

  #ifndef CCODE
    #include <Rcpp.h>
  #endif

#endif