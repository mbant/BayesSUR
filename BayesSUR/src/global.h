#ifndef GLOBAL
#define GLOBAL

  #include <random>
  #include <vector>
  
  // to get std::beta
  #define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1


  #ifdef _OPENMP
    #include <omp.h>
  #endif

  #ifndef CCODE
    #include <RcppArmadillo.h>
  #else
    #include <armadillo>
  #endif

#endif
