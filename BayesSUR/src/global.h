#ifndef GLOBAL
#define GLOBAL

#ifdef _OPENMP
  #include <omp.h>
#endif

#ifndef CCODE
  #include <Rcpp.h>
#endif

#include <vector>
#include <random>
#include <armadillo> // include this globally
#include <iostream>

std::ostream& console_out();
std::ostream& err_out();

#endif
