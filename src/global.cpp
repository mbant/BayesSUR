#include "global.h"

#include <omp.h>
#include <vector>
#include <random>

omp_lock_t RNGlock; // GLOBAL NAMESPACE CAUSE I NEED IT NO MATTER WHERE
std::vector<std::mt19937_64> rng;
