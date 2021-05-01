[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/BayesSUR)](https://cran.r-project.org/package=BayesSUR)

# BayesSUR

This repository contains a [new and improved](https://github.com/mbant/BayesSUR/blob/master/BayesSUR/inst/doc/BayesSUR.pdf) [R]() package for high-dimensional multivariate Bayesian variable and covariance selection in linear regression, started as an interface to the [Bayesian SSUR](https://github.com/mbant/Bayesian_SSUR) C++-only, UNIX-specific, code.

## R package instructions
See `test.R` file for usage; see the package vignette [`BayesSUR.pdf`](https://github.com/mbant/BayesSUR/blob/master/BayesSUR/inst/doc/BayesSUR.pdf) for more information.

## Installation

Install the current release with
``` r
install.packages("BayesSUR")
```

Install the current development version with
``` r
library("devtools")
devtools::install_github("mbant/BayesSUR/BayesSUR")
```

## Citation

Zhi Zhao, Marco Banterle, Leonardo Bottolo, Sylvia Richardson, Alex Lewin, Manuela Zucknick (2021). BayesSUR: An R package for high-dimensional multivariate Bayesian variable and covariance selection in linear regression. **Journal of Statistical Software**, to appear. arXiv: [2104.14008](https://arxiv.org/abs/2104.14008).
