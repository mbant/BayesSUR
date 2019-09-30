# BayesSUR

This repository contains a [new and improved](https://github.com/mbant/BayesSUR/blob/master/BayesSUR/vignettes/vignettes.pdf) [R]() package for high-dimensional multivariate Bayesian variable and covariance selection in linear regression, started as an interface to the [Bayesian SSUR](https://github.com/mbant/Bayesian_SSUR) C++-only, UNIX-specific, code.

## R package instructions
See `test.R` file for usage; see the package vignette `BayesSUR.pdf` (submitted to the Journal of Statistical Software) for more information.

## Update

### New in version `BayesSUR_1.0-1.tar.gz` (30 September 2019):

The function `plotCPO()` allows different cut-off for each outcome.

### Version `BayesSUR_1.0-0.tar.gz` (27 September 2019):

Improved help files.

### Version `BayesSUR_0.1.24.tar.gz` (27 September 2019):

Fixed bugs and improved help files.

### Version `BayesSUR_0.1.23.tar.gz` (27 September 2019):

Fixed bugs and improved help files.

### Version `BayesSUR_0.1.22.tar.gz` (26 September 2019):

Fixed bugs for the posterior predictives of the HRR models.

### Version `BayesSUR_0.1.21.tar.gz` (23 September 2019):

1) Fixed the CPO, elpd (lpd and WAIC) of the HRR models.
2) Modified the function `summary()`.
3) Added the functions `print()`, `coef()`, `fitted()`, `predict()`.

### Version `BayesSUR_0.1.20.tar.gz` (22 September 2019):

Removed last warning about strand std::cou.

### Version `BayesSUR_0.1.19.tar.gz` (27 August 2019):

Fixing the warnings from the `devtools::check("BayesSUR")`......

### Version `BayesSUR_0.1.18.tar.gz` (20 August 2019):

1) The function plotCPO() has been updated to check outliers overall response variables alternatively.

2) The model prediction accuracy is measured by the function elpd() (expected log pointwise predictive density), which has two methods to estimate elpd, i.e., the LOO (leave-one-out cross- validation) and WAIC (widely applicable information criterion
).

### Version `BayesSUR_0.1.17.tar.gz` (13 August 2019):

1) The function plotCPO() uses posterior predictive (i.e., conditional predictive ordinate) to check outliers.

2) The function plotEstimator() can extract the CPOs matrix through argument `estimator="CPO"`.


### Version `BayesSUR_0.1.16.tar.gz` (30 July 2019):

1) The function plotEstimator() can choose to plot beta_hat and gamma_hat with the labeled axes rather than numbers only, though the arguments name.responses and name.predictors.


2) The function plotEstimator() prints the legend of gamma_hat with fixed bar range [0,1].

3) The function plotManhattan() can show axis labels by giving the argument axis.label and can also mark the corresponding responses in the first Manhattan-like plot (mPIP).

4) Some "hard coding" in the function plotMCMCdiag() has been fixed.

5) The function runSUR() allows to give the data.matrix and data.frame type of dataset.

6) The arguments outFilePath and tmpFolder in the runSUR() have been fixed hopefully so that they can automatically identify the given relative or absolute directory path no matter in WinOS or MacOS.


Remaining issues:

1) The runSUR() doesn’t print the running iterations information in the RGui of WinOS, but it works in the RStudio of WinOS, RStudio of MacOS and R Console of MacOS.

### Version `BayesSUR_0.1.13.tar.gz` (30 June 2019)# BayesSUR

### Libraries

[Armadillo](http://arma.sourceforge.net/) and its respective [Rcpp interface](https://github.com/RcppCore/RcppArmadillo)

[PugiXML](http://pugixml.org/)

[Boost](www.boost.org) and [BoostHeaders](https://github.com/eddelbuettel/bh)
