# BayesSUR

This repository contains a [new and improved](https://github.com/mbant/BayesSUR/blob/master/BayesSUR/vignettes/vignettes.pdf) [R]() package for high-dimensional multivariate Bayesian variable and covariance selection in linear regression, started as an interface to the [Bayesian SSUR](https://github.com/mbant/Bayesian_SSUR) C++-only, UNIX-specific, code.

## R package instructions
See `test.R` file for usage; see the package vignette `BayesSUR.pdf` for more information.

## Installation

To install the current release, use
``` r
install.packages("BayesSUR")
```

To install the current development version, use
``` r
library("devtools")
devtools::install_github("mbant/BayesSUR")
```

## Update

### New in verion `BayesSUR_1.1-2.tar.gz` (18 April 2020):

Fixed the issue of the hyperparameter w_j initialization in the HRR model.

### New in verion `BayesSUR_1.1-1.tar.gz` (16 April 2020):

1) Use of `Rcpp::RNGScope scope` for the reproducibility of results.
2) Generic function for plotting of R objects.

### Version `BayesSUR_1.1-0.tar.gz` (02 April 2020):

Random effects have been included in the SUR mdoels.

### Version `BayesSUR_1.0-4.tar.gz` (17 December 2019):

1) An argument `maxThreads` was added in the main function `BayesSUR()` to allow users to specify the maximum number of threads for parallelisation.

2) One issue in `plotEstimator()` was fixed, which didn't print the estimated Gamma matrix correctly when using random effects in previous version.

3) More warning information can be printed by function `plotMCMCdiag()` if the user specifies very small MCMC iterations.

4) The name of dataset `example_GDSC_target.rda` was changed to `targetGene.rda`. The corresponding help file `example_GDSC_target.R` was also changed to `targetGene.R`. 

### Version `BayesSUR_1.0-3.tar.gz` (08 December 2019):

1) In the file `DESCRIPTION`, Waldir was added as a contributor, RoxygenNote version was increased and `gRbase` was removed from the list of suggested packages.

2) Function `plotMCMCdiag()` was modified to allow any more than one MCMC iteration.

3) Deprecated awk `-e` flag in file `configure.ac` was removed.

4) The vignette-generation procedure was changed. Instead of having an Rnw file calling the pdf, now there is a `.pdf.asis` file containing the proper vignette engine LaTeX calls. This change puts the vignette-generation procedure in accordance with the [R.rsp package manual](https://cran.r-project.org/web/packages/R.rsp/vignettes/R_packages-Static_PDF_and_HTML_vignettes.pdf).

5) The example code in `example_GDSC.R`, `example_GDSC_targets.R` and `example_eQTL.R` were wrapped in `\dontrun{}`, because `\donttest{}` was yielding errors on Ubuntu 18.04 (Linux 5.0.0) and Manjaro (5.3.12).

6) Code in the files mentioned on the previous point was reformatted in order to reduce line width to less than 100 characters.

7) `.gitignore` was modified to ignore Rcheck files, tarballs and log files.

### version `BayesSUR_1.0-2.tar.gz` (29 October 2019):

Fixed log(2) and similar issues in the file `junction_tree.cpp` for the solaris-x86 check. 

### Version `BayesSUR_1.0-1.tar.gz` (23 October 2019):

1) Updated the vignettes, especailly adding an appendix for the elpd.
2) Fixed the issues from the clang compling with AddressSanitizer.
3) Added an argument `output_model_visit` in the main function `BayesSUR()` to print all visited models (indices of nonzero gamma and response graph) of MCMC iterations.

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
