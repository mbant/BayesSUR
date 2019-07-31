# BayesSUR

This is WIP repository for an R package interface to the [Bayesian SSUR](github.com/mbant/Bayesian_SSUR) C++ code.
Please use only the C++ code for the moment.

## R package instructions
see `test.R` file for usage; see the package vignette `BayesSUR.pdf` (submitted to the Journal of Statistical Software) for more information.

## Update of the version `BayesSUR_0.1.16.tar.gz`:

# the function plotEstimator() can choose to plot beta_hat and gamma_hat with the labeled axes rather than numbers only, though the arguments name.responses and name.predictors;

# the function plotEstimator() prints the legend of gamma_hat with fixed bar range [0,1];

3) the function plotManhattan() can show axis labels by giving the argument axis.label and can also mark the corresponding responses in the first Manhattan-like plot (mPIP);

4) some "hard coding" in the function plotMCMCdiag() has been fixed;

5) the function runSUR() allows to give the data.matrix and data.frame type of dataset;

6) the arguments outFilePath and tmpFolder in the runSUR() have been fixed hopefully so that they can automatically identify the given relative or absolute directory path no matter in WinOS or MacOS.


Remaining issues:

1) The runSUR() doesn’t print the running iterations information in the RGui of WinOS, but it works in the RStudio of WinOS, RStudio of MacOS and R Console of MacOS;

2) ….
