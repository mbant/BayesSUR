# BayesSUR 2.2-1

* Fixed a bug about `SUR_Chain.cpp` introduced in version 2.2-0
* Fixed argument mrfG index issue in function `BayesSUR()`

# BayesSUR 2.2-0

* Fixed a bug about the use of temperature parameter in `HRR_Chain.cpp`, also minor update in `SUR_Chain.cpp`
* Fixed random effects' sampler `HRR_Chain::stepW0()`
* Added argument `beta.type` in function `plotEstimator()` to plot MPM coefficient estimates
* Fixed argument `mrfG` index issue in function `BayesSUR()`

# BayesSUR 2.1-7

* Increased the threshold of predictor dimension to 100000 (previously 5000) when pre-computing XtX
* Minor changes of the parameter updates of the HRR models, especially in function `step()`

# BayesSUR 2.1-6

* Fixed issues of knitr VignetteBuilder and README

# BayesSUR 2.1-5

* Updated citation, added README and one more vignette

# BayesSUR 2.1-4

* Fixed a documentation issue by adding `@alias BayesSUR-package`

# BayesSUR 2.1-3

* Updated C++11 specification to current default of C++17
* Cleaned R scripts

# BayesSUR 2.1-2

* Fixed gcc-UBSAN error `reference binding to null pointer` in SUR_Chain.cpp:3818

# BayesSUR 2.1-1

* Fixed R session aborted fatal error

# BayesSUR 2.1-0

* Fixed an issue with `omp_set_num_threads` that did not work

# BayesSUR 2.0-1

* Code update

# BayesSUR 2.0-0

* Major update, introducing random effects for mandatory variables without variable selection.
