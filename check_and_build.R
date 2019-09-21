devtools::check("BayesSUR/",cran = TRUE,build_args=c('--compact-vignettes=both'))
devtools::build("BayesSUR/",args=c('--compact-vignettes=both'))
