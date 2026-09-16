This repository contains the data and scripts necessary to reproduce all analyses for:

Paper title and authors

Preprint here: Link

To reproduce analyses, an additional folder "data" has to be downloaded from Zenodo (LINK), and placed into the same working directory as the remaining folders.

Software requirements

The analysis code was written in R (v4.5.1.) Statistical models are fit using the Stan MCMC engine (v2.32.2) and cmdstanr (v0.9.0.9000) packages, which require a C++ compiler. 
Installation instructions are available at https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started and https://mc-stan.org/cmdstanr/articles/cmdstanr.html. 
See also the Stan user guide at https://mc-stan.org/users/documentation.

All other required R packages and their versions can be found in "functions/used-packages.R"
