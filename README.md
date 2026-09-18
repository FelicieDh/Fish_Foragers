# Code repository

This repository contains the scripts necessary to reproduce all analyses for:

Dhellemmes F., Schakowski A., Chirkov V., Kavelaars M.M., Korzilius F.F.A., Deffner D., Romanczuk P., Kortet R., Pihlasvaara P., Kurvers R.H.J.M. (Preprint 2026)
Tracking human foragers and their prey reveals adaptive predator–prey dynamics

Preprint here: Link

To reproduce analyses, an additional folder "data" has to be downloaded from Zenodo (https://doi.org/10.5281/zenodo.22798712), and placed into the same working directory as the remaining folders.

Your working directory should contain the following folders to reproduce the analysis:


```
└── Fish_foragers/
    ├── 0_descriptive_analysis/
    ├── 1_step_selection/
    ├── 2_ecological_validation_step_selection/
    ├── 3_spot_leaving/
    ├── 4_ecological_validation_spot_leaving/
    ├── 5_angler_effects/
    ├── 6_recovery/
    ├── 7_revisit_analysis/
    ├── data/
    └── functions/
```

# Software requirements

The analysis code was written in R (v4.5.1.) Statistical models are fit using the Stan MCMC engine (v2.32.2) and cmdstanr (v0.9.0.9000) packages, which require a C++ compiler. 
Installation instructions are available at https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started and https://mc-stan.org/cmdstanr/articles/cmdstanr.html. 
See also the Stan user guide at https://mc-stan.org/users/documentation.

All other required R packages and their versions can be found in "functions/used-packages.R"
