################################################################################
#
# Title: 3. Run full model
#
# Author F Dhellemmes
#
# last update: 03/08/2026 (DD/MM/YYYY)
#
################################################################################

library(cmdstanr)
library(posterior)
library(bayesplot)
library(data.table)

# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################

stan_data<-readRDS(here::here("data","3_spot_leaving", "stan_spot_leaving.RDS"))

stan_file_full <- here::here("3_spot_leaving","stan","spot_leaving_full.stan")
mod_full <- cmdstan_model(stan_file_full, cpp_options = list(stan_threads = TRUE))

fit_full <- mod_full$sample(
  data = stan_data,
  seed = 42,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  parallel_chains = 4,
  threads_per_chain = round(32/4)
)

fit_full$save_object(file =  here::here("data","3_spot_leaving","model_fit","spot_leaving_full_fit.RDS"))

fit_full<-readRDS(here::here("data","3_spot_leaving","model_fit","spot_leaving_full_fit.RDS"))

params_full <- as_draws_df(fit_full)

beta_fit_full_draws <- subset_draws(params_full, variable = "beta") 

color_scheme_set("teal")#
bayesplot_theme_set(theme_default())
library(ggplot2)

dplyr::bind_cols(summarise_draws(beta_fit_full_draws),c("Int catch 0",
                                                     "Int catch 1",
                                                     "Fish",
                                                     "time since event",
                                                     "distance sonar",
                                                     "social",
                                                     "success",
                                                     "loss"))

bayesplot::mcmc_areas(params_full,regex_pars = "beta\\[[1-8]\\]", prob = 0.8)+ scale_y_discrete(
  labels = c("beta[1]" = "Int catch 0",
             "beta[2]" = "Int catch 1",
             "beta[3]" = "Fish",
             "beta[4]" = "Time since event",
             "beta[5]" = "distance sonar",
             "beta[6]" = "social",
             "beta[7]" = "success",
             "beta[8]" = "loss"))


################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################
