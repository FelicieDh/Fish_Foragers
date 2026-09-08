################################################################################
#
# Title: 4. Run consecutive models
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

stan_file_nosocial <- here::here("3_spot_leaving","stan","spot_leaving_nosocial.stan")
mod_nosocial <- cmdstan_model(stan_file_nosocial, cpp_options = list(stan_threads = TRUE))

stan_file_nopersonal <- here::here("3_spot_leaving","stan","spot_leaving_nopersonal.stan")
mod_nopersonal <- cmdstan_model(stan_file_nopersonal, cpp_options = list(stan_threads = TRUE))

stan_file_nocatch <- here::here("3_spot_leaving","stan","spot_leaving_nocatch.stan")
mod_nocatch <- cmdstan_model(stan_file_nocatch, cpp_options = list(stan_threads = TRUE))


fit_nosocial <- mod_nosocial$sample(
  data = stan_data,
  seed = 42,
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  parallel_chains = 4,
  threads_per_chain = round(32/4)
)

fit_nosocial$save_object(file =  here::here("data","3_spot_leaving","model_fit","spot_leaving_nosocial_fit.RDS"))

fit_nopersonal <- mod_nopersonal$sample(
  data = stan_data,
  seed = 42,
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  parallel_chains = 4,
  threads_per_chain = round(32/4)
)

fit_nopersonal$save_object(file =  here::here("data","3_spot_leaving","model_fit","spot_leaving_nopersonal_fit.RDS"))

fit_nocatch <- mod_nocatch$sample(
  data = stan_data,
  seed = 42,
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  parallel_chains = 4,
  threads_per_chain = round(32/4)
)

fit_nocatch$save_object(file =  here::here("data","3_spot_leaving","model_fit","spot_leaving_nocatch_fit.RDS"))


params_nosocial <- as_draws_df(fit_nosocial)
beta_fit_nosocial_draws <- subset_draws(params_nosocial, variable = "beta")

params_nopersonal <- as_draws_df(fit_nopersonal)
beta_fit_nopersonal_draws <- subset_draws(params_nopersonal, variable = "beta")

params_nocatch <- as_draws_df(fit_nocatch)
beta_fit_nocatch_draws <- subset_draws(params_nocatch, variable = "beta")


color_scheme_set("teal")#
bayesplot_theme_set(theme_default())


dplyr::bind_cols(summarise_draws(beta_fit_nosocial_draws),c("Int catch 0",
                                                     "Int catch 1",
                                                     "Fish",
                                                     "time since event",
                                                     "distance sonar",
                                                     "success",
                                                     "loss"
))

bayesplot::mcmc_areas(params_nosocial,regex_pars = "beta\\[[1-7]\\]", prob = 0.8)+ scale_y_discrete(
  labels = c("beta[1]" = "Int catch 0",
             "beta[2]" = "Int catch 1",
             "beta[3]" = "Fish",
             "beta[4]" = "Time since event",
             "beta[5]" = "distance sonar",
             "beta[6]" = "success",
             "beta[7]" = "loss"))




dplyr::bind_cols(summarise_draws(beta_fit_nopersonal_draws),c("Int catch 0",
                                                            "Int catch 1",
                                                            "Fish",
                                                            "time since event",
                                                            "distance sonar"
))

bayesplot::mcmc_areas(params_nopersonal,regex_pars = "beta\\[[1-5]\\]", prob = 0.8)+ scale_y_discrete(
  labels = c("beta[1]" = "Int catch 0",
             "beta[2]" = "Int catch 1",
             "beta[3]" = "Fish",
             "beta[4]" = "Time since event",
             "beta[5]" = "distance sonar"))

dplyr::bind_cols(summarise_draws(beta_fit_nocatch_draws),c("Int",
                                                              "Fish",
                                                              "distance sonar"
))


bayesplot::mcmc_areas(params_nopersonal,regex_pars = "beta\\[[1-3]\\]", prob = 0.8)+ scale_y_discrete(
  labels = c("beta[1]" = "Int",
             "beta[2]" =  "Fish",
             "beta[3]" = "distance sonar"))

################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
##################