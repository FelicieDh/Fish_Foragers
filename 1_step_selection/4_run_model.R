################################################################################
#
# Title: 4. Running step selection
#
# Author F Dhellemmes
#
# last update: 31/07/2026 (DD/MM/YYYY)
#
################################################################################


library(cmdstanr)
library(posterior)
library(bayesplot)
library(data.table)
library(ggplot2)

# set the seed

set.seed(42)

################################################################################
################################################################################
################################################################################


stan_data<-readRDS(here::here("data","1_step_selection","stan_step_selection.RDS"))

stan_file <- here::here("1_step_selection","stan","step_selection.stan")


# Compile the model
step_selection_model <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

step_selection_fit <- step_selection_model$sample(
  data = stan_data,
  seed = 42,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  parallel_chains = 4,
  threads_per_chain = round(32/4)#,
  #init = inits_fun
)


params <- as_draws_df(step_selection_fit)


beta_draws <- subset_draws(params,variable = "beta")  
bandwidth_draws <- subset_draws(params,variable = "bandwidth_log")  


dplyr::bind_cols(summarise_draws(beta_draws),c("social",
                                                    "success",
                                                    "loss",
                                                    "distance_sonar",
                                                    "locality",
                                                    "locality*time"))


dplyr::bind_cols(summarise_draws(bandwidth_draws),c("social",
                                                         "success_space",
                                                         "success_time",
                                                         "loss_space",
                                                         "loss_time"))





step_selection_fit$save_object(file =  here::here("data","1_step_selection","model_fit","step_selection_fit.RDS"))

#step_selection_fit<-readRDS(here::here("data","1_step_selection","model_fit","step_selection_fit.RDS"))




################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################











