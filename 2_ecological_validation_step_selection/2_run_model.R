################################################################################
#
# Title: 2. Run models
#
# Authors: Dhellemmes F. 
#
# Last updated: 03/08/2026 (DD/MM/YYYY)
#
################################################################################

library(cmdstanr)
library(posterior)
library(bayesplot)
library(data.table)
library(ggplot2)
library(here)
library(brms)

# set the seed

set.seed(42)

################################################################################


data<-readRDS(here("data","2_ecological_validation_step_selection","Spotchoice_features.rds"))

metadata<-fread(here::here("data","metadata_anon.csv"))



data <- do.call(rbind, data)

data[,fish_star := fish_3/max(fish_3)] 

N <- length(data$fish_star)
data$fish_star <- (data$fish_star * (N - 1) + 0.5) / N ##beta distribution


data[,time_std:=(minute-mean(minute))/sd(minute)] ##std
data[,dist_sonar_std:=(dist_sonar-mean(dist_sonar))/sd(dist_sonar)] ##std


data[, `:=`(
  success_density_std =
    (success_density - mean(success_density, na.rm = TRUE)) /
    sd(success_density, na.rm = TRUE),
  
  loss_density_std =
    (loss_density - mean(loss_density, na.rm = TRUE)) /
    sd(loss_density, na.rm = TRUE),
  
  social_density_std =
    (social_density - mean(social_density, na.rm = TRUE)) /
    sd(social_density, na.rm = TRUE)
), by = compID_num] #std


sd(data$dist_sonar)

ev_ss_fit <- brm(
  fish_star ~ dist_sonar_std + time_std + social_density_std + success_density_std + loss_density_std + (time_std|compID),
  data = data,
  family = Beta(link = "logit"),
  iter = 3000,
  warmup = 1000,
  chains = 4,
  cores = 4)

#summary(ev_ss_fit)


saveRDS(ev_ss_fit, file = here::here("data","2_ecological_validation_step_selection","model_fit","ev_ss_fit.RDS"))


################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################