################################################################################
#
# Title: 1. Preparing data and running model
#
# Author F Dhellemmes
#
# last update: 05/08/2026 (DD/MM/YYYY)
#
################################################################################

library(here)
library(data.table)
library(terra)
library(dplyr)
library(tidybayes)
library(ggplot2)
library(posterior)
library(bayesplot)
library(ggplot2)
library(brms)

source(here("functions","library.R"))


# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################


stan_data<-readRDS(here::here("data","3_spot_leaving", "stan_spot_leaving.RDS"))

data<-fread(here::here("data","2511_timebin_data.csv"))

model_data<-as.data.table(data.frame(catch = stan_data[["catch"]],
                       fish_3 = stan_data[["fish3"]],
                       time_since_event = stan_data[["time_since_event"]],
                       compID = stan_data[["compID"]],
                       dist_sonar = stan_data[["dist_sonar"]],
                       social_density = stan_data[["social_density"]],
                       success_density = stan_data[["success_density"]],
                       loss_density = stan_data[["loss_density"]]))
             
model_data$spot_num<-data$spot_num
model_data$time_bins<-data$time_bins
model_data$compIDwatch<-data$compIDwatch
model_data[,min_bin:=floor(time_bins/6)+1]

#################################################
##Binning model_data

model_data_bin<-model_data[, .(compID = max(compID),
               bi_catch = max(catch),
               time_since = min(time_since_event),
               avg_social = mean(social_density),
               avg_success = mean(success_density),
               avg_loss = mean(loss_density),
               avg_fish = mean(fish_3),
               avg_sonar = mean(dist_sonar)), by = .(min_bin, spot_num, compIDwatch)]

## Beta transformation
model_data_bin[,fish_star := avg_fish/max(avg_fish,na.rm=T)]
N <- length(model_data_bin$fish_star)
model_data_bin$fish_star <- (model_data_bin$fish_star * (N - 1) + 0.5) / N

model_data_bin[,time_since := floor(time_since/6)+1] #per minute

model_data_bin[, `:=`(
  success_density_std =
    (avg_success - mean(avg_success, na.rm = TRUE)) /
    sd(avg_success, na.rm = TRUE),
  
  loss_density_std =
    (avg_loss - mean(avg_loss, na.rm = TRUE)) /
    sd(avg_loss, na.rm = TRUE),
  
  social_density_std =
    (avg_social - mean(avg_social, na.rm = TRUE)) /
    sd(avg_social, na.rm = TRUE),
  
  avg_sonar_std =
    (avg_sonar - mean(avg_sonar, na.rm = TRUE)) /
    sd(avg_sonar, na.rm = TRUE)
), by = compID]




##Run the model
fit <- brm(
  fish_star ~ bi_catch + time_since + social_density_std + success_density_std + loss_density_std + avg_sonar_std +(1|compID),
  data = model_data_bin,
  family = Beta(link = "logit"),
  iter = 3000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 42)

summary(fit)

conditional_effects(fit)


saveRDS(fit, file = here::here("data","4_ecological_validation_spot_leaving","model_fit","model_validation_fit.RDS"))


################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################


