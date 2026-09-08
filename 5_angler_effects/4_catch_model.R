################################################################################
#
# Title: 4. run the catch model
#
# Author F Dhellemmes
#
# last update: 06/08/2026 (DD/MM/YYYY)
#
################################################################################


library(here)
library(data.table)
library(ggplot2)
library(terra)
library(kableExtra)
library(brms)
library(tidybayes)

source(here("functions","library.R"))

sourceCpp(here("functions","fast_dist.cpp"))
sourceCpp(here("functions","ptinpoly.cpp"))


# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################

real<-readRDS(here::here("data","5_angler_effects","RealSpots.rds"))
real[,counter:=seq(1,.N,by=1), by=anglingID]

real[,time_at_spot := max(counter), by=anglingID]

model_catch<-brm(bi_catch~counter+fish_3+dist_sonar+time_at_spot+(1|compID)+(1|anglingID),
                 data = real[success==1,],
                 family = bernoulli,
                 iter = 3000,
                 warmup = 1000,
                 chains = 4,
                 cores = 4,
                 seed=42)

# round(summary(model_catch)$fixed,2)
# conditional_effects(model_catch)


saveRDS(model_catch, here::here("data","5_angler_effects", "model_fit", "model_catch.RDS"))


# model_catch<-readRDS(here::here("data","5_angler_effects", "model_fit", "model_catch.RDS"))
################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################
