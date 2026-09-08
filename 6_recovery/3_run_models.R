################################################################################
#
# Title: 3. run model
#
# Author F Dhellemmes
#
# last update: 08/08/2026 (DD/MM/YYYY)
#
################################################################################


library(data.table)
library(terra)
library(here)
library(brms)
library(tidybayes)
library(ggplot2)

# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################


grid_data15<-readRDS(here::here("data","6_recovery","grid_data_15.RDS"))
grid15<-readRDS(here::here("data", "6_recovery", "voronoi_grid_15.rds"))
grid15 <- lapply(grid15, terra::unwrap)

grid_data15[,cell_name:=paste0(compID,"_",cell_id)]
grid_data15_events<-grid_data15[!is.na(event_name),]


summary_grid<-grid_data15_events[,.(compID=max(compID),
                                    minute=min(minute),
                                    cell_id=max(cell_id),
                                    angler_presence=max(angler_presence),
                                    angler_nr=max(angler_nr),
                                    catch_bi=max(catch_bi),
                                    catch_cum=sum(catch_abs),
                                    event_nr=min(event_nr),
                                    time_away=min(time_away),
                                    successful_event=min(successful_event),
                                    cumulative_time=min(cumulative_time),
                                    fish=mean(mean_fish),
                                    event_time=max(event_time),
                                    cell_name=max(cell_name)),
                                 by=event_name]

grid_data15_successful<-grid_data15_events[
  , if (any(catch_bi > 0)) .SD, 
  by = cell_name
] ##only succesful cells



mod_full<-brm(catch_bi~event_time+event_nr+time_away+angler_nr+mean_fish+(1|cell_name), grid_data15_successful, 
              family = bernoulli,
              iter = 3000,
              warmup = 1000,
              chains = 4,
              cores = 4,
              seed=42)

#summary(mod_full)

fixef(mod_full)

conditional_effects(mod_full)

saveRDS(mod_full, here::here("data","6_recovery","model_fit","catch_recovery_model_fit.RDS"))




################################################################################
############ Grid 20 - sensitivity #############################################
################################################################################


grid_data20<-readRDS(here::here("data","6_recovery","grid_data_20.RDS"))
grid20<-readRDS(here::here("data", "6_recovery", "voronoi_grid_20.rds"))
grid20 <- lapply(grid20, terra::unwrap)

grid_data20[,cell_name:=paste0(compID,"_",cell_id)]
grid_data20_events<-grid_data20[!is.na(event_name),]


summary_grid<-grid_data20_events[,.(compID=max(compID),
                                    minute=min(minute),
                                    cell_id=max(cell_id),
                                    angler_presence=max(angler_presence),
                                    angler_nr=max(angler_nr),
                                    catch_bi=max(catch_bi),
                                    catch_cum=sum(catch_abs),
                                    event_nr=min(event_nr),
                                    time_away=min(time_away),
                                    successful_event=min(successful_event),
                                    cumulative_time=min(cumulative_time),
                                    fish=mean(mean_fish),
                                    event_time=max(event_time),
                                    cell_name=max(cell_name)),
                                 by=event_name]

grid_data20_successful<-grid_data20_events[
  , if (any(catch_bi > 0)) .SD,
  by = cell_name
] ##only succesful cells



mod_20<-brm(catch_bi~event_time+event_nr+time_away+angler_nr+mean_fish+(1|cell_name), grid_data20_successful,
              family = bernoulli,
              iter = 3000,
              warmup = 1000,
              chains = 4,
              cores = 4,
              seed=42)

#summary(mod_full)

fixef(mod_20)

saveRDS(mod_20, here::here("data","6_recovery","model_fit","catch_recovery_model_fit_20.RDS"))




################################################################################
############ Grid 10 - sensitivity #############################################
################################################################################


grid_data10<-readRDS(here::here("data","6_recovery","grid_data_10.RDS"))
grid10<-readRDS(here::here("data", "6_recovery", "voronoi_grid_10.rds"))
grid10 <- lapply(grid10, terra::unwrap)

grid_data10[,cell_name:=paste0(compID,"_",cell_id)]
grid_data10_events<-grid_data10[!is.na(event_name),]


summary_grid<-grid_data10_events[,.(compID=max(compID),
                                    minute=min(minute),
                                    cell_id=max(cell_id),
                                    angler_presence=max(angler_presence),
                                    angler_nr=max(angler_nr),
                                    catch_bi=max(catch_bi),
                                    catch_cum=sum(catch_abs),
                                    event_nr=min(event_nr),
                                    time_away=min(time_away),
                                    successful_event=min(successful_event),
                                    cumulative_time=min(cumulative_time),
                                    fish=mean(mean_fish),
                                    event_time=max(event_time),
                                    cell_name=max(cell_name)),
                                 by=event_name]

grid_data10_successful<-grid_data10_events[
  , if (any(catch_bi > 0)) .SD,
  by = cell_name
] ##only succesful cells



mod_10<-brm(catch_bi~event_time+event_nr+time_away+angler_nr+mean_fish+(1|cell_name), grid_data10_successful,
            family = bernoulli,
            iter = 3000,
            warmup = 1000,
            chains = 4,
            cores = 4,
            seed=42)

#summary(mod_full)

fixef(mod_10)

saveRDS(mod_10, here::here("data","6_recovery","model_fit","catch_recovery_model_fit_10.RDS"))



################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################

