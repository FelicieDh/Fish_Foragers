################################################################################
#
# Title: 1. Calculate features for every chosen spot
#
# Authors: Dhellemmes F. 
#
# Last updated: 03/08/2026 (DD/MM/YYYY)
#
################################################################################

library(here)
library(data.table)
library(terra)
library(dplyr)
library(tidybayes)
library(posterior)
library(ggplot2)

source(here("functions","library.R"))

# source helper functions
sourceCpp(here("functions","fast_dist.cpp"))
sourceCpp(here("functions","ptinpoly.cpp"))

# set the seed

set.seed(42)

################################################################################

## All data imports

# Import Spot Choices
model_data = fread(here("data","2511_angling_events.csv"))

metadata<-fread(here::here("data","metadata_anon.csv"))

model_data[,compID_num:=metadata[match(model_data$compID,metadata$compID),"compID_num"]]
model_data[,partID_num:=metadata[match(model_data$compIDwatch,metadata$compIDwatch),"partID"]]

model_data<-model_data[,c(2,4,5,6,7,9,11,14,15,16,20,21,22)]

model_data[,minute:=round(posix_since_start/60)]

data_list <- split(model_data, seq(nrow(model_data)))

## Import metadata
metadata<-fread(here::here("data","metadata_anon.csv"))

## Import fish presence
file_list <- readRDS(here::here("data","2512_rasters_file_paths.rds"))

fish_presence <- list()
base_path <- here::here("data")

for (comp in names(file_list)) {
  full_paths <- file.path(base_path, file_list[[comp]])
  
  # Load rasters from these full paths
  fish_presence[[comp]] <- lapply(full_paths, terra::rast)
}


## Open arena and sonar data
arena<-fread(here::here("data","2408_ArenaExtents.csv"))
sonars<-rbind(arena[,c(6,7)],
              arena[,c(8,9)],
              arena[,c(10,11)],use.names=FALSE)



## Import model results and calculate lambdas
fit <- readRDS(here::here("1_step_selection","model_fit","step_selection_fit.RDS"))

model_fitted <- as_draws_df(fit)

tmp_main <- spread_draws(model_fitted, bandwidth_log[feature])

main_lambdas <- aggregate(
  bandwidth_log ~ feature,
  data = tmp_main,
  FUN = mean)

names(main_lambdas)[2] <- "mean_bw"

tmp_comp <- spread_draws(model_fitted, v_compID[compID, feature])

# filter feature %in% 7:11 (so all the bandwidths)
tmp_comp <- tmp_comp[tmp_comp$feature %in% 7:11, ]

# mutate(feature = feature - 6) # back to 1
tmp_comp$feature <- tmp_comp$feature - 6

# group_by(feature, compID) %>% summarize(mean)
comp_lambdas <- aggregate(
  v_compID ~ feature + compID,
  data = tmp_comp,
  FUN = mean)

names(comp_lambdas)[3] <- "mean_bw_offset"

lambdas <- merge(
  main_lambdas,
  comp_lambdas,
  by = "feature",
  all.x = TRUE)

lambdas$comp_bw <- exp(
  lambdas$mean_bw_offset + lambdas$mean_bw)

lambdas$compID_org<-metadata[match(lambdas$compID,metadata$compID_num), compID]


saveRDS(lambdas, here::here("data","1_step_selection", "model_fit", "fitted_bw.RDS")) #save the bandwidth for later


## Import social matrix
social_distance_matrix<-readRDS(here("data","1_step_selection","social_distance_matrix.rds"))
columns_social <- readRDS(here::here("data","1_step_selection","columns_social.rds"))

## Import success matrices
success_spatial_matrix<-readRDS(here("data","1_step_selection","success_spatial_matrix.rds"))
success_time_matrix<-readRDS(here("data","1_step_selection","success_time_matrix.rds"))

## Import success matrices
loss_spatial_matrix<-readRDS(here("data","1_step_selection","loss_spatial_matrix.rds"))
loss_time_matrix<-readRDS(here("data","1_step_selection","loss_time_matrix.rds"))

## Import columns 
columns_success <- readRDS(here::here("data","1_step_selection","columns_success.rds"))
columns_loss <- readRDS(here::here("data","1_step_selection","columns_loss.rds"))



###############################################################################
###############################################################################
###############################################################################
### Fish presence
###############################################################################
###############################################################################
###############################################################################


extract_fish_feature <- function(sub, fish_presence){
  
  compID_tmp = as.character(sub[1, "compID"])
  coord_tmp = sub[, c("E_utm", "N_utm")]
  min_tmp = sub[1,minute]
  
  cat(paste("Processing compID = ", compID_tmp,
            ", minute = ", min_tmp,"\r"))
  
  
  r <- raster(fish_presence[[compID_tmp]][[min_tmp]])
  
  sub[,"fish_3"] <-  raster::extract(r,coord_tmp, buffer = 3, fun=mean) 
  
  # Free memory
  rm(r)       # remove raster object
  gc(verbose=FALSE)  # trigger garbage collection quietly
  
  return(sub)
}

data_list <- lapply(data_list, extract_fish_feature, fish_presence = fish_presence)


# #check results
# i=750
# i=75
# plot(fish_presence[[data_list[[i]]$compID]][[data_list[[i]]$minute]])
# text(data_list[[i]]$E_utm, data_list[[i]]$N_utm, round(data_list[[i]]$fish_3,2), cex=0.6)
# 
# #this looks great



###############################################################################
###############################################################################
###############################################################################
### distance to the sonar
###############################################################################
###############################################################################
###############################################################################

extract_sonar_feature <- function(sub, sonars){
  
  result = matrix(apply(fastPdist2(as.matrix(sub[1,c("E_utm", "N_utm")]), as.matrix(sonars[,c(1,2)])), 1, FUN = min), ncol=1)
  
  sub$dist_sonar = result
  
  return(sub)
}

data_list <- lapply(data_list,  extract_sonar_feature, sonars)


# #check results
# i = 760
# plot(fish_density[[data_list[[i]]$compID]][[data_list[[i]]$minute]])
# points(sonars$xpink, sonars$ypink, col="hotpink", pch=19)
# text(data_list[[i]]$E_utm, data_list[[i]]$N_utm, round(data_list[[i]]$dist_sonar,2), cex=.6, col="white")
# 
# #this looks great


################################################################################
###############################################################################
###############################################################################
### social_kernel
###############################################################################
###############################################################################
###############################################################################

data_list<- Map(cbind, data_list, index = seq_along(data_list))

sub<-data_list[[1]]

social_kernel <- function(sub, social_distances, columns_social, lambdas){
  
  # subset spot identifier
  index_tmp = sub[1, index]
  compID_tmp = sub[1, compID]
  
  # compute distances to target
  sub$social_density <- sum(exp( - (social_distances[[index_tmp]][1,1:columns_social[index_tmp]] ^ 2) /
                                   (2*(lambdas$comp_bw[lambdas$compID_org == compID_tmp &
                                    lambdas$feature == 1]^ 2)))  * (1.0 / (1 + exp(-20 * (social_distances[[index_tmp]][1,1:columns_social[index_tmp]] - 5)))))
  return(sub)
  
}

data_list<-lapply(data_list,  social_kernel, social_distance_matrix, columns_social, lambdas)


gps_data<-fread(here("data","2604_gps_data.csv"))
gps_data$min<-round(gps_data$posix_since_start/60,0)

# i=670
# data_list[[i]]
# plot(fish_density[[data_list[[i]]$compID]][[data_list[[i]]$minute]])
# points(gps_data[compID==data_list[[i]]$compID & min == data_list[[i]]$minute, c("E_utm", "N_utm")], col=scales::alpha("aquamarine4",0.4), pch=19)
# points(data_list[[i]][, c("E_utm", "N_utm")], col="red", pch=19)
# 
# #Looks good



################################################################################
###############################################################################
###############################################################################
### success and loss kernels
###############################################################################
###############################################################################
###############################################################################


success_loss_kernel_time <- function(sub, success_spatial_matrix, success_time_matrix, loss_spatial_matrix, loss_time_matrix, columns_success, columns_loss, lambdas){
  
  # subset spot identifier
  compID_tmp = sub[1, compID]
  index_tmp = sub[1, index]
  
  if (success_spatial_matrix[[index_tmp]][1,1]<0) {
    sub$success_density <- NA
  } else {
    # compute distances to target
    sub$success_density <- sum(exp( - (success_spatial_matrix[[index_tmp]][1,1:columns_success[index_tmp]] ^ 2) /
                                      (2*(lambdas$comp_bw[lambdas$compID_org == compID_tmp &
                                                            lambdas$feature == 2]^ 2))) *
                                 exp( - (success_time_matrix[[index_tmp]][1,1:columns_success[index_tmp]] ^ 2) /
                                        (2*(lambdas$comp_bw[lambdas$compID_org == compID_tmp &
                                                              lambdas$feature == 3]^ 2))))
  }
  
  if (loss_spatial_matrix[[index_tmp]][1,1]<0) {
    sub$loss_density <- NA
  } else {
    sub$loss_density <- sum(exp( - (loss_spatial_matrix[[index_tmp]][1,1:columns_loss[index_tmp]] ^ 2) /
                                   (2*(lambdas$comp_bw[lambdas$compID_org == compID_tmp &
                                                         lambdas$feature == 4]^ 2))) *
                              exp( - (loss_time_matrix[[index_tmp]][1,1:columns_loss[index_tmp]] ^ 2) /
                                     (2*(lambdas$comp_bw[lambdas$compID_org == compID_tmp &
                                                           lambdas$feature == 5]^ 2))) )
  }
  return(sub)
  
}


data_list<-lapply(data_list, success_loss_kernel_time, success_spatial_matrix, success_time_matrix, loss_spatial_matrix, loss_time_matrix, columns_success, columns_loss, lambdas)



saveRDS(data_list, here("data","2_ecological_validation_step_selection","Spotchoice_features.rds"))


################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################
