################################################################################
#
# Title: 3. Create stan data
#
# Author F Dhellemmes
#
# last update: 31/07/2026 (DD/MM/YYYY)
#
################################################################################

library(here)
library(data.table)

# set the seed

set.seed(42)

################################################################################
################################################################################
################################################################################


## load data

model_data<-fread(here("data","2511_angling_events.csv"))

## load identifier

identifiers<-readRDS(here("data","1_step_selection","identifiers.rds"))
Feature_matrix <-  readRDS(here("data","1_step_selection","Feature_matrix.rds"))

## load features

social_distance_matrix<-readRDS(here("data","1_step_selection","social_distance_matrix.rds"))

success_spatial_matrix<-readRDS(here("data","1_step_selection","success_spatial_matrix.rds"))
success_time_matrix<-readRDS(here("data","1_step_selection","success_time_matrix.rds"))

loss_spatial_matrix<-readRDS(here("data","1_step_selection","loss_spatial_matrix.rds"))
loss_time_matrix<-readRDS(here("data","1_step_selection","loss_time_matrix.rds"))

locality_matrix<-readRDS(here("data","1_step_selection","locality_matrix.rds"))
locality_matrix <- lapply(locality_matrix, as.vector)

sonar_distance_matrix<-readRDS(here("data","1_step_selection","sonar_distance_matrix.rds"))
sonar_distance_matrix <- lapply(sonar_distance_matrix, as.vector)


columns_social <- readRDS(here::here("data","1_step_selection","columns_social.rds"))
columns_success <- readRDS(here::here("data","1_step_selection","columns_success.rds"))
columns_loss <- readRDS(here::here("data","1_step_selection","columns_loss.rds"))


## Merging data in one megadataset

data<-do.call(rbind,Feature_matrix)

metadata<-fread(here::here("data","metadata_anon.csv"))

data$tripID_num<-metadata[match(data$compIDwatch,metadata$compIDwatch),"tripID_num"] ##add trip ID


data$time_since_start<-(data$posix_since_start)/(2*60*60) #make time between 0 and 1

data[,spot_ID:=paste0(compIDwatch, spot_num)]

indices<-data[!duplicated(data$spot_ID),]

data[, spot_rank := seq_len(.N), by = spot_ID]

#head(data)


stan_data <- list(
  N = length(unique(data$spot_ID)), #number of spot choices
  N_spots = max(data$spot_rank), #number of alternative options
  time_since_start = round(indices$time_since_start, digits=2),
  
  locality = lapply(locality_matrix, function(x) x / 120), ## Make it between 0 an 1, with 120m the maximum possible distance in the arena
  
  sonar_distance = lapply(sonar_distance_matrix, function(x) x / 36.65666), ##make dist_sonar be between 0 and 1
  
  columns_success = columns_success,
  max_columns_success = max(columns_success),
  spatial_success = success_spatial_matrix,
  time_success = success_time_matrix,
  
  columns_loss = columns_success,
  max_columns_loss = max(columns_loss),
  spatial_loss = loss_spatial_matrix, 
  time_loss = loss_time_matrix, 
  
  columns_social=columns_social,
  max_columns_social = max(columns_social),
  distance_social = social_distance_matrix,
  
  N_tripID = length(unique(data$tripID)),
  tripID = data[chosen==1,tripID_num],
  N_compID = length(unique(data$compID)),
  compID = data[chosen==1,compID_num],
  N_partID = length(unique(data$partID_num)),
  partID = data[chosen==1,partID_num]
)



saveRDS(stan_data, here::here("data","1_step_selection","stan_step_selection.RDS"))


################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################