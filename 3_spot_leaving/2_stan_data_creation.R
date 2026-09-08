################################################################################
#
# Title: 2. Create stan data
#
# Author F Dhellemmes
#
# last update: 03/08/2026 (DD/MM/YYYY)
#
################################################################################


library(data.table)

# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################

## Open data
data<-fread(here::here("data","2511_timebin_data.csv"))

## Open social matrix
social_distance_matrix <- readRDS(here::here("data","3_spot_leaving","social_distance_matrix.rds"))
social_distance_matrix <- lapply(social_distance_matrix, function(x) as.vector(x))

## Open success and loss matrix
success_spatial_matrix <- readRDS(here::here("data","3_spot_leaving","success_spatial_matrix.rds"))
success_time_matrix <- readRDS(here::here("data","3_spot_leaving","success_time_matrix.rds"))
loss_spatial_matrix <- readRDS(here::here("data","3_spot_leaving","loss_spatial_matrix.rds"))
loss_time_matrix <- readRDS(here::here("data","3_spot_leaving","loss_time_matrix.rds"))

success_spatial_matrix <- lapply(success_spatial_matrix, function(x) as.vector(x))
success_time_matrix <- lapply(success_time_matrix, function(x) as.vector(x))
loss_spatial_matrix <- lapply(loss_spatial_matrix, function(x) as.vector(x))
loss_time_matrix <- lapply(loss_time_matrix, function(x) as.vector(x))

## Open columns
columns_social <- readRDS(here::here("data","3_spot_leaving","columns_social.rds"))
columns_success <- readRDS(here::here("data","3_spot_leaving","columns_success.rds"))
columns_loss <- readRDS(here::here("data","3_spot_leaving","columns_loss.rds"))
names(columns_social) <- NULL
names(columns_success) <- NULL
names(columns_loss) <- NULL

## Open locality and sonar
locality <- readRDS(here::here("data","3_spot_leaving","locality_matrix.rds"))
sonar_distance <- readRDS(here::here("data","3_spot_leaving","sonar_distance_matrix.rds"))
names(sonar_distance) <- NULL

locality <- lapply(locality, function(x) as.vector(x))
sonar_distance <- lapply(sonar_distance, function(x) as.vector(x))


## Open  fish
fish_matrix <- readRDS(here::here("data","3_spot_leaving","fish_feature.rds"))
fish_matrix <- lapply(fish_matrix, function(x) as.vector(x))
names(fish_matrix) <- NULL
fish3<-unlist(fish_matrix)

data<-cbind(data,fish3)

trip_bounds <- data[, .(
  first_row = .I[1],   # index of first row in this group
  last_row  = .I[.N]   # index of last row in this group
), by = .(tripID_num)]




################################################################################
#### Social Success and Loss features calculation using the step selection model
################################################################################

##open the bandwidths
bw<-readRDS(here::here("1_step_selection", "model_fit", "fitted_bw.RDS"))

data_list <- split(data, seq(nrow(data)))
data_list<- Map(cbind, data_list, index = seq_along(data_list))

social_density<-list()

social_kernel <- function(sub, social_distances, columns_social, lambdas){
  
  # subset spot identifier
  index_tmp = sub[1, index]
  compID_tmp = sub[1, compID]
  
  
  # compute distances to target
  temp <- sum(exp( - (social_distances[[index_tmp]][1:columns_social[index_tmp]] ^ 2) /
                     (2*(lambdas$comp_bw[lambdas$compID_org == compID_tmp &
                                           lambdas$feature == 1]^ 2)))  * (1.0 / (1 + exp(-20 * (social_distances[[index_tmp]][1:columns_social[index_tmp]] - 5)))))
  return(temp)
  
}

social_density<-lapply(data_list,  social_kernel, social_distance_matrix, columns_social, bw)


success_density<-list()
loss_density<-list()

success_kernel_time <- function(sub, success_spatial_matrix, success_time_matrix, columns_success, lambdas){
  
  
  # subset spot identifier
  compID_tmp = sub[1, compID]
  index_tmp = sub[1, index]
  
  cat(paste(index_tmp, "\r"))
  if (success_spatial_matrix[[index_tmp]][1]<0) {
    temp_success <- 0
  } else {
    # compute distances to target
    temp_success <- sum(exp( - (as.numeric(success_spatial_matrix[[index_tmp]][1:columns_success[index_tmp]]) ^ 2) /
                               (2*(lambdas$comp_bw[lambdas$compID_org == compID_tmp &
                                                     lambdas$feature == 2]^ 2))) *
                          exp( - (as.numeric(success_time_matrix[[index_tmp]][1:columns_success[index_tmp]]) ^ 2) /
                                 (2*(lambdas$comp_bw[lambdas$compID_org == compID_tmp &
                                                       lambdas$feature == 3]^ 2))))
  }
  
  
  return(temp_success)
  
}

success_density<-lapply(data_list,  success_kernel_time, success_spatial_matrix, success_time_matrix, columns_success, bw)




loss_kernel_time <- function(sub, loss_spatial_matrix, loss_time_matrix, columns_loss, lambdas){
  
  # subset spot identifier
  compID_tmp = sub[1, compID]
  index_tmp = sub[1, index]
  cat(paste(index_tmp, "\r"))
  
  if (loss_spatial_matrix[[index_tmp]][1]<0) {
    temp_loss <- 0
  } else {
    temp_loss <- sum(exp( - (as.numeric(loss_spatial_matrix[[index_tmp]][1:columns_loss[index_tmp]]) ^ 2) /
                            (2*(lambdas$comp_bw[lambdas$compID_org == compID_tmp &
                                                  lambdas$feature == 4]^ 2))) *
                       exp( - (as.numeric(loss_time_matrix[[index_tmp]][1:columns_loss[index_tmp]]) ^ 2) /
                              (2*(lambdas$comp_bw[lambdas$compID_org == compID_tmp &
                                                    lambdas$feature == 5]^ 2))) )
  }
  return(temp_loss)
  
}

loss_density<-lapply(data_list,  loss_kernel_time, loss_spatial_matrix, loss_time_matrix, columns_loss, bw)




# Prepare data (replace these with your actual data vectors)
stan_data <- list(
  N = nrow(data), #number of spot choices
  leave = data$leave,
  
  catch = data$catch,
  time_since_event = data$time_since_event,
  cumulative_catch = data$cum_catch,
  
  Ntrips = length(unique(data$tripID_num)),
  tripID = data$tripID_num,
  slice_trips = 1:length(unique(data$tripID_num)),
  
  Ncomp = length(unique(data$compID_num)),
  compID = data$compID_num,
  Nids = length(unique(data$partID_num)),
  partID = data$partID_num,
  
  trip_start = trip_bounds[,first_row],
  trip_end = trip_bounds[,last_row],
  
  dist_sonar = unlist(sonar_distance)/max(unlist(sonar_distance)), ##make dist_edge be between 0 and 1
  fish3 = data$fish3,
  dist_last_spot = unlist(locality),
  
  social_density = as.numeric(do.call(rbind, social_density)),
  success_density = as.numeric(do.call(rbind, success_density)),
  loss_density = as.numeric(do.call(rbind, loss_density))
  
)

saveRDS(stan_data, here::here("data","3_spot_leaving", "stan_spot_leaving.RDS"))




################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################
