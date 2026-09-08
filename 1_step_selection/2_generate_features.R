################################################################################
#
# Title: 2. Create feature matrices for spot selection
#
# Author F Dhellemmes
#
# last update: 31/07/2026 (DD/MM/YYYY)
#
################################################################################

library(data.table)
library(terra)
library(here)

source(here("functions","library.R"))

# source helper functions
sourceCpp(here("functions","fast_dist.cpp"))


# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################



#####################################################
# Pad matrices so they all have the same dimension 
#####################################################

pad_matrix <- function(matrix, new_rows, new_columns, filler){
  
  if ( new_rows > nrow(matrix) & new_columns > ncol(matrix)){
    
    # add columns
    column_diff = new_columns - ncol(matrix)
    columns = replicate(n = column_diff, expr = rep(filler, nrow(matrix)))
    matrix = cbind(matrix, columns)
    
    # add rows
    row_diff = new_rows - nrow(matrix)
    rows = t(replicate(n = row_diff, expr = rep(filler, ncol(matrix))))
    matrix = rbind(matrix, rows)
    
  } else if (new_rows > nrow(matrix) & new_columns == ncol(matrix)){
    
    # add rows
    row_diff = new_rows - nrow(matrix)
    rows = t(replicate(n = row_diff, expr = rep(filler, ncol(matrix))))
    matrix = rbind(matrix, rows)
    
  } else if (new_rows == nrow(matrix) & new_columns > ncol(matrix) & new_rows > 1){
    
    # add columns
    column_diff = new_columns - ncol(matrix)
    columns = replicate(n = column_diff, expr = rep(filler, nrow(matrix)))
    matrix = cbind(matrix, columns)
    
  } else if (new_rows == 1){
    
    ## add columns
    column_diff = new_columns - ncol(matrix)
    columns = replicate(n = column_diff, expr = rep(filler, 1))
    matrix = c(matrix, columns)
    
  } else {
    # do nothing
    #print("New dimensions are old dimensions.")
  }
  
  return(matrix)
}

################################################################################
################################################################################
################################################################################

#time_option=0
#####################################################
# Extract social distances overall
#####################################################
extract_social <- function(sub, gps_data, time_option){
  
  #sub = Feature_matrix[[identifiers$spot_id != 1]]
  #sub = sub[[45]]
  
  # define cpp function for cluster export
  sourceCpp("functions/fast_dist.cpp")
  
  # subset spot identifier
  compID_tmp = sub[1, "compID"]
  id_tmp = sub[1, "compIDwatch"]
  coord_tmp = sub[, c("E_utm", "N_utm")]
  time_tmp = sub[1, "posix_time"]
  
  # subset coordinates of other participants
  coordinates <- gps_data[compID == as.character(compID_tmp) & 
                            compIDwatch != as.character(id_tmp) &
                            posix_time == as.numeric(time_tmp) - time_option &
                            !is.na(E_utm),c("compIDwatch","E_utm","N_utm")]
  
 
  # compute distances to target
  tmp_distances = fastPdist2(as.matrix(coord_tmp), as.matrix(coordinates[, 2:3]))
  
  return(tmp_distances)
  
}

################################################################################
################################################################################
################################################################################

#sub<-Feature_matrix[[26]]
#option=0
#########################################
# Compute distance to successes and losses
########################################

extract_spot_distance <- function(sub, model_data, option){ #option either 0 or 1
  
  # define cpp function for cluster export
  #sourceCpp("functions/fast_dist.cpp")
  
  # subset spot identifier
  compID_tmp = sub[1, "compID"]
  id_tmp = sub[1, "compIDwatch"]
  coord_tmp = sub[, c("E_utm", "N_utm")]
  spot_tmp = sub[1, "spot_num"]
  time_tmp = sub[1, "posix_time"]
  
  
  # get coordinates of last spots
  spots <- model_data[compID == as.character(compID_tmp) &
                        compIDwatch == as.character(id_tmp) &
                        spot_num < as.numeric(spot_tmp) &
                        catch_at_spot == option, c("spot_num", "E_utm", "N_utm", "posix_stop")]
  
  spots[,"spot_since"]<- round(as.numeric(spot_tmp)-spots$spot_num)
  
  spots<-spots[order(spots$spot_since),]
  
  if(nrow(spots) == 0){
    space = matrix(rep(-99, nrow(sub)),ncol=1)
    time = matrix(rep(-99, nrow(sub)),ncol=1)
  } else {
    space = matrix(fastPdist2(as.matrix(coord_tmp), as.matrix(spots[, c("E_utm", "N_utm")]))[, 1:min(nrow(spots),50)], nrow=nrow(sub)) 
    time = matrix(rep(spots$spot_since[1:min(nrow(spots),50)],each=nrow(sub)), nrow=nrow(sub))
    
    # sometimes missing, if distances very close
    space[is.na(space)]<-.00001
    
  }
  
  return(list(
    space  = space, 
    time = time))
  
}

################################################################################
################################################################################
################################################################################

###################################
# extract distance to last spot
###################################

extract_locality_feature <- function(sub, model_data){
  
  # define cpp function for cluster export
  #sourceCpp("functions/fast_dist.cpp")
  
  # subset spot identifier
  compID_tmp = sub[1, "compID"]
  id_tmp = sub[1, "compIDwatch"]
  coord_tmp = sub[, c("E_utm", "N_utm")]
  spot_tmp = sub[1, "spot_num"]
  
  
  # get coordinates of last spots
  spots <- model_data[compID == as.character(compID_tmp) &
                        compIDwatch == as.character(id_tmp) &
                        spot_num < as.numeric(spot_tmp), c("spot_num", "E_utm", "N_utm")]
  
  
  if(nrow(spots) == 0){
    result = matrix(rep(0, nrow(sub)), ncol=1)
  } else {
    result = fastPdist2(as.matrix(coord_tmp), as.matrix(spots[nrow(spots), c("E_utm", "N_utm")]))
  }
  
  # sometimes missing, if distances very close
  result[is.na(result)]<-0
  
  return(result)
  
}

################################################################################
################################################################################
################################################################################


###################################
# extract distance to sonar
###################################

extract_sonar_feature <- function(sub, sonars){
  
  result = matrix(apply(fastPdist2(as.matrix(sub[,c("E_utm", "N_utm")]), as.matrix(sonars[,c(1,2)])), 1, FUN = min), ncol=1)
  
  return(result)
}
################################################################################
################################################################################
################################################################################



################################################################################
################################################################################
################################################################################
# LOAD DATA
################################################################################
################################################################################
################################################################################

model_data<-fread(here("data","2511_angling_events.csv"))

metadata<-fread(here::here("data","metadata_anon.csv"))

gps_data<-fread(here("data","2604_gps_data.csv"))

Feature_matrix <-  readRDS(here("data", "1_step_selection","Feature_matrix.rds"))

arena<-fread(here::here("data","2408_ArenaExtents.csv"))



model_data[, start_min := round(posix_since_start/60)]
model_data[, stop_min := round((posix_since_start+posix_stop-posix_start)/60)]

model_data[,compID_num:=metadata[match(model_data$compID,metadata$compID),"compID_num"]]
model_data[,partID_num:=metadata[match(model_data$compIDwatch,metadata$compIDwatch),"partID"]]

gps_data[,compID_num:=metadata[match(gps_data$compID,metadata$compID),"compID_num"]]
gps_data[,partID_num:=metadata[match(gps_data$compIDwatch,metadata$compIDwatch),"partID"]]


sonars<-rbind(arena[,c(6,7)],
              arena[,c(8,9)],
              arena[,c(10,11)],use.names=FALSE)



################################################################################
################################################################################
################################################################################

################################################################################
# extract social distances
################################################################################

social_distance_matrix = lapply(Feature_matrix, function(x) extract_social(x, gps_data, time_option = 0))

columns_social<-sapply(social_distance_matrix, ncol)
social_distance_matrix<-lapply(social_distance_matrix, function(x) {pad_matrix(x, max(sapply(social_distance_matrix, nrow)) ,max(sapply(social_distance_matrix, ncol)), -99)})


saveRDS(columns_social, file = here::here("data","1_step_selection","columns_social.rds"))
saveRDS(social_distance_matrix, file = here::here("data","1_step_selection","social_distance_matrix.rds"))

# ##check
# polygons<-readRDS(here("data","arena_polygons.rds"))
# polygons <- lapply(polygons, terra::unwrap)
# plot(polygons[["som15am"]])
# points(Feature_matrix[[1]][1,c("E_utm",  "N_utm")], pch=19, col="red")
# points(gps_data[compID == as.character(Feature_matrix[[1]][1,c("compID")]) &
#                   compIDwatch != as.character(Feature_matrix[[1]][1,c("compIDwatch")]) &
#                   posix_time == as.numeric(Feature_matrix[[1]][1,c("posix_time")]) - 0 &
#                   !is.na(E_utm),c("E_utm","N_utm")], pch=19, col="aquamarine")
# 
# social_distance_matrix[[1]] ##looks right

rm(social_distance_matrix)
gc()
gc()


################################################################################
################################################################################
################################################################################


###############################################################################
# extract success loss
###############################################################################


success_distance_matrix = lapply(Feature_matrix, function(x) extract_spot_distance(x, model_data, option = 1))

loss_distance_matrix = lapply(Feature_matrix, function(x) extract_spot_distance(x, model_data, option = 0))

success_spatial_matrix  <- lapply(success_distance_matrix, `[[`, "space")
success_time_matrix <- lapply(success_distance_matrix, `[[`, "time")
loss_spatial_matrix  <- lapply(loss_distance_matrix, `[[`, "space")
loss_time_matrix <- lapply(loss_distance_matrix, `[[`, "time")



##check
# polygons<-readRDS(here("data","arena_polygons.rds"))
# polygons <- lapply(polygons, terra::unwrap)
# plot(polygons[[as.character(Feature_matrix[[26]][1,c("compID")])]])
# points(Feature_matrix[[26]][1,c("E_utm",  "N_utm")], pch=19, col="blue")
# points(model_data[compID == as.character(Feature_matrix[[26]][1,c("compID")]) &
#                   compIDwatch == as.character(Feature_matrix[[26]][1,c("compIDwatch")]) &
#                   spot_num < as.numeric(Feature_matrix[[26]][1,c("spot_num")]) &
#                   catch_at_spot == 0,c("E_utm","N_utm")], pch=19, col="red")
# points(model_data[compID == as.character(Feature_matrix[[26]][1,c("compID")]) &
#                     compIDwatch == as.character(Feature_matrix[[26]][1,c("compIDwatch")]) &
#                     spot_num < as.numeric(Feature_matrix[[26]][1,c("spot_num")]) &
#                     catch_at_spot == 1,c("E_utm","N_utm")], pch=19, col="green")
# 
# success_distance_matrix[[26]] ##looks right
# loss_distance_matrix[[26]] ##loooks good



##pad the matrices
columns_success<-sapply(success_spatial_matrix, ncol)
columns_loss<-sapply(loss_spatial_matrix, ncol)

success_spatial_matrix<-lapply(success_spatial_matrix, function(x) {pad_matrix(x, max(sapply(success_spatial_matrix, nrow)) ,max(sapply(success_spatial_matrix, ncol)), -99)})
success_time_matrix<-lapply(success_time_matrix, function(x) {pad_matrix(x, max(sapply(success_time_matrix, nrow)) ,max(sapply(success_time_matrix, ncol)), -99)})

loss_spatial_matrix<-lapply(loss_spatial_matrix, function(x) {pad_matrix(x, max(sapply(loss_spatial_matrix, nrow)) ,max(sapply(loss_spatial_matrix, ncol)), -99)})
loss_time_matrix<-lapply(loss_time_matrix, function(x) {pad_matrix(x, max(sapply(loss_time_matrix, nrow)) ,max(sapply(loss_time_matrix, ncol)), -99)})

saveRDS(columns_success, file = here("data","1_step_selection","columns_success.rds"))
saveRDS(columns_loss, file = here("data","1_step_selection","columns_loss.rds"))

saveRDS(success_spatial_matrix, file = here("data","1_step_selection","success_spatial_matrix.rds"))
saveRDS(loss_spatial_matrix, file = here("data","1_step_selection","loss_spatial_matrix.rds"))

saveRDS(success_time_matrix, file = here("data","1_step_selection","success_time_matrix.rds"))
saveRDS(loss_time_matrix, file = here("data","1_step_selection","loss_time_matrix.rds"))


rm(success_distance_matrix)
rm(loss_distance_matrix)

gc()
gc()

################################################################################
################################################################################
################################################################################

################################################################################
# extract locality feature
################################################################################
locality_matrix = list()

for (i in 1:length(Feature_matrix)){
  locality_matrix[[i]] <- extract_locality_feature(Feature_matrix[[i]], model_data)
}


# ##check
# polygons<-readRDS(here("data","arena_polygons.rds"))
# polygons <- lapply(polygons, terra::unwrap)
# plot(polygons[[as.character(Feature_matrix[[3]][1,c("compID")])]])
# points(Feature_matrix[[3]][1,c("E_utm",  "N_utm")], pch=19, col="blue")
# points(model_data[compID == as.character(Feature_matrix[[3]][1,c("compID")]) &
#                     compIDwatch == as.character(Feature_matrix[[2]][1,c("compIDwatch")]) &
#                     spot_num == as.numeric(Feature_matrix[[3]][1,c("spot_num")])-1, c( "E_utm", "N_utm")], pch=19, col="red")
# 
# 
# locality_matrix[[3]] ##looks right
# 


#locality_matrix = map(Feature_matrix, function(x) extract_locality_feature(x, model_data), .progress = TRUE)
saveRDS(locality_matrix, file = here("data","1_step_selection","locality_matrix.rds"))

rm(locality_matrix)
gc()


################################################################################
################################################################################
################################################################################

################################################################################
# extract edge feature
################################################################################


sonar_distance_matrix = lapply(Feature_matrix, function(x) extract_sonar_feature(x, sonars))

#hist(unlist(sonar_distance_matrix))
# #check
# polygons<-readRDS(here("data","arena_polygons.rds"))
# polygons <- lapply(polygons, terra::unwrap)
# plot(polygons[[as.character(Feature_matrix[[7]][1,c("compID")])]])
# points(Feature_matrix[[7]][1,c("E_utm",  "N_utm")], pch=19, col="blue")
# 
# sonar_distance_matrix[[7]] ##looks right


saveRDS(sonar_distance_matrix, file = here("data","1_step_selection","sonar_distance_matrix.rds"))


################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################