################################################################################
#
# Title: 1. Create features for spot leaving
#
# Author F Dhellemmes
#
# last update: 03/08/2026 (DD/MM/YYYY)
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
    columns = matrix(replicate(n = column_diff, expr = rep(filler, 1)), nrow=1)
    matrix = cbind(matrix, columns)
    
  } else {
    # do nothing
    #print("New dimensions are old dimensions.")
  }
  
  return(matrix)
}


#####################################################
# Extract social distances
#####################################################
extract_social <- function(sub, gps_data, time_option){
  
  #sub = matrix[[45]]
  
  # define cpp function for cluster export
  sourceCpp("functions/fast_dist.cpp")
  
  # subset spot identifier
  compID_tmp = sub[1, "compID"]
  id_tmp = sub[1, "watch_id"]
  coord_tmp = sub[, c("E_utm", "N_utm")]
  time_tmp = sub[1, "posix_time"]
  
  # subset coordinates of other participants
  coordinates <- gps_data[compID == as.character(compID_tmp) & 
                            watch_id != as.numeric(id_tmp) &
                            posix_time == as.numeric(time_tmp) - time_option &
                            !is.na(E_utm),c("watch_id","E_utm","N_utm")]
  
  
  # compute distances to target
  tmp_distances = fastPdist2(as.matrix(coord_tmp), as.matrix(coordinates[, 2:3]))
  
  return(tmp_distances)
  
}



#########################################
# compute distance to successes and losses
########################################

#sub= matrix[[101]]
# option =0
extract_spot_distance <- function(sub, model_data, option){ #option either 0 or 1
  
  # define cpp function for cluster export
  #sourceCpp("functions/fast_dist.cpp")
  
  # subset spot identifier
  compID_tmp = sub[1, "compID"]
  id_tmp = sub[1, "watch_id"]
  coord_tmp = sub[, c("E_utm", "N_utm")]
  spot_tmp = sub[1, "spot_num"]
  time_tmp = sub[1, "posix_time"]
  
  
  
  # get coordinates of last spots
  if(option==0){
    spots <- model_data[compID == as.character(compID_tmp) &
                          watch_id == as.numeric(id_tmp) &
                          spot_num < as.numeric(spot_tmp) &
                          tot_catch == option, c("catch","spot_num", "E_utm", "N_utm", "max_time")]
  } else {
    spots <- model_data[compID == as.character(compID_tmp) &
                          watch_id == as.numeric(id_tmp) &
                          spot_num < as.numeric(spot_tmp) &
                          tot_catch >= option, c("catch","spot_num", "E_utm", "N_utm", "max_time")] 
  }
  
  spots<-spots[order(spots$max_time,decreasing=T),]
  
  if(nrow(spots) == 0){
    space = matrix(-99, ncol=1) ##the last success or loss is outside the arena
    time = matrix(-99, ncol=1) ##the last success or loss is outside the arena
  } else {
    if (length(unique(spots$spot_num))>1){
      space = matrix(fastPdist2(as.matrix(coord_tmp), as.matrix(spots[!duplicated(spots$spot_num), c("E_utm", "N_utm")]))[1:min(length(unique(spots$spot_num)),50)], nrow=1)
      time = matrix(round((as.numeric(time_tmp)-spots[!duplicated(spots$spot_num), max_time][1:min(length(unique(spots$spot_num)),50)])/60)+1, nrow=1)
      
    } else {
      space = matrix(fastPdist2(as.matrix(coord_tmp), as.matrix(spots[1, c("E_utm", "N_utm")]))[1:min(length(unique(spots$spot_num)),50)], nrow=1) 
      time = matrix(round((as.numeric(time_tmp)-spots$max_time[1:min(length(unique(spots$spot_num)),50)])/60)+1, nrow=1)
    }
    
    # sometimes missing, if distances very close
    space[is.na(space)]<-.00001
    
  }
  
  return(list(
    space  = space, 
    time = time))
  
  
}

###################################
# extract distance to last spot
###################################


extract_locality_feature <- function(sub, model_data){
  
  # define cpp function for cluster export
  #sourceCpp("functions/fast_dist.cpp")
  
  # subset spot identifier
  compID_tmp = sub[1, "compID"]
  id_tmp = sub[1, "watch_id"]
  coord_tmp = sub[, c("E_utm", "N_utm")]
  spot_tmp = sub[1, "spot_num"]
  
  
  # get coordinates of last spots
  spots <- model_data[compID == as.character(compID_tmp) &
                        watch_id == as.numeric(id_tmp) &
                        spot_num < as.numeric(spot_tmp), c("spot_num", "E_utm", "N_utm")]
  
  
  if(nrow(spots) == 0){
    result = rep(-99, nrow(sub))
  } else {
    result = fastPdist2(as.matrix(coord_tmp), as.matrix(spots[nrow(spots), c("E_utm", "N_utm")]))
  }
  
  # sometimes missing, if distances very close
  result[is.na(result)]<-0
  
  return(result)
  
}



###################################
# extract distance to edge
###################################

extract_sonar_feature <- function(sub, sonars){
  
  result = apply(fastPdist2(as.matrix(sub[,c("E_utm", "N_utm")]), as.matrix(sonars[,c(1,2)])), 1, FUN = min)
  
  return(result)
}



###################################
# extract fish presence
###################################

#sub<-matrix[[29593]]
extract_fish_feature <- function(sub,fish_presence){
  
  compID_tmp = as.character(sub[1, "compID"])
  id_tmp = as.numeric(sub[1, "watch_id"])
  coord_tmp = sub[, c("E_utm", "N_utm")]
  spot_tmp = as.numeric(sub[1, "spot_num"])
  min_tmp = as.numeric(sub[1,"minute"])
  
  cat(paste("Processing compID = ", compID_tmp,
            ", watch = ", id_tmp,
            ", minute = ", min_tmp,"\r"))
  
  
  if(min_tmp<=length(fish_presence[[compID_tmp]])){
    
    fish_3 <- do.call(cbind,lapply(lapply(fish_presence[[compID_tmp]][min_tmp],
                                            function (x) lapply(raster::extract(raster(x),coord_tmp, buffer=3), function(y) mean(y))), 
                                     function(x) do.call(rbind, x)))
  } else {
    fish_3 <- do.call(cbind,lapply(lapply(fish_presence[[compID_tmp]][length(fish_presence[[compID_tmp]])],
                                          function (x) lapply(raster::extract(raster(x),coord_tmp, buffer=3), function(y) mean(y))), 
                                   function(x) do.call(rbind, x)))
  }
  result<-fish_3
  
  return(result)
}



##############################################################################
# LOAD STUFF
##############################################################################
model_data<-fread(here("data","2511_timebin_data.csv"))

matrix<-split(model_data, 1:nrow(model_data))

gps_data<-fread(here("data","2604_gps_data.csv"))

arena<-fread(here::here("data","2408_ArenaExtents.csv"))

sonars<-rbind(arena[,c(6,7)],
              arena[,c(8,9)],
              arena[,c(10,11)],use.names=FALSE)

file_list <- readRDS(here::here("data","2512_rasters_file_paths.rds"))

fish_presence <- list()
base_path <- here::here("data")

for (comp in names(file_list)) {
  full_paths <- file.path(base_path, file_list[[comp]])
  
  # Load rasters from these full paths
  fish_presence[[comp]] <- lapply(full_paths, terra::rast)
}




################################################################################
# extract social distances
################################################################################

social_distance_matrix = lapply(matrix, function(x) extract_social(x, gps_data, time_option = 0))

columns_social<-sapply(social_distance_matrix, ncol)

social_distance_matrix<-lapply(social_distance_matrix, function(x) {pad_matrix(x, 1 ,max(columns_social), -99)})

social_distance_matrix<- lapply(social_distance_matrix, function(x) {
  x <- if (is.list(x)) unlist(x) else x
  matrix(x, nrow = 1)
})



saveRDS(social_distance_matrix, file = here("data","3_spot_leaving","social_distance_matrix.rds"))
saveRDS(columns_social, file = here::here("data","3_spot_leaving","columns_social.rds"))


# ##check
# polygons<-readRDS(here("data","arena_polygons.rds"))
# polygons <- lapply(polygons, terra::unwrap)
# plot(polygons[["som15am"]])
# points(matrix[[1]][1,c("E_utm",  "N_utm")], pch=19, col="red")
# points(gps_data[compID == as.character(matrix[[1]][1,c("compID")]) &
#                   compIDwatch != as.character(matrix[[1]][1,c("compIDwatch")]) &
#                   posix_time == as.numeric(matrix[[1]][1,c("posix_time")]) - 0 &
#                   !is.na(E_utm),c("E_utm","N_utm")], pch=19, col="aquamarine")
# 
# social_distance_matrix[[1]] ##looks right

rm(social_distance_matrix)
gc()
gc()


###############################################################################
# extract gain loss
###############################################################################

success_distance_matrix = lapply(matrix, function(x) extract_spot_distance(x, model_data, option = 1))

loss_distance_matrix = lapply(matrix, function(x) extract_spot_distance(x, model_data, option = 0))

repetitions<-rbindlist(matrix)

success_spatial_matrix  <- lapply(success_distance_matrix, `[[`, "space")
success_time_matrix <- lapply(success_distance_matrix, `[[`, "time")
loss_spatial_matrix  <- lapply(loss_distance_matrix, `[[`, "space")
loss_time_matrix <- lapply(loss_distance_matrix, `[[`, "time")

repetitions[, idx := .I] ##this is used to make all distances within an angling spot the same in every 10sec bins
repetitions[, {
  success_spatial_matrix[idx] <<- rep(list(success_spatial_matrix[[idx[1]]]), .N)
  NULL  # return nothing useful, just keep data.table happy
}, by = .(compIDwatch, spot_num)] ##works

repetitions[, {
  loss_spatial_matrix[idx] <<- rep(list(loss_spatial_matrix[[idx[1]]]), .N)
  NULL  # return nothing useful, just keep data.table happy
}, by = .(compIDwatch, spot_num)] ##works


columns_success<-sapply(success_spatial_matrix, ncol)

columns_loss<-sapply(loss_spatial_matrix, ncol)

success_spatial_matrix<-lapply(success_spatial_matrix, function(x) {pad_matrix(x, 1 ,max(columns_success), -99)})
loss_spatial_matrix<-lapply(loss_spatial_matrix, function(x) {pad_matrix(x, 1 ,max(columns_loss), -99)})

success_time_matrix<-lapply(success_time_matrix, function(x) {pad_matrix(x, 1 ,max(columns_success), -99)})
loss_time_matrix<-lapply(loss_time_matrix, function(x) {pad_matrix(x, 1 ,max(columns_loss), -99)})


# ##check
# polygons<-readRDS(here("data","arena_polygons.rds"))
# polygons <- lapply(polygons, terra::unwrap)
# plot(polygons[[as.character(matrix[[101]][1,c("compID")])]])
# points(matrix[[101]][1,c("E_utm",  "N_utm")], pch=19, col="blue")
# points(model_data[compID == as.character(matrix[[101]][1,c("compID")]) &
#                   compIDwatch == as.character(matrix[[101]][1,c("compIDwatch")]) &
#                   spot_num < as.numeric(matrix[[101]][1,c("spot_num")]) &
#                   catch == 0,c("E_utm","N_utm")], pch=19, col="red")
# points(model_data[compID == as.character(matrix[[101]][1,c("compID")]) &
#                     compIDwatch == as.character(matrix[[101]][1,c("compIDwatch")]) &
#                     spot_num < as.numeric(matrix[[101]][1,c("spot_num")]) &
#                     catch == 1,c("E_utm","N_utm")], pch=19, col="green")
# 
# unique(model_data[compID == as.character(matrix[[101]][1,c("compID")]) &
#              compIDwatch == as.character(matrix[[101]][1,c("compIDwatch")]) &
#              spot_num < as.numeric(matrix[[101]][1,c("spot_num")]) &
#              tot_catch == 0,anglingID])
# unique(model_data[compID == as.character(matrix[[101]][1,c("compID")]) &
#                     compIDwatch == as.character(matrix[[101]][1,c("compIDwatch")]) &
#                     spot_num < as.numeric(matrix[[101]][1,c("spot_num")]) &
#                     tot_catch >= 1,anglingID])
# 
# success_distance_matrix[[101]] ##looks right
# loss_distance_matrix[[101]] ##loooks good

saveRDS(success_spatial_matrix, file = here("data","3_spot_leaving","success_spatial_matrix.rds"))
saveRDS(success_time_matrix, file = here("data","3_spot_leaving","success_time_matrix.rds"))
saveRDS(columns_success, file = here("data","3_spot_leaving","columns_success.rds"))

saveRDS(loss_spatial_matrix, file = here("data","3_spot_leaving","loss_spatial_matrix.rds"))
saveRDS(loss_time_matrix, file = here("data","3_spot_leaving","loss_time_matrix.rds"))
saveRDS(columns_loss, file = here("data","3_spot_leaving","columns_loss.rds"))

rm(success_spatial_matrix)
rm(loss_spatial_matrix)
rm(success_time_matrix)
rm(loss_time_matrix)
gc()
gc()

################################################################################
# extract locality feature
################################################################################
locality_matrix = list()

for (i in 1:length(matrix)){
  locality_matrix[[i]] <- extract_locality_feature(matrix[[i]], model_data)
}

# #check
# polygons<-readRDS(here("data","arena_polygons.rds"))
# polygons <- lapply(polygons, terra::unwrap)
# plot(polygons[[as.character(matrix[[100]][1,c("compID")])]])
# points(matrix[[100]][1,c("E_utm",  "N_utm")], pch=19, col="blue")
# points(model_data[compID == as.character(matrix[[100]][1,c("compID")]) &
#                     compIDwatch == as.character(matrix[[2]][1,c("compIDwatch")]) &
#                     spot_num == as.numeric(matrix[[100]][1,c("spot_num")])-1, c( "E_utm", "N_utm")], pch=19, col="red")
# 
# 
# locality_matrix[[100]] ##looks right

saveRDS(locality_matrix, file = here("data","3_spot_leaving","locality_matrix.rds"))

rm(locality_matrix)
gc()


################################################################################
# extract edge feature
################################################################################

sonar_distance_matrix = lapply(matrix, function(x) extract_sonar_feature(x, sonars))

# hist(unlist(sonar_distance_matrix))
# ###check
# polygons<-readRDS(here("data","arena_polygons.rds"))
# polygons <- lapply(polygons, terra::unwrap)
# plot(polygons[[as.character(matrix[[7]][1,c("compID")])]])
# points(matrix[[7]][1,c("E_utm",  "N_utm")], pch=19, col="blue")
# 
# sonar_distance_matrix[[7]] ##looks right

saveRDS(sonar_distance_matrix, file = here("data","3_spot_leaving","sonar_distance_matrix.rds"))

rm(sonar_distance_matrix)
gc()

################################################################################
# extract fish feature
################################################################################


fish_feature = lapply(matrix, function(x) extract_fish_feature(x, fish_presence)) #this runs overnight


# # ##check
# terra::plot(fish_presence[[as.character(matrix[[3050]][1,c("compID")])]][[as.numeric(matrix[[3050]][1,c("minute")])]])
# points(matrix[[3050]][1,c("E_utm", "N_utm")], pch=19, col="red")
# fish_feature[[3050]] ##looks right
# 
# terra::plot(fish_presence[[as.character(matrix[[2050]][1,c("compID")])]][[as.numeric(matrix[[2050]][1,c("minute")])]])
# points(matrix[[2050]][1,c("E_utm", "N_utm")], pch=19, col="red")
# fish_feature[[2050]] ##looks right

saveRDS(fish_feature, file = here("data","3_spot_leaving","fish_feature.rds"))


################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################


