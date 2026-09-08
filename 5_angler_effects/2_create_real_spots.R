################################################################################
#
# Title: 2. Preparing the dataset of real spots
#
# Author F Dhellemmes
#
# last update: 06/08/2026 (DD/MM/YYYY)
#
################################################################################

library(here)
library(data.table)
library(terra)
source(here("functions","library.R"))

# source helper functions
sourceCpp(here("functions","fast_dist.cpp"))
sourceCpp(here("functions","ptinpoly.cpp"))

# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################

##open GPS data
fulldata<-fread(here("data", "2604_gps_data.csv"))

##Open fish data
file_list <- readRDS(here::here("data","2512_rasters_file_paths.rds"))

if(!exists("fish_presence")){
  fish_presence <- list()
  base_path <- here::here("data")
  for (comp in names(file_list)) {
    full_paths <- file.path(base_path, file_list[[comp]])
    
    # Load rasters from these full paths
    fish_presence[[comp]] <- lapply(full_paths, terra::rast)
  }}


min_data<-fulldata[,.(compID = max(compID, na.rm=TRUE),
                      compIDwatch = max(compIDwatch, na.rm=TRUE),
                      watch_id = max(watch_id, na.rm=TRUE),
                      real_time = min(real_time),
                      posix_time = min(posix_time),
                      posix_since_start = min(posix_since_start),
                      E_utm = mean(E_utm),
                      N_utm = mean(N_utm),
                      cum_catch_at_spot = max(cum_catch_at_spot, na.rm=TRUE)#,
), by=list(anglingID, min)]

min_data[,min:=min+1]
min_data[,success := as.integer(any(cum_catch_at_spot > 0)), by=anglingID] ##is the spot successful

# Create the 'catch' column
min_data[, catch := c(0, diff(cum_catch_at_spot)), by = anglingID]
min_data[is.na(anglingID), catch:=NA]
min_data[, bi_catch := as.integer(catch > 0)] #catch per minute

##subset by anglingID
real_data<-min_data[-which(min_data$anglingID==""),] ##remove moments where peope aren't angling
real_data[, N_utm := mean(N_utm), by = anglingID]
real_data[, E_utm := mean(E_utm), by = anglingID]

real_data$minute<-real_data$min+1
real_data[, fish_3 := {
  cat(sprintf("Processing compID = %s, minute = %s (%d pts)\r", 
              compID[1], minute[1], .N))
  
  if(minute[1]<121){
  r <- raster::raster(fish_presence[[as.character(compID[1])]][[minute[1]]])
  } else {
    r <- raster::raster(fish_presence[[as.character(compID[1])]][[120]])
    
  }
  
  vals <- raster::extract(r, cbind(E_utm, N_utm), buffer = 3, fun = mean)
  
  rm(r)
  vals
}, by = .(compID, min)]


## now we fill the social data
real_data$comp_in_10m<-as.numeric(NA)
real_data$nearest_social<-as.numeric(NA)
real_data$nearest_social_capt<-as.numeric(NA)
#real_data$nearest_social_capt<-NA


for (i in 1:nrow(real_data)){
  
  if(real_data[compID==real_data[i,compID] & min ==  real_data[i,min] & watch_id != real_data[i,watch_id], .N] ==0) next
  
  cat(i/8704*100,"\r")
  
  
  real_data[i,nearest_social := apply(fastPdist2(as.matrix(real_data[i,c("E_utm","N_utm")]),
                                                 as.matrix(real_data[compID==real_data[i,compID] & min == real_data[i,min] & watch_id != real_data[i,watch_id], c("E_utm", "N_utm")])
  ), 1, FUN = min, na.rm = TRUE)]
  
  
  real_data[i,comp_in_10m := apply(fastPdist2(as.matrix(real_data[i,c("E_utm","N_utm")]),
                                              as.matrix(real_data[compID==real_data[i,compID] & min == real_data[i,min] & watch_id != real_data[i,watch_id], c("E_utm", "N_utm")])
  ), 1, function(x) length(which(x<=10)))]
  
  if(real_data[compID==real_data[i,compID] & min == real_data[i,min] & watch_id != real_data[i,watch_id] & cum_catch_at_spot>0, .N] ==0) {
    real_data[i,nearest_social_capt := NA]
  } else {
    real_data[i,nearest_social_capt:=apply(fastPdist2(as.matrix(real_data[i,c("E_utm","N_utm")]),
                                                      as.matrix(real_data[compID==real_data[i,compID] & min == real_data[i,min] & watch_id != real_data[i,watch_id] & cum_catch_at_spot>0, c("E_utm", "N_utm")])
    ), 1, FUN = min, na.rm = TRUE)]
  }
  
}


# #check
# i=2013
# plot(fish_presence[[real_data[i,compID]]][[real_data[i,min]]])
# #plot(polygons[[real_data[i,compID]]], add=T)
# points(real_data[i,c("E_utm","N_utm")], pch=19, col="purple")
# points(real_data[compID==real_data[i,compID] & min == real_data[i,min] & watch_id != real_data[i,watch_id], c("E_utm", "N_utm")], pch=19, col="royalblue")
# points(real_data[compID==real_data[i,compID] & min == real_data[i,min] & watch_id != real_data[i,watch_id] & cum_catch_at_spot>0, c("E_utm", "N_utm")], pch=19, col="blue")
# real_data[i,]
# ##it works!!

## Distance to sonar
arena<-fread(here::here("data","2408_ArenaExtents.csv"))

sonars<-rbind(arena[,c(6,7)],
              arena[,c(8,9)],
              arena[,c(10,11)],use.names=FALSE)

i=1

for (i in 1:nrow(real_data)){
  
  real_data[i,"dist_sonar"]<-apply(fastPdist2(as.matrix(real_data[i,c("E_utm", "N_utm")]), as.matrix(sonars[,c(1,2)])), 1, FUN = min)
  
}



saveRDS(real_data, file = here("data","5_angler_effects","RealSpots.rds"))


################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################
