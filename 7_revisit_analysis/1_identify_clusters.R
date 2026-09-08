################################################################################
#
# Title: 1. Identify revisits
#
# Author F Dhellemmes
#
# last update: 09/08/2026 (DD/MM/YYYY)
#
################################################################################

library(data.table)
library(geosphere)  
library(dbscan)     
library(ggplot2)
library(here)
library(terra)
library(brglm2)
library(terra)
library(brms)

source(here("functions","library.R"))


# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################


data <- fread(here("data","2511_angling_events.csv"))

# Import Lake Data
polygons<-readRDS(here("data","arena_polygons.rds"))
polygons <- lapply(polygons, terra::unwrap)

# Threshold
EPS_M <- 10     # meters "same spot" diameter
MIN_PTS <- 1   # 1 = every spot becomes at least its own cluster


data[,start_min:=floor(posix_since_start/60)+1]
data[,stop_min:=floor((posix_since_start+time_at_spot_s)/60)-1]
data[,success:=ifelse(max_catch_at_spot>0,1,0)]
data[, prev_success := shift(success, 1L), by = compIDwatch]

#max(data$stop_min)

file_list <- readRDS(here::here("data","2512_rasters_file_paths.rds"))

if(!exists("fish_presence")){
  fish_presence <- list()
  base_path <- here::here("data")
  for (comp in names(file_list)) {
    full_paths <- file.path(base_path, file_list[[comp]])
    
    # Load rasters from these full paths
    fish_presence[[comp]] <- lapply(full_paths, terra::rast)
  }}

i <- 0
n <- sum(lengths(fish_presence))
if(!exists("fish_presence_raster")){
  fish_presence_raster <- lapply(fish_presence, function(sublist) {
    lapply(sublist, function(r) {
      i <<- i + 1
      cat("\rConversion", round(i/n*100), "%")
      flush.console()
      raster::raster(r)})  # Converts only the first layer
  })}


# Fish presence when anglers arrive and leave
for (i in 1:nrow(data)){
  data[i, fish_start:=mean(raster::extract(fish_presence_raster[[data[i,compID]]][[data[i,start_min]]], data[i,c("E_utm","N_utm")], buffer=3)[[1]])]
  data[i, fish_stop:=mean(raster::extract(fish_presence_raster[[data[i,compID]]][[data[i,stop_min]]], data[i,c("E_utm","N_utm")], buffer=3)[[1]])]
}

data[, prev_fish := shift(fish_stop, 1L), by = compIDwatch]



# Generate the dataset
data[,is_revisit := FALSE]
data[,will_revisit := FALSE]
data[,success_in_cluster := 0]
data[, final_cluster_id := 0]   


for (i in 1:51){
  data[,paste0("new",i):=as.numeric(NA)]
} ##prepare columns for each new cluster
data_it<-data[1,]

for(j in 1:length(unique(data$compIDwatch))){
  it<-data[compIDwatch==unique(data$compIDwatch)[j],]
  for (i in 2:nrow(it)){
    coords <- it[1:i,.(E_utm, N_utm)]
    clust<-dbscan(coords, eps = EPS_M, minPts = MIN_PTS)$cluster
    if(any(clust[1:i-1] == clust[i] & clust[i-1] != clust[i])){
      it[i,is_revisit := TRUE]
      it[which(clust[1:i-1] == clust[i]), will_revisit := TRUE]
    }
    it[1:i,paste0("new",i):=clust]
    
    it[i, success_in_cluster:=ifelse(any(it[1:i,success]==1),1,0)]
    
  }
  it[, final_cluster_id:= clust]   
  data_it<-rbind(data_it, it)
}


data_it<-data_it[-1,]



saveRDS(data_it, here::here("data","7_revisit_analysis","revisits_10.RDS"))

################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################
