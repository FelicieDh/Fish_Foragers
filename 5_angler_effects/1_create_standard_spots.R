################################################################################
#
# Title: 1. Preparing the 100 control spots
#
# Author F Dhellemmes
#
# last update: 05/08/2026 (DD/MM/YYYY)
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


# Import Lake Data
polygons<-readRDS(here("data","arena_polygons.rds"))

polygons <- lapply(polygons, terra::unwrap)


## Create the sampling spots
i=1
comp_list<-list()

for (i in 1:length(polygons)){
  lake = polygons[[i]]
  
  r <- rast(lake, resolution = 10)   # 10 m x 10 m grid
  values(r) <- 1
  r <- mask(r, lake)
  pts <- as.points(r)
  
  simulated_points = terra::intersect(pts, lake)
  simulated_points = geom(simulated_points )
  simulated_points = as.data.frame(simulated_points)
  # simulate randomly from lake
  
  simulated_points = as.data.frame(simulated_points)
  
  simulated_points$nearest_social<-NA
  
  sim_list <- rep(list(simulated_points), 121)
  
  comp_list[[i]]<-sim_list
  names(comp_list)[i]<-names(polygons)[i]
  
}


##load the fish data

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


###############################################################################
### Extract fish presence for each spot on the comp_list
###############################################################################

for (i in 1:length(comp_list)){
  for (j in 1:length(comp_list[[i]])){
    if(j > length(fish_presence[[names(comp_list)[i]]]))next
    cat(paste(names(comp_list)[i],j),"\r")
    comp_list[[i]][[j]]$fish_3<-rowMeans(do.call(rbind,raster::extract(fish_presence_raster[[names(comp_list)[i]]][[j]],  comp_list[[i]][[j]][,3:4], buffer=3)))
  }}


#############################################################################
# Now the social variables
#############################################################################

fulldata<-fread(here("data", "2604_gps_data.csv"))

min_data<-fulldata[,.(compID = max(compID, na.rm=TRUE),
                      compIDwatch = max(compIDwatch, na.rm=TRUE),
                      watch_id = max(watch_id, na.rm=TRUE),
                      real_time = min(real_time),
                      posix_time = min(posix_time),
                      posix_since_start = min(posix_since_start),
                      E_utm = mean(E_utm),
                      N_utm = mean(N_utm),
                      cum_catch_at_spot = max(cum_catch_at_spot, na.rm=TRUE)
), by=list(anglingID, min)]

min_data[,min:=min+1]
min_data[,success := as.integer(any(cum_catch_at_spot > 0)), by=anglingID]

# Create the 'catch' column
min_data[, catch := c(0, diff(cum_catch_at_spot)), by = anglingID]
min_data[is.na(anglingID), catch:=NA]
min_data[, bi_catch := as.integer(catch > 0)]


## We create 3 columns, one with the nearest competitor, one with (0,1) is there a competitor within 10m, one with the nearest successful compeptitor
for (i in 1:length(comp_list)){
  for (j in 1:length(comp_list[[i]])){
    
    if(min_data[compID==names(comp_list)[i] & min == j & !is.na(anglingID), .N] ==0) next
    cat(paste(names(comp_list)[i],j),"\r")
    
    comp_list[[i]][[j]]$nearest_social<-apply(fastPdist2(as.matrix(comp_list[[i]][[j]][, 3:4]),
                                                         as.matrix(min_data[compID==names(comp_list)[i] & min == j & !is.na(anglingID), c("E_utm", "N_utm")])
    ), 1, FUN = min, na.rm = TRUE)
    
    
    comp_list[[i]][[j]]$comp_in_10m<-apply(fastPdist2(as.matrix(comp_list[[i]][[j]][, 3:4]),
                                                      as.matrix(min_data[compID==names(comp_list)[i] & min == j & !is.na(anglingID), c("E_utm", "N_utm")])
    ), 1, function(x) length(which(x<=10)))
    
    if(min_data[compID==names(comp_list)[i] & min == j  & !is.na(anglingID) & cum_catch_at_spot>0, .N] ==0) {
      comp_list[[i]][[j]]$nearest_social_capt<-NA
    } else {
      comp_list[[i]][[j]]$nearest_social_capt<-apply(fastPdist2(as.matrix(comp_list[[i]][[j]][, 3:4]),
                                                                as.matrix(min_data[compID==names(comp_list)[i] & min == j & !is.na(anglingID) & cum_catch_at_spot>0, c("E_utm", "N_utm")])
      ), 1, FUN = min, na.rm = TRUE)}
    
  }}



#############################################################################
# Then the sonar distances
#############################################################################


arena<-fread(here::here("data","2408_ArenaExtents.csv"))

sonars<-rbind(arena[,c(6,7)],
              arena[,c(8,9)],
              arena[,c(10,11)],use.names=FALSE)



extract_sonar_feature <- function(sub, sonars){
  
  sub[,"dist_sonar"]<-matrix(apply(fastPdist2(as.matrix(sub[,c("x", "y")]), as.matrix(sonars[,c(1,2)])), 1, FUN = min), ncol=1)
  
  return(sub)
}

comp_list <- lapply(comp_list, function(inner_list) {
  lapply(inner_list, function(df) {
    extract_sonar_feature(df, sonars)
  })
})


saveRDS(comp_list, file = here("data","5_angler_effects","StandardSpots.rds"))

################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################

