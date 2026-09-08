################################################################################
#
# Title: 1. Preparing the voronoi grid
#
# Author F Dhellemmes
#
# last update: 08/08/2026 (DD/MM/YYYY)
#
################################################################################



library(data.table)
library(terra)
library(here)
library(sf)

# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################



arena<-fread(here::here("data","2408_ArenaExtents.csv"))
arena$compID<-paste0(tolower(substr(arena$compID, 1,3)), substr(arena$compID,nchar(arena$compID)-3,nchar(arena$compID)))

sonars<-rbind(arena[,c(1,6,7)],
              arena[,c(1,8,9)],
              arena[,c(1,10,11)],use.names=FALSE)

cellsizes<-c(10,15,20)

for (i in 1:length(cellsizes)){
  
  grid_list<-list()
  
  for (j in 1:length(unique(sonars$compID))){
    
    ind<-which(sonars$compID==unique(sonars$compID)[j])
    
    
    pt1 <- terra::vect(
      sonars[ind[1], c("xpink","ypink")],c("xpink","ypink"),
      crs="+proj=utm +zone=35 +datum=WGS84 +units=m")
    outer1 <- buffer(pt1, 35)
    
    pt2 <- terra::vect(
      sonars[ind[2], c("xpink","ypink")],c("xpink","ypink"),
      crs="+proj=utm +zone=35 +datum=WGS84 +units=m")
    outer2 <- buffer(pt2, 35)
    
    pt3 <- terra::vect(
      sonars[ind[3], c("xpink","ypink")],c("xpink","ypink"),
      crs="+proj=utm +zone=35 +datum=WGS84 +units=m")
    outer3 <- buffer(pt3, 35)
    
    all_outer<-rbind(outer1, outer2, outer3)
    merged <- terra::aggregate(all_outer)
    
    merged_sf <- sf::st_as_sf(merged)
    
    res <- cellsizes[i]  # approximate spacing between Voronoi seeds
    pts <- sf::st_make_grid(merged_sf, cellsize = res, what = "centers")
    pts <- sf::st_intersection(pts, merged_sf)  # only points inside the polygon
    
    # convert to sf points
    pts_sf <- sf::st_sf(geometry = pts)
    
    # compute Voronoi polygons
    voronoi_sf <- sf::st_voronoi(do.call(sf::st_union, list(pts_sf)))
    voronoi_sf <- sf::st_collection_extract(voronoi_sf, "POLYGON")
    
    # clip Voronoi cells to merged polygon
    voronoi_sf <- sf::st_intersection(voronoi_sf, merged_sf)
    
    # convert back to terra
    voronoi_terra <- terra::vect(voronoi_sf)
    voronoi_single <- terra::as.polygons(voronoi_terra)
    voronoi_single$cell_id <- 1:length(voronoi_single)
    voronoi_single$area_m2 <- expanse(voronoi_single, unit = "m") 
    
    grid_list<-as.list(append(grid_list,  assign(unique(sonars$compID)[j], list(voronoi_single)))) ##add the UD to my fish list
    names(grid_list)[length(grid_list)] = unique(sonars$compID)[j] #name the list accordi
    # grid_list[[j]]<-voronoi_terra
    # names(grid_list)[j]<-unique(sonars$compID)[j]
  }
  
  grid_list_wrap<- lapply(grid_list, terra::wrap)
  saveRDS(grid_list_wrap, here::here("data", "6_recovery", paste0("voronoi_grid_",cellsizes[i],".rds")))
  rm(grid_list)
  gc()
}

# grid15<-readRDS(here::here("data", "6_recovery", "voronoi_grid_10.rds"))
# grid15<-readRDS(here::here("data", "6_recovery", "voronoi_grid_15.rds"))
# grid20<-readRDS(here::here("data", "6_recovery", "voronoi_grid_20.rds"))



################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################

