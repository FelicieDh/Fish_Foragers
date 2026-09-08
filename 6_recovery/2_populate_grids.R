################################################################################
#
# Title: 2. populate grids
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

grid15<-readRDS(here::here("data", "6_recovery", "voronoi_grid_15.rds"))
grid20<-readRDS(here::here("data", "6_recovery", "voronoi_grid_20.rds"))

grid15 <- lapply(grid15, terra::unwrap)

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

min_data[, catch := c(0, diff(cum_catch_at_spot)), by = anglingID]
min_data[is.na(anglingID), catch:=NA]
min_data[, bi_catch := as.integer(catch > 0)]

min_data<-min_data[-which(min_data$anglingID==""),]


min_data$cell_id<-as.numeric(NA)
min_data$cell_area_m2<-as.numeric(NA)

### First i don't do the fish density. I'll only add it for active cells to save processing time

for (i in 1:length(grid15)){ #measure gridcells
  
  for (j in 1:121){
    pts <- terra::vect(as.matrix(min_data[compID==names(grid15)[i] & min==j,c("E_utm", "N_utm")]),  crs=crs(grid15[[i]]))
    
    min_data[compID==names(grid15)[i] & min==j, c("cell_id", "cell_area_m2") := terra::extract(grid15[[i]], pts)[,2:3]]
    
  }
}


##make the storage dataset
for (i in 1:length(grid15)){
  
  temp_data<-cbind(rep(names(grid15)[i],nrow(as.data.frame(grid15[[i]]))*120),rep(1:120, each=nrow(as.data.frame(grid15[[i]]))), 
                   do.call(rbind, replicate(120, as.data.frame(grid15[[i]]), simplify = FALSE)))
  
  colnames(temp_data)[1:2]<-c("compID","minute")
  
  if (i == 1){
    grid_data<-temp_data
  } else {
    grid_data<-rbind(grid_data, temp_data)
  }
  
}



grid_data$angler_presence<-as.numeric(NA)
grid_data$angler_nb<-as.numeric(NA)
grid_data$catch_bi<-as.numeric(NA)
grid_data$catch_abs<-as.numeric(NA)

grid_data<-as.data.table(grid_data)

##populate the dataset
#i=104
for (i in 1:nrow(grid_data)){
  if (min_data[compID==grid_data[i,compID] &
               min==grid_data[i,minute] & 
               cell_id==grid_data[i,cell_id],.N]==0)next
  
  grid_data[i,angler_presence:=1]
  grid_data[i,angler_nr:=length(unique(min_data[compID==grid_data[i,compID] &
                                                  min==grid_data[i,minute] & 
                                                  cell_id==grid_data[i,cell_id],watch_id]))]
  grid_data[i,catch_bi:=as.numeric(any(min_data[compID==grid_data[i,compID] &
                                                  min==grid_data[i,minute] & 
                                                  cell_id==grid_data[i,cell_id],catch]>0))]
  grid_data[i,catch_abs:=sum(min_data[compID==grid_data[i,compID] &
                                        min==grid_data[i,minute] & 
                                        cell_id==grid_data[i,cell_id],catch])]
  
  cat(i/nrow(grid_data)*100,"\r")
  
}

grid_data[catch_abs<0, catch_abs:=0]

# #test
# min_data[compID=="kuo19am" &
#            min==68 & 
#            cell_id==37,] ##here 5 catch among 2 anglers

#### Now I must separate the events with rleid

grid_data[, event_nr := {
  x <- rleid(angler_presence)
  x[is.na(angler_presence)] <- NA
  match(x, unique(na.omit(x)))
}, by = .(compID, cell_id)]

grid_data[!is.na(angler_presence), event_name:=paste0(compID,"_",cell_id,"_",event_nr)]



##add time since last event and whether the event is successful
grid_data$time_away<-as.numeric(NA)


for (i in 1:length(unique(grid_data$event_name))){
  if(grid_data[compID == grid_data[event_name==unique(grid_data$event_name)[i], compID][1] &
               cell_id == grid_data[event_name==unique(grid_data$event_name)[i], cell_id][1] &
               minute < grid_data[event_name==unique(grid_data$event_name)[i], minute][1] &
               !is.na(event_name),.N]==0) {
    grid_data[event_name==unique(grid_data$event_name)[i], time_away := 3*60]
  } else {
    grid_data[event_name==unique(grid_data$event_name)[i], time_away := min(minute) - max(grid_data[compID == grid_data[event_name==unique(grid_data$event_name)[i], compID][1] &
                                                                                                      cell_id == grid_data[event_name==unique(grid_data$event_name)[i], cell_id][1] &
                                                                                                      minute < grid_data[event_name==unique(grid_data$event_name)[i], minute][1] &
                                                                                                      !is.na(event_name),minute])]
  }
}


grid_data[,successful_event:=ifelse(any(catch_bi>0),1,0), by=event_name]

grid_data[!is.na(event_name),cumulative_time :=seq_len(.N), by=.(compID, cell_id)]

grid_data[!is.na(event_name),event_time :=seq_len(.N), by=.(event_name)]

##Add the fish density

file_list <- readRDS(here::here("data","2512_rasters_file_paths.rds"))

fish_presence <- list()
base_path <- here::here("data")

for (comp in names(file_list)) {
  full_paths <- file.path(base_path, file_list[[comp]])
  
  # Load rasters from these full paths
  fish_presence[[comp]] <- lapply(full_paths, terra::rast)
}

grid_data$mean_fish<-as.numeric(NA)
grid_data$min_fish<-as.numeric(NA)
grid_data$max_fish<-as.numeric(NA)

for (i in 1:nrow(grid_data)){
  
  if(is.na(grid_data[i,event_name])) next
  
  grid_data[i,mean_fish := terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]], 
                                      grid15[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = mean, na.rm = TRUE)[,2]]
  
  grid_data[i,min_fish := terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]], 
                                     grid15[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = min, na.rm = TRUE)[,2]]
  
  grid_data[i,max_fish := terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]], 
                                     grid15[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = max, na.rm = TRUE)[,2]]
  cat(i/nrow(grid_data)*100,"\r")
  
}


# ##test
# i=134 #i=500
# 
# par(mfrow=c(1,1))
# plot(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]])
# plot(grid15[[grid_data[i,compID]]][grid_data[i,cell_id]], border="white", add=T)                                                               
# 
# terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]], 
#         grid15[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = mean, na.rm = TRUE)[,2]
# ##great


saveRDS(grid_data, here::here("data","6_recovery","grid_data_15.RDS"))



################################################################################
################################################################################
################################################################################


###################################################################################
################################## 20m ############################################
###################################################################################

# 
# grid20<-readRDS(here::here("data", "6_recovery", "voronoi_grid_20.rds"))
# 
# grid20 <- lapply(grid20, terra::unwrap)
# 
# fulldata<-fread(here("data", "2604_gps_data.csv"))
# 
# min_data<-fulldata[,.(compID = max(compID, na.rm=TRUE),
#                       compIDwatch = max(compIDwatch, na.rm=TRUE),
#                       watch_id = max(watch_id, na.rm=TRUE),
#                       real_time = min(real_time),
#                       posix_time = min(posix_time),
#                       posix_since_start = min(posix_since_start),
#                       E_utm = mean(E_utm),
#                       N_utm = mean(N_utm),
#                       cum_catch_at_spot = max(cum_catch_at_spot, na.rm=TRUE)
# ), by=list(anglingID, min)]
# 
# min_data[,min:=min+1]
# min_data[,success := as.integer(any(cum_catch_at_spot > 0)), by=anglingID]
# 
# min_data[, catch := c(0, diff(cum_catch_at_spot)), by = anglingID]
# min_data[is.na(anglingID), catch:=NA]
# min_data[, bi_catch := as.integer(catch > 0)]
# 
# min_data<-min_data[-which(min_data$anglingID==""),]
# 
# 
# min_data$cell_id<-as.numeric(NA)
# min_data$cell_area_m2<-as.numeric(NA)
# 
# ### First i donÄt do the fish density. I'll only add it for active cells to save processing time
# 
# for (i in 1:length(grid20)){
# 
#   for (j in 1:121){
#     pts <- terra::vect(as.matrix(min_data[compID==names(grid20)[i] & min==j,c("E_utm", "N_utm")]),  crs=crs(grid20[[i]]))
# 
#     min_data[compID==names(grid20)[i] & min==j, c("cell_id", "cell_area_m2") := terra::extract(grid20[[i]], pts)[,2:3]]
# 
#   }
# }
# 
# 
# ##make the storage dataset
# for (i in 1:length(grid20)){
# 
#   temp_data<-cbind(rep(names(grid20)[i],nrow(as.data.frame(grid20[[i]]))*120),rep(1:120, each=nrow(as.data.frame(grid20[[i]]))),
#                    do.call(rbind, replicate(120, as.data.frame(grid20[[i]]), simplify = FALSE)))
# 
#   colnames(temp_data)[1:2]<-c("compID","minute")
# 
#   if (i == 1){
#     grid_data<-temp_data
#   } else {
#     grid_data<-rbind(grid_data, temp_data)
#   }
# 
# }
# 
# 
# 
# grid_data$angler_presence<-as.numeric(NA)
# grid_data$angler_nb<-as.numeric(NA)
# grid_data$catch_bi<-as.numeric(NA)
# grid_data$catch_abs<-as.numeric(NA)
# 
# grid_data<-as.data.table(grid_data)
# ##populate the dataset
# 
# i=104
# for (i in 1:nrow(grid_data)){
#   if (min_data[compID==grid_data[i,compID] &
#                min==grid_data[i,minute] &
#                cell_id==grid_data[i,cell_id],.N]==0)next
# 
#   grid_data[i,angler_presence:=1]
#   grid_data[i,angler_nr:=length(unique(min_data[compID==grid_data[i,compID] &
#                                                   min==grid_data[i,minute] &
#                                                   cell_id==grid_data[i,cell_id],watch_id]))]
#   grid_data[i,catch_bi:=as.numeric(any(min_data[compID==grid_data[i,compID] &
#                                                   min==grid_data[i,minute] &
#                                                   cell_id==grid_data[i,cell_id],catch]>0))]
#   grid_data[i,catch_abs:=sum(min_data[compID==grid_data[i,compID] &
#                                         min==grid_data[i,minute] &
#                                         cell_id==grid_data[i,cell_id],catch])]
# 
#   cat(i/nrow(grid_data)*100,"\r")
# 
# }
# 
# grid_data[catch_abs<0, catch_abs:=0]
# #test
# min_data[compID=="kuo19am" &
#            min==68 &
#            cell_id==22,] ##here the 5 catch among 2 anglers are split into two cells (15 and 22)!
# 
# #### Now I must separate the events with rleid
# 
# grid_data[, event_nr := {
#   x <- rleid(angler_presence)
#   x[is.na(angler_presence)] <- NA
#   match(x, unique(na.omit(x)))
# }, by = .(compID, cell_id)]
# 
# grid_data[!is.na(angler_presence), event_name:=paste0(compID,"_",cell_id,"_",event_nr)]
# 
# 
# 
# ##add time since last event and whether the event is successful
# grid_data$time_away<-as.numeric(NA)
# 
# 
# for (i in 1:length(unique(grid_data$event_name))){
#   if(grid_data[compID == grid_data[event_name==unique(grid_data$event_name)[i], compID][1] &
#                cell_id == grid_data[event_name==unique(grid_data$event_name)[i], cell_id][1] &
#                minute < grid_data[event_name==unique(grid_data$event_name)[i], minute][1] &
#                !is.na(event_name),.N]==0) {
#     grid_data[event_name==unique(grid_data$event_name)[i], time_away := 3*60]
#   } else {
#     grid_data[event_name==unique(grid_data$event_name)[i], time_away := min(minute) - max(grid_data[compID == grid_data[event_name==unique(grid_data$event_name)[i], compID][1] &
#                                                                                                       cell_id == grid_data[event_name==unique(grid_data$event_name)[i], cell_id][1] &
#                                                                                                       minute < grid_data[event_name==unique(grid_data$event_name)[i], minute][1] &
#                                                                                                       !is.na(event_name),minute])]
#   }
# }
# 
# 
# grid_data[,successful_event:=ifelse(any(catch_bi>0),1,0), by=event_name]
# 
# grid_data[!is.na(event_name),cumulative_time :=seq_len(.N), by=.(compID, cell_id)]
# 
# grid_data[!is.na(event_name),event_time :=seq_len(.N), by=.(event_name)]
# 
# ##Add the fish density
# 
# file_list <- readRDS(here::here("data","2512_rasters_file_paths.rds"))
# 
# fish_presence <- list()
# base_path <- here::here("data")
# 
# for (comp in names(file_list)) {
#   full_paths <- file.path(base_path, file_list[[comp]])
# 
#   # Load rasters from these full paths
#   fish_presence[[comp]] <- lapply(full_paths, terra::rast)
# }
# 
# grid_data$mean_fish<-as.numeric(NA)
# grid_data$min_fish<-as.numeric(NA)
# grid_data$max_fish<-as.numeric(NA)
# 
# for (i in 1:nrow(grid_data)){
# 
#   if(is.na(grid_data[i,event_name])) next
# 
#   grid_data[i,mean_fish := terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]],
#                                       grid20[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = mean, na.rm = TRUE)[,2]]
# 
#   grid_data[i,min_fish := terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]],
#                                      grid20[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = min, na.rm = TRUE)[,2]]
# 
#   grid_data[i,max_fish := terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]],
#                                      grid20[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = max, na.rm = TRUE)[,2]]
#   cat(i/nrow(grid_data)*100,"\r")
# 
# }
# 
# 
# ##test
# i=134 #i=500
# 
# plot(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]])
# plot(grid20[[grid_data[i,compID]]][grid_data[i,cell_id]], border="white", add=T)
# 
# terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]],
#         grid20[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = mean, na.rm = TRUE)[,2]
# ##great
# 
# 
# saveRDS(grid_data, here::here("data","6_recovery","grid_data_20.RDS"))






###################################################################################
################################## 10m ############################################
###################################################################################

# 
# grid10<-readRDS(here::here("data", "6_recovery", "voronoi_grid_10.rds"))
# 
# grid10 <- lapply(grid10, terra::unwrap)
# 
# fulldata<-fread(here("data", "2604_gps_data.csv"))
# 
# min_data<-fulldata[,.(compID = max(compID, na.rm=TRUE),
#                       compIDwatch = max(compIDwatch, na.rm=TRUE),
#                       watch_id = max(watch_id, na.rm=TRUE),
#                       real_time = min(real_time),
#                       posix_time = min(posix_time),
#                       posix_since_start = min(posix_since_start),
#                       E_utm = mean(E_utm),
#                       N_utm = mean(N_utm),
#                       cum_catch_at_spot = max(cum_catch_at_spot, na.rm=TRUE)
# ), by=list(anglingID, min)]
# 
# min_data[,min:=min+1]
# min_data[,success := as.integer(any(cum_catch_at_spot > 0)), by=anglingID]
# 
# min_data[, catch := c(0, diff(cum_catch_at_spot)), by = anglingID]
# min_data[is.na(anglingID), catch:=NA]
# min_data[, bi_catch := as.integer(catch > 0)]
# 
# min_data<-min_data[-which(min_data$anglingID==""),]
# 
# 
# min_data$cell_id<-as.numeric(NA)
# min_data$cell_area_m2<-as.numeric(NA)
# 
# ### First i donÄt do the fish density. I'll only add it for active cells to save processing time
# 
# for (i in 1:length(grid10)){
# 
#   for (j in 1:121){
#     pts <- terra::vect(as.matrix(min_data[compID==names(grid10)[i] & min==j,c("E_utm", "N_utm")]),  crs=crs(grid10[[i]]))
# 
#     min_data[compID==names(grid10)[i] & min==j, c("cell_id", "cell_area_m2") := terra::extract(grid10[[i]], pts)[,2:3]]
# 
#   }
# }
# 
# 
# ##make the storage dataset
# for (i in 1:length(grid10)){
# 
#   temp_data<-cbind(rep(names(grid10)[i],nrow(as.data.frame(grid10[[i]]))*120),rep(1:120, each=nrow(as.data.frame(grid10[[i]]))),
#                    do.call(rbind, replicate(120, as.data.frame(grid10[[i]]), simplify = FALSE)))
# 
#   colnames(temp_data)[1:2]<-c("compID","minute")
# 
#   if (i == 1){
#     grid_data<-temp_data
#   } else {
#     grid_data<-rbind(grid_data, temp_data)
#   }
# 
# }
# 
# 
# 
# grid_data$angler_presence<-as.numeric(NA)
# grid_data$angler_nb<-as.numeric(NA)
# grid_data$catch_bi<-as.numeric(NA)
# grid_data$catch_abs<-as.numeric(NA)
# 
# grid_data<-as.data.table(grid_data)
# ##populate the dataset
# 
# i=104
# for (i in 1:nrow(grid_data)){
#   if (min_data[compID==grid_data[i,compID] &
#                min==grid_data[i,minute] &
#                cell_id==grid_data[i,cell_id],.N]==0)next
# 
#   grid_data[i,angler_presence:=1]
#   grid_data[i,angler_nr:=length(unique(min_data[compID==grid_data[i,compID] &
#                                                   min==grid_data[i,minute] &
#                                                   cell_id==grid_data[i,cell_id],watch_id]))]
#   grid_data[i,catch_bi:=as.numeric(any(min_data[compID==grid_data[i,compID] &
#                                                   min==grid_data[i,minute] &
#                                                   cell_id==grid_data[i,cell_id],catch]>0))]
#   grid_data[i,catch_abs:=sum(min_data[compID==grid_data[i,compID] &
#                                         min==grid_data[i,minute] &
#                                         cell_id==grid_data[i,cell_id],catch])]
# 
#   cat(i/nrow(grid_data)*100,"\r")
# 
# }
# 
# grid_data[catch_abs<0, catch_abs:=0]
# #test
# min_data[compID=="kuo19am" &
#            min==68 &
#            cell_id==22,] ##here the 5 catch among 2 anglers are split into two cells (15 and 22)!
# 
# #### Now I must separate the events with rleid
# 
# grid_data[, event_nr := {
#   x <- rleid(angler_presence)
#   x[is.na(angler_presence)] <- NA
#   match(x, unique(na.omit(x)))
# }, by = .(compID, cell_id)]
# 
# grid_data[!is.na(angler_presence), event_name:=paste0(compID,"_",cell_id,"_",event_nr)]
# 
# 
# 
# ##add time since last event and whether the event is successful
# grid_data$time_away<-as.numeric(NA)
# 
# 
# for (i in 1:length(unique(grid_data$event_name))){
#   if(grid_data[compID == grid_data[event_name==unique(grid_data$event_name)[i], compID][1] &
#                cell_id == grid_data[event_name==unique(grid_data$event_name)[i], cell_id][1] &
#                minute < grid_data[event_name==unique(grid_data$event_name)[i], minute][1] &
#                !is.na(event_name),.N]==0) {
#     grid_data[event_name==unique(grid_data$event_name)[i], time_away := 3*60]
#   } else {
#     grid_data[event_name==unique(grid_data$event_name)[i], time_away := min(minute) - max(grid_data[compID == grid_data[event_name==unique(grid_data$event_name)[i], compID][1] &
#                                                                                                       cell_id == grid_data[event_name==unique(grid_data$event_name)[i], cell_id][1] &
#                                                                                                       minute < grid_data[event_name==unique(grid_data$event_name)[i], minute][1] &
#                                                                                                       !is.na(event_name),minute])]
#   }
# }
# 
# 
# grid_data[,successful_event:=ifelse(any(catch_bi>0),1,0), by=event_name]
# 
# grid_data[!is.na(event_name),cumulative_time :=seq_len(.N), by=.(compID, cell_id)]
# 
# grid_data[!is.na(event_name),event_time :=seq_len(.N), by=.(event_name)]
# 
# ##Add the fish density
# 
# file_list <- readRDS(here::here("data","2512_rasters_file_paths.rds"))
# 
# fish_presence <- list()
# base_path <- here::here("data")
# 
# for (comp in names(file_list)) {
#   full_paths <- file.path(base_path, file_list[[comp]])
# 
#   # Load rasters from these full paths
#   fish_presence[[comp]] <- lapply(full_paths, terra::rast)
# }
# 
# grid_data$mean_fish<-as.numeric(NA)
# grid_data$min_fish<-as.numeric(NA)
# grid_data$max_fish<-as.numeric(NA)
# 
# for (i in 1:nrow(grid_data)){
# 
#   if(is.na(grid_data[i,event_name])) next
# 
#   grid_data[i,mean_fish := terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]],
#                                       grid10[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = mean, na.rm = TRUE)[,2]]
# 
#   grid_data[i,min_fish := terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]],
#                                      grid10[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = min, na.rm = TRUE)[,2]]
# 
#   grid_data[i,max_fish := terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]],
#                                      grid10[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = max, na.rm = TRUE)[,2]]
#   cat(i/nrow(grid_data)*100,"\r")
# 
# }
# 
# 
# ##test
# i=134 #i=500
# 
# plot(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]])
# plot(grid10[[grid_data[i,compID]]][grid_data[i,cell_id]], border="white", add=T)
# 
# terra::extract(fish_presence[[grid_data[i,compID]]][[grid_data[i,minute]]],
#         grid10[[grid_data[i,compID]]][grid_data[i,cell_id]], fun = mean, na.rm = TRUE)[,2]
# ##great
# 
# 
# saveRDS(grid_data, here::here("data","6_recovery","grid_data_10.RDS"))


################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################
