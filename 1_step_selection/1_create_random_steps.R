################################################################################
#
# Title: 1. Simulate available angling spots step selection
#
# Authors: Dhellemmes F. 
#
# Last updated: 31/07/2026 (DD/MM/YYYY)
#
################################################################################

library(here)
library(data.table)
library(terra)

# source helper functions
source(here("functions","library.R"))
sourceCpp(here("functions","fast_dist.cpp"))
sourceCpp(here("functions","ptinpoly.cpp"))


# set the seed

set.seed(42)

################################################################################





# Import arena boudaries
polygons<-readRDS(here("data","arena_polygons.rds"))

polygons <- lapply(polygons, terra::unwrap)


# Import Spot Choices
model_data = fread(here("data","2511_angling_events.csv"))

metadata<-fread(here::here("data","metadata_anon.csv"))


# Add participant ID and Competition ID to model_data
model_data[,compID_num:=metadata[match(model_data$compID,metadata$compID),"compID_num"]]
model_data[,partID_num:=metadata[match(model_data$compIDwatch,metadata$compIDwatch),"partID"]]




# Import GPS Data
gps_data <- fread(here("data","2604_gps_data.csv"))


# Add participant ID and Competition ID to gps_data
gps_data[,compID_num:=metadata[match(gps_data$compID,metadata$compID),"compID_num"]]
gps_data[,partID_num:=metadata[match(gps_data$compIDwatch,metadata$compIDwatch),"partID"]]



# how many spots to analyze
length(unique(paste(model_data$compIDwatch[model_data$spot_num > 1], model_data$spot_num[model_data$spot_num > 1])))



# expand model data by N_spots
N_spots = 30

# simulate coordinates for new spots
T = nrow(model_data)

model_data$chosen = 1
model_data$total_spot_id = 1:nrow(model_data)


# compute starting point of each competition (where anglers started walking from. All were in a group)
starting_point = data.table(compID = unique(gps_data$compID))

for (i in 1:nrow(starting_point)){
  starting_point[i, x := mean(gps_data[compID == starting_point[i,compID] &
                                         posix_time == min(gps_data[compID == starting_point[i,compID],posix_time]),E_utm])]
  starting_point[i, y := mean(gps_data[compID == starting_point[i,compID] &
                                         posix_time == min(gps_data[compID == starting_point[i,compID],posix_time]),N_utm])]
}




# preallocate features
Feature_matrix = list()

# start simulation
# initialize step length distribution
steps = model_data$step[model_data$step > -99]
steps = ifelse(steps == 0, .1, steps)
angles = model_data$angle_d[model_data$angle_d>-99]




# start sim
time_0 <- Sys.time()
#t=2
for (t in 1:T){#1:T T = 
  
  # get spot data
  chosen_spot = model_data[t,]
  
  # subset all choices of ID
  choices_id = model_data[model_data$compIDwatch == chosen_spot$compIDwatch, ]
  first_choice = min(choices_id$spot_num)
  
  # lake surface as polygon
  lake = polygons[[which(names(polygons)==chosen_spot$compID)]]
  
  social_spot<-gps_data[compID==chosen_spot$compID & 
                           posix_time == chosen_spot$posix_time
                         & watch_id != chosen_spot$watch_id,] #we make sure simulated spot don't fall on top of another angler
  
  # min_dist<-rnorm(1, mean=5, sd=.5)
  # min_dist<-max(5,min_dist)
  #  min_dist<-5
  # simulate coordinates of alternative spots
  if (t == 1 || chosen_spot$compIDwatch != model_data[t-1, "compIDwatch"]){
    
    # simulate randomly from lake
    simulated_points = spatSample(lake, N_spots)
    simulated_points = geom(simulated_points)
    
    simulated_x = simulated_points[,3]
    simulated_y = simulated_points[,4]
    
  } else {
    
    # compute heading angle:
    if(chosen_spot$spot_num == 2){
      
      heading = heading_angle(x = as.matrix(starting_point[compID == model_data[t, "compID"], c("x","y")]), 
                              y = as.matrix(model_data[t-1,c("E_utm", "N_utm")]))
      
    } else {
      
      heading = heading_angle(x = as.matrix(model_data[t-2, c("E_utm", "N_utm")]), 
                              y = as.matrix(model_data[t-1,c("E_utm", "N_utm")]))
      
    }
    
    # if simulated spots are not in competition area
    radius = sample(steps[-which(is.na(steps))], N_spots, replace = T)
    
    # simulate turning angle from empirical turning angle distribution
    turning_angle = sample(angles[-which(is.na(angles))], N_spots, replace = F)
    
    # add to coordinates of last spot
    simulated_x = pull(model_data[t-1, "E_utm"]) + radius * cos(heading + turning_angle)
    simulated_y = pull(model_data[t-1, "N_utm"]) + radius * sin(heading + turning_angle)
    
    repeat {
      
      simulated <- cbind(simulated_x, simulated_y)
      
      # dist <- fastPdist2(
      #   as.matrix(simulated),
      #   as.matrix(social_spot[, 7:8])
      # )
      
      #bad_dist <- apply(dist < min_dist, 1, any)
      
      inside <- relate(
        lake,
        terra::vect(simulated, crs="+proj=utm +zone=35 +datum=WGS84 +units=m"),
        "intersects"
      ) # we make sure the new spot falls inside the arena boudaries
      
      bad <- which( !inside) #bad_dist |
      
      if (length(bad) == 0) break
      
      rad <- sample(steps[!is.na(steps)], length(bad), replace = TRUE)
      tur <- sample(angles[!is.na(angles)], length(bad), replace = FALSE)
      
      simulated_x[bad] <- pull(model_data[t-1, "E_utm"]) +
        rad * cos(heading + tur)
      
      simulated_y[bad] <- pull(model_data[t-1, "N_utm"]) +
        rad * sin(heading + tur)
    }
    
  }
  
  # create features object
  features = data.frame(E_utm = simulated_x,
                        N_utm = simulated_y,
                        spot_num = model_data[t, "spot_num"],
                        posix_since_start = model_data[t, "posix_since_start"],
                        posix_time = model_data[t, "posix_time"],
                        compID = model_data[t, "compID"], 
                        compID_num = model_data[t, "compID_num"], 
                        compIDwatch = model_data[t, "compIDwatch"], 
                        partID_num = model_data[t, "partID_num"],
                        chosen = 0)
  
  # chosen spot in first row
  features<-rbind(chosen_spot[,c("E_utm", "N_utm", "spot_num","posix_since_start", "posix_time", "compID", "compID_num", "compIDwatch", "partID_num", "chosen")],features)
  
  # add to feature matrix
  Feature_matrix[[t]] <- features
  
  # progress
  cat(paste(t,"/", T, sep = ""),"\r")
}

time_1 <- Sys.time()

time_1-time_0 #55 seconds


# save results
saveRDS(Feature_matrix, file = here("data","1_step_selection","Feature_matrix.rds"))

# create identifier list
identifiers = data.frame(compID = vector(),
                         compID_num = vector(),
                         compIDwatch = vector(),
                         spot_num = vector(),
                         partID_num = vector(),
                         t = vector(),
                         posix_since_start = vector(),
                         posix_time = vector())

for (t in 1:length(Feature_matrix)){
  compID = Feature_matrix[[t]]$compID[1]
  compID_num = Feature_matrix[[t]]$compID_num[1]
  compIDwatch = Feature_matrix[[t]]$compIDwatch[1]
  spot_num = Feature_matrix[[t]]$spot_num[1]
  partID_num = Feature_matrix[[t]]$partID_num[1]
  posix_since_start = Feature_matrix[[t]]$posix_since_start[1]
  posix_time = Feature_matrix[[t]]$posix_time[1]
  identifiers[t,] <- c(compID,compID_num,compIDwatch, spot_num, partID_num, t,posix_since_start,posix_time)
}


saveRDS(identifiers, file = here("data","1_step_selection","identifiers.rds"))

rm(list = ls())

################################################################################
# END
################################################################################