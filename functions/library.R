## Library

## If a package is installed, it will be loaded. If any 
## are not, the missing package(s) will be installed 
## from CRAN and then loaded.

## First specify the packages of interest
packages = c("readr",     #read data
             "tidyverse", #data processing
             "lubridate", #date transformations
             "eeptools",  #date transformations
             "ggplot2",   #visualizations
             "reshape2",  #data transformations
             "geosphere",
             "Rcpp",
             "devtools",
             "loo",
             "abind",
             "terra",
             "momentuHMM",
             "data.table",
             "igraph",
             "tidybayes",
             "concaveman",
             "future",
             "furrr",
             "purrr",
             "ggforce",
             #"ggalt",
             "forcats",
             "readxl",
             "rstan",
             "brms",
             "cmdstanr",
             "RColorBrewer",
             "raster",
             "sf",
             "sp",
             "fasterize",
             "R.utils",
             "gtools",
             "fields"
)

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

############################# helper functions #################################

# compute mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# inverse logit
inv_logit <- function(x){
  exp(x)/(1+exp(x))
}

# compute x and y coordinates#Function: https://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example# +datum=WGS84
  res <- spTransform(xy, CRS(paste("+proj=utm +datum=WGS84 +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}
library(sp)
# This function takes as a .csv file that was converted using the GARMIN SDK Tool, 
# and extracts the variable of interest (distance, heart rate, gps position, ...)
# It returns a two column dataframe. The first column contains the variable of interest,
# the second column contains the respective timestamp.


extract_variable_from_csv <- function(df, # output file of .fit to .csv converter
                                      variable #variable of interest
){
  
  # get all rows with data on variable
  positions = which(df == variable, arr.ind=TRUE)
  
  # get all rows that contain timestamps
  positions_of_timestamps = which(df == "timestamp", arr.ind=TRUE)
  
  # now select entries to the right of time and variable entries and keep row index
  positions[,2] = positions[,2] + 1
  positions_of_timestamps[,2] = positions_of_timestamps[,2] + 1
  
  variable_data = data.frame(row_index = vector(), variable_data = vector())
  for (i in 1:length(positions[,1])){
    variable_data[i,] = c(positions[i,1],as.numeric(df[positions[i,1], positions[i,2]]))
  }
  
  time_data = data.frame(row_index = vector(), time_data = vector())
  for (i in 1:length(positions_of_timestamps[,1])){
    time_data[i,] = c(positions_of_timestamps[i,1],as.numeric(df[positions_of_timestamps[i,1], positions_of_timestamps[i,2]]))
  }
  
  matched_data = full_join(variable_data, time_data, by = "row_index")
  
  matched_data[order(matched_data$row_index),]
  
  colnames(matched_data) <- c("row_index",paste(variable), "timestamp")
  return(matched_data)
}

# heading angles
heading_angle <- function(x, y){
  v <- c(y[1]-x[1],y[2]-x[2])
  angle <- atan2(v[2],v[1])
  
  while(angle<=(-pi))
    angle <- angle + 2*pi
  while(angle>pi)
    angle <- angle - 2*pi
  
  return(angle)
}


rm(list = c("package.check", "packages"))

# https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/


# key linking lake id 2022 and 2023
key = cbind(y2023 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
            y2022 = c(7, 10, 6, 9, 8, 1, 2, 3, 5, 4))


