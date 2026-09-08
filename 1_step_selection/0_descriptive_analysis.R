################################################################################
#
# Title: 0. Step selection descriptive
#
# Author F Dhellemmes
#
# last update: 05/08/2026 (DD/MM/YYYY)
#
################################################################################

library(here)
library(data.table)
library(ggplot2)
library(terra)
library(brms)

# source helper functions
sourceCpp(here("functions","fast_dist.cpp"))
sourceCpp(here("functions","ptinpoly.cpp"))

# set the seed
set.seed(42)

################################################################################
################################################################################
################################################################################


spot_data<-fread(here("data","2511_angling_events.csv"))

metadata<-fread(here::here("data","metadata_anon.csv"))

# Add participant ID and Competition ID to spot_data
spot_data[,compID_num:=metadata[match(spot_data$compID,metadata$compID),"compID_num"]]
spot_data[,partID_num:=metadata[match(spot_data$compIDwatch,metadata$compIDwatch),"partID"]]
spot_data[,minute:=round(posix_since_start/60)]


# Import Lake Data
polygons<-readRDS(here("data","arena_polygons.rds"))

polygons <- lapply(polygons, terra::unwrap)
spot_data[i,"compID"]


##creating the dataset
i=1
spot_data<-spot_data[,c(1:6,23)]
spot_data$data_type<-"real"
pts<-spatSample(polygons[[spot_data[i,compID]]], size = 5)
coords <- crds(pts)
control_data<-rbind(spot_data[i,],spot_data[i,],spot_data[i,],spot_data[i,],spot_data[i,])
control_data$data_type<-"control"
control_data$N_utm<-coords[,2]
control_data$E_utm<-coords[,1]

# filling the dataset
for (i in 2:nrow(spot_data)){
  pts<-spatSample(polygons[[spot_data[i,compID]]], size = 5)
  coords <- crds(pts)
  control_add<-rbind(spot_data[i,],spot_data[i,],spot_data[i,],spot_data[i,],spot_data[i,])
  control_add$data_type<-"control"
  control_add$N_utm<-coords[,2]
  control_add$E_utm<-coords[,1]
  
  control_data<-rbind(control_data, control_add)
}

#merge datasets
spot_data<-rbind(spot_data, control_data)

## Load fish presence data
file_list <- readRDS(here::here("data","2512_rasters_file_paths.rds"))

fish_presence <- list()
base_path <- here::here("data")

for (comp in names(file_list)) {
  full_paths <- file.path(base_path, file_list[[comp]])
  
  # Load rasters from these full paths
  fish_presence[[comp]] <- lapply(full_paths, terra::rast)
}


##add fish presence
spot_data[, fish_3 := {
  cat(sprintf("Processing compID = %s, minute = %s (%d pts)\r", 
              compID[1], minute[1], .N))
  
  r <- raster::raster(fish_presence[[as.character(compID[1])]][[minute[1]]])
  
  vals <- raster::extract(r, cbind(E_utm, N_utm), buffer = 3, fun = mean)
  
  rm(r)
  vals
}, by = .(compID, minute)]


### Model the fish presence at real and control spots

m1<-brm(fish_3~data_type+(1|compID),
        data=spot_data,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores=4,
        seed=42)

summary(m1)


# #png(here::here("figures","supp","S19.png"), res=300, units="in", width=4.5, height=4.5)
# #svg(here::here("figures","supp","S19.png"), width=4.5, height=4.5)
# par(mar=c(4,4,3,1))
# vioplot::vioplot(fish_3~data_type, data=spot_data, col=scales::alpha(c("#0097A7", "magenta3"), c(0.3,0.3) ),
#                  border=c("#0097A7", "magenta3"), rectCol=c("#0097A7", "magenta3"), lineCol=c("#0097A7", "magenta3"),
#                  las=2 , xaxt="n", xlab="", ylab="Fish presence", ylim=c(0,1.2))
# 
# text(1,1.18,"Estimate = 0.58 \n[0.5, 0.65]", cex=0.8)
# text(2,1.18,"Estimate = 0.18 \n[0.58, 0.76]", cex=0.8)
# axis(1, at=c(1,2), c("Control spots","Chosen spots"))
# dev.off()


################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################
