################################################################################
#
# Title: 3. Create control dataset and run models
#
# Author F Dhellemmes
#
# last update: 06/08/2026 (DD/MM/YYYY)
#
################################################################################


library(here)
library(data.table)
library(ggplot2)
library(terra)
library(kableExtra)
library(brms)
library(tidybayes)

source(here("functions","library.R"))

sourceCpp(here("functions","fast_dist.cpp"))
sourceCpp(here("functions","ptinpoly.cpp"))


# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################


standard<-readRDS(here::here("data","5_angler_effects","StandardSpots.rds"))
real<-readRDS(here::here("data","5_angler_effects","RealSpots.rds"))

real[,counter:=seq(1,.N,by=1), by=anglingID]

fulldata<-fread(here("data", "2604_gps_data.csv"))

for (i in 1:length(unique(real$anglingID))){
  real[anglingID==unique(real$anglingID)[i], posix_start := min(fulldata[anglingID==unique(real$anglingID)[i],posix_time])]
  real[anglingID==unique(real$anglingID)[i], posix_stop := max(fulldata[anglingID==unique(real$anglingID)[i],posix_time])]
}

real[,time_at_spot_s:=posix_stop-posix_start]

real[,data_type:="real"]

stdspot <- setNames(
  lapply(standard, function(inner) rbindlist(inner, idcol = "source_id", fill = TRUE)),
  names(standard)
) ##make the standard spots into a list

stdspot <- lapply(stdspot, function(dt) {
  dt[, part := as.character(part)]
  return(dt)
}) ##make the standard spots into a list

#colnames(real)

real$spot_id<-0



alldata<-real[1,c(1:3,26,25,9:10,12,20,21,17,18,15)] #create a receiving dataset

threshold<-0.02

i=1


for (i in 1:nrow(real)){
  
  alldata<-rbind(alldata, real[i,c(1:3,26,25,9:10,12,20,21,17,18,15)])
  
  if (real[i,counter]>1) next
  
  spotIDs<-stdspot[[real[i,compID]]][source_id==real[i,min] & fish_3>real[i,fish_3]-threshold & fish_3<real[i,fish_3]+threshold & comp_in_10m==0,geom]
  n <- min(5, length(spotIDs))
  spotIDs <- sample(spotIDs, n)
  
  add<-stdspot[[real[i,compID]]][source_id %in% real[anglingID==real[i,anglingID],min] & geom %in% spotIDs,]
  add$data_type<-"control"
  add$anglingID<-real[i,anglingID]
  add$compID<-real[i,compID]
  add$counter<-add$source_id-min(real[anglingID==real[i,anglingID],min])+1
  add$success<-real[i,success]
  
  add<-add[,c(13,1,14,2,12,4,5,16,9,15,10,7,8)] #we match add to alldata
  colnames(add)<-colnames(alldata)
  
  alldata<-rbind(alldata, add) #and bind
  
}

alldata<-alldata[-1,]

alldata[,uniqueID:=paste0(anglingID,"_", spot_id)] #create a unique id

alldata[  , delta_fish := fish_3 - fish_3[counter == 1][1], by = uniqueID] #calculate deltafish

alldata[, fish_t0 := first(fish_3), by = uniqueID] #inform fish t0


summary_dt <- alldata[, .N, by = .(data_type, min, anglingID)]
summary(summary_dt[which(summary_dt$data_type=="control"),"N"]) ##how many controls do we have per real data



model_loss<-brm(delta_fish~counter*fish_t0*data_type+dist_sonar+(1|compID)+(1|anglingID),
                data = alldata[success==0,],
                family = gaussian(),
                iter = 3000,
                warmup = 1000,
                chains = 4,
                cores = 4,
                seed = 42,
                control = list(max_treedepth = 15))

saveRDS(model_loss, here::here("data","5_angler_effects", "model_fit","model_loss.RDS"))

model_success<-brm(delta_fish~counter*fish_t0*data_type+dist_sonar+(1|compID)+(1|anglingID),
                   data = alldata[success==1,],
                   family = gaussian(),
                   iter = 3000,
                   warmup = 1000,
                   chains = 4,
                   cores = 4,
                   seed = 42,
                   control = list(max_treedepth = 15))

saveRDS(model_success, here::here("data","5_angler_effects", "model_fit","model_success.RDS")) #(X10 makes the delta density X 10)


# model_loss<-readRDS(here::here("data","5_angler_effects", "model_fit","model_loss.RDS"))
# model_success<-readRDS(here::here("data","5_angler_effects", "model_fit","model_success.RDS"))

# # loss model
# fixef(model_loss)
# summary(model_loss)
# conditional_effects(model_loss)
# 
# 
# # success model
# summary(model_success)
# conditional_effects(model_success)
# fixef(model_success)

################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################
