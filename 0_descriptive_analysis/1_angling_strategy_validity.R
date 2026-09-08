################################################################################
#
# Title: 1. Angling strategy validity
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




## Open arena and sonar data
arena<-fread(here::here("data","2408_ArenaExtents.csv"))
sonars<-rbind(arena[,c(6,7)],
              arena[,c(8,9)],
              arena[,c(10,11)],use.names=FALSE)


##calculate dist_sonar
dmat <- fastPdist2(as.matrix(spot_data[, .(E_utm, N_utm)]), as.matrix(sonars[, 1:2]))
spot_data[, dist_sonar := apply(dmat, 1, min)]



################################################################################
######## Model 1: Test if fish presence is higher at successful spots
################################################################################

spot_data[,success:=as.factor(catch_at_spot)]

m1<-brm(fish_3~success+dist_sonar+(1|compID),
        data=spot_data,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores=4,
        seed=42)

summary(m1)

#### Figure
# png(here::here("figures","supp","S16.png"), res=300, units="in", width=4.5, height=4.5)
# svg(here::here("figures","supp","S16.svg"), width=4.5, height=4.5)
# par(mar=c(4,4,3,1))
# vioplot::vioplot(fish_3~success, data=spot_data, col=scales::alpha(c( "magenta3","#0097A7"), c(0.3,0.3) ), border=c( "magenta3","#0097A7"), rectCol=c( "magenta3","#0097A7"), lineCol=c( "magenta3","#0097A7"), las=2 , xaxt="n", xlab="", ylab="Fish presence")
# axis(1, at=c(1,2), c("No","Yes"))
# mtext("Spot successful ?", 1, 2)
# mtext(paste("Successful spot estimate = 0.12 [0.09, 0.15] \n with intercept taken as unsuccessful spot"), 3, line=1, cex=.75)
# dev.off()






################################################################################
#### Model 2: fish under angler vs in the whole arena per minute
################################################################################

fish_min<-readRDS(here::here("data", "0_descriptive_analysis", "fish_minute.RDS"))

mod_df_1<-fish_min[,c(1,2,3)]
mod_df_2<-fish_min[,c(1,2,4)]

colnames(mod_df_1)[3]<-"fish"
colnames(mod_df_2)[3]<-"fish"

mod_df_1$group<-"arena"
mod_df_2$group<-"anglers"

mod_df<-rbind(mod_df_1, mod_df_2)

m2<-brm(fish~group+minute +(minute|compID), 
        data=mod_df,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores=4,
        seed=42)

summary(m2)
print(summary(m2), digits = 5)

#saveRDS(m1,here::here("0_descriptive_analysis","fish_min_model.RDS"))
model_fish <- readRDS(here::here("0_descriptive_analysis","fish_min_model.RDS"))

fitval<-conditional_effects(model_fish)
fitval_2<-conditional_effects(model_fish, conditions = make_conditions(model_fish, "group"))


##Figure model output
# #png(here::here("figures","supp","S20.png"), res=300, units="in", width=6, height=4.5)
# #svg(here::here("figures","supp","S20.svg"), width=6, height=4.5)
# plot(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="arena")],fitval_2[[2]]$estimate__[which(fitval_2[[2]]$group=="arena")], type = "l", 
#      lwd=3, ylab="", las= 1, ylim=c(min(fitval_2[[2]]$lower__)-0.1,
#                            max(fitval_2[[2]]$upper__)+0.1),col=scales::alpha("#0097A7",1),  xlab="")
# polygon(c(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="arena")], rev(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="arena")])),
#         c(fitval_2[[2]]$upper__[which(fitval_2[[2]]$group=="arena")], rev(fitval_2[[2]]$lower__[which(fitval_2[[2]]$group=="arena")])),
#         col = scales::alpha("#0097A7",0.2), lty = 0)
# lines(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="arena")],fitval_2[[2]]$estimate__[which(fitval_2[[2]]$group=="arena")], type = "l", lwd=3, col="#0097A7")
# lines(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="arena")],fitval_2[[2]]$lower__[which(fitval_2[[2]]$group=="arena")], type = "l", lwd=2, lty=2,col="#0097A7")
# lines(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="arena")],fitval_2[[2]]$upper__[which(fitval_2[[2]]$group=="arena")], type = "l", lwd=2, lty=2,col="#0097A7")
# 
# points(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="anglers")],fitval_2[[2]]$estimate__[which(fitval_2[[2]]$group=="anglers")], type = "l", 
#      lwd=3, ylab="", las= 1, col=scales::alpha("magenta3",1),  xlab="")
# polygon(c(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="anglers")], rev(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="anglers")])),
#         c(fitval_2[[2]]$upper__[which(fitval_2[[2]]$group=="anglers")], rev(fitval_2[[2]]$lower__[which(fitval_2[[2]]$group=="anglers")])),
#         col = scales::alpha("magenta3",0.2), lty = 0)
# lines(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="anglers")],fitval_2[[2]]$estimate__[which(fitval_2[[2]]$group=="anglers")], type = "l", lwd=3, col="magenta3")
# lines(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="anglers")],fitval_2[[2]]$lower__[which(fitval_2[[2]]$group=="anglers")], type = "l", lwd=2, lty=2,col="magenta3")
# lines(fitval_2[[2]]$minute[which(fitval_2[[2]]$group=="anglers")],fitval_2[[2]]$upper__[which(fitval_2[[2]]$group=="anglers")], type = "l", lwd=2, lty=2,col="magenta3")
# mtext(side=2, line=2.5,"Fish presence")
# mtext(side=1, line=2.5,"Minute")
# legend("topright",
#        legend = c("Estimated fish presence", "Estimated fish presence under anglers"),
#        col = c("#0097A7", "magenta3"),
#        lwd = 2,
#        #fill = c(scales::alpha("#0097A7", 0.3), scales::alpha("deeppink",  0.3)),
#        border = NA,
#        bg=NA, bty="n")
# 
# dev.off()


# Aggregate mean and SE per minute
agg_mean <- aggregate(mean_density ~ minute, data = fish_min, FUN = mean)
agg_mean$sd <- aggregate(mean_density ~ minute, data = fish_min, 
                         FUN = function(x) sd(x))$mean_density

agg_angler <- aggregate(angler_density ~ minute, data = fish_min, FUN = mean)
agg_angler$sd <- aggregate(angler_density ~ minute, data = fish_min, 
                           FUN = function(x) sd(x))$angler_density





# #### Figure 1 F
# #png(here::here("figures","Fig1_f.png"), res=300, units="in", width=6, height=4.5)
# #svg(here::here("figures","Fig1_f.svg"), width=6, height=4.5)
# par(mfrow = c(1,1))
# 
# ylim_range <- range(c(agg_mean$mean_density - agg_mean$sd,
#                       agg_mean$mean_density + agg_mean$sd,
#                       agg_angler$angler_density - agg_angler$sd,
#                       agg_angler$angler_density + agg_angler$sd), na.rm = TRUE)
# 
# plot(agg_mean$minute, agg_mean$mean_density,
#      type = "n",
#      ylim = ylim_range,
#      xlab = "Minute", ylab = "Presence probability", las=1)
# 
# # Green ribbon for mean_density
# polygon(c(agg_mean$minute, rev(agg_mean$minute)),
#         c(agg_mean$mean_density + agg_mean$sd,
#           rev(agg_mean$mean_density - agg_mean$sd)),
#         col = scales::alpha("#0097A7",  0.3), border = NA)
# 
# # Pink ribbon for angler_density
# polygon(c(agg_angler$minute, rev(agg_angler$minute)),
#         c(agg_angler$angler_density + agg_angler$sd,
#           rev(agg_angler$angler_density - agg_angler$sd)),
#         col = scales::alpha("magenta3",  0.3), border = NA)
# 
# # Lines
# lines(agg_mean$minute, agg_mean$mean_density, col = "#0097A7", lwd = 2)
# lines(agg_angler$minute, agg_angler$angler_density, col = "magenta3", lwd = 2)
# 
# # Legend
# legend("bottomright",
#        legend = c("Average fish presence", "Fish presence under anglers"),
#        col = c("#0097A7", "magenta3"),
#        lwd = 2,
#        #fill = c(scales::alpha("#0097A7", 0.3), scales::alpha("deeppink",  0.3)),
#        border = NA,
#        bg=NA, bty="n")
# 
# dev.off()



# ## Figure per lake
# 
# fish_min$lake<-substr(fish_min$compID, 1,3)
# 
# # Aggregate mean and SE per minute and lake
# agg_mean <- aggregate(mean_density ~ minute + lake, data = fish_min, FUN = mean)
# agg_mean$sd <- aggregate(mean_density ~ minute + lake, data = fish_min, 
#                          FUN = function(x) sd(x))$mean_density
# 
# agg_angler <- aggregate(angler_density ~ minute + lake, data = fish_min, FUN = mean)
# agg_angler$sd <- aggregate(angler_density ~ minute + lake, data = fish_min, 
#                            FUN = function(x) sd(x))$angler_density
# 
#
# #png(here::here("figures","supp","S13.png"), res=300, units="in", width=8, height=4)
# #svg(here::here("figures","supp","S13.svg"), width=8, height=4)
# 
# par(mfrow=c(1,3), mar=c(4.5,4,2,1))
# 
# ylim_range <- range(c(agg_mean$mean_density - agg_mean$sd,
#                       agg_mean$mean_density + agg_mean$sd,
#                       agg_angler$angler_density - agg_angler$sd,
#                       agg_angler$angler_density + agg_angler$sd), na.rm = TRUE)
# 
# plot(agg_mean$minute, agg_mean$mean_density,
#      type = "n",
#      ylim = ylim_range,
#      xlab = "Minute", ylab = "Fish presence",
#      main = "Kuorinka", las=1)
# 
# # Green ribbon for mean_density
# polygon(c(agg_mean$minute[which(agg_mean$lake=="kuo")], rev(agg_mean$minute[which(agg_mean$lake=="kuo")])),
#         c(agg_mean$mean_density[which(agg_mean$lake=="kuo")] + agg_mean$sd[which(agg_mean$lake=="kuo")],
#           rev(agg_mean$mean_density[which(agg_mean$lake=="kuo")] - agg_mean$sd[which(agg_mean$lake=="kuo")])),
#         col = scales::alpha("#0097A7",  0.3), border = NA)
# 
# # Pink ribbon for angler_density
# polygon(c(agg_angler$minute[which(agg_angler$lake=="kuo")], rev(agg_angler$minute[which(agg_angler$lake=="kuo")])),
#         c(agg_angler$angler_density[which(agg_angler$lake=="kuo")] + agg_angler$sd[which(agg_angler$lake=="kuo")],
#           rev(agg_angler$angler_density[which(agg_angler$lake=="kuo")] - agg_angler$sd[which(agg_angler$lake=="kuo")])),
#         col = scales::alpha("magenta3",  0.3), border = NA)
# 
# # Lines
# lines(agg_mean$minute[which(agg_mean$lake=="kuo")], agg_mean$mean_density[which(agg_mean$lake=="kuo")], col = "#0097A7", lwd = 2)
# lines(agg_angler$minute[which(agg_angler$lake=="kuo")], agg_angler$angler_density[which(agg_angler$lake=="kuo")], col = "magenta3", lwd = 2)
# 
# 
# 
# 
# plot(agg_mean$minute, agg_mean$mean_density,
#      type = "n",
#      ylim = ylim_range,
#      xlab = "Minute", ylab = "Fish presence",
#      main = "Sompalampi", las=1)
# 
# # Green ribbon for mean_density
# polygon(c(agg_mean$minute[which(agg_mean$lake=="som")], rev(agg_mean$minute[which(agg_mean$lake=="som")])),
#         c(agg_mean$mean_density[which(agg_mean$lake=="som")] + agg_mean$sd[which(agg_mean$lake=="som")],
#           rev(agg_mean$mean_density[which(agg_mean$lake=="som")] - agg_mean$sd[which(agg_mean$lake=="som")])),
#         col = scales::alpha("#0097A7",  0.3), border = NA)
# 
# # Pink ribbon for angler_density
# polygon(c(agg_angler$minute[which(agg_angler$lake=="som")], rev(agg_angler$minute[which(agg_angler$lake=="som")])),
#         c(agg_angler$angler_density[which(agg_angler$lake=="som")] + agg_angler$sd[which(agg_angler$lake=="som")],
#           rev(agg_angler$angler_density[which(agg_angler$lake=="som")] - agg_angler$sd[which(agg_angler$lake=="som")])),
#         col = scales::alpha("magenta3",  0.3), border = NA)
# 
# # Lines
# lines(agg_mean$minute[which(agg_mean$lake=="som")], agg_mean$mean_density[which(agg_mean$lake=="som")], col = "#0097A7", lwd = 2)
# lines(agg_angler$minute[which(agg_angler$lake=="som")], agg_angler$angler_density[which(agg_angler$lake=="som")], col = "magenta3", lwd = 2)
# 
# 
# 
# 
# 
# 
# 
# plot(agg_mean$minute, agg_mean$mean_density,
#      type = "n",
#      ylim = ylim_range,
#      xlab = "Minute", ylab = "Fish presence",
#      main = "Pitkänen", las=1)
# 
# # Green ribbon for mean_density
# polygon(c(agg_mean$minute[which(agg_mean$lake=="pit")], rev(agg_mean$minute[which(agg_mean$lake=="pit")])),
#         c(agg_mean$mean_density[which(agg_mean$lake=="pit")] + agg_mean$sd[which(agg_mean$lake=="pit")],
#           rev(agg_mean$mean_density[which(agg_mean$lake=="pit")] - agg_mean$sd[which(agg_mean$lake=="pit")])),
#         col = scales::alpha("#0097A7",  0.3), border = NA)
# 
# # Pink ribbon for angler_density
# polygon(c(agg_angler$minute[which(agg_angler$lake=="pit")], rev(agg_angler$minute[which(agg_angler$lake=="pit")])),
#         c(agg_angler$angler_density[which(agg_angler$lake=="pit")] + agg_angler$sd[which(agg_angler$lake=="pit")],
#           rev(agg_angler$angler_density[which(agg_angler$lake=="pit")] - agg_angler$sd[which(agg_angler$lake=="pit")])),
#         col = scales::alpha("magenta3",  0.3), border = NA)
# 
# # Lines
# lines(agg_mean$minute[which(agg_mean$lake=="pit")], agg_mean$mean_density[which(agg_mean$lake=="pit")], col = "#0097A7", lwd = 2)
# lines(agg_angler$minute[which(agg_angler$lake=="pit")], agg_angler$angler_density[which(agg_angler$lake=="pit")], col = "magenta3", lwd = 2)
# 
# 
# # Legend
# legend("bottomright",
#        legend = c("Average fish \n presence", "Fish under anglers"),
#        col = c("#0097A7", "magenta3"),
#        lwd = 2,
#        bty = "n",
#        #fill = c(scales::alpha("cyan4", 0.3), scales::alpha("magenta3",  0.3)),
#        border = NA)
# dev.off()


################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################

