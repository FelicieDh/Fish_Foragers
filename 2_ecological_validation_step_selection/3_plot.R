################################################################################
#
# Title: 3. Make plot
#
# Author F Dhellemmes
#
# last update: 03/08/2026 (DD/MM/YYYY)
#
################################################################################


library(here)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(data.table)
library(ggplot2)
library(terra)
library(kableExtra)
library(brms)
library(lme4)
library(dplyr)

source(here::here("functions","prep_model_results.R"))
source(here::here("functions","plot_model_results.R"))



# set the seed

set.seed(42)

################################################################################
################################################################################
################################################################################

#### load data
model_data <- fread(here("data","2511_angling_events.csv"))


#### Load step selection
step_selection_fit<-readRDS(here::here("data","1_step_selection","model_fit","step_selection_fit.RDS"))
params <- as_draws_df(step_selection_fit)
params$`beta[4]`<-params$`beta[4]`*(0-1) ## we invert sonar so that positive values mean *attraction to sonar*
beta_draws <- subset_draws(params,variable = "beta")  
bandwidth_draws <- subset_draws(params,variable = "bandwidth_log")  

#### Load validation
validation_fit<-readRDS(here::here("data","2_ecological_validation_step_selection","model_fit","ev_ss_fit.RDS"))
data_features<-readRDS(here("data","2_ecological_validation_step_selection","Spotchoice_features.rds"))
data_features <- do.call(rbind, data_features)


validation_fit[["fit"]]@sim[["samples"]][[1]][["b_dist_sonar_std"]]<-validation_fit[["fit"]]@sim[["samples"]][[1]][["b_dist_sonar_std"]]*(0-1) ## we invert sonar so that large values mean *close to sonar*
validation_fit[["fit"]]@sim[["samples"]][[2]][["b_dist_sonar_std"]]<-validation_fit[["fit"]]@sim[["samples"]][[2]][["b_dist_sonar_std"]]*(0-1) ## we invert sonar so that large values mean *close to sonar*
validation_fit[["fit"]]@sim[["samples"]][[3]][["b_dist_sonar_std"]]<-validation_fit[["fit"]]@sim[["samples"]][[3]][["b_dist_sonar_std"]]*(0-1) ## we invert sonar so that large values mean *close to sonar*
validation_fit[["fit"]]@sim[["samples"]][[4]][["b_dist_sonar_std"]]<-validation_fit[["fit"]]@sim[["samples"]][[4]][["b_dist_sonar_std"]]*(0-1) ## we invert sonar so that large values mean *close to sonar*


#### Standardizing the step selection model results 
data_features[,time_since_start := posix_since_start/(2*60*60)]
data_features[,time_since_start := time_since_start/max(time_since_start), by=(compID)]

data_features[,locality:= model_data$step]


trip_sd<-data_features[,.(compID = compID_num,
                     success = sd(success_density, na.rm=T),
                     loss = sd(loss_density, na.rm=T),
                     social = sd(social_density, na.rm=T),
                     sonar = sd(dist_sonar),
                     time = sd(time_since_start),
                     locality = sd(locality, na.rm=T)),by=(compIDwatch)]

trip_sd[,timeloc := time*locality]

comp_sd<-trip_sd[,.(success = mean(success),
                    loss = mean(loss),
                    social = mean(social),
                    sonar = mean(sonar),
                    locality = mean(locality),
                    timeloc = mean(timeloc)),by=(compID)]


# Build a single global SD per feature, matching beta[1:6] order: social, success, loss, sonar, locality, timeloc
feature_sd_ss <- c(
  social   = sd(data_features$social_density, na.rm = TRUE),
  success  = sd(data_features$success_density, na.rm = TRUE),
  loss     = sd(data_features$loss_density, na.rm = TRUE),
  sonar    = sd(data_features$dist_sonar),
  locality = sd(data_features$locality, na.rm = TRUE),
  timeloc  = sd(data_features$time_since_start * data_features$locality, na.rm = TRUE)
)

beta_sd_ss <- setNames(feature_sd_ss, as.character(1:6))
names(beta_sd_ss)   # should be "1" "2" "3" "4" "5" "6"

beta_effects<-prep_model_results(step_selection_fit, beta, feature, compID, 1:6)

bandwidth_effects<-prep_model_results(step_selection_fit, bandwidth_log, feature, compID, 7:11)

beta_std<-as.data.table(beta_effects)

long_sd<-as.data.table(melt(comp_sd, id.vars = "compID"))

feature_nr<-data.table(number=1:6, feature=c("social","success","loss",
                                             "sonar","locality","timeloc"))
for (i in 1:nrow(long_sd)){
  beta_std[feature == feature_nr[feature==long_sd[i,variable], number] &
             group == long_sd[i,compID],  `:=`(
               mean_effect_offset = mean_effect_offset*long_sd[i,value],
               lower_effect_offset = lower_effect_offset*long_sd[i,value],
               upper_effect_offset = upper_effect_offset*long_sd[i,value])]
}

for (f in unique(beta_std$feature)) {
  if (as.character(f) %in% names(beta_sd_ss)) {
    beta_std[beta_std$feature == f, "mean_effect"] <-
      beta_std[beta_std$feature == f, "mean_effect"][1] * beta_sd_ss[[as.character(f)]]
  }
}
beta_std<-as.data.frame(beta_std)

#beta_std[which(beta_std$feature==4),"mean_effect"]<-beta_std[which(beta_std$feature==4),"mean_effect"]*(0-1)




#### preparing the bandwidths plots 


#space
x      <- seq(-50, 50, length.out = 400)
dist <- seq(0, 35, length.out=60)
centres <- c(-20, 10, 30)

bandwidths_spat <- as.numeric(c(exp(round(summarise_draws(bandwidth_draws)[2:5],2)[1,1]), exp(round(summarise_draws(bandwidth_draws)[2:5],2)[2,1]), exp(round(summarise_draws(bandwidth_draws)[2:5],2)[4,1])))

bandwidths_spat_ci<-as.numeric(c(exp(round(summarise_draws(bandwidth_draws)[6:7],2)[1,1]), exp(round(summarise_draws(bandwidth_draws)[6:7],2)[1,2]),
                                 exp(round(summarise_draws(bandwidth_draws)[6:7],2)[2,1]), exp(round(summarise_draws(bandwidth_draws)[6:7],2)[2,2]),
                                 exp(round(summarise_draws(bandwidth_draws)[6:7],2)[4,1]),exp(round(summarise_draws(bandwidth_draws)[6:7],2)[4,2])))

labels     <- c("Social", "Loss", "Success")
colors_spat     <- c("darkorange", "magenta3", "#0097A7")  # orange, blue, purple

kernel_data_spat <- expand.grid(d = dist, bandwidth = bandwidths_spat) %>%
  mutate(kernel = exp(-d^2 / (2 * bandwidth^2))) %>%
  mutate(bandwidth = factor(bandwidth, levels = bandwidths_spat,
                            labels = c("Social", "Loss", "Success")))

kernel_data_spat_ci <- expand.grid(d = dist, bandwidth = bandwidths_spat_ci) %>%
  mutate(kernel = exp(-d^2 / (2 * bandwidth^2))) %>%
  mutate(bandwidth = factor(bandwidth, levels = bandwidths_spat_ci,
                            labels = c("Social_l", "Social_u",
                                       "Loss_l", "Loss_u",
                                       "Success_l", "Success_u")))

####
#time


bandwidths <- as.numeric(c(exp(round(summarise_draws(bandwidth_draws)[2:5],2)[3,1]),exp(round(summarise_draws(bandwidth_draws)[2:5],2)[5,1])))  # small, medium, large

kernel_data <- expand.grid(d = dist, bandwidth = bandwidths) %>%
  mutate(kernel = exp(-d^2 / (2 * bandwidth^2))) %>%
  mutate(bandwidth = factor(bandwidth, labels = c("Loss", "Success")))

bandwidths_ci <- as.numeric(c(exp(round(summarise_draws(bandwidth_draws)[6:7],2)[3,1]),exp(round(summarise_draws(bandwidth_draws)[6:7],2)[3,2]),
                              exp(round(summarise_draws(bandwidth_draws)[6:7],2)[5,1]),exp(round(summarise_draws(bandwidth_draws)[6:7],2)[5,2])))  # small, medium, large

bandwidths_ci^2

kernel_data_ci <- expand.grid(d = dist, bandwidth = bandwidths_ci) %>%
  mutate(kernel = exp(-dist^2 / (2 * bandwidth^2))) %>%
  mutate(bandwidth = factor(bandwidth, labels = c("Loss_l","Loss_u", "Success_l", "Success_u")))



#### prepping the validation model results 
fitval<-conditional_effects(validation_fit)
fitval_draws<-tidybayes::spread_draws(validation_fit, c(b_intercept,b_dist_sonar_std,b_time_std, b_social_density_std, b_success_density_std, b_loss_density_std))

fitval_draws<-cbind(fitval_draws[,4:8], fitval_draws[,1:3])

var<-fitval_draws[,c(2,5,4,3)]

densities <-lapply(var, density)
cis <- lapply(var, quantile, probs = c(0.025, 0.975))
cis2 <- lapply(var, quantile, probs = c(0.17, 0.83))

medians <- lapply(var, median)

col<-"#0097A7"
y_pos <- c(1:length(var))   # vertical positions
colors <- rep(col,length(var))#
x_lab<-"Standardized posterior estimate \n of the effect on fish presence"
cex<-1

xlim <- range(sapply(densities, function(d) range(d$x)))















# ####
# #png(here::here("figures","Fig2_step_selection.png"), res=300, units="in", width=7.5, height=5.5)
# svg(here::here("figures","Fig2_step_selection.svg"), width=7.5, height=5.5, bg="transparent", pointsize = 12)
# 
# layout(matrix(c(1, 3, 2, 4),
#               nrow = 2,
#               ncol = 2), widths=c(3,3), height=c(3,2))
# 
# par(mar=c(5,5,2,1), oma=c(0,0,0,0))
# 
# plotfunc(params, beta_std, 3:1, "Standardized posterior estimate\nof the spot selection model","#0097A7", c("Loss","Success","Social"), xlim=c(-2,2.5), cex=0.6, beta_sd = beta_sd_ss)#"Sonar",
# text(2.4,3.4,"A", font=2)
# mtext("Spatial features", 3,  line=0.6, cex=.7, font=2)
# 
# par(mar=c(5,5,2,1))
# plot(NA,
#      xlim = xlim,
#      ylim = c(min(y_pos)-0.7, max(y_pos)+0.5),
#      yaxt = "n",
#      xlab = x_lab,
#      ylab = "",
#      cex=cex)
# mtext("Ecological validity", 3,  line=0.6, cex=.7, font=2)
# for (i in 1:length(var)) {
#   d <- densities[[i]]
#   ci <- cis[[i]]
#   ci2 <- cis2[[i]]
#   
#   # Scale density height
#   y <- d$y / max(d$y) * 0.35 + y_pos[i]
#   
#   # Density line
#   lines(d$x, y, col = colors[i], lwd = 1)
#   lines(c(min(d$x), max(d$x)), y[1:2], col = colors[i], lwd = 1)
#   lines(c(ci[1], ci[2]), y[1:2], col = colors[i], lwd = 5, lend=1)
#   #lines(c(ci2[1], ci2[2]), y[1:2], col = colors[i], lwd = 6, lend=1)
#   
#   points(medians[[i]], y_pos[i], pch=19, col=colors[i], cex=1.6)
#   # CI shading
#   idx <- d$x >= ci[1] & d$x <= ci[2]
#   idx2 <- d$x >= ci2[1] & d$x <= ci2[2]
#   
#   polygon(
#     c(d$x, rev(d$x)),
#     c(y, rep(y_pos[i], length(d$x))),
#     col = adjustcolor(colors[i], alpha.f = 0.2),
#     border = NA
#   )
# }
# text(0.41,4.4,"B", font=2)
# axis(2,
#      at = y_pos,
#      labels = c("Time", "Loss", "Success", "Social"),#"Sonar",
#      las = 1, cex=cex)
# 
# abline(v=0, las=2, lty=2, col="darkgrey")
# 
# # --- Single overlapping plot ---
# par(mar=c(5,5,2,1))
# plot(kernel~d, data=kernel_data_spat[which(kernel_data_spat$bandwidth=="Social"),], type="l",
#      col="darkorange", lwd=3, las=1, ylab="Spatial kernel",  xlab="", ylim=c(0,1))
# mtext("Meters", 1, 2.1, cex=.8)
# mtext("Spatial kernel", 3,  line=0.6, cex=.7, font=2)
# 
# polygon(c(seq(0, 35, length.out=60), rev(seq(0, 35, length.out=60))),
#         c(kernel_data_spat_ci[which(kernel_data_spat_ci$bandwidth=="Social_u"), "kernel"],
#           rev(kernel_data_spat_ci[which(kernel_data_spat_ci$bandwidth=="Social_l"), "kernel"])),
#         col = scales::alpha("darkorange",  0.3), border = NA)
# points(kernel~d, data=kernel_data_spat[which(kernel_data_spat$bandwidth=="Loss"),], type="l",
#        col="#0097A7", lwd=3)
# polygon(c(seq(0, 35, length.out=60), rev(seq(0, 35, length.out=60))),
#         c(kernel_data_spat_ci[which(kernel_data_spat_ci$bandwidth=="Loss_u"), "kernel"],
#           rev(kernel_data_spat_ci[which(kernel_data_spat_ci$bandwidth=="Loss_l"), "kernel"])),
#         col = scales::alpha("#0097A7",  0.3), border = NA)
# points(kernel~d, data=kernel_data_spat[which(kernel_data_spat$bandwidth=="Success"),], type="l",
#        col="magenta3", lwd=3)
# polygon(c(seq(0, 35, length.out=60), rev(seq(0, 35, length.out=60))),
#         c(kernel_data_spat_ci[which(kernel_data_spat_ci$bandwidth=="Success_u"), "kernel"],
#           rev(kernel_data_spat_ci[which(kernel_data_spat_ci$bandwidth=="Success_l"), "kernel"])),
#         col = scales::alpha("magenta3",  0.3), border = NA)
# 
# legend(20,1.1, c("Social","Success", "Loss"), col=c("darkorange","magenta3","#0097A7"), bty = "n", lwd=2, lty=1)
# 
# text(35,.95,"C", font=2)
# 
# par(mar=c(5,5,2,1))
# plot(kernel~d, data=kernel_data[which(kernel_data$bandwidth=="Loss"),], type="l",
#      col="#0097A7", lwd=3, las=1, ylab="Time kernel",  xlab="")
# mtext("Spot number", 1, 2.1, cex=.8)
# mtext("Temporal kernel", 3,  line=0.6, cex=.7, font=2)
# 
# polygon(c(seq(0, 35, length.out=60), rev(seq(0, 35, length.out=60))),
#         c(kernel_data_ci[which(kernel_data_ci$bandwidth=="Loss_u"), "kernel"],
#           rev(kernel_data_ci[which(kernel_data_ci$bandwidth=="Loss_l"), "kernel"])),
#         col = scales::alpha("#0097A7",  0.3), border = NA)
# points(kernel~d, data=kernel_data[which(kernel_data$bandwidth=="Success"),], type="l",
#        col="magenta3", lwd=3)
# polygon(c(seq(0, 35, length.out=60), rev(seq(0, 35, length.out=60))),
#         c(kernel_data_ci[which(kernel_data_ci$bandwidth=="Success_u"), "kernel"],
#           rev(kernel_data_ci[which(kernel_data_ci$bandwidth=="Success_l"), "kernel"])),
#         col = scales::alpha("magenta3",  0.3), border = NA)
# 
# #legend(18,1.1, c("Success", "Loss"), col=c("#0097A7","magenta3"), bty = "n", lwd=2, lty=1)
# text(35,.95,"D", font=2)
# 
# dev.off()
# # layout(1)


################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################
