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

fit_full<-readRDS(here::here("data","3_spot_leaving","model_fit","spot_leaving_full_fit.RDS"))

params_full <- as_draws_df(fit_full)

full_effects<-prep_model_results(fit_full, beta, feature, comp, 1:8)

validation_fit <-readRDS(file = here::here("data","4_ecological_validation_spot_leaving","model_fit","model_validation_fit.RDS"))




stan_data<-readRDS(here::here("data","3_spot_leaving", "stan_spot_leaving.RDS"))


model_data<-as.data.table(data.frame(catch = stan_data[["catch"]],
                                     fish_3 = stan_data[["fish3"]],
                                     time_since_event = stan_data[["time_since_event"]],
                                     compID = stan_data[["compID"]],
                                     partID = stan_data[["partID"]],
                                     dist_sonar = stan_data[["dist_sonar"]],
                                     social_density = stan_data[["social_density"]],
                                     success_density = stan_data[["success_density"]],
                                     loss_density = stan_data[["loss_density"]]))



trip_sd<-model_data[,.(compID = compID,
                     success = sd(success_density),
                     loss = sd(loss_density),
                     social = sd(social_density),
                     fish = sd(fish_3),
                     sonar = sd(dist_sonar)#,
                     #time_since_event = sd(time_since_event)
                     ),by=(partID)]


comp_sd<-trip_sd[,.(success = mean(success),
                    loss = mean(loss),
                    social = mean(social),
                    sonar = mean(sonar),
                    fish = mean(fish)#,
                    #time_since_event = mean(time_since_event)
                    ),by=(compID)]

feature_sd <- c(
  fish3            = sd(stan_data[["fish3"]]),
  #time_since_event = sd(stan_data[["time_since_event"]]),
  dist_sonar       = sd(stan_data[["dist_sonar"]]),
  social_density   = sd(stan_data[["social_density"]]),
  success_density  = sd(stan_data[["success_density"]]),
  loss_density     = sd(stan_data[["loss_density"]])
)

# Map feature index -> SD, matching your beta[3..8] slope order
beta_sd <- setNames(
  c(feature_sd["fish3"], #feature_sd["time_since_event"], 
    feature_sd["dist_sonar"],
    feature_sd["social_density"], feature_sd["success_density"], feature_sd["loss_density"]),
  as.character(c(3,5:8))
)

full_effects_std<-as.data.table(full_effects)

long_sd<-as.data.table(melt(comp_sd, id.vars = "compID"))

feature_nr<-data.table(number=c(3,5:8), feature=c("fish","sonar","social","success","loss"))

i=1


for (i in 1:nrow(long_sd)){
  full_effects_std[feature == feature_nr[feature==long_sd[i,variable], number] &
                     group == long_sd[i,compID],  `:=`(
                       mean_effect_offset = mean_effect_offset*long_sd[i,value],
                       lower_effect_offset = lower_effect_offset*long_sd[i,value],
                       upper_effect_offset = upper_effect_offset*long_sd[i,value])]
}

for (f in unique(full_effects_std$feature)) {
  if (as.character(f) %in% names(beta_sd)) {
    full_effects_std[full_effects_std$feature == f, "mean_effect"] <-
      full_effects_std[full_effects_std$feature == f, "mean_effect"][1] * beta_sd[[as.character(f)]]
  }
}


full_effects_std[feature ==4, mean_effect := mean_effect*6]
full_effects_std[feature ==4, mean_effect_offset := mean_effect_offset*6]
full_effects_std[feature ==4, lower_effect_offset := lower_effect_offset*6]
full_effects_std[feature ==4, upper_effect_offset := upper_effect_offset*6] #put time in minutes

params_full_std<-params_full
params_full_std[,5]<-params_full_std[,5]*6


full_effects_std<-as.data.frame(full_effects_std)





# #png(here::here("figures","Fig3_spot_leaving_figure_rev.png"), res=300, units="in", width=8, height=4)
# svg(here::here("figures","Fig3_spot_leaving_figure.svg"), width = 180 / 25.4,
#     height = 100 / 25.4, bg="transparent", pointsize=12)
# 
# par(mfrow=c(1,2),mar=c(5,6,2,1))
# 
# plotfunc(params_full_std, full_effects_std, c(3,2,4,8,7,6), "Standardized posterior estimates \nof the spot leaving model","#0097A7", xlim=rev(c(-3,1.5)) ,
#          c("Fish presence","Fish catch > 0", "Time since\nevent (min)", "Loss","Success","Social" ),beta_sd = beta_sd, cex=0.9)
# text(0-2.8,6.4,"A", font=2)
# mtext("Information features", 3, .6, cex=.8, font=2)
# 
# fitval_draws<-tidybayes::spread_draws(validation_fit, c(b_Intercept,b_bi_catch,b_time_since, b_social_density_std, b_success_density_std, b_loss_density_std, b_avg_sonar_std))
# #fitval_draws<-tidybayes::spread_draws(validation_fit, c(b_Intercept,b_catch,b_time_since_event_std, b_social_density_std, b_success_density_std, b_loss_density_std, b_sonar_std))
# 
# fitval_draws<-cbind(fitval_draws[,4:10], fitval_draws[,1:3])
# 
# fitval_random_draws<-tidybayes::spread_draws(validation_fit, r_compID[condition,term])
# 
# col<-"#0097A7"
# 
# var<-fitval_draws[,c(2,3,6,5,4)]
# 
#   densities <-lapply(var, density)
#   cis <- lapply(var, quantile, probs = c(0.025, 0.975))
#   cis2 <- lapply(var, quantile, probs = c(0.17, 0.83))
# 
#   medians <- lapply(var, median)
# 
#   y_pos <- c(1:(length(var)+1))   # vertical positions
#   colors <- rep(col,length(var))#
#   x_lab<-"Standardized estimate \n of the effect on fish presence"
#   cex<-0.9
# 
#    xlim <- range(sapply(densities, function(d) range(d$x)))#rev()
# 
#    #par(mar=c(5,6,2,1))
#     plot(NA,
#        xlim = xlim,
#        ylim = c(min(y_pos)-0.7, max(y_pos)+0.5),
#        yaxt = "n",
#        xlab = x_lab,
#        ylab = "",
#        cex.lab=cex,
#        cex.axis =cex)
#     mtext("Ecological validation", 3, .6, cex=.8, font=2)
#      for (i in 1:length(var)) {
#     d <- densities[[i]]
#     ci <- cis[[i]]
#     ci2 <- cis2[[i]]
# 
#     # Scale density height
#     y <- d$y / max(d$y) * 0.35 + (y_pos[i]+1)
# 
#     # Density line
#     lines(d$x, y, col = colors[i], lwd = 1)
#     lines(c(min(d$x), max(d$x)), y[1:2], col = colors[i], lwd = 1)
#     lines(c(ci[1], ci[2]), y[1:2], col = colors[i], lwd = 5, lend=1)
#     #lines(c(ci2[1], ci2[2]), y[1:2], col = colors[i], lwd = 6, lend=1)
# 
#     points(medians[[i]], (y_pos[i]+1), pch=19, col=colors[i], cex=1)
#     # CI shading
#     idx <- d$x >= ci[1] & d$x <= ci[2]
#     idx2 <- d$x >= ci2[1] & d$x <= ci2[2]
# 
#     polygon(
#       c(d$x, rev(d$x)),
#       c(y, rep((y_pos[i]+1), length(d$x))),
#       col = adjustcolor(colors[i], alpha.f = 0.2),
#       border = NA
#     )
#      }
# 
#      axis(2,
#        at = y_pos[2:length(y_pos)],
#        labels = c("Fish catch > 0", "Time since\nevent (min)", "Loss", "Success","Social"),
#        las = 1, cex.axis=cex)
# 
#   abline(v=0, las=2, lty=2, col="darkgrey")
#   text(0.5,6.4,"B", font=2)
# 
#   dev.off()
#   
#   
  ################################################################################
  ################################################################################
  ################################################################################
  
  
  
  rm(list = ls())
  
  ################################################################################
  # END
  ################################################################################
  