################################################################################
#
# Title: 6. Supplementary figure
#
# Author F Dhellemmes
#
# last update: 06/08/2026 (DD/MM/YYYY)
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

source(here("functions","library.R"))

# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################


model_success<-readRDS(here::here("data","5_angler_effects", "model_fit", "model_success.RDS"))
model_loss<-readRDS(here::here("data","5_angler_effects", "model_fit","model_loss.RDS"))
model_catch<-readRDS(here::here("data","5_angler_effects", "model_fit", "model_catch.RDS"))

# variables(model_success)
# variables(model_catch)

success_draws<-tidybayes::spread_draws(model_success, c(b_Intercept,                        
                                                        b_counter,                         
                                                        b_fish_t0,
                                                        b_data_typereal,                  
                                                        b_dist_sonar,
                                                        `b_counter:fish_t0`,             
                                                        `b_counter:data_typereal`,
                                                        `b_fish_t0:data_typereal`,        
                                                        `b_counter:fish_t0:data_typereal`
))


loss_draws<-tidybayes::spread_draws(model_loss, c(b_Intercept,                        
                                                  b_counter,                         
                                                  b_fish_t0,
                                                  b_data_typereal,                  
                                                  b_dist_sonar,
                                                  `b_counter:fish_t0`,             
                                                  `b_counter:data_typereal`,
                                                  `b_fish_t0:data_typereal`,        
                                                  `b_counter:fish_t0:data_typereal`
))


catch_draws<-tidybayes::spread_draws(model_catch, c(b_Intercept,                      
                                                    b_counter,                        
                                                    b_fish_3,                       
                                                    b_dist_sonar,                     
                                                    b_time_at_spot       
))
params_catch <- as_draws_df(model_catch)





col<-"#0097A7"

var_succ<-success_draws[,c(5:12)]

densities_succ <-lapply(var_succ, density)
cis_succ <- lapply(var_succ, quantile, probs = c(0.025, 0.975))
cis2_succ <- lapply(var_succ, quantile, probs = c(0.17, 0.83))

medians_succ <- lapply(var_succ, median)

y_pos_succ <- c(1:length(var_succ))   # vertical positions
colors_succ <- rep(col,length(var_succ))#
x_lab<-"Posterior estimates"
cex<-1

xlim_succ <- range(sapply(densities_succ, function(d) range(d$x)))






col<-"#0097A7"

var_loss<-loss_draws[,c(5:12)]

densities_loss <-lapply(var_loss, density)
cis_loss <- lapply(var_loss, quantile, probs = c(0.025, 0.975))
cis2_loss <- lapply(var_loss, quantile, probs = c(0.17, 0.83))

medians_loss <- lapply(var_loss, median)

y_pos_loss <- c(1:length(var_loss))   # vertical positions
colors_loss <- rep(col,length(var_loss))#
x_lab<-"Posterior estimates"
cex<-1

xlim_loss <- range(sapply(densities_loss, function(d) range(d$x)))




var_catch<-catch_draws[,c(5:7)]

densities_catch <-lapply(var_catch, density)
cis_catch <- lapply(var_catch, quantile, probs = c(0.025, 0.975))
cis2_catch <- lapply(var_catch, quantile, probs = c(0.17, 0.83))

medians_catch <- lapply(var_catch, median)

y_pos_catch <- c(1:length(var_catch))   # vertical positions
colors_catch <- rep(col,length(var_catch))#
x_lab<-"Posterior estimates"
cex<-1

xlim_catch <- range(sapply(densities_catch, function(d) range(d$x)))










# 
# 
# #svg(here::here("figures", "supp","Fig4_estimates.svg"),  width=3.5, height=8)
# png(here::here("figures", "supp","Fig4_estimates.png"), res=300, units="in", width=3.5, height=9)
# 
# layout(matrix(c(1,1,2,2,3), ncol=1))
# 
# par(mar=c(5,8,3,1))
# plot(NA,
#      xlim = xlim_loss,
#      ylim = c(min(y_pos_loss)-0.7, max(y_pos_loss)+0.5),
#      yaxt = "n",
#      xlab = x_lab,
#      ylab = "",
#      cex=cex)
# 
# for (i in 1:length(var_loss)) {
#   d <- densities_loss[[i]]
#   ci <- cis_loss[[i]]
#   ci2 <- cis2_loss[[i]]
#   
#   # Scale density height
#   y <- d$y / max(d$y) * 0.35 + y_pos_loss[i]
#   
#   # Density line
#   lines(d$x, y, col = colors_loss[i], lwd = 1)
#   lines(c(min(d$x), max(d$x)), y[1:2], col = colors_loss[i], lwd = 1)
#   lines(c(ci[1], ci[2]), y[1:2], col = colors_loss[i], lwd = 5, lend=1)
#   #lines(c(ci2[1], ci2[2]), y[1:2], col = colors[i], lwd = 6, lend=1)
#   
#   points(medians_loss[[i]], y_pos_loss[i], pch=19, col=colors_loss[i], cex=1.6)
#   # CI shading
#   idx <- d$x >= ci[1] & d$x <= ci[2]
#   idx2 <- d$x >= ci2[1] & d$x <= ci2[2]
#   
#   polygon(
#     c(d$x, rev(d$x)),
#     c(y, rep(y_pos_loss[i], length(d$x))),
#     col = adjustcolor(colors_loss[i], alpha.f = 0.2),
#     border = NA
#   )
# }
# text(.0115,8.5,"A",font=2)
# mtext("Fish dynamics at \nunsuccessful spots", 3, line=.6, cex=.7, font=2)
# 
# axis(2,
#      at = y_pos_loss,
#      labels = c("Time (min)", "Starting\nfish presence", "Spot type =\n real","Distance \n to sonar","Time:Starting fish",
#                 "Time:\nSpot type=real", "Starting fish:\n Spot type=real","Time:Starting fish:\nSpot type=real"),
#      las = 1, cex=cex)
# 
# 
# abline(v=0, las=2, lty=2, lwd=2, col="darkgrey")
# 
# 
# plot(NA,
#      xlim = xlim_succ,
#      ylim = c(min(y_pos_succ)-0.7, max(y_pos_succ)+0.5),
#      yaxt = "n",
#      xlab = x_lab,
#      ylab = "",
#      cex=cex)
# 
# for (i in 1:length(var_succ)) {
#   d <- densities_succ[[i]]
#   ci <- cis_succ[[i]]
#   ci2 <- cis2_succ[[i]]
#   
#   # Scale density height
#   y <- d$y / max(d$y) * 0.35 + y_pos_succ[i]
#   
#   # Density line
#   lines(d$x, y, col = colors_succ[i], lwd = 1)
#   lines(c(min(d$x), max(d$x)), y[1:2], col = colors_succ[i], lwd = 1)
#   lines(c(ci[1], ci[2]), y[1:2], col = colors_succ[i], lwd = 5, lend=1)
#   #lines(c(ci2[1], ci2[2]), y[1:2], col = colors[i], lwd = 6, lend=1)
#   
#   points(medians_succ[[i]], y_pos_succ[i], pch=19, col=colors_succ[i], cex=1.6)
#   # CI shading
#   idx <- d$x >= ci[1] & d$x <= ci[2]
#   idx2 <- d$x >= ci2[1] & d$x <= ci2[2]
#   
#   polygon(
#     c(d$x, rev(d$x)),
#     c(y, rep(y_pos_succ[i], length(d$x))),
#     col = adjustcolor(colors_succ[i], alpha.f = 0.2),
#     border = NA
#   )
# }
# 
# text(.058,8.5,"B",font=2)
# mtext("Fish dynamics at \nsuccessful spots", 3, line=.6, cex=.7, font=2)
# 
# axis(2,
#      at = y_pos_succ,
#      labels = c("Time (min)", "Starting\nfish presence", "Spot type =\n real","Distance \n to sonar","Time:Starting fish",
#                 "Time:\nSpot type=real", "Starting fish:\n Spot type=real","Time:Starting fish:\nSpot type=real"),
#      las = 1, cex=cex)
# 
# 
# abline(v=0, las=2, lty=2, lwd=2, col="darkgrey")
# 
# 
# 
# 
# plot(NA,
#      xlim = xlim_catch,
#      ylim = c(min(y_pos_catch)-0.7, max(y_pos_catch)+0.5),
#      yaxt = "n",
#      xlab = x_lab,
#      ylab = "",
#      cex=cex)
# 
# for (i in 1:length(var_catch)) {
#   d <- densities_catch[[i]]
#   ci <- cis_catch[[i]]
#   ci2 <- cis2_catch[[i]]
#   
#   # Scale density height
#   y <- d$y / max(d$y) * 0.35 + y_pos_catch[i]
#   
#   # Density line
#   lines(d$x, y, col = colors_catch[i], lwd = 1)
#   lines(c(min(d$x), max(d$x)), y[1:2], col = colors_catch[i], lwd = 1)
#   lines(c(ci[1], ci[2]), y[1:2], col = colors_catch[i], lwd = 5, lend=1)
#   #lines(c(ci2[1], ci2[2]), y[1:2], col = colors[i], lwd = 6, lend=1)
#   
#   points(medians_catch[[i]], y_pos_catch[i], pch=19, col=colors_catch[i], cex=1.6)
#   # CI shading
#   idx <- d$x >= ci[1] & d$x <= ci[2]
#   idx2 <- d$x >= ci2[1] & d$x <= ci2[2]
#   
#   polygon(
#     c(d$x, rev(d$x)),
#     c(y, rep(y_pos_catch[i], length(d$x))),
#     col = adjustcolor(colors_catch[i], alpha.f = 0.2),
#     border = NA
#   )
# }
# 
# text(1.45,3,"C",font=2)
# mtext("Catch probability model", 3, line=.6, cex=.7, font=2)
# 
# 
# axis(2,
#      at = y_pos_catch,
#      labels = c("Time (min)", "Fish presence", "Distance \n to sonar"),
#      las = 1, cex=cex)
# 
# 
# abline(v=0, las=2, lty=2, lwd=2, col="darkgrey")
# 
# dev.off()


################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################
