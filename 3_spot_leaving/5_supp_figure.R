################################################################################
#
# Title: 5. Make consecutive figures
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

source(here::here("functions","prep_model_results.R"))
source(here::here("functions","plot_model_results.R"))

fit_nocatch <- readRDS(here::here("data","3_spot_leaving","model_fit","spot_leaving_nocatch_fit.RDS"))
params_nocatch <- as_draws_df(fit_nocatch)
rm(fit_nocatch)

fit_nopersonal <- readRDS(here::here("data","3_spot_leaving","model_fit","spot_leaving_nopersonal_fit.RDS"))
params_nopersonal <- as_draws_df(fit_nopersonal)
rm(fit_nopersonal)

fit_nosocial <- readRDS(here::here("data","3_spot_leaving","model_fit","spot_leaving_nosocial_fit.RDS"))
params_nosocial <- as_draws_df(fit_nosocial)
rm(fit_nosocial)

fit_full <- readRDS(here::here("data","3_spot_leaving","model_fit","spot_leaving_full_fit.RDS"))
params_full <- as_draws_df(fit_full)
rm(fit_full)



params_full$`beta[5]`<-params_full$`beta[5]`*(0-1) ##distance to sonar as proximity to sonar
params_nosocial$`beta[5]`<-params_nosocial$`beta[5]`*(0-1) ##distance to sonar as proximity to sonar
params_nopersonal$`beta[5]`<-params_nopersonal$`beta[5]`*(0-1) ##distance to sonar as proximity to sonar
params_nocatch$`beta[3]`<-params_nocatch$`beta[3]`*(0-1) ##distance to sonar as proximity to sonar


params_nopersonal$`beta[4]`<-params_nopersonal$`beta[4]`*6 #time in minutes
params_nosocial$`beta[4]`<-params_nosocial$`beta[4]`*6 #time in minutes
params_full$`beta[4]`<-params_full$`beta[4]`*6 #time in minutes


fit_nocatch_thin<-posterior::thin_draws(params_nocatch, thin=10) #thin everything
fit_nopersonal_thin<-posterior::thin_draws(params_nopersonal, thin=10)
fit_nosocial_thin<-posterior::thin_draws(params_nosocial, thin=10)
fit_full_thin<-posterior::thin_draws(params_full, thin=10)


nocatch_effects<-prep_model_results(fit_nocatch_thin, beta, feature, comp, 1:3)
nopersonal_effects<-prep_model_results(fit_nopersonal_thin, beta, feature, comp, 1:5)
nosocial_effects<-prep_model_results(fit_nosocial_thin, beta, feature, comp, 1:6)
full_effects<-prep_model_results(fit_full_thin, beta, feature, comp, 1:8)



##

var2<-cbind(fit_nocatch_thin[,c(3)],fit_nopersonal_thin[,c(4)],fit_nosocial_thin[,c(4)],fit_full_thin[,c(4)])


rand_effects2<-rbind(nocatch_effects[nocatch_effects$feature==2,],
                     nopersonal_effects[nopersonal_effects$feature==3,],
                     nosocial_effects[nosocial_effects$feature==3,],
                     full_effects[full_effects$feature==3,])

rand_effects2$feature<-c(rep(1,12), rep(2,12), rep(3,12),rep(4,12))

col<-"#0097A7"
densities2 <-lapply(var2, density)
cis2 <- lapply(var2, quantile, probs = c(0.025, 0.975))
cis22 <- lapply(var2, quantile, probs = c(0.17, 0.83))
medians2 <- lapply(var2, median)

all_labels2<-c("- catch",
               "- success and loss",
               "- social",
               "Full model")

#y_pos2 <- c(length(var2):1)   # vertical positions
y_pos2 <- c(1:length(var2)) 
colors2 <- c(rep("#0097A7",length(var2)))#
x_lab2<-"Posterior estimate for fish presence across \n  consecutive spot leaving models (log-odds)"

cex<-1

xlim2 <- range(sapply(densities2, function(d) range(d$x)))





# #png(here::here("figures","supp","S18.png"), res=300, units="in", width=5, height=4)
# #svg(here::here("figures","supp","S18.svg"), width=5, height=4)
# 
# par(mar=c(5,8,1,1))
# plot(NA,
#      xlim = xlim2,
#      ylim = c(min(y_pos2)-0.7, max(y_pos2)+0.7),
#      yaxt = "n",
#      xlab = x_lab2,
#      ylab = "",
#      cex=cex)
# 
# for (i in length(var2):1) {
#   d <- densities2[[i]]
#   ci <- cis2[[i]]
#   ci2 <- cis22[[i]]
#   
#   # Scale density height
#   y <- d$y / max(d$y) * 0.35 + y_pos2[i]
#   
#   # Density line
#   lines(d$x, y, col = colors2[i], lwd = 1)
#   lines(c(min(d$x), max(d$x)), y[1:2], col = colors2[i], lwd = 1)
#   lines(c(ci[1], ci[2]), y[1:2], col = colors2[i], lwd = 5, lend=1)
#   #lines(c(ci2[1], ci2[2]), y[1:2], col = colors[i], lwd = 6, lend=1)
#   
#   points(medians2[[i]], y_pos2[i], pch=19, col=colors2[i], cex=1.6)
#   # CI shading
#   idx <- d$x >= ci[1] & d$x <= ci[2]
#   idx2 <- d$x >= ci2[1] & d$x <= ci2[2]
#   
#   polygon(
#     c(d$x, rev(d$x)),
#     c(y, rep(y_pos2[i], length(d$x))),
#     col = adjustcolor(colors2[i], alpha.f = 0.2),
#     border = NA
#   )
#   
#   
#   rand_y<-y_pos2[i]-(rand_effects2[which(rand_effects2$feature == y_pos2[i]),"group"]/20)
#   
#   points((rand_effects2[which(rand_effects2$feature == y_pos2[i]),"mean_effect"]+
#             rand_effects2[which(rand_effects2$feature == y_pos2[i]),"mean_effect_offset"]), rand_y, pch=19, cex=1, col=scales::alpha("darkgrey",.5))
#   
#   for (j in 1:length(rand_effects2[which(rand_effects2$feature == y_pos2[i]),"mean_effect"])){
#     segments(
#       x0 = (rand_effects2[which(rand_effects2$feature == y_pos2[i]),"mean_effect"][j]+
#               rand_effects2[which(rand_effects2$feature == y_pos2[i]),"lower_effect_offset"][j]),
#       y0 = rand_y[j],
#       x1 = (rand_effects2[which(rand_effects2$feature == y_pos2[i]),"mean_effect"][j]+
#               rand_effects2[which(rand_effects2$feature == y_pos2[i]),"upper_effect_offset"][j]),
#       y1 = rand_y[j],
#       col=scales::alpha("darkgrey",.5),
#       lwd = 1
#     )}
#   
#   text(medians2[[i]], y_pos2[i]+0.5, paste0(round(medians2[[i]],2), "[",round(ci[1],2), ", ",round(ci[2],2), "]"), cex=.8)
# }
# 
# axis(2,
#      at = y_pos2,
#      labels = all_labels2,
#      las = 1, cex=cex)
# 
# abline(v=0, las=2, lty=2, col="darkgrey")
# 
# dev.off()

################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################

