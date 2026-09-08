################################################################################
#
# Title: 3. Make figure
#
# Author F Dhellemmes
#
# last update: 09/08/2026 (DD/MM/YYYY)
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
library(dbscan)

source(here::here("functions","prep_model_results.R"))
source(here::here("functions","plot_model_results.R"))


# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################


data <- fread(here("data","2511_angling_events.csv"))
data <- readRDS(here::here("data","7_revisit_analysis","revisits_10.RDS"))
# Import Lake Data
polygons<-readRDS(here("data","arena_polygons.rds"))
polygons <- lapply(polygons, terra::unwrap)

my_colours <- c(
  "#FF8C00", "#0097A7", "#CD00CD",
  "#E05A00", "#004D5B", "#6A1B9A",
  "#FF6347", "#00BCD4", "#E040FB",
  "#FF4500", "#006064", "#7B00D4",
  "#FFB347", "#26C6DA", "#BA68C8",
  "#FFA040", "#80DEEA", "#CE93D8",
  "#D4500A", "#00B0C8", "#A000A0",
  "#C0392B", "#00838F", "#880E4F",
  "#E57373", "#4DB6AC", "#F48FB1",
  "#FF8A65", "#80CBC4", "#F06292",
  "#FFCC80", "#B2EBF2"
)

for (i in 1:nrow(data)){
  data[i, col:=my_colours[data[i,final_cluster_id]]]
}

head(data)




mod_future_state<-readRDS(here::here("data","7_revisit_analysis", "model_fit", "mod_future_state_10.RDS"))
mod_future_fish<-readRDS(here::here("data","7_revisit_analysis", "model_fit", "mod_future_fish_10.RDS"))
mod_revisit_state<-readRDS(here::here("data","7_revisit_analysis", "model_fit", "mod_revisit_state_10.RDS"))
mod_revisit_fish<-readRDS(here::here("data","7_revisit_analysis", "model_fit", "mod_revisit_fish_10.RDS"))


nd_state <- data.frame(
  success  = c(0, 1),
  stop_min = mean(data$stop_min),
  compID   = NA  # NA to ignore random effect (population-level)
)

# Get posterior draws of predicted probabilities
draws_state <- posterior_epred(mod_future_state, newdata = nd_state, 
                               allow_new_levels = TRUE, re_formula = NA)
# draws is a matrix: rows = posterior samples, cols = [success=0, success=1]

# Summarise
p_mean_state <- apply(draws_state, 2, mean)
p_lo_state   <- apply(draws_state, 2, quantile, 0.025)
p_hi_state   <- apply(draws_state, 2, quantile, 0.975)

# --- Plot ---
x <- c(0, 1)


nd_fish <- data.frame(
  fish_stop = seq(min(data$fish_stop), max(data$fish_stop), length.out = 100),
  stop_min  = mean(data$stop_min),
  compID    = NA
)

draws_fish <- posterior_epred(mod_future_fish, newdata = nd_fish,
                              allow_new_levels = TRUE, re_formula = NA)

p_mean_fish <- apply(draws_fish, 2, mean)
p_lo_fish   <- apply(draws_fish, 2, quantile, 0.025)
p_hi_fish   <- apply(draws_fish, 2, quantile, 0.975)





data[,step_1:= data.table::shift(step, 1L), by=compIDwatch]
# New data: prev_success = 0 and 1, continuous vars fixed at mean
nd_prev_success <- data.frame(
  prev_success = as.integer(c(0, 1)),
  start_min    = as.integer(floor(mean(data$start_min))),
  step_1       = mean(data$step_1, na.rm=T),
  spot_num     = as.integer(floor(mean(data$spot_num))),
  compID       = NA
)

draws_prev_success <- posterior_epred(mod_revisit_state, newdata = nd_prev_success,
                                      allow_new_levels = TRUE, re_formula = NA)

p_mean_prev_success <- apply(draws_prev_success, 2, mean)
p_lo_prev_success   <- apply(draws_prev_success, 2, quantile, 0.025)
p_hi_prev_success   <- apply(draws_prev_success, 2, quantile, 0.975)

x_prev_success <- c(0, 1)




nd_prev_fish <- data.frame(
  prev_fish = seq(min(data$prev_fish, na.rm=T), max(data$prev_fish, na.rm=T), length.out = 100),
  start_min    = as.integer(floor(mean(data$start_min))),
  step_1       = mean(data$step_1, na.rm=T),
  spot_num     = as.integer(floor(mean(data$spot_num))),
  compID       = NA
)

draws_prev_fish <- posterior_epred(mod_revisit_fish, newdata = nd_prev_fish,
                                   allow_new_levels = TRUE, re_formula = NA)

p_mean_prev_fish <- apply(draws_prev_fish, 2, mean)
p_lo_prev_fish   <- apply(draws_prev_fish, 2, quantile, 0.025)
p_hi_prev_fish   <- apply(draws_prev_fish, 2, quantile, 0.975)




i=46 #53, 46

cluster_centres <- data[compIDwatch==unique(data$compIDwatch)[i], 
                        .(E_utm = mean(E_utm), N_utm = mean(N_utm), col=col[1]), 
                        by = final_cluster_id]
par(mfrow=c(1,1),mar=c(1,0,1,0))
plot(polygons[[data[compIDwatch==unique(data$compIDwatch)[i], compID][1]]], las=1, axes=FALSE)
points(data[compIDwatch==unique(data$compIDwatch)[i], c("E_utm", "N_utm")], col=scales::alpha(data[compIDwatch==unique(data$compIDwatch)[i], col],.6), pch=19)
lines(terra::buffer(terra::vect(data[compIDwatch==unique(data$compIDwatch)[i], c("E_utm", "N_utm")],geom = c("E_utm", "N_utm")),5),
      col=scales::alpha(data[compIDwatch==unique(data$compIDwatch)[i], col],.6),
      lwd=2)
ext <- terra::ext(polygons[[data[compIDwatch==unique(data$compIDwatch)[i], compID][1]]])
x0 <- ext$xmin + (ext$xmax - ext$xmin) * 0.05
y0 <- ext$ymin + (ext$ymax - ext$ymin) * 0.05
lines(c(x0, x0 + 10), c(y0, y0), lwd=2)
text(x0 + 5, y0 + (ext$ymax - ext$ymin) * 0.04, "10 m", cex=1, adj=0.5)
#text(ext$xmax-5,ext$ymax-5,"A", font=2)
text(ext$xmax, ext$ymax, "A", cex=1, font=2, adj=1)
text(cluster_centres$E_utm[1:nrow(unique(data[compIDwatch==unique(data$compIDwatch)[i], c("final_cluster_id")]))]-15, 
     cluster_centres$N_utm[1:nrow(unique(data[compIDwatch==unique(data$compIDwatch)[i], c("final_cluster_id")]))], 
     cluster_centres$final_cluster_id[1:nrow(unique(data[compIDwatch==unique(data$compIDwatch)[i], c("final_cluster_id")]))], 
     col=cluster_centres$col[1:nrow(unique(data[compIDwatch==unique(data$compIDwatch)[i], c("final_cluster_id")]))])




# #png(here::here("figures","Fig5_revisits.png"), res=300, units="in", width=7, height=6)
# svg(here::here("figures","Fig5_revisits_simple.svg"), width=7, height=6, pointsize=12, bg="transparent")
# 
# 
# layout(
#   matrix(c(1, 2, 3,
#            
# 
#            4, 5, 6), ncol = 3, byrow = TRUE),
#   widths  = c(1.5,2,2),
#   heights = c(2, 2)
# )
# 
# i=46
# par(mar=c(1,0,1,0))
# plot.new()
# 
# par(mar=c(5,5,3,1))
# 
# 
# plot(x, p_mean_state,
#      ylim = c(0, 0.8), xlim = c(-0.5, 1.5),
#      pch = 19, cex = 1.5, xaxt = "n",
#      xlab = "",
#      ylab = "Predicted probability to revisit cluster",
#      #,main = "Posterior predicted probability by success"
#      col=c(scales::alpha("#0097A7",1), scales::alpha("magenta3",1)), las=1 )
# mtext("Was spot successful ?", 1, 2.1, cex=.7)
# axis(1, at = c(0, 1), labels = c("No", "Yes"))
# mtext("Probability of future revisits \n given current success", 3, .6, cex=.7, font = 2)
# 
# arrows(x, p_lo_state, x, p_hi_state,
#        angle = 90, code = 3, length = 0.05, lwd = 2,col=c(scales::alpha("#0097A7",1), scales::alpha("magenta3",1)))
# text(1.45,.8, "B", font=2)
# # Density mountains
# for (i in 1:2) {
#   d <- density(draws_state[, i])
# 
#   # Scale density width for visual clarity
#   d_scaled <- d$y / max(d$y) * 0.35
# 
#   # Draw filled polygon (horizontal density, centered on x[i])
#   polygon(c(x[i] + d_scaled, x[i] - rev(d_scaled)),
#           c(d$x, rev(d$x)),
#           col = c(scales::alpha("#0097A7",0.2), scales::alpha("magenta3",0.2))[i],
#           border = NA)
# }
# 
# # Redraw points on top
# points(x, p_mean_state, pch = 19, cex = 1.5,col=c(scales::alpha("#0097A7",1), scales::alpha("magenta3",1)))
# 
# plot(nd_fish$fish_stop, p_mean_fish,
#      type = "l", lwd = 2,
#      ylim = c(0, .8),
#      xlab = "Fish presence when leaving \n the spot", ylab = "Predicted probability to revisit cluster",
#      col=scales::alpha("#0097A7",1), las=1)
# 
# mtext("Probability of future revisits \n given fish presence", 3, .6, cex=.7, font = 2)
# 
# polygon(c(nd_fish$fish_stop, rev(nd_fish$fish_stop)),
#         c(p_lo_fish, rev(p_hi_fish)),
#         col = scales::alpha("#0097A7",0.2),
#         border = NA)
# 
# lines(nd_fish$fish_stop, p_mean_fish, lwd = 2, col = scales::alpha("#0097A7",1))
# text(1.05,.8, "C", font=2)
# 
# 
# i=8
# par(mar=c(1,0,1,0))
# plot(polygons[[data[compIDwatch==unique(data$compIDwatch)[i], compID][1]]], las=1, axes=FALSE)
# points(data[compIDwatch==unique(data$compIDwatch)[i], c("E_utm", "N_utm")], col=scales::alpha(data[compIDwatch==unique(data$compIDwatch)[i], col],.6), pch=19)
# lines(terra::buffer(terra::vect(data[compIDwatch==unique(data$compIDwatch)[i], c("E_utm", "N_utm")],geom = c("E_utm", "N_utm")),5),
#       col=scales::alpha(data[compIDwatch==unique(data$compIDwatch)[i], col],.6),
#       lwd=2)
# ext <- terra::ext(polygons[[data[compIDwatch==unique(data$compIDwatch)[i], compID][1]]])
# x0 <- ext$xmin + (ext$xmax - ext$xmin) * 0.05
# y0 <- ext$ymin + (ext$ymax - ext$ymin) * 0.05
# lines(c(x0, x0 + 10), c(y0, y0), lwd=2)
# text(x0 + 5, y0 + (ext$ymax - ext$ymin) * 0.04, "10 m", cex=1, adj=0.5)
# text(ext$xmax-5,ext$ymax-5,"D", font=2)
# #text(ext$xmax, ext$ymax, "A", cex=1, font=2, adj=1)
# text(cluster_centres$E_utm-15, cluster_centres$N_utm, cluster_centres$final_cluster_id, col=cluster_centres$col)
# 
# 
# par(mar=c(5,5,3,1))
# plot(x_prev_success, p_mean_prev_success,
#      ylim = c(0, .25), xlim = c(-0.5, 1.5),
#      pch = 19, cex = 1.5, xaxt = "n",
#      xlab = "",
#      ylab = "Predicted probability of next spot choice \n being a revisit",
#      #main = "Posterior predicted probability by previous success",
#      col=c(scales::alpha("#0097A7",1), scales::alpha("magenta3",1)), las=1)
# mtext("Was spot successful ?", 1, 2.1, cex=.7)
# axis(1, at = c(0, 1), labels = c("No", "Yes"))
# mtext("Probability of the next spot being \n a revisit given current success", 3, .6, cex=.7, font = 2)
# 
# 
# arrows(x_prev_success, p_lo_prev_success, x_prev_success, p_hi_prev_success,
#        angle = 90, code = 3, length = 0.05, lwd = 2, col=c(scales::alpha("#0097A7",1), scales::alpha("magenta3",1)))
# 
# for (i in 1:2) {
#   d <- density(draws_prev_success[, i])
#   d_scaled <- d$y / max(d$y) * 0.35
#   polygon(c(x_prev_success[i] + d_scaled, x_prev_success[i] - rev(d_scaled)),
#           c(d$x, rev(d$x)),
#           col = c(scales::alpha("#0097A7",0.2), scales::alpha("magenta3",0.2))[i],
#           border = NA)
# }
# 
# points(x_prev_success, p_mean_prev_success, pch = 19, cex = 1.5, col=c(scales::alpha("#0097A7",1), scales::alpha("magenta3",1)))
# text(1.45,.25, "E", font=2)
# 
# 
# 
# 
# 
# plot(nd_prev_fish$prev_fish, p_mean_prev_fish,
#      type = "l", lwd = 2,
#      ylim = c(0, .25),
#      xlab = "Fish presence when leaving \n the spot", ylab = "Predicted probability of next spot choice \n being a revisit",
#      col=scales::alpha("darkgrey",1), las=1)
# #main = "Posterior predicted probability by fish stop")
# mtext("Probability of next spot being \n a revisit given fish presence", 3, .6, cex=.7, font = 2)
# 
# polygon(c(nd_prev_fish$prev_fish, rev(nd_prev_fish$prev_fish)),
#         c(p_lo_prev_fish, rev(p_hi_prev_fish)),
#         col = scales::alpha("darkgrey",0.2),
#         border = NA)
# 
# lines(nd_prev_fish$prev_fish, p_mean_prev_fish, lwd = 2, col = scales::alpha("darkgrey",1))
# text(1.05,.25, "F", font=2)
# 
# 
# dev.off()



################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################
