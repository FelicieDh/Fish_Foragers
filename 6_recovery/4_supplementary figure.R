################################################################################
#
# Title: 4. Supplementary figure
#
# Author F Dhellemmes
#
# last update: 08/08/2026 (DD/MM/YYYY)
#
################################################################################


library(data.table)
library(terra)
library(here)
library(brms)
library(tidybayes)
library(ggplot2)

# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################


model_recovery<-readRDS(here::here("data","6_recovery","model_fit","catch_recovery_model_fit.RDS"))


grid_data15<-readRDS(here::here("data","6_recovery","grid_data_15.RDS"))
grid15<-readRDS(here::here("data", "6_recovery", "voronoi_grid_15.rds"))
grid15 <- lapply(grid15, terra::unwrap)

grid_data15[,cell_name:=paste0(compID,"_",cell_id)]
grid_data15_events<-grid_data15[!is.na(event_name),]


summary_grid<-grid_data15_events[,.(compID=max(compID),
                                    minute=min(minute),
                                    cell_id=max(cell_id),
                                    angler_presence=max(angler_presence),
                                    angler_nr=max(angler_nr),
                                    catch_bi=max(catch_bi),
                                    catch_cum=sum(catch_abs),
                                    event_nr=min(event_nr),
                                    time_away=min(time_away),
                                    successful_event=min(successful_event),
                                    cumulative_time=min(cumulative_time),
                                    fish=mean(mean_fish),
                                    event_time=max(event_time),
                                    cell_name=max(cell_name)),
                                 by=event_name]

grid_data15_successful<-grid_data15_events[
  , if (any(catch_bi > 0)) .SD, 
  by = cell_name
] ##only succesful cells



recovery_draws<-tidybayes::spread_draws(model_recovery, c(b_Intercept,b_event_time,b_event_nr, b_time_away, b_angler_nr, b_mean_fish))

#fitval_draws<-cbind(fitval_draws[,4:8], fitval_draws[,1:3])

recovery_random_draws<-tidybayes::spread_draws(model_recovery, r_cell_name[condition,term])

sd_table<-grid_data15_successful[,.(sd_event_time = sd(event_time),
                                    sd_event_nr = sd(event_nr),
                                    sd_time_away = sd(time_away),
                                    sd_angler_nr = sd(angler_nr),
                                    sd_mean_fish = sd(mean_fish))]


recovery_draws$b_event_time_std<-recovery_draws$b_event_time*sd_table$sd_event_time
recovery_draws$b_event_nr_std<-recovery_draws$b_event_nr*sd_table$sd_event_nr
recovery_draws$b_time_away_std<-recovery_draws$b_time_away*sd_table$sd_time_away
recovery_draws$b_angler_nr_std<-recovery_draws$b_angler_nr*sd_table$sd_angler_nr
recovery_draws$b_mean_fish_std<-recovery_draws$b_mean_fish*sd_table$sd_mean_fish


col<-"#0097A7"

var<-recovery_draws[,c(10:14)]

densities <-lapply(var, density)
cis <- lapply(var, quantile, probs = c(0.025, 0.975))
cis2 <- lapply(var, quantile, probs = c(0.17, 0.83))

medians <- lapply(var, median)

y_pos <- c(1:length(var))   # vertical positions
colors <- rep(col,length(var))#
x_lab<-"Standardized estimate \n of the effect on catch probability"
cex<-1

xlim <- range(sapply(densities, function(d) range(d$x)))

# png(here::here("figures","supp", "Fig_revisit_catch.png"), units="in", res=300, width=4.5, height=4.5)
# par(mar=c(5,8,1,1))
# plot(NA,
#      xlim = xlim,
#      ylim = c(min(y_pos)-0.7, max(y_pos)+0.5),
#      yaxt = "n",
#      xlab = x_lab,
#      ylab = "",
#      cex=cex)
# 
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
# 
# axis(2,
#      at = y_pos,
#      labels = c("Event time", "Event #","Time away","Number \n of anglers","Fish presence " ),
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
