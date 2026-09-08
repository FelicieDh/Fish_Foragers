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

spot_data[,time_at_spot_m:=time_at_spot_s/60]
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

head(spot_data)

m1<-brm(time_at_spot_m~fish_3+dist_sonar+(1|compID),
        data=spot_data,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores=4,
        seed=42)

summary(m1)


post <- as_draws_df(m1)

# Value of dist_sonar to condition on
d0 <- mean(spot_data$dist_sonar, na.rm = TRUE)

# Sequence of x values
x <- seq(min(spot_data$fish_3),
         max(spot_data$fish_3),
         length.out = 200)

# Posterior fitted values
y_post <- sapply(x, function(xx)
  post$b_Intercept +
    post$b_fish_3 * xx +
    post$b_dist_sonar * d0)

# Posterior summaries
fit <- apply(y_post, 2, mean)
lwr <- apply(y_post, 2, quantile, 0.025)
upr <- apply(y_post, 2, quantile, 0.975)

# #png(here::here("figures","supp","S11_2.png"), res=300, units="in", width=4.5, height=4)
# #svg(here::here("figures","supp","S11_2.png"), width=4.5, height=4.5)
# par(mar=c(4,4,3,1))
# plot(time_at_spot_m~fish_3, data=spot_data, col=scales::alpha("#0097A7",.3),
#                   pch=19,
#                  las=2 ,  xlab="Fish presence", ylab="Time spent at spot (min)", ylim=c(0,20))
# polygon(c(x, rev(x)),
#         c(lwr, rev(upr)),
#         col = adjustcolor("#0097A7", alpha.f = 0.3),
#         border = NA)
# 
# # Posterior mean regression line
# lines(x, fit, lwd = 2, col = "#0097A7")
# text(0.5,19,"Estimate = 2.9 [2.06, 3.72]", cex=0.8)
# dev.off()


################################################################################
################################################################################
################################################################################



rm(list = ls())

################################################################################
# END
################################################################################
