################################################################################
#
# Title: 2. Run models
#
# Author F Dhellemmes
#
# last update: 09/08/2026 (DD/MM/YYYY)
#
################################################################################

library(data.table)
library(geosphere)  
library(dbscan)     
library(ggplot2)
library(here)
library(terra)
library(brglm2)
library(terra)
library(brms)

source(here("functions","library.R"))


# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################


data<-readRDS(here::here("data","7_revisit_analysis","revisits_10.RDS"))

data$will_revisit <- as.integer(data$will_revisit)
data$is_revisit <- as.integer(data$is_revisit)

cat(sprintf("Revisits detected: %d / %d data (%.1f%%)\n",
            data[is_revisit==T,.N], nrow(data),
            100 * data[is_revisit==T,.N]/nrow(data)))


# png(here::here("figures","supp","S17.png"), res=300, units="in", width=4.5, height=4.5)
# par(mar=c(5,5,1,1))
# hist(table(data$is_revisit, data$compIDwatch)[2,], 
#      main="", xlab="Number of revisits per angler", las=1,
#      col=scales::alpha("#0097A7",.2), border="#0097A7")
# text(7,14, paste("Median =", round(median(table(data$is_revisit, data$compIDwatch)[2,], na.rm=T),1), "\n +-", round(sd(table(data$is_revisit, data$compIDwatch)[2,], na.rm=T),1), "SD"))
# dev.off()

### Does success predict how likely a spot is to be revisited later ?

mod_future_state<-brm(will_revisit~success+stop_min+(1|compID), 
                      data=data, 
                      family="Bernoulli",
                      iter = 3000,
                      warmup = 1000,
                      chains = 4,
                      cores = 4,
                      seed = 42)

summary(mod_future_state)

saveRDS(mod_future_state, here::here("data","7_revisit_analysis", "model_fit", "mod_future_state_10.RDS"))



### Does fish presence predict how likely a spot is to be revisited later?

mod_future_fish<-brm(will_revisit~fish_stop+stop_min+(1|compID), 
                     data=data, 
                     family="Bernoulli",
                     iter = 3000,
                     warmup = 1000,
                     chains = 4,
                     cores = 4,
                     seed = 42)

summary(mod_future_fish)

saveRDS(mod_future_fish, here::here("data","7_revisit_analysis", "model_fit", "mod_future_fish_10.RDS"))




mod_revisit_state<-brm(is_revisit~prev_success+start_min+spot_num+(1|compID),
                       data=data, 
                       family="Bernoulli",
                       iter = 3000,
                       warmup = 1000,
                       chains = 4,
                       cores = 4,
                       seed = 42)

summary(mod_revisit_state)
saveRDS(mod_revisit_state, here::here("data","7_revisit_analysis", "model_fit", "mod_revisit_state_10.RDS"))



mod_revisit_fish<-brm(is_revisit~prev_fish+start_min+spot_num+(1|compID),
                      data=data, 
                      family="Bernoulli",
                      iter = 3000,
                      warmup = 1000,
                      chains = 4,
                      cores = 4,
                      seed = 42)

summary(mod_revisit_fish)
saveRDS(mod_revisit_fish, here::here("data","7_revisit_analysis", "model_fit", "mod_revisit_fish_10.RDS"))



# cluster_summary <- data_it[, .(
#   n_spots   = .N,
#   n_success = sum(success)
# ), by = .(compID, compIDwatch, final_cluster)]
# 
# setnames(cluster_summary, "final_cluster", "cluster_id")
# 
# plot(n_success~n_spots, data=cluster_summary)

################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################
