################################################################################
#
# Title: 5. Make the figure
#
# Author F Dhellemmes
#
# last update: 06/08/2026 (DD/MM/YYYY)
#
################################################################################

library(here)
library(data.table)
library(ggplot2)
library(brms)

source(here("functions","library.R"))

# set the seed

set.seed(42)


################################################################################
################################################################################
################################################################################

real<-readRDS(here::here("data","5_angler_effects","RealSpots.rds"))


model_success<-readRDS(here::here("data","5_angler_effects", "model_fit", "model_success.RDS"))
model_loss<-readRDS(here::here("data","5_angler_effects", "model_fit","model_loss.RDS"))



## make model predictions
newdata_successloss<- expand_grid(data_type = c("real", "control"),
                              fish_t0 = seq(0,1,by=0.1),
                              counter = seq(1,31,by=5),
                              dist_sonar=15,
                              compID = unique(real$compID),
                              anglingID =NA)



# loss_pred <- predicted_draws(model_loss,newdata = newdata_successloss,re_formula = ~(1|compID), # Keeps 'block', sets 'actor' to 0
#                              allow_new_levels = TRUE)
# saveRDS(loss_pred, here::here("data","5_angler_effects", "model_fit","loss_pred.RDS"))
# success_pred <- predicted_draws(model_success,newdata = newdata_successloss,re_formula = ~(1|compID), # Keeps 'block', sets 'actor' to 0
#                                 allow_new_levels = TRUE)
# saveRDS(success_pred, here::here("data","5_angler_effects", "model_fit","success_pred.RDS"))


loss_pred<-readRDS(here::here("data","5_angler_effects", "model_fit","loss_pred.RDS"))
success_pred<-readRDS(here::here("data","5_angler_effects", "model_fit","success_pred.RDS"))

## Compute differences between real and control
loss_pred <- as.data.table(loss_pred)
success_pred <- as.data.table(success_pred)

# Key columns identifying a unique prediction: covariates + which posterior draw
key_cols <- c("fish_t0", "counter", "dist_sonar", "compID", "anglingID", ".draw")

loss_pred_real    <- loss_pred[data_type == "real"]
loss_pred_control <- loss_pred[data_type == "control"]

success_pred_real    <- success_pred[data_type == "real"]
success_pred_control <- success_pred[data_type == "control"]

# Merge control predictions onto real rows by matching key_cols exactly
loss_pred_real <- merge(
  loss_pred_real,
  loss_pred_control[, c(..key_cols, ".prediction"), with = FALSE],
  by = key_cols,
  suffixes = c("", "_control"),
  all.x = TRUE
)
setnames(loss_pred_real, ".prediction_control", "control_pred")

success_pred_real <- merge(
  success_pred_real,
  success_pred_control[, c(..key_cols, ".prediction"), with = FALSE],
  by = key_cols,
  suffixes = c("", "_control"),
  all.x = TRUE
)
setnames(success_pred_real, ".prediction_control", "control_pred")

# Sanity check: every real row should have found a match
stopifnot(
  loss_pred_real[, sum(is.na(control_pred))] == 0,
  success_pred_real[, sum(is.na(control_pred))] == 0
)

loss_pred_real[, predictive_diff := .prediction - control_pred]
success_pred_real[, predictive_diff := .prediction - control_pred]

loss_avg_pred <- loss_pred_real[
  , .(avg_pred = mean(predictive_diff)),
  by = .(fish_t0, counter)
]

success_avg_pred <- success_pred_real[
  , .(avg_pred = mean(predictive_diff)),
  by = .(fish_t0, counter)
]

## same for catch

model_catch<-readRDS(here::here("data","5_angler_effects", "model_fit", "model_catch.RDS"))

newdata_catch<- tidyr::expand_grid(fish_3 = seq(0,1,by=0.1),
                                   counter = seq(1,101,by=5),
                                   dist_sonar=15,
                                   time_at_spot= 5,
                                   compID = unique(real$compID),
                                   anglingID=NA)



catch_pred <- predicted_draws(model_catch,newdata = newdata_catch, allow_new_levels = TRUE)


#saveRDS(catch_pred, here::here("data","5_angler_effects", "model_fit","success_pred.RDS"))


p <- posterior_epred(
  model_catch,
  newdata = newdata_catch,
  re_formula     = ~ (1 | compID),  # keep anglingID RE, drop compID
  allow_new_levels = FALSE
)

newdata_catch$post_prob<-colMeans(p)

newdata_catch<-as.data.table(newdata_catch)
catch_avg_pred <- newdata_catch[
  , .(avg_post_prob = mean(post_prob)),
  by = .(fish_3, counter)
]

# Build newdata at density_t0 = 0.5 across observed counter range
counter_seq <- seq(min(success_avg_pred$counter), max(success_avg_pred$counter), length.out = 30)

newdat <- data.frame(
  counter     = counter_seq,
  fish_t0  = 0.5,
  dist_sonar  = mean(model_success$data$dist_sonar),   # held at mean; adjust if you average differently
  data_type   = "real"                                  # or whichever level "avg_pred" was built from
)

# Posterior draws (rows = draws, cols = newdata rows)
epred_succ <- posterior_epred(model_success, newdata = newdat, re_formula = NA)

# Summarize: median + 95% CI, then shift by baseline 0.5
med_succ    <- 0.5 + apply(epred_succ, 2, median)
lower_succ  <- 0.5 + apply(epred_succ, 2, quantile, probs = 0.025)
upper_succ  <- 0.5 + apply(epred_succ, 2, quantile, probs = 0.975)



# Posterior draws (rows = draws, cols = newdata rows)
epred_loss <- posterior_epred(model_loss, newdata = newdat, re_formula = NA)

# Summarize: median + 95% CI, then shift by baseline 0.5
med_loss    <- 0.5 + apply(epred_loss, 2, median)
lower_loss  <- 0.5 + apply(epred_loss, 2, quantile, probs = 0.025)
upper_loss  <- 0.5 + apply(epred_loss, 2, quantile, probs = 0.975)



newdat_catch <- data.frame(
  fish_3    = 0.5,
  counter      = counter_seq,
  dist_sonar   = 15,
  time_at_spot = 5
)

epred_catch <- posterior_epred(model_catch, newdata = newdat_catch, re_formula = NA)

med_catch   <- apply(epred_catch, 2, median)
lower_catch <- apply(epred_catch, 2, quantile, probs = 0.025)
upper_catch <- apply(epred_catch, 2, quantile, probs = 0.975)



# svg(here::here("figures","Fig_4_predictions.svg"), width=7, height=3.5, pointsize=12, bg="transparent")
# #png(here::here("figures","Fig_4_predictions.png"), width=7, units="in",res=300, height=3.5)
# 
# par(mfrow=c(1,2),mar = c(4, 5, 2, 1))
# plot(counter_seq, med_succ, type = "n",
#      ylim = range(0.4, 0.6), las=1,
#      xlab = "", ylab = "Predicted fish presence")
# mtext("Time (minutes)",1,2)
# polygon(c(counter_seq, rev(counter_seq)),
#         c(upper_succ, rev(lower_succ)),
#         col = scales::alpha("magenta3", 0.2), border = NA)
# 
# lines(counter_seq, med_succ, col = "magenta3", lwd = 2)
# 
# polygon(c(counter_seq, rev(counter_seq)),
#         c(upper_loss, rev(lower_loss)),
#         col = scales::alpha("#0097A7", 0.2), border = NA)
# 
# lines(counter_seq, med_loss, col = "#0097A7", lwd = 2)
# abline(h = 0.5, lty = 2, lwd=2, col = "darkgrey")
# 
# legend(0,.61, c("Successful", "Unsuccessful"), col=c("magenta3","#0097A7"), bty = "n", lwd=2, lty=1)
# text(30,0.595,"A", font=2)
# 
# plot(counter_seq, med_catch, type = "n",
#      ylim = range(0,0.4),
#      xlab = "", ylab = "Predicted catch probability", las=1)
# mtext("Time (minutes)",1,2)
# polygon(c(counter_seq, rev(counter_seq)),
#         c(upper_catch, rev(lower_catch)),
#         col = scales::alpha("magenta3", 0.2), border = NA)
# 
# lines(counter_seq, med_catch, col = "magenta3", lwd = 2)
# text(30,0.39,"B", font=2)
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
