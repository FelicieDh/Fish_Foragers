################################################################################
#
# Title: 5. supp figures
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


data2<-readRDS(here::here("data","7_revisit_analysis","sensitivity","revisits_2.RDS"))
data4<-readRDS(here::here("data","7_revisit_analysis","sensitivity","revisits_4.RDS"))
data6<-readRDS(here::here("data","7_revisit_analysis","sensitivity","revisits_6.RDS"))
data8<-readRDS(here::here("data","7_revisit_analysis","sensitivity","revisits_8.RDS"))
data10<-readRDS(here::here("data","7_revisit_analysis","revisits_10.RDS"))
data12<-readRDS(here::here("data","7_revisit_analysis","sensitivity","revisits_12.RDS"))
data14<-readRDS(here::here("data","7_revisit_analysis","sensitivity","revisits_14.RDS"))


data2$will_revisit <- as.integer(data2$will_revisit)
data2$is_revisit <- as.integer(data2$is_revisit)
data4$will_revisit <- as.integer(data4$will_revisit)
data4$is_revisit <- as.integer(data4$is_revisit)
data6$will_revisit <- as.integer(data6$will_revisit)
data6$is_revisit <- as.integer(data6$is_revisit)
data8$will_revisit <- as.integer(data8$will_revisit)
data8$is_revisit <- as.integer(data8$is_revisit)
data10$will_revisit <- as.integer(data10$will_revisit)
data10$is_revisit <- as.integer(data10$is_revisit)
data12$will_revisit <- as.integer(data12$will_revisit)
data12$is_revisit <- as.integer(data12$is_revisit)
data14$will_revisit <- as.integer(data14$will_revisit)
data14$is_revisit <- as.integer(data14$is_revisit)


tab<-data.frame(threshold=seq(2,14,by=2),
                counts=c(data2[is_revisit==T,.N],
                         data4[is_revisit==T,.N],
                         data6[is_revisit==T,.N],
                         data8[is_revisit==T,.N],
                         data10[is_revisit==T,.N],
                         data12[is_revisit==T,.N],
                         data14[is_revisit==T,.N]))


par(mar=c(5,5,1,1))
plot(counts~threshold, data=tab, pch=19, col="#0097A7", las=1)
lines(counts~threshold, data=tab, col="#0097A7")




mod_future_state_2<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_state_2.RDS"))
mod_future_state_4<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_state_4.RDS"))
mod_future_state_6<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_state_6.RDS"))
mod_future_state_8<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_state_8.RDS"))
mod_future_state_10<-readRDS(here::here("data","7_revisit_analysis", "model_fit", "mod_future_state_10.RDS"))
mod_future_state_12<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_state_12.RDS"))
mod_future_state_14<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_state_14.RDS"))

mod_future_fish_2<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_fish_2.RDS"))
mod_future_fish_4<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_fish_4.RDS"))
mod_future_fish_6<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_fish_6.RDS"))
mod_future_fish_8<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_fish_8.RDS"))
mod_future_fish_10<-readRDS(here::here("data","7_revisit_analysis", "model_fit", "mod_future_fish_10.RDS"))
mod_future_fish_12<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_fish_12.RDS"))
mod_future_fish_14<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_future_fish_14.RDS"))

mod_revisit_state_2<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_state_2.RDS"))
mod_revisit_state_4<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_state_4.RDS"))
mod_revisit_state_6<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_state_6.RDS"))
mod_revisit_state_8<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_state_8.RDS"))
mod_revisit_state_10<-readRDS(here::here("data","7_revisit_analysis", "model_fit", "mod_revisit_state_10.RDS"))
mod_revisit_state_12<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_state_12.RDS"))
mod_revisit_state_14<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_state_14.RDS"))

mod_revisit_fish_2<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_fish_2.RDS"))
mod_revisit_fish_4<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_fish_4.RDS"))
mod_revisit_fish_6<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_fish_6.RDS"))
mod_revisit_fish_8<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_fish_8.RDS"))
mod_revisit_fish_10<-readRDS(here::here("data","7_revisit_analysis", "model_fit", "mod_revisit_fish_10.RDS"))
mod_revisit_fish_12<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_fish_12.RDS"))
mod_revisit_fish_14<-readRDS(here::here("data","7_revisit_analysis", "model_fit","sensitivity", "mod_revisit_fish_14.RDS"))




sum_future_state<-data.frame(threshold=seq(2,14,by=2),
                             est = c(fixef(mod_future_state_2)[2,1],
                                     fixef(mod_future_state_4)[2,1],
                                     fixef(mod_future_state_6)[2,1],
                                     fixef(mod_future_state_8)[2,1],
                                     fixef(mod_future_state_10)[2,1],
                                     fixef(mod_future_state_12)[2,1],
                                     fixef(mod_future_state_14)[2,1]),
                             low=c(fixef(mod_future_state_2)[2,3],
                                   fixef(mod_future_state_4)[2,3],
                                   fixef(mod_future_state_6)[2,3],
                                   fixef(mod_future_state_8)[2,3],
                                   fixef(mod_future_state_10)[2,3],
                                   fixef(mod_future_state_12)[2,3],
                                   fixef(mod_future_state_14)[2,3]),
                             high=c(fixef(mod_future_state_2)[2,4],
                                    fixef(mod_future_state_4)[2,4],
                                    fixef(mod_future_state_6)[2,4],
                                    fixef(mod_future_state_8)[2,4],
                                    fixef(mod_future_state_10)[2,4],
                                    fixef(mod_future_state_12)[2,4],
                                    fixef(mod_future_state_14)[2,4]))

sum_future_fish<-data.frame(threshold=seq(2,14,by=2),
                             est = c(fixef(mod_future_fish_2)[2,1],
                                     fixef(mod_future_fish_4)[2,1],
                                     fixef(mod_future_fish_6)[2,1],
                                     fixef(mod_future_fish_8)[2,1],
                                     fixef(mod_future_fish_10)[2,1],
                                     fixef(mod_future_fish_12)[2,1],
                                     fixef(mod_future_fish_14)[2,1]),
                             low=c(fixef(mod_future_fish_2)[2,3],
                                   fixef(mod_future_fish_4)[2,3],
                                   fixef(mod_future_fish_6)[2,3],
                                   fixef(mod_future_fish_8)[2,3],
                                   fixef(mod_future_fish_10)[2,3],
                                   fixef(mod_future_fish_12)[2,3],
                                   fixef(mod_future_fish_14)[2,3]),
                             high=c(fixef(mod_future_fish_2)[2,4],
                                    fixef(mod_future_fish_4)[2,4],
                                    fixef(mod_future_fish_6)[2,4],
                                    fixef(mod_future_fish_8)[2,4],
                                    fixef(mod_future_fish_10)[2,4],
                                    fixef(mod_future_fish_12)[2,4],
                                    fixef(mod_future_fish_14)[2,4]))

sum_revisit_state<-data.frame(threshold=seq(2,14,by=2),
                             est = c(fixef(mod_revisit_state_2)[2,1],
                                     fixef(mod_revisit_state_4)[2,1],
                                     fixef(mod_revisit_state_6)[2,1],
                                     fixef(mod_revisit_state_8)[2,1],
                                     fixef(mod_revisit_state_10)[2,1],
                                     fixef(mod_revisit_state_12)[2,1],
                                     fixef(mod_revisit_state_14)[2,1]),
                             low=c(fixef(mod_revisit_state_2)[2,3],
                                   fixef(mod_revisit_state_4)[2,3],
                                   fixef(mod_revisit_state_6)[2,3],
                                   fixef(mod_revisit_state_8)[2,3],
                                   fixef(mod_revisit_state_10)[2,3],
                                   fixef(mod_revisit_state_12)[2,3],
                                   fixef(mod_revisit_state_14)[2,3]),
                             high=c(fixef(mod_revisit_state_2)[2,4],
                                    fixef(mod_revisit_state_4)[2,4],
                                    fixef(mod_revisit_state_6)[2,4],
                                    fixef(mod_revisit_state_8)[2,4],
                                    fixef(mod_revisit_state_10)[2,4],
                                    fixef(mod_revisit_state_12)[2,4],
                                    fixef(mod_revisit_state_14)[2,4]))

sum_revisit_fish<-data.frame(threshold=seq(2,14,by=2),
                            est = c(fixef(mod_revisit_fish_2)[2,1],
                                    fixef(mod_revisit_fish_4)[2,1],
                                    fixef(mod_revisit_fish_6)[2,1],
                                    fixef(mod_revisit_fish_8)[2,1],
                                    fixef(mod_revisit_fish_10)[2,1],
                                    fixef(mod_revisit_fish_12)[2,1],
                                    fixef(mod_revisit_fish_14)[2,1]),
                            low=c(fixef(mod_revisit_fish_2)[2,3],
                                  fixef(mod_revisit_fish_4)[2,3],
                                  fixef(mod_revisit_fish_6)[2,3],
                                  fixef(mod_revisit_fish_8)[2,3],
                                  fixef(mod_revisit_fish_10)[2,3],
                                  fixef(mod_revisit_fish_12)[2,3],
                                  fixef(mod_revisit_fish_14)[2,3]),
                            high=c(fixef(mod_revisit_fish_2)[2,4],
                                   fixef(mod_revisit_fish_4)[2,4],
                                   fixef(mod_revisit_fish_6)[2,4],
                                   fixef(mod_revisit_fish_8)[2,4],
                                   fixef(mod_revisit_fish_10)[2,4],
                                   fixef(mod_revisit_fish_12)[2,4],
                                   fixef(mod_revisit_fish_14)[2,4]))






# png(here::here("figures","supp","Supp_revisits.png"), res=300, units="in", width=6, height=6)
# 
# par(mfrow=c(2,2), mar=c(5,5,1,1))
# plot(threshold~est, data=sum_future_state, pch=19, col="#0097A7", ylab="Cluster epsilon (m)", xlab="Estimate",
#      las=1, xlim=c(min(sum_future_state$low-0.1), max(sum_future_state$high+0.1)), main="")
# for (i in 1:nrow(sum_future_state)){
#   segments(x0=sum_future_state[i,"low"],
#            y0 = sum_future_state[i,"threshold"],
#            x1=sum_future_state[i,"high"], col="#0097A7")
# }
# abline(v=0, col="darkgrey", lty=2)
# text(1.9,13.5,"A",font=2)
# 
# plot(threshold~est, data=sum_future_fish, pch=19, col="#0097A7", ylab="Cluster epsilon (m)", xlab="Estimate",
#      las=1, xlim=c(min(sum_future_fish$low-0.1), max(sum_future_fish$high+0.1)), main="")
# for (i in 1:nrow(sum_future_fish)){
#   segments(x0=sum_future_fish[i,"low"],
#            y0 = sum_future_fish[i,"threshold"],
#            x1=sum_future_fish[i,"high"], col="#0097A7")
# }
# abline(v=0, col="darkgrey", lty=2)
# text(1.8,13.5,"B",font=2)
# 
# plot(threshold~est, data=sum_revisit_state, pch=19, col="#0097A7", ylab="Cluster epsilon (m)", xlab="Estimate",
#      las=1, xlim=c(min(sum_revisit_state$low-0.1), max(sum_revisit_state$high+0.1)), main="")
# for (i in 1:nrow(sum_revisit_state)){
#   segments(x0=sum_revisit_state[i,"low"],
#            y0 = sum_revisit_state[i,"threshold"],
#            x1=sum_revisit_state[i,"high"], col="#0097A7")
# }
# abline(v=0, col="darkgrey", lty=2)
# text(-0.05,13.5,"C",font=2)
# 
# plot(threshold~est, data=sum_revisit_fish, pch=19, col="#0097A7", ylab="Cluster epsilon (m)", xlab="Estimate",
#      las=1, xlim=c(min(sum_revisit_fish$low-0.1), max(sum_revisit_fish$high+0.1)), main="")
# for (i in 1:nrow(sum_revisit_fish)){
#   segments(x0=sum_revisit_fish[i,"low"],
#            y0 = sum_revisit_fish[i,"threshold"],
#            x1=sum_revisit_fish[i,"high"], col="#0097A7")
# }
# abline(v=0, col="darkgrey", lty=2)
# text(0.35,13.5,"D",font=2)
# 
# dev.off()

################################################################################
################################################################################
################################################################################


rm(list = ls())

################################################################################
# END
################################################################################



