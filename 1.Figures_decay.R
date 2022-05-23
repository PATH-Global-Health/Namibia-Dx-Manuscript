#UCSF Namibia Dataset Visualization
#Authors: William Sheahan and Hannah Slater based on previous work by Emily Reichert
#Date Updated: May 20th, 2022


##--- biphasic exp decay ---##
#From Reichert et al. 
#https://github.com/PATH-Global-Health/Reichert_Mali_antigens/blob/master/mali_tbv_hrp2_decay.Rmd
#test for optimal switch point posttreatment, minimizing AIC and BIC
library(tidyverse)
library(lme4)
library(nlme)
library(arm)
library(lmerTest)
library(purrr)

remove(list = ls())
options(scipen = 10000000)
#Read in cleared infection subset analysis file
cleared_infection <- read_csv("../Namibia_Analysis_WS_21/Data/Output_Data/cleared_infection_subset_2022-02-03.csv")

#set max time based on cleaned infection subset (should be 105)
days <- seq(0, max(cleared_infection$days_since_tx, na.rm = TRUE), 1)

#prepare dataset for biphasic decay modeling (remove missing data ~ 10 obs)
Namibia_lme <- cleared_infection %>% 
  dplyr::select(individual_lme = pa_id, HRP2_lme = HRP2_pg_ml, Time_lme = days_since_tx) %>%
  filter(!is.na(Time_lme)) %>%
  filter(!is.na(HRP2_lme))

#Initialize AIC/BIC for determining switch time
aic_values <- rep(NA, 105)
bic_values <- rep(NA, 105)

for(i in 1:105){
  switch_time <- i
  lmefit_biphase_Namibia <- lme4::lmer(log10(HRP2_lme) ~ Time_lme + ifelse(Time_lme <= (switch_time), 0, Time_lme-(switch_time)) + (1|individual_lme), REML = F, data = Namibia_lme)
  aic_values[i] <- AIC(lmefit_biphase_Namibia)
  bic_values[i] <- BIC(lmefit_biphase_Namibia)
}

#find point where AIC/BIC are minimized (day 13)
min(aic_values)
min(bic_values)
#confirm that estimated minimum date is correct = matches value from min(aic_values)
aic_values[13]
bic_values[13]

# AIC/BIC smallest for day 13 and biphasic decay
switch_time <- 13
lmefit_biphase_Namibia <- lme4::lmer(log10(HRP2_lme) ~ Time_lme + ifelse(Time_lme <= (switch_time), 0, Time_lme-(switch_time)) + (1|individual_lme), REML = F, data = Namibia_lme)

#Pullout Coefficients for random and fixed effects
coef(lmefit_biphase_Namibia) 
lme4::fixef(lmefit_biphase_Namibia)
decay2 = lme4::fixef(lmefit_biphase_Namibia)[3] + lme4::fixef(lmefit_biphase_Namibia)[2]
lme4::ranef(lmefit_biphase_Namibia)

# Using simulate to estimate the posterior distribution of intercept and slopes
lmefit_biphase_Namibia_sim <- arm::sim(lmefit_biphase_Namibia, n.sim = 10000)

# 95% CI of paras based on sim()
lmefit_biphase_Namibia_sim_interceptCI <- as.numeric(quantile(lme4::fixef(lmefit_biphase_Namibia_sim)[,1], c(.025, .975)))
lmefit_biphase_Namibia_sim_slope1CI <- as.numeric(quantile(lme4::fixef(lmefit_biphase_Namibia_sim)[,2], c(.025, .975)))
lmefit_biphase_Namibia_sim_slope2CI <- as.numeric(quantile((lme4::fixef(lmefit_biphase_Namibia_sim)[,3]), c(.025, .975)))
# predict the average value over all posttreatment days using lmer() fitted paras
lmepred_biphase_Namibia_pop <- lme4::fixef(lmefit_biphase_Namibia)["(Intercept)"] + lme4::fixef(lmefit_biphase_Namibia)["Time_lme"] * days +
  lme4::fixef(lmefit_biphase_Namibia)[3] * ifelse(days <= (switch_time), 0, days- switch_time)

# predict inds value for each day using lmer() fitted paras
lmepred_biphase_Namibia_inds <- matrix(NA, nrow=length(unique(Namibia_lme$individual_lme)), ncol=length(days))
for(i in 1:length(unique(Namibia_lme$individual_lme))){
  lmepred_biphase_Namibia_inds[i, ] <- coef(lmefit_biphase_Namibia)$individual_lme[i,1] + coef(lmefit_biphase_Namibia)$individual_lme[i,2] * days +
    coef(lmefit_biphase_Namibia)$individual_lme[i,3] * ifelse(days <= (switch_time), 0, days-(switch_time))
}

### Figure X -----
# plot on log10-scale
#tiff("Images/final/biphase_hrp2_log_HS.tiff", width = 6, height = 5, units = "in", res = 300, compression = "lzw")
plot.new()
par(mfrow=c(1,1), mar = c(4,7,3,1.5))
plot(days, lmepred_biphase_Namibia_pop, tck=.01, type="l", 
     ylim=range(c(log10(Namibia_lme$HRP2_lme), lmepred_biphase_Namibia_pop, 
                  lmepred_biphase_Namibia_inds)), ylab="", 
     xlab = "Days Post Treatment", main = "Biphasic HRP2 decay, log10 scale", lwd=2, font=2, font.lab=2,
     yaxt = "n")
for(i in 1:length(unique(Namibia_lme$individual_lme))){
  lines(days, lmepred_biphase_Namibia_inds[i, ], col="pink")
}
points(Namibia_lme$Time_lme, log10(Namibia_lme$HRP2_lme), pch=20, col = alpha("firebrick3",0.3))
lines(days, lmepred_biphase_Namibia_pop, col = "white", lwd=4)
lines(days, lmepred_biphase_Namibia_pop, lwd=3)
axis(2, seq(-5, 10, by = 2), 10^seq(-5, 10, by = 2), las=1)
mtext("HRP2 concentration (pg/ml)", side = 2, line = 5, cex = 1.4, font = 2)

dev.off()


############################# LDH Version -------------------------
Namibia_lme2 <- cleared_infection %>% 
  dplyr::select(individual_lme = pa_id, pLDH_lme = LDH_Pan_pg_ml, Time_lme = days_since_tx) %>%
  filter(!is.na(Time_lme)) %>%
  filter(!is.na(pLDH_lme))

#Initialize AIC/BIC for determining switch time
aic_values_ldh <- rep(NA, 105)
bic_values_ldh <- rep(NA, 105)

for(i in 1:105){
  switch_time2 <- i
  lmefit_biphase_Namibia2 <- lme4::lmer(log10(pLDH_lme) ~ Time_lme + ifelse(Time_lme <= (switch_time2), 0, Time_lme-(switch_time2)) + (1|individual_lme), REML = F, data = Namibia_lme2)
  aic_values_ldh[i] <- AIC(lmefit_biphase_Namibia2)
  bic_values_ldh[i] <- BIC(lmefit_biphase_Namibia2)
}

#find point where AIC/BIC are minimized
min(aic_values_ldh)
min(bic_values_ldh)
#confirm that the switch point has the minimized aic/bic 
aic_values_ldh[8]
bic_values_ldh[8]

# AIC/BIC smallest for (8 for ldh biphasic)
switch_time_ldh <- 8
lmefit_biphase_Namibia2 <- lme4::lmer(log10(pLDH_lme) ~ Time_lme + ifelse(Time_lme <= (switch_time_ldh), 0, Time_lme-(switch_time_ldh)) + (1|individual_lme), REML = F, data = Namibia_lme2)

#Pullout Coefficients for random and fixed effects
coef(lmefit_biphase_Namibia2)
lme4::fixef(lmefit_biphase_Namibia2)
decay2 = lme4::fixef(lmefit_biphase_Namibia2)[3] + lme4::fixef(lmefit_biphase_Namibia2)[2]
lme4::ranef(lmefit_biphase_Namibia2)

# Using simulate to estimate the posterior distribution of intercept and slopes
lmefit_biphase_Namibia_sim2 <- arm::sim(lmefit_biphase_Namibia2, n.sim = 10000)
par(mfrow=c(2,2))

# 95% CI of paras based on sim()
lmefit_biphase_Namibia_sim2_interceptCI <- as.numeric(quantile(lme4::fixef(lmefit_biphase_Namibia_sim2)[,1], c(.025, .975)))
lmefit_biphase_Namibia_sim2_slope1CI <- as.numeric(quantile(lme4::fixef(lmefit_biphase_Namibia_sim2)[,2], c(.025, .975)))
lmefit_biphase_Namibia_sim2_slope2CI <- as.numeric(quantile((lme4::fixef(lmefit_biphase_Namibia_sim2)[,3]), c(.025, .975)))
# predict the average value over all posttreatment days using lmer() fitted paras
lmepred_biphase_Namibia_pop2 <- lme4::fixef(lmefit_biphase_Namibia2)["(Intercept)"] + lme4::fixef(lmefit_biphase_Namibia2)["Time_lme"] * days +
  lme4::fixef(lmefit_biphase_Namibia2)[3] * ifelse(days <= (switch_time_ldh), 0, days- switch_time_ldh)

# predict inds value for each day using lmer() fitted paras
lmepred_biphase_Namibia_inds2 <- matrix(NA, nrow=length(unique(Namibia_lme2$individual_lme)), ncol=length(days))
for(i in 1:length(unique(Namibia_lme2$individual_lme))){
  lmepred_biphase_Namibia_inds2[i, ] <- coef(lmefit_biphase_Namibia2)$individual_lme[i,1] + coef(lmefit_biphase_Namibia2)$individual_lme[i,2] * days +
    coef(lmefit_biphase_Namibia2)$individual_lme[i,3] * ifelse(days <= (switch_time_ldh), 0, days-(switch_time_ldh))
}

# plot on log10-scale
#tiff("Images/final/biphase_ldh_log.tiff", width = 6, height = 5, units = "in", res = 300, compression = "lzw")

plot.new()
par(mfrow=c(1,1), mar = c(4,7,3,1.5))
plot(days, lmepred_biphase_Namibia_pop2, tck=.01, type="l", 
     ylim=range(c(log10(Namibia_lme2$pLDH_lme), lmepred_biphase_Namibia_pop2, lmepred_biphase_Namibia_inds2)), 
     ylab="", yaxt = "n", 
     xlab = "Days Post Treatment", main = "Biphasic pLDH decay, log10 scale", lwd=2, font=2, font.lab=2)
for(i in 1:length(unique(Namibia_lme2$individual_lme))){
  lines(days, lmepred_biphase_Namibia_inds2[i, ], col="pink")
}
points(Namibia_lme2$Time_lme, log10(Namibia_lme2$pLDH_lme), pch=20, col = alpha("firebrick3", 0.3))
lines(days, lmepred_biphase_Namibia_pop2, col = "white", lwd=4)
lines(days, lmepred_biphase_Namibia_pop2, lwd=3)
axis(2, seq(0, 8, by = 2), 10^seq(0, 8, by = 2), las=1)
mtext("pLDH concentration (pg/ml)", side = 2, line = 5, cex = 1.4, font = 2)

dev.off()
