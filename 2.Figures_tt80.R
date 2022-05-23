#UCSF Namibia Dataset Visualization
#Authors: William Sheahan and Hannah Slater based on previous work by Emily Reichert
#Date Updated: May 20th, 2022


##--- biphasic exp decay ---##
#From Reichert et al. 
#https://github.com/PATH-Global-Health/Reichert_Mali_antigens/blob/master/mali_tbv_hrp2_decay.Rmd
#test for optimal switch point posttreatment, minimizing AIC and BIC
#then calculate time to 80 pg/ml and half-life of each antigen
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


##### Figure X -----
#histogram of time to 80 pg/ml HRP2
time_to_80 <- rep(NA, length(unique(Namibia_lme$individual_lme)))

for(i in 1:length(unique(Namibia_lme$individual_lme))){
  time_to_80[i] = which.min(abs(10^lmepred_biphase_Namibia_inds[i, ] - 80))
}

time_to_80_plotting <- as.data.frame(time_to_80) 

brx <- pretty(range(time_to_80_plotting$time_to_80), 
              n = nclass.Sturges(time_to_80_plotting$time_to_80),min.n = 1)

tt80_hrp2 <- ggplot(time_to_80_plotting, aes(time_to_80)) +
  geom_histogram(color="darkgray",fill="goldenrod2", breaks=brx) + 
  scale_x_continuous("Days Since Treatment") + 
  labs(title = "Histogram of Time to Reach 80 pg/ml HRP2") + 
  theme_classic()

#half-life calc biphasic hrp2
a <- 10^max(lmepred_biphase_Namibia_pop)/2
b <- which.min(abs(10^lmepred_biphase_Namibia_pop - a))
days[b]

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

#calculate time to 80 pg/ml LDH
time_to_80_ldh <- rep(NA, length(unique(Namibia_lme$individual_lme)))

for(i in 1:length(unique(Namibia_lme$individual_lme))){
  time_to_80_ldh[i] = which.min(abs(10^lmepred_biphase_Namibia_inds2[i, ] - 80))
}

time_to_80_plotting_ldh <- as.data.frame(time_to_80_ldh) 

brx_ldh <- pretty(range(time_to_80_plotting_ldh$time_to_80_ldh), 
                  n = nclass.Sturges(time_to_80_plotting_ldh$time_to_80_ldh),min.n = 1)

tt80_pldh <- ggplot(time_to_80_plotting_ldh, aes(time_to_80_ldh)) +
  geom_histogram(color="darkgray",fill="firebrick3", breaks=brx_ldh) + 
  scale_x_continuous("Days Since Treatment") + 
  labs(title = "Histogram of Time to Reach 80 pg/ml pLDH") +
  theme_classic()

#half-life calc biphasic pLDH
a2 <- 10^max(lmepred_biphase_Namibia_pop2)/2
b2 <- which.min(abs(10^lmepred_biphase_Namibia_pop2 - a2))
days[b2]

#Arrange two figures side by side
ggarrange(tt80_hrp2, tt80_pldh)
