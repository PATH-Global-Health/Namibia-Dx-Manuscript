#UCSF Namibia Dataset Visualization
#Authors: William Sheahan and Hannah Slater
#Date Updated: May 20th, 2022

##---Time to Negativity Analysis---##
remove(list = ls())
options(scipen = 10000000)
#Packages
library(dplyr)
library(fields)
library(tidyverse)
library(raster)
library(gridExtra)
library(survival)
library(survminer)
library(ggfortify)
library(GGally)

#Read in full analysis file with probability surface data
Namibia_analysis <- read_csv("Data/Output_Data/Namibia_full_flags_analysis.csv")


##### Survival Curves -----
#dataframe for fk90 survival object
fk90_cleared <- Namibia_analysis %>%
  filter(!is.na(days_since_tx)) %>%
  group_by(pa_id) %>%
  mutate(time_to_first_neg = min(days_since_tx[Abbot_FK90_flag == 1], na.rm = TRUE)) %>%
  filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
  filter(days_since_tx == max(days_since_tx)) %>%
  dplyr::select(pa_id, label, days_since_tx, time_to_first_neg, urdt_surv, prob_Abbot_FK90, Abbot_FK90_flag)

#dataframe for fk90 survival object using LDH flag
fk90_cleared_ldh <- Namibia_analysis %>%
  filter(!is.na(days_since_tx)) %>%
  group_by(pa_id) %>%
  mutate(time_to_first_neg = min(days_since_tx[LDH_Abbot90_flag == 1], na.rm = TRUE)) %>%
  filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
  filter(days_since_tx == max(days_since_tx)) %>%
  dplyr::select(pa_id, label, time_point_days, days_since_tx, time_to_first_neg, urdt_surv, prob_Abbot_FK90, Abbot_FK90_flag, LDH_Abbot90_flag)

#dataframe for fk90 survival object using HRP2 flag
fk90_cleared_hrp2 <- Namibia_analysis %>%
  filter(!is.na(days_since_tx)) %>%
  group_by(pa_id) %>%
  mutate(time_to_first_neg = min(days_since_tx[HRP2_Abbot90_flag == 1], na.rm = TRUE)) %>%
  filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
  filter(days_since_tx == max(days_since_tx)) %>%
  dplyr::select(pa_id, label, time_point_days, days_since_tx, time_to_first_neg, urdt_surv, prob_Abbot_FK90, Abbot_FK90_flag, HRP2_Abbot90_flag)

#dataframe for rapigen PF survival object
rapigen_PF_cleared <- Namibia_analysis %>%
  filter(!is.na(days_since_tx)) %>%
  group_by(pa_id) %>%
  mutate(time_to_first_neg = min(days_since_tx[rapigen_PF_flag == 1], na.rm = TRUE)) %>%
  filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
  filter(days_since_tx == max(days_since_tx)) %>%
  dplyr::select(pa_id, label, time_point_days, days_since_tx, time_to_first_neg, urdt_surv, prob_rapigen_PF, rapigen_PF_flag)

#dataframe for rapigen PF survival object using LDH flag
rapigen_PF_cleared_ldh <- Namibia_analysis %>%
  filter(!is.na(days_since_tx)) %>%
  group_by(pa_id) %>%
  mutate(time_to_first_neg = min(days_since_tx[LDH_rapipf_flag == 1], na.rm = TRUE)) %>%
  filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
  filter(days_since_tx == max(days_since_tx)) %>%
  dplyr::select(pa_id, label, time_point_days, days_since_tx, time_to_first_neg, urdt_surv, prob_rapigen_PF, rapigen_PF_flag, LDH_rapipf_flag)

#dataframe for rapigen PF survival object using HRP2 flag
rapigen_PF_cleared_hrp2 <- Namibia_analysis %>%
  filter(!is.na(days_since_tx)) %>%
  group_by(pa_id) %>%
  mutate(time_to_first_neg = min(days_since_tx[HRP2_rapipf_flag == 1], na.rm = TRUE)) %>%
  filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
  filter(days_since_tx == max(days_since_tx)) %>%
  dplyr::select(pa_id, label, time_point_days, days_since_tx, time_to_first_neg, urdt_surv, prob_rapigen_PF, rapigen_PF_flag, HRP2_rapipf_flag)

#Survival Fits for curve generation
fit_fk90_cleared <- survival::survfit(Surv(days_since_tx, Abbot_FK90_flag) ~ 1, data = fk90_cleared)
fit_fk90_cleared_ldh <- survival::survfit(Surv(days_since_tx, LDH_Abbot90_flag) ~ 1, data = fk90_cleared_ldh)
fit_fk90_cleared_hrp2 <- survival::survfit(Surv(days_since_tx, HRP2_Abbot90_flag) ~ 1, data = fk90_cleared_hrp2)
fit_rapigen_PF_cleared <- survival::survfit(Surv(days_since_tx, rapigen_PF_flag) ~ 1, data = rapigen_PF_cleared)
fit_rapigen_PF_cleared_ldh <- survival::survfit(Surv(days_since_tx, LDH_rapipf_flag) ~ 1, data = rapigen_PF_cleared_ldh)
fit_rapigen_PF_cleared_hrp2 <- survival::survfit(Surv(days_since_tx, HRP2_rapipf_flag) ~ 1, data = rapigen_PF_cleared_hrp2)

#Pull out individual time to negativity for each iteration to include on plot/tables
fk90_full <- surv_median(fit_fk90_cleared, combine = TRUE)
fk90_ldh <- surv_median(fit_fk90_cleared_ldh, combine = TRUE)
fk90_hrp2 <- surv_median(fit_fk90_cleared_hrp2, combine = TRUE)
rapipf_full <- surv_median(fit_rapigen_PF_cleared, combine = TRUE)
rapipf_ldh <- surv_median(fit_rapigen_PF_cleared_ldh, combine = TRUE)
rapipf_hrp2 <- surv_median(fit_rapigen_PF_cleared_hrp2, combine = TRUE)

##### Combined Plots -----
#Combine multiple surv objects/curves onto one plot with Dummy Dataframes
demo_data_fk90 <- rbind((cbind(fk90_cleared, type = "Abbot FK90 Full")),
                        (cbind(fk90_cleared_ldh, type = "Abbot FK90 LDH Only")),
                        (cbind(fk90_cleared_hrp2, type = "Abbot FK90 HRP2 Only")))

demo_data_rapipf <- rbind((cbind(rapigen_PF_cleared, type = "rapigenPF Full")),
                          (cbind(rapigen_PF_cleared_ldh, type = "rapigenPF LDH Only")),
                          (cbind(rapigen_PF_cleared_hrp2, type = "rapigenPFPV HRP2 Only")))

#Combined model fits from dummy data
all.fit.fk90 <- survfit(Surv(days_since_tx) ~ type, data = demo_data_fk90)
all.fit.rapipf <- survfit(Surv(days_since_tx) ~ type, data = demo_data_rapipf)


#Abbot fk90 Combined, LDH only, and HRP2 only
surv_median_fk90 <- as.vector(summary(all.fit.fk90)$table[, "median"])
df_fk90 <- data.frame(x1 = surv_median_fk90, x2 = surv_median_fk90,
                      y1 = rep(0, length(surv_median_fk90)), y2 = rep(0.5, length(surv_median_fk90)))

#Rapigen pfpv Combined, LDH only, and HRP2 only
surv_median_rapipf <- as.vector(summary(all.fit.rapipf)$table[, "median"])
df_rapipf <- data.frame(x1 = surv_median_rapipf, x2 = surv_median_rapipf,
                        y1 = rep(0, length(surv_median_rapipf)), y2 = rep(0.5, length(surv_median_rapipf)))


#Plot Abbot FK90 lines comparing combined vs HRP2 and LDH only lines
fk90_combined <- ggsurv(all.fit.fk90, plot.cens = FALSE, 
                        surv.col = c("dodgerblue", "goldenrod2", "firebrick3"), size.est = 1.0) +
  guides(linetype = "none") +
  scale_color_manual(labels = c("Abbot FK90 Full", "Abbot FK90 HRP2 Only", "Abbot FK90 LDH Only"), 
                     values = c("dodgerblue", "goldenrod2", "firebrick3"),
                     "Diagnostic Type") +
  labs(title = "Combined Survival Curve for Abbot FK90 Test", subtitle = "Participants Who Cleared Infection",
       x = "Days Post Treatment", y = "Probability of Positive Malaria Test") +
  theme_classic() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3))

#add lines for median ttn's from surv_median dataframes
fk90_final <- fk90_combined + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_segment(aes(x = 7, xend = 7, y = 0, yend = 0.5), linetype = "dashed", color = "firebrick3", size = 1.0) +
  geom_segment(aes(x = 43, xend = 43, y = 0, yend = 0.5), linetype = "dashed", color = "goldenrod2", size = 1.0) +
  geom_segment(aes(x = 49, xend = 49, y = 0, yend = 0.5), linetype = "dashed", color = "dodgerblue", size = 1.0)

#Plot rapigen pfpv lines comparing combined vs HRP2 and LDH only lines
rapipf_combined <- ggsurv(all.fit.rapipf, plot.cens = FALSE, 
                          surv.col = c("dodgerblue", "goldenrod2", "firebrick3"), size.est = 1.0) +
  guides(linetype = "none") +
  scale_color_manual(labels = c("Rapigen PF Full", "Rapigen PF HRP2 Only", "Rapigen PF LDH Only"), 
                     values = c("dodgerblue", "goldenrod2", "firebrick3"),
                     "Diagnostic Type") +
  labs(title = "Combined Survival Curve for Rapigen PF Test", subtitle = "Participants Who Cleared Infection",
       x = "Days Post Treatment", y = "Probability of Positive Malaria Test") +
  theme_classic() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3))

#add lines for median ttn's from surv_median dataframes
#slight adjustment to line for HRP2/LDH lines so they are both visible despite having the same true value
rapipf_final <- rapipf_combined + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_segment(aes(x = 9, xend = 9, y = 0, yend = 0.5), linetype = "dashed", color = "firebrick3", size = 1.0) +
  geom_segment(aes(x = 55.75, xend = 55.75, y = 0, yend = 0.5), linetype = "dashed", color = "goldenrod2", size = 1.0) +
  geom_segment(aes(x = 56, xend = 56, y = 0, yend = 0.5), linetype = "dashed", color = "dodgerblue", size = 1.0)

#arrange two figures side-by-side
ggarrange(fk90_final, rapipf_final)
