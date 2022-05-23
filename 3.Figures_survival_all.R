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
#dataframe for HS-RDT survival object
hs_rdt_cleared <- Namibia_analysis %>%
  filter(!is.na(days_since_tx)) %>%
  group_by(pa_id) %>%
  mutate(left_study = max(days_since_tx)) %>%
  filter(left_study == days_since_tx) %>%
  dplyr::select(pa_id, label, left_study, days_since_tx, urdt_surv)

#dataframe for standard rdt survival object
std_rdt_cleared <- Namibia_analysis %>%
  filter(!is.na(days_since_tx)) %>%
  group_by(pa_id) %>%
  mutate(time_to_first_neg = min(days_since_tx[std_rdt_surv == 1], na.rm = TRUE)) %>%
  filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
  filter(days_since_tx == max(days_since_tx)) %>%
  dplyr::select(pa_id, label, days_since_tx, std_rdt_surv)

#dataframe for fk90 survival object
fk90_cleared <- Namibia_analysis %>%
  filter(!is.na(days_since_tx)) %>%
  group_by(pa_id) %>%
  mutate(time_to_first_neg = min(days_since_tx[Abbot_FK90_flag == 1], na.rm = TRUE)) %>%
  filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
  filter(days_since_tx == max(days_since_tx)) %>%
  dplyr::select(pa_id, label, days_since_tx, time_to_first_neg, urdt_surv, prob_Abbot_FK90, Abbot_FK90_flag)

#dataframe for rapigen PF survival object
rapigen_PF_cleared <- Namibia_analysis %>%
  filter(!is.na(days_since_tx)) %>%
  group_by(pa_id) %>%
  mutate(time_to_first_neg = min(days_since_tx[rapigen_PF_flag == 1], na.rm = TRUE)) %>%
  filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
  filter(days_since_tx == max(days_since_tx)) %>%
  dplyr::select(pa_id, label, time_point_days, days_since_tx, time_to_first_neg, urdt_surv, prob_rapigen_PF, rapigen_PF_flag)

#dataframe for rapigen PF/PV survival object
rapigen_PFPV_cleared <- Namibia_analysis %>%
  filter(!is.na(days_since_tx)) %>%
  group_by(pa_id) %>%
  mutate(time_to_first_neg = min(days_since_tx[rapigen_PFPV_flag == 1], na.rm = TRUE)) %>%
  filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
  filter(days_since_tx == max(days_since_tx)) %>%
  dplyr::select(pa_id, label, time_point_days, days_since_tx, time_to_first_neg, urdt_surv, prob_rapigen_PFPV, rapigen_PFPV_flag)

#Survival Fits for curve generation
fit_hs_rdt_cleared <- survival::survfit(Surv(days_since_tx, urdt_surv) ~ 1, data = hs_rdt_cleared)
fit_std_rdt_cleared <- survival::survfit(Surv(days_since_tx, std_rdt_surv) ~ 1, data = std_rdt_cleared)                   
fit_fk90_cleared <- survival::survfit(Surv(days_since_tx, Abbot_FK90_flag) ~ 1, data = fk90_cleared)
fit_rapigen_PF_cleared <- survival::survfit(Surv(days_since_tx, rapigen_PF_flag) ~ 1, data = rapigen_PF_cleared)
fit_rapigen_PFPV_cleared <- survival::survfit(Surv(days_since_tx, rapigen_PFPV_flag) ~ 1, data = rapigen_PFPV_cleared)

#Pull out individual time to negativity for each iteration to include on plot/tables
urdt_median_ttn <- surv_median(fit_hs_rdt_cleared, combine = TRUE)
std_rdt_median_ttn <- surv_median(fit_std_rdt_cleared, combine = TRUE)
fk90_full <- surv_median(fit_fk90_cleared, combine = TRUE)
rapipf_full <- surv_median(fit_rapigen_PF_cleared, combine = TRUE)
rapipfpv_full <- surv_median(fit_rapigen_PFPV_cleared, combine = TRUE)

##### Combined Plots -----
#Combine multiple surv objects/curves onto one plot with Dummy Dataframes
demo_data <- rbind((cbind(hs_rdt_cleared, type = "uRDT")), (cbind(fk90_cleared, type = "fk90")),
                   (cbind(rapigen_PF_cleared, type = "rapigenPF")), (cbind(rapigen_PFPV_cleared, type = "rapigenPFPV")))

#Combined model fits from dummy data
all.fit <- survfit(Surv(days_since_tx) ~ type, data = demo_data)

#Combined LDH and HRP2 lines for uRDT, Abbot, and Rapigen
surv_median_all <- as.vector(summary(all.fit)$table[, "median"])
df_all <- data.frame(x1 = surv_median_all, x2 = surv_median_all,
                     y1 = rep(0, length(surv_median_all)), y2 = rep(0.5, length(surv_median_all)))


#Plotting the combined model fits
#Plot combined pLDH and HRP2 lines for HS-RDT, Abbot, and Rapigen
all_full <- ggsurv(all.fit, plot.cens = FALSE, size.est = 0.75, 
                   surv.col = c("dodgerblue", "goldenrod2", "firebrick3", "darkolivegreen")) +
  guides(linetype = "none") +
  scale_color_manual(labels = c("HS-RDT", "Abbot FK90", "Rapigen pf", "Rapigen pf/pv"), 
                     values = c("dodgerblue", "goldenrod2", "firebrick3", "darkolivegreen"),
                     "Diagnostic Type") +
  labs(title = "Combined Survival Curve for Abbot, Rapigen, and HS-RDT Tests", subtitle = "Participants Who Cleared Infection",
       x = "Days Post Treatment", y = "Probability of Positive Malaria Test") +
  theme_classic()

#add lines for median ttn's from surv_median dataframes
all_full_final <- all_full + 
  geom_segment(aes(x = 0, xend = 68, y = 0.5, yend = 0.5), linetype = "dashed") +
  geom_segment(aes(x = 10, xend = 10, y = 0, yend = 0.5), linetype = "dashed", color = "darkolivegreen") +
  geom_segment(aes(x = 49, xend = 49, y = 0, yend = 0.5), linetype = "dashed", color = "goldenrod2") +
  geom_segment(aes(x = 56, xend = 56, y = 0, yend = 0.5), linetype = "dashed", color = "firebrick3") +
  geom_segment(aes(x = 68, xend = 68, y = 0, yend = 0.5), linetype = "dashed", color = "dodgerblue")
