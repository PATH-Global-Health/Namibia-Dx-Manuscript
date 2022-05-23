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


#Survival Fits for curve generation
fit_hs_rdt_cleared <- survival::survfit(Surv(days_since_tx, urdt_surv) ~ 1, data = hs_rdt_cleared)
fit_std_rdt_cleared <- survival::survfit(Surv(days_since_tx, std_rdt_surv) ~ 1, data = std_rdt_cleared)                   

#Pull out individual time to negativity for each iteration to include on plot/tables
urdt_median_ttn <- surv_median(fit_hs_rdt_cleared, combine = TRUE)
std_rdt_median_ttn <- surv_median(fit_std_rdt_cleared, combine = TRUE)

##### Combined Plots -----
#Combine multiple surv objects/curves onto one plot with Dummy Dataframes
demo_data_std <- rbind((cbind(hs_rdt_cleared, type = "uRDT")), (cbind(std_rdt_cleared, type = "standard rdt")))

#Combined model fits from dummy data
all.fit.std <- survfit(Surv(days_since_tx) ~ type, data = demo_data_std)

#Combined LDH and HRP2 lines for uRDT and standard RDT
surv_median_std <- as.vector(summary(all.fit.std)$table[, "median"])
df_std <- data.frame(x1 = surv_median_std, x2 = surv_median_std,
                     y1 = rep(0, length(surv_median_std)), y2 = rep(0.5, length(surv_median_std)))


#Plot standard RDT vs uRDT, only available as both pLDH and HRP2 combined lines
std_urdt_full <- ggsurv(all.fit.std, plot.cens = FALSE, size.est = 0.75, 
                        surv.col = c("dodgerblue", "mediumvioletred")) +
  guides(linetype = "none") +
  scale_color_manual(labels = c("HS-RDT", "Standard RDT"), 
                     values = c("dodgerblue", "mediumvioletred"),
                     "Diagnostic Type") +
  labs(title = "Combined Survival Curve for Standard and HS-RDT Tests", subtitle = "Participants Who Cleared Infection",
       x = "Days Post Treatment", y = "Probability of Positive Malaria Test") +
  theme_classic()

#add lines for median ttn's from surv_median dataframes
std_urdt_full_final <- std_urdt_full + 
  geom_segment(aes(x = 0, xend = 68, y = 0.5, yend = 0.5), linetype = "dashed") +
  geom_segment(aes(x = 42, xend = 42, y = 0, yend = 0.5), linetype = "dashed", color = "mediumvioletred") +
  geom_segment(aes(x = 68, xend = 68, y = 0, yend = 0.5), linetype = "dashed", color = "dodgerblue")
