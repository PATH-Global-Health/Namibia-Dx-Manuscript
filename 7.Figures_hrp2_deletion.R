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

#Read-in RDS if of previously-run simulations from Setup_hrp2_deletion.R
combined_median_all <- readRDS("Data/Output_Data/combined_ttn_1k_2-100percent.rds")

#Prepare for plotting
combined_result <- combined_median_all %>% 
  group_by(strata) %>%
  summarise(mean_mTTN = round(mean(median, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(strata = factor(strata, levels = c("Abbot90-2 Percent", "Abbot90-10 Percent",
                                            "Abbot90-25 Percent", "Abbot90-50 Percent",
                                            "Abbot90-75 Percent", "Abbot90-100 Percent",
                                            "Rapigen pf-2 Percent", "Rapigen pf-10 Percent",
                                            "Rapigen pf-25 Percent", "Rapigen pf-50 Percent",
                                            "Rapigen pf-75 Percent", "Rapigen pf-100 Percent"))) %>%
  arrange(strata, desc = FALSE) %>%
  separate(col = strata, into = c("Brand", "Deletion Prevalence"), sep = "-", remove = TRUE) %>%
  pivot_wider(names_from = "Deletion Prevalence", values_from = mean_mTTN)  

#Plot combined ttn estimates as boxplot
plot_all <- ggplot(combined_median_all) +
  geom_boxplot(aes(x = id, y = median, fill = id)) +
  labs(x = "Rapid Diagnostic Test", y = "Median Time to Negativity (days post-treatment)",
       title = "Simulated Median Time to Negativity by Prevalence of HRP2 Deletion", fill = "RDT") +
  facet_wrap(. ~ deletion2, nrow = 1)
