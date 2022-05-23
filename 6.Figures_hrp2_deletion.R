#UCSF Namibia Dataset Visualization
#Authors: William Sheahan and Hannah Slater
#Date Updated: May 20th, 2022

##---HRP2/3 Deletion---##
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

##### Bootstrapping for HRP2/3 Deletion -----
#set the seed
set.seed(59)

#make dataset for bootstrapping from analysis dataset
new_data <- Namibia_analysis %>%
  dplyr::select(pa_id, days_since_tx, Abbot_FK90_flag, LDH_Abbot90_flag, HRP2_Abbot90_flag,
                rapigen_PF_flag, LDH_rapipf_flag, HRP2_rapipf_flag,
                rapigen_PFPV_flag, LDH_rapipfpv_flag)

#set pre-loop parameters
participants <- unique(new_data$pa_id)
hrp2_deletion_prev <- c(0.02, 0.1, 0.25, 0.5, 0.75, 1.0)
n_iterations <- 1000

#Rapigen PFPV not included because there is no HRP2 line, all values are the same regardless of deletion prevalence
#Be warned that this step can take quite a while
#If you want to test the output first, try reducing n_iterations to 10
for (j in 1:length(hrp2_deletion_prev)) {
  
  fk90_surv_object = list()
  rapipf_surv_object = list()
  
  for (i in 1:n_iterations){
    
    sample_index = sample(1:length(participants), 
                          round(length(participants)*hrp2_deletion_prev[j]), replace = FALSE)
    
    sample_participants = participants[sample_index]
    
    new_data_sim = new_data
    
    #If participant is hrp2 deleted, use ldh only flag for survival, otherwise combined
    new_data_sim$sim_fk90_flag = case_when(new_data_sim$pa_id %in% sample_participants & new_data_sim$LDH_Abbot90_flag == 0 ~ 0,
                                           new_data_sim$pa_id %in% sample_participants & new_data_sim$LDH_Abbot90_flag == 1 ~ 1,
                                           new_data_sim$pa_id %not_in% sample_participants ~ new_data_sim$Abbot_FK90_flag)
    new_data_sim$sim_rapipf_flag = case_when(new_data_sim$pa_id %in% sample_participants & new_data_sim$LDH_rapipf_flag == 0 ~ 0,
                                             new_data_sim$pa_id %in% sample_participants & new_data_sim$LDH_rapipf_flag == 1 ~ 1,
                                             new_data_sim$pa_id %not_in% sample_participants ~ new_data_sim$rapigen_PF_flag)
    
    fk90_unique_surv_data = new_data_sim %>%
      filter(!is.na(days_since_tx)) %>%
      group_by(pa_id) %>%
      mutate(time_to_first_neg = min(days_since_tx[sim_fk90_flag == 1], na.rm = TRUE)) %>%
      filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
      filter(days_since_tx == max(days_since_tx))
    
    rapipf_unique_surv_data = new_data_sim %>%
      filter(!is.na(days_since_tx)) %>%
      group_by(pa_id) %>%
      mutate(time_to_first_neg = min(days_since_tx[sim_rapipf_flag == 1], na.rm = TRUE)) %>%
      filter(time_to_first_neg == days_since_tx | is.infinite(time_to_first_neg)) %>%
      filter(days_since_tx == max(days_since_tx))
    
    
    fk90_deletion_model <- survival::survfit(Surv(days_since_tx, sim_fk90_flag) ~ 1, data = fk90_unique_surv_data)
    rapipf_deletion_model <- survival::survfit(Surv(days_since_tx, sim_rapipf_flag) ~ 1, data = rapipf_unique_surv_data)
    
    fk90_surv_object[[i]] = fk90_deletion_model 
    rapipf_surv_object[[i]] = rapipf_deletion_model 
    
  }
  saveRDS(fk90_surv_object, file = paste0("Data/Output_Data/Bootstrapping/fk90_surv_object_j_", j, ".RDS"))
  saveRDS(rapipf_surv_object, file = paste0("Data/Output_Data/Bootstrapping/rapipf_surv_object_j_", j, ".RDS"))
  
  gc()
  
}

#Read in the output RDS files
fk90_prev1 <- readRDS("Data/Output_Data/Bootstrapping/fk90_surv_object_j_1.RDS")
fk90_prev2 <- readRDS("Data/Output_Data/Bootstrapping/fk90_surv_object_j_2.RDS")
fk90_prev3 <- readRDS("Data/Output_Data/Bootstrapping/fk90_surv_object_j_3.RDS")
fk90_prev4 <- readRDS("Data/Output_Data/Bootstrapping/fk90_surv_object_j_4.RDS")
fk90_prev5 <- readRDS("Data/Output_Data/Bootstrapping/fk90_surv_object_j_5.RDS")
fk90_prev6 <- readRDS("Data/Output_Data/Bootstrapping/fk90_surv_object_j_6.RDS")

rapipf_prev1 <- readRDS("Data/Output_Data/Bootstrapping/rapipf_surv_object_j_1.RDS")
rapipf_prev2 <- readRDS("Data/Output_Data/Bootstrapping/rapipf_surv_object_j_2.RDS")
rapipf_prev3 <- readRDS("Data/Output_Data/Bootstrapping/rapipf_surv_object_j_3.RDS")
rapipf_prev4 <- readRDS("Data/Output_Data/Bootstrapping/rapipf_surv_object_j_4.RDS")
rapipf_prev5 <- readRDS("Data/Output_Data/Bootstrapping/rapipf_surv_object_j_5.RDS")
rapipf_prev6 <- readRDS("Data/Output_Data/Bootstrapping/rapipf_surv_object_j_6.RDS")

#Extract median time to negativity estimates from surv objects
Abbot90_median_ttn1 <- surv_median(fk90_prev1, combine = TRUE)
Abbot90_median_ttn1$id <- "Abbot90"
Abbot90_median_ttn1$deletion <- "2 Percent"

Abbot90_median_ttn2 <- surv_median(fk90_prev2, combine = TRUE)
Abbot90_median_ttn2$id <- "Abbot90"
Abbot90_median_ttn2$deletion <- "10 Percent"

Abbot90_median_ttn3 <- surv_median(fk90_prev3, combine = TRUE)
Abbot90_median_ttn3$id <- "Abbot90"
Abbot90_median_ttn3$deletion <- "25 Percent"

Abbot90_median_ttn4 <- surv_median(fk90_prev4, combine = TRUE)
Abbot90_median_ttn4$id <- "Abbot90"
Abbot90_median_ttn4$deletion <- "50 Percent"

Abbot90_median_ttn5 <- surv_median(fk90_prev5, combine = TRUE)
Abbot90_median_ttn5$id <- "Abbot90"
Abbot90_median_ttn5$deletion <- "75 Percent"

Abbot90_median_ttn6 <- surv_median(fk90_prev6, combine = TRUE)
Abbot90_median_ttn6$id <- "Abbot90"
Abbot90_median_ttn6$deletion <- "100 Percent"

rapipf_median_ttn1 <- surv_median(rapipf_prev1, combine = TRUE)
rapipf_median_ttn1$id <- "Rapigen pf"
rapipf_median_ttn1$deletion <- "2 Percent"

rapipf_median_ttn2 <- surv_median(rapipf_prev2, combine = TRUE)
rapipf_median_ttn2$id <- "Rapigen pf"
rapipf_median_ttn2$deletion <- "10 Percent"

rapipf_median_ttn3 <- surv_median(rapipf_prev3, combine = TRUE)
rapipf_median_ttn3$id <- "Rapigen pf"
rapipf_median_ttn3$deletion <- "25 Percent"

rapipf_median_ttn4 <- surv_median(rapipf_prev4, combine = TRUE)
rapipf_median_ttn4$id <- "Rapigen pf"
rapipf_median_ttn4$deletion <- "50 Percent"

rapipf_median_ttn5 <- surv_median(rapipf_prev5, combine = TRUE)
rapipf_median_ttn5$id <- "Rapigen pf"
rapipf_median_ttn5$deletion <- "75 Percent"

rapipf_median_ttn6 <- surv_median(rapipf_prev6, combine = TRUE)
rapipf_median_ttn6$id <- "Rapigen pf"
rapipf_median_ttn6$deletion <- "100 Percent"

#Combine disparate ttn estimate files
combined_median_ttn1 <- rbind(Abbot90_median_ttn1, rapipf_median_ttn1)
combined_median_ttn2 <- rbind(Abbot90_median_ttn2, rapipf_median_ttn2)
combined_median_ttn3 <- rbind(Abbot90_median_ttn3, rapipf_median_ttn3)
combined_median_ttn4 <- rbind(Abbot90_median_ttn4, rapipf_median_ttn4)
combined_median_ttn5 <- rbind(Abbot90_median_ttn5, rapipf_median_ttn5)
combined_median_ttn6 <- rbind(Abbot90_median_ttn6, rapipf_median_ttn6)

combined_median_all <- rbind(combined_median_ttn1, combined_median_ttn2,
                             combined_median_ttn3, combined_median_ttn4,
                             combined_median_ttn5, combined_median_ttn6)

combined_median_all <- combined_median_all %>%
  mutate(deletion2 = if_else(deletion == "2 Percent", .02, 0),
         deletion2 = if_else(deletion == "10 Percent", .10, deletion2),
         deletion2 = if_else(deletion == "25 Percent", .25, deletion2),
         deletion2 = if_else(deletion == "50 Percent", .50, deletion2),
         deletion2 = if_else(deletion == "75 Percent", .75, deletion2),
         deletion2 = if_else(deletion == "100 Percent", 1.00, deletion2))

combined_median_all <- combined_median_all %>%
  unite('strata', c(id, deletion), remove = FALSE, sep = "-")


#Read-in RDS if of previously-run simulations from Setup_hrp2_deletion.R
#saveRDS(combined_median_all, file = "Data/Output_Data/combined_ttn_1k_2-100percent.rds")
#combined_median_all <- readRDS("Data/Output_Data/combined_ttn_1k_2-100percent.rds")

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
