library(dplyr)
library(tidyverse)
library(huxtable)

source("Functions/aux_tables.R")
lab <- c("Model","Estimate","Bias","Coverage")
lab_trans <- c("Transition: 1->2","Transition: 1->3","Transition: 2->3")

####################################################################
#################### pop-based study 3/6 years ##################### 
####################################################################
load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xb_3000.RData")
load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xb_10000.RData")

table_xb_3000 <- rbind(rep("n=3000",4),fix_table(pm_df_est_xb_3000))
table_xb_10000 <- rbind(rep("n=10000",4),fix_table(pm_df_est_xb_10000))
table_xb <- cbind(table_xb_3000,
                  table_xb_10000 %>% dplyr::select(-model))

table_xb 


# %>% 
#   as_hux(add_colnames = F) %>% 
#   merge_cells(row=1) %>% 
#   merge_cells(row=2) %>% 
#   merge_cells(row=3) %>%
#   merge_cells(row=10) %>% 
#   merge_cells(row=17) %>%
#   merge_cells(row=24) %>%
#   merge_cells(row=25) %>% 
#   merge_cells(row=32) %>% 
#   merge_cells(row=39) %>% 
#   merge_cells(row=46) %>%
#   merge_cells(row=47) %>% 
#   merge_cells(row=54) %>%
#   merge_cells(row=61)%>% 
#   to_latex() %>% 
#   write("Tables/Simulation_result_table_xb.txt")


####################################################################
#################### irregular observations ######################## 
####################################################################


load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xc_3000.RData")
load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xc_10000.RData")

table_xc_3000 <- rbind(rep("n=3000",4),fix_table(pm_df_est_xc_3000))
table_xc_10000 <- rbind(rep("n=10000",4),fix_table(pm_df_est_xc_10000))
table_xc <- cbind(table_xc_3000,
                  table_xc_10000 %>% dplyr::select(-model))

table_xc

save(table_xb,table_xc,file = c("Tables/Simulation_result_tables.RData"))
