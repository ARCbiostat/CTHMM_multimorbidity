library(dplyr)
library(tidyverse)
library(huxtable)

source("Functions/aux_tables.R")
lab <- c("Model","Estimate","Bias","S.E. Bias","Coverage")

####################################################################
#################### pop-based study 3/6 years ##################### 
####################################################################
load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xb_3000.RData")
load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xb_10000.RData")

table_xb_3000 <-fix_table(pm_df_est_xb_3000)
table_xb_10000 <- fix_table(pm_df_est_xb_10000)

####################################################################
#################### irregular observations ######################## 
####################################################################


load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xc_3000.RData")
load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xc_10000.RData")

table_xc_3000 <- fix_table(pm_df_est_xc_3000)
table_xc_10000 <- fix_table(pm_df_est_xc_10000)

####################################################################
#################### CREATE LATEX TABLES ######################## 
####################################################################

# TRANSITION 1

table_t1 <- rbind(
                  c("firstcol",rep("PS",8)),
                  c("firstcol",rep("n=3000",4),rep("n=10000",4)),
                  c(lab,lab[-1]),
                  cbind(table_xb_3000 %>% filter(trans==1) %>% dplyr::select(-trans),
                  table_xb_10000 %>% filter(trans==1) %>% dplyr::select(-trans,-model)),
                  c("firstcol",rep("IV",8)),
                  c("firstcol",rep("n=3000",4),rep("n=10000",4)),
                  c(lab,lab[-1]),
                  cbind(table_xc_3000 %>% filter(trans==1) %>% dplyr::select(-trans),
                        table_xc_10000 %>% filter(trans==1) %>% dplyr::select(-trans,-model)))



table_t1%>%
  as_hux(add_colnames = F) %>%
  merge_cells(row=1,col=c(2:9)) %>%
  merge_cells(row=2,col=c(2:5)) %>%
  merge_cells(row=2,col=c(6:9)) %>%
  merge_cells(row=4) %>%
  merge_cells(row=10) %>%
  merge_cells(row=16) %>%
  merge_cells(row=22,col=c(2:9)) %>%
  merge_cells(row=23,col=c(2:5)) %>%
  merge_cells(row=23,col=c(6:9)) %>%
  merge_cells(row=25) %>% 
  merge_cells(row=31) %>% 
  merge_cells(row=37) %>% 
  to_latex() %>%
  write("Tables/Simulation_result_table_trans1.txt")


# TRANSITION 2

table_t2 <- rbind(
  c("firstcol",rep("PS",8)),
  c("firstcol",rep("n=3000",4),rep("n=10000",4)),
  c(lab,lab[-1]),
  cbind(table_xb_3000 %>% filter(trans==2) %>% dplyr::select(-trans),
        table_xb_10000 %>% filter(trans==2) %>% dplyr::select(-trans,-model)),
  c("firstcol",rep("IV",8)),
  c("firstcol",rep("n=3000",4),rep("n=10000",4)),
  c(lab,lab[-1]),
  cbind(table_xc_3000 %>% filter(trans==2) %>% dplyr::select(-trans),
        table_xc_10000 %>% filter(trans==2) %>% dplyr::select(-trans,-model)))



table_t2%>%
  as_hux(add_colnames = F) %>%
  merge_cells(row=1,col=c(2:9)) %>%
  merge_cells(row=2,col=c(2:5)) %>%
  merge_cells(row=2,col=c(6:9)) %>%
  merge_cells(row=4) %>%
  merge_cells(row=10) %>%
  merge_cells(row=16) %>%
  merge_cells(row=22,col=c(2:9)) %>%
  merge_cells(row=23,col=c(2:5)) %>%
  merge_cells(row=23,col=c(6:9)) %>%
  merge_cells(row=25) %>% 
  merge_cells(row=31) %>% 
  merge_cells(row=37) %>% 
  to_latex() %>%
  write("Tables/Simulation_result_table_trans2.txt")

# TRANSITION 3

table_t3 <- rbind(
  c("firstcol",rep("PS",8)),
  c("firstcol",rep("n=3000",4),rep("n=10000",4)),
  c(lab,lab[-1]),
  cbind(table_xb_3000 %>% filter(trans==3) %>% dplyr::select(-trans),
        table_xb_10000 %>% filter(trans==3) %>% dplyr::select(-trans,-model)),
  c("firstcol",rep("IV",8)),
  c("firstcol",rep("n=3000",4),rep("n=10000",4)),
  c(lab,lab[-1]),
  cbind(table_xc_3000 %>% filter(trans==3) %>% dplyr::select(-trans),
        table_xc_10000 %>% filter(trans==3) %>% dplyr::select(-trans,-model)))



table_t3%>%
  as_hux(add_colnames = F) %>%
  merge_cells(row=1,col=c(2:9)) %>%
  merge_cells(row=2,col=c(2:5)) %>%
  merge_cells(row=2,col=c(6:9)) %>%
  merge_cells(row=4) %>%
  merge_cells(row=10) %>%
  merge_cells(row=16) %>%
  merge_cells(row=22,col=c(2:9)) %>%
  merge_cells(row=23,col=c(2:5)) %>%
  merge_cells(row=23,col=c(6:9)) %>%
  merge_cells(row=25) %>% 
  merge_cells(row=31) %>% 
  merge_cells(row=37) %>% 
  to_latex() %>%
  write("Tables/Simulation_result_table_trans3.txt")




save(table_t1,table_t2,table_t3,file = c("Tables/Simulation_result_tables.RData"))

