library(dplyr)
library(tidyverse)
library(huxtable)
library(ggplot2)

source("Functions/aux_tables.R")
lab <- c("Model","Estimate","Bias","S.E. Bias","Coverage")

load("results/results_estimates/pm_df_est_xb_3000.RData")
load("results/results_estimates/pm_df_est_xb_10000.RData")

table_xb_3000 <-fix_table(pm_df_est_xb_3000)
table_xb_10000 <- fix_table(pm_df_est_xb_10000)

gg