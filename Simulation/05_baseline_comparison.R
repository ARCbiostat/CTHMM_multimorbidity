library(tidyverse)

##################################################################
################# Evaluate baseline estimates ##################### 
source("Functions//performance_measures.R")
t_vals <- seq(60, 90, length.out = 100) 


####################################################################
#################### pop-based study 3/6 years ##################### 
####################################################################
pm_df_3000<- get(load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xb_3000.RData"))
pm_df_10000<- get(load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xb_10000.RData"))

base_meas_3000 <- eval_cum_hazard(pm_df_3000,true_param,t_vals, nsim= 3000)
base_meas_10000 <- eval_cum_hazard(pm_df_10000,true_param,t_vals, nsim= 10000)

base_meas_3000[, sapply(base_meas_3000, is.numeric)] <- 
  lapply(base_meas_3000[, sapply(base_meas_3000, is.numeric)], function(x) {
    ifelse(abs(x) < 0.001, formatC(x, format = "e", digits = 3), round(x, 4))
  })
base_meas_10000[, sapply(base_meas_10000, is.numeric)] <- 
  lapply(base_meas_10000[, sapply(base_meas_10000, is.numeric)], function(x) {
    ifelse(abs(x) < 0.001, formatC(x, format = "e", digits = 3), round(x, 4))
  })


library(xtable)
latex_table <- xtable(base_meas_3000)
print(latex_table, type = "latex", include.rownames = FALSE)
latex_table <- xtable(base_meas_10000)
print(latex_table, type = "latex", include.rownames = FALSE)


####################################################################
#################### Irregular visits ##################### 
####################################################################
pm_df_3000<- get(load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xc_3000.RData"))
pm_df_10000<- get(load("~/Documents/CTHMM_multimorbidity/results/results_estimates/pm_df_est_xc_10000.RData"))

base_meas_3000 <- eval_cum_hazard(pm_df_3000,true_param,t_vals, nsim= 3000)
base_meas_10000 <- eval_cum_hazard(pm_df_10000,true_param,t_vals, nsim= 10000)

base_meas_3000[, sapply(base_meas_3000, is.numeric)] <- 
  lapply(base_meas_3000[, sapply(base_meas_3000, is.numeric)], function(x) {
    ifelse(abs(x) < 0.001, formatC(x, format = "e", digits = 3), round(x, 4))
  })
base_meas_10000[, sapply(base_meas_10000, is.numeric)] <- 
  lapply(base_meas_10000[, sapply(base_meas_10000, is.numeric)], function(x) {
    ifelse(abs(x) < 0.001, formatC(x, format = "e", digits = 3), round(x, 4))
  })


library(xtable)
latex_table <- xtable(base_meas_3000)
print(latex_table, type = "latex", include.rownames = FALSE)
latex_table <- xtable(base_meas_10000)
print(latex_table, type = "latex", include.rownames = FALSE)