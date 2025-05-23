library(magrittr)
library(tidyverse)
library(gtsummary)
library(ggplot2)
library(ggprism)
library(ggpubr)
library(poLCA)
library(sBIC)
library(dplyr)
library(future.apply)
library(parallel)
library(msm)
library(nhm)
library(flexsurv)
library(flexsurvcure)
library(ggplot2)
library(dplyr)
library(survminer)
library(tidyverse)


###########################################################################################
#################### Extract estimates from fitted model #################
###########################################################################################

source("Functions//aux_extract_estimates.R")

# folder where estimates will be saved
result_folder <-  "results/results_estimates"

if (!dir.exists(result_folder)) {
  dir.create(result_folder)
  message("Folder 'results_estimates' created.")
} else {
  message("Folder 'result_estimate' already exists.")
}

## load true parameters
load("Data Simulation/LCA object.RData")
true_param_obj <- LCA_obj
true_param_obj$trans_mod
ntrans <- length(true_param_obj$trans_mod)
true_param <- data.frame(
  trans = character(),
  rate = numeric(),
  shape = numeric(),
  beta_cov1 = numeric(),
  beta_cov2 = numeric(),
  beta_cov3 = numeric()
)
trans_label <- c("1->2", "1->3", "2->3")
for (i in 1:ntrans) {
  curr <- true_param_obj$trans_mod[[i]]
  true_param <- rbind(
    true_param,
    data.frame(
      trans = trans_label[i],
      rate = curr$res.t["rate", "est"],
      shape = curr$res.t["shape", "est"],
      beta_cov1 = curr$res.t["dm_sexFemales", "est"],
      beta_cov2 = curr$res.t["educ_el", "est"],
      beta_cov3 = 0
    )
  )
}
true_param


models <- c("TIMM", "TIHMM", "ApproxTIMM", "ApproxTIHMM")


####################################################################
#################### benchmark                 #####################
####################################################################
#3000
model_est_bench_3000 <- load_and_extract_est(
  results_path = paste0("results/results_20250522_101420_benchmark_3000", "/"),
  nsim = 3000,
  N = 100,
  st = "benchmark",
  models = "benchmark_model"
)

model_est_bench_10000 <- load_and_extract_est(
  results_path = paste0("results/results_20250522_101850_benchmark_10000", "/"),
  nsim = 10000,
  N = 100,
  st = "benchmark",
  models = "benchmark_model"
)

####################################################################
#################### pop-based study 3/6 years #####################
####################################################################

## 3000
model_est_xb_3000 <- load_and_extract_est(
  results_path = paste0("Simulation/results/results_20250516_xb_3000", "/"),
  nsim = 3000,
  N = 100,
  st = "xb",
  models = models
)
## 10000
model_est_xb_10000 <- load_and_extract_est(
  results_path = paste0("Simulation/results/results_20250517_081436_xb_10000", "/"),
  nsim = 10000,
  N = 100,
  st = "xb",
  models = models
)
####################################################################
#################### irregular observations ########################
####################################################################

### 3000
model_est_xc_3000 <- load_and_extract_est(
  results_path = paste0("Simulation/results/results_20250517_141958_xc_3000", "/"),
  nsim = 3000,
  N = 100,
  st = "xc",
  models = models
)
### 10000
model_est_xc_10000 <- load_and_extract_est(
  results_path = paste0("Simulation/results/results_20250517_190836_xc_10000", "/"),
  nsim = 10000,
  N = 100,
  st = "xc",
  models = models
)


########################################
######## Compute Performance measures ###########
########################################
# compute performance measures: bias, empirical standard error, MSE, Coverage, bias-eliminated coverage
source("Functions//performance_measures.R")

model_est_bench_3000 <- get(
  load(
    "~/Documents/CTHMM_multimorbidity/results/results_estimates/model_est_benchmark_3000.RData"
  )
)
model_est_bench_10000 <- get(
  load(
    "~/Documents/CTHMM_multimorbidity/results/results_estimates/model_est_benchmark_10000.RData"
  )
)
model_est_xb_3000 <- get(
  load(
    "~/Documents/CTHMM_multimorbidity/results/results_estimates/model_est_xb_3000.RData"
  )
)
model_est_xb_10000 <- get(
  load(
    "~/Documents/CTHMM_multimorbidity/results/results_estimates/model_est_xb_10000.RData"
  )
)
model_est_xc_3000 <- get(
  load(
    "~/Documents/CTHMM_multimorbidity/results/results_estimates/model_est_xc_3000.RData"
  )
)
model_est_xc_10000 <- get(
  load(
    "~/Documents/CTHMM_multimorbidity/results/results_estimates/model_est_xc_10000.RData"
  )
)

################## bechmark ##############
# 3000
pm_df_est_bench_3000 <- compute_pm(model_est_bench_3000$param_df, true_param)
# 10000
pm_df_est_bench_10000 <- compute_pm(model_est_bench_10000$param_df, true_param)
####################################################################
#################### pop-based study 3/6 years #####################
####################################################################

## 3000
pm_df_est_xb_3000 <- compute_pm(model_est_xb_3000$param_df, true_param)
pm_df_est_xb_3000 <-  pm_df_est_xb_3000 %>% filter(model != "benchmark_model")
pm_df_est_xb_3000 <- rbind(pm_df_est_xb_3000, pm_df_est_bench_3000)
save(pm_df_est_xb_3000,
     file = paste0(result_folder, "/pm_df_est_xb_3000.RData"))

## 10000
pm_df_est_xb_10000 <- compute_pm(model_est_xb_10000$param_df, true_param)
pm_df_est_xb_10000 <-  pm_df_est_xb_10000 %>% filter(model != "benchmark_model")
pm_df_est_xb_10000 <- rbind(pm_df_est_xb_10000, pm_df_est_bench_10000)
save(pm_df_est_xb_10000,
     file = paste0(result_folder, "/pm_df_est_xb_10000.RData"))
####################################################################
#################### irregular observations ########################
####################################################################

## 3000
pm_df_est_xc_3000 <- compute_pm(model_est_xc_3000$param_df, true_param)
pm_df_est_xc_3000 <-  pm_df_est_xc_3000 %>% filter(model != "benchmark_model")
pm_df_est_xc_3000 <- rbind(pm_df_est_xc_3000, pm_df_est_bench_3000)
save(pm_df_est_xc_3000,
     file = paste0(result_folder, "/pm_df_est_xc_3000.RData"))

## 10000
pm_df_est_xc_10000 <- compute_pm(model_est_xc_10000$param_df, true_param)
pm_df_est_xc_10000 <-  pm_df_est_xc_10000 %>% filter(model != "benchmark_model")
pm_df_est_xc_10000 <- rbind(pm_df_est_xc_10000, pm_df_est_bench_10000)
save(pm_df_est_xc_10000,
     file = paste0(result_folder, "/pm_df_est_xc_10000.RData"))
