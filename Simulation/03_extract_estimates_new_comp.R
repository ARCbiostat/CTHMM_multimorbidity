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
source("Functions/performance_measures.R")

# folder where estimates will be saved
result_folder <-  "results/results_estimates_comp"

if (!dir.exists(result_folder)) {
  dir.create(result_folder)
  message("Folder 'results_estimates_comp' created.")
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
#################### pop-based study 3/6 years #####################
####################################################################

## 3000
model_est_xb_3000 <- load_and_extract_est(
  results_path = paste0("Simulation/results/results_20250902_130945_xb_under_3000", "/"),
  nsim = 3000,
  N = 100,
  st = "xb_under",
  models = models
)

####################################################################
#################### irregular observations ########################
####################################################################

### 3000
# model_est_xc_3000 <- load_and_extract_est(
#   results_path = paste0("Simulation/results/results_20250517_141958_xc_3000", "/"),
#   nsim = 3000,
#   N = 100,
#   st = "xc",
#   models = models
# )

########################################
######## Compute Performance measures ###########
########################################
# compute performance measures: bias, empirical standard error, MSE, Coverage, bias-eliminated coverage
source("Functions//performance_measures.R")


model_est_xb_3000_under <- get(
  load(
    "~/Documents/CTHMM_multimorbidity/results/results_estimates_comp/model_est_xb_under_3000.RData"
  )
)

# model_est_xc_3000 <- get(
#   load(
#     "~/Documents/CTHMM_multimorbidity/results/results_estimates_comp/model_est_xc_3000.RData"
#   )
# )


## 3000
pm_df_est_xb_3000_under <- compute_pm(model_est_xb_3000_under$param_df, true_param)
pm_df_est_xb_3000_under <-  pm_df_est_xb_3000_under %>% filter(model != "benchmark_model")
save(pm_df_est_xb_3000_under,
     file = paste0(result_folder, "/pm_df_est_xb_3000_under.RData"))


pm_df_est_xb_3000_under$type <- "under"
load("~/Documents/CTHMM_multimorbidity/results/results_estimates/model_est_xb_3000.RData")
pm_df_est_xb_3000$type <- "normal"

pm_df_3000 <- rbind(pm_df_est_xb_3000,pm_df_est_xb_3000_under)
pm_df_3000_comp <- pm_df_3000 %>% group_by(model,trans,parameter) %>% 
  filter(parameter%in%c("beta_cov1","beta_cov2","beta_cov3")) %>% 
  arrange(type) %>% 
  transmute(diff_bias=bias-lag(bias),
            diff_coverage=coverage-lag(coverage)) %>% 
    drop_na()


####################################################################
#################### irregular observations ########################
####################################################################

## 3000
pm_df_est_xc_3000 <- compute_pm(model_est_xc_3000$param_df, true_param)
pm_df_est_xc_3000 <-  pm_df_est_xc_3000 %>% filter(model != "benchmark_model")
pm_df_est_xc_3000 <- rbind(pm_df_est_xc_3000, pm_df_est_bench_3000)
save(pm_df_est_xc_3000,
     file = paste0(result_folder, "/pm_df_est_xc_3000.RData"))


