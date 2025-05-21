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
gt <-load(paste0("Data SNAC-K/scenario_B.RData"))
true_param_obj <- get(gt)
true_param_obj$trans_mod
ntrans <- length(true_param_obj$trans_mod)
true_param<- data.frame(trans=character(),
                        rate=numeric(),
                        shape=numeric(),
                        beta_cov1=numeric(),
                        beta_cov2= numeric(),
                        beta_cov3= numeric())
trans_label <- c("1->2","1->3","2->3")
for (i in 1:ntrans){
  curr <- true_param_obj$trans_mod[[i]]
  true_param <-rbind(true_param, data.frame(trans=trans_label[i], 
                                            rate= curr$res.t["rate", "est"],
                                            shape= curr$res.t["shape", "est"],
                                            beta_cov1= curr$res.t["cov1", "est"],
                                            beta_cov2= curr$res.t["cov2", "est"],
                                            beta_cov3= curr$res.t["cov3", "est"]
  ))
}
true_param


models<-c("TIMM","TIHMM","ApproxTIMM", "ApproxTIHMM","benchmark_model")

 
####################################################################
#################### pop-based study 3/6 years ##################### 
####################################################################

## 3000 
model_est_xb_3000 <- load_and_extract_est(results_path =paste0("Simulation/results/results_20250516_xb_3000","/"), nsim=3000,
                                           N=100, st="xb",models=models)
## 10000
model_est_xb_10000 <- load_and_extract_est(results_path =paste0("Simulation/results/results_20250517_081436_xb_10000","/"), nsim=10000,
                                           N=100, st="xb",models=models)
####################################################################
#################### irregular observations ######################## 
####################################################################

### 3000
model_est_xc_3000 <- load_and_extract_est(results_path =paste0("Simulation/results/results_20250517_141958_xc_3000","/"), nsim=10000,
                                           N=100, st="xc",models=models)
### 10000
model_est_xc_10000 <- load_and_extract_est(results_path =paste0("Simulation/results/results_20250517_190836_xc_10000","/"), nsim=10000,
                                           N=100, st="xc",models=models)
## load desired models
nsim<-10000
N <-100 # total number of datasets
st <- "xc"

# folder from where
results_path <- paste0("Simulation/results/results_20250517_190836_xc_10000","/")

load_model <- FALSE
if (load_model){
  pt <-load(paste0("results/results_estimates/model_est_",st,"_",nsim,".RData"))
  param <- get(pt)
  param_df <- param$param_df
  conv_time <- param$time_df
}else{
  param_df <-data.frame(dataset_id = integer(), 
                        trans=character(),
                        rate=numeric(),
                        shape=numeric(),
                        rate_se=numeric(),
                        shape_se=numeric(),
                        beta_cov1 = numeric(),
                        beta_cov1_se = numeric(),
                        beta_cov2 = numeric(),
                        beta_cov2_se = numeric(),
                        beta_cov3 = numeric(),
                        beta_cov3_se = numeric()
  )
  
  # extract rate, shape, their se and the transition name and store it in param for each dataset
  conv_time <- data.frame(dataset_id = integer(), 
                          model_name = character(),
                          time = numeric()
  )
  
  for(model_name in models){
    print(model_name)
    for(i in 1:N){
      model_file <- paste0(results_path, model_name, "_", nsim, "_", i, ".RData")
      model_n <- tryCatch({
        load(model_file)
      }, error = function(e) {
        message(paste("File not found or cannot load:", model_file))
        return(NULL) 
      })
      
      if (is.null(model_n)) {
        next
      }
      model_n <-load(paste0(results_path,model_name,"_",nsim,"_",i,".RData"))
      model <- get(model_n)
      conv_time <- rbind(conv_time, data.frame(dataset_id = i, model_name = model_name, time=model$time))
      #print(param_df)
      param_df <- rbind(param_df,get_estimates(model, model_name,i))
    }
  }
}
param_df
conv_time


model_est <-list(
  param_df = param_df,
  time_df = conv_time
)

if (!load_model){
  save(model_est, file =paste0(result_folder,"/model_est_",st,"_", nsim,".RData"))
}


########################################
######## Compute Performance measures ###########
########################################
# compute performance measures: bias, empirical standard error, MSE, Coverage, bias-eliminated coverage
source("Functions//performance_measures.R")
pm_df <- compute_pm(param_df,true_param) 

if (!load_model){
  save(pm_df, file =paste0(result_folder,"/pm_df_",scenario,"_", nsim,".RData"))
}