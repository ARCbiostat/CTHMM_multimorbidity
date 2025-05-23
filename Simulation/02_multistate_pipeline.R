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

#functions LCA
source("Functions//misc_matrix.R")
source("Functions//Apply_LCA.R")


#functions for multistate
source("Functions//msm.R")
source("Functions//nhm.R")
source("Functions//flexsurv.R")

# auxiliary functions to run analysis on simulated data
source("Functions/aux_sim_run.R")

# load LCA object 
load("Data Simulation/LCA object.RData")

####################################################################
#################### BENCHMARK                ##################### 
####################################################################

#3000
run_analysis_benchmark(nsim = 3000,LCA_obj = LCA_obj)

#10000
run_analysis_benchmark(nsim = 10000,LCA_obj = LCA_obj)
####################################################################
#################### pop-based study 3/6 years ##################### 
####################################################################

# 3000
run_analysis(study_type = "xb",nsim = 3000,LCA_obj = LCA_obj)
# 10 000
run_analysis(study_type = "xb",nsim = 10000,LCA_obj = LCA_obj)

####################################################################
#################### irregular observations ######################## 
####################################################################

# 3000
run_analysis(study_type = "xc",nsim = 3000,LCA_obj = LCA_obj)
# 10 000
run_analysis(study_type = "xc",nsim = 10000,LCA_obj = LCA_obj, avoid= c(1:22, 31:50, 71:100))





