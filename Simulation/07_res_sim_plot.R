library(dplyr)
library(tidyverse)
library(huxtable)
library(ggplot2)
library(ggprism)

source("Functions/aux_tables.R")
lab <- c("Model","Estimate","Bias","S.E. Bias","Coverage")

load("results/results_estimates/pm_df_est_xb_3000.RData")
load("results/results_estimates/pm_df_est_xb_10000.RData")

pm_df_est_xb_3000 %<>%  filter(parameter%in%c("beta_cov1","beta_cov2","beta_cov3"))  
pm_df_est_xb_10000 %<>%  filter(parameter%in%c("beta_cov1","beta_cov2","beta_cov3"))  

ggplot(subset(pm_df_est_xb_3000,model!="benchmark_model"))+
  geom_errorbar(aes(x=model,y=bias,ymin=bias-2*se_bias,ymax=bias+2*se_bias))+
  geom_point(aes(x=model,y=bias))+
  geom_hline(aes(yintercept=0),col="grey60",linetype="dashed")+
  theme_prism()+
  facet_grid(parameter~trans,scales="free")+
  scale_x_discrete("")+
  scale_y_continuous("Bias")+
  theme(axis.text.x  = element_text(angle=90))



ggplot(subset(pm_df_est_xb_10000,model!="benchmark_model"))+
  geom_errorbar(aes(x=model,y=bias,ymin=bias-2*se_bias,ymax=bias+2*se_bias))+
  geom_point(aes(x=model,y=bias))+
  geom_hline(aes(yintercept=0),col="grey60",linetype="dashed")+
  theme_prism()+
  facet_grid(parameter~trans,scales="free")+
  scale_x_discrete("")+
  scale_y_continuous("Bias")+
  theme(axis.text.x  = element_text(angle=90))
