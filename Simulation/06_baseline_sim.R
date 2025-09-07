library(parallel)
source("~/Documents/CTHMM_multimorbidity/Functions/run_sim_baseline.R")

load(paste0("Data Simulation/schema_xb_3000_all.RData"))
datasets <- split(data_mm, data_mm$dataset_id)
rm(data_mm)


res <-mclapply(1:100, function(s) run_sim_baseline2(datasets[[s]],0.05),  mc.cores = detectCores()-10)


res_table <- do.call("rbind",res)
table(res_table[,1])
table(res_table[,2])
