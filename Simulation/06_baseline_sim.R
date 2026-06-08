library(parallel)
source("~/Documents/CTHMM_multimorbidity/Functions/run_sim_baseline.R")

load(paste0("Data Simulation/schema_xb_3000_all.RData"))
datasets <- split(data_mm, data_mm$dataset_id)
rm(data_mm)


res <-mclapply(1:100, function(s) run_sim_baseline2(datasets[[s]],0.05)
,  mc.cores = detectCores()-10)


res_table <- do.call("rbind",res)
table(res_table[,1])
table(res_table[,2])

res2 <-mclapply(1:100, function(s) run_sim_baseline2(datasets[[s]],0.08,classes = 2)
               ,  mc.cores = detectCores()-10)


res_table2 <- do.call("rbind",res2)
table(res_table2[,1])
summary(res_table2[,2])

load(paste0("Data Simulation/schema_xb_10000_all.RData"))
datasets <- split(data_mm, data_mm$dataset_id)
rm(data_mm)


res3 <-mclapply(1:100, function(s) run_sim_baseline2(datasets[[s]],0.05,classes=2)
               ,  mc.cores = detectCores()-10)

res_table3 <- do.call("rbind",res3)
table(res_table3[,1])
summary(res_table3[,2])
