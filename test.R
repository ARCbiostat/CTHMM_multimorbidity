load("Data Simulation/schema_xb_3000_all.RData")
data_mm <- as.data.frame(data_mm)
debug(apply_LCA)
apply_LCA(data_mm,LCA_obj)



test <- data_mm %>% drop_na()
