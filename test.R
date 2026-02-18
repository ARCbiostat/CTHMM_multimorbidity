load("Data Simulation/schema_xb_3000_all.RData")
data_mm <- as.data.frame(data_mm)
debug(apply_LCA)
test <- apply_LCA(data_mm,LCA_obj)
table(test$MP_sim,test$MP)

