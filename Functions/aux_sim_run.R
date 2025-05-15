# auxiliary functions:
create_unique_folder <- function(prefix = "results",study_type) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  folder_name <- paste0(prefix, "_", timestamp,"_",study_type)
  return(folder_name)
}
run_multistate <- function(pop,result_folder,LCA_obj,nsim){
  
  pop_ms_death_1 <- pop_ms %>%
    dplyr::filter(dht == 1) %>%
    group_by(dataset_id, subject_id) %>%
    group_split() %>%
    map_df(~ bind_rows(
      .x,
      slice_tail(.x, n = 1) %>% mutate(MP = dim(LCA_obj$tmat)[1],MP_sim=dim(LCA_obj$tmat)[1], age = age_exit)
    )) %>%
    ungroup()
  
  pop_ms_death_0 <- pop_ms %>%
    dplyr::filter(dht == 0)
  
  pop_ms <- bind_rows(pop_ms_death_0, pop_ms_death_1)
  
  pop_ms <- pop_ms %>%
    group_by(dataset_id,subject_id) %>%
    filter(n() > 1) %>%
    ungroup()
  
  pop_ms <- pop_ms %>%
    group_by(subject_id) %>%
    mutate(flag = cumsum(MP == 2),
           # Adjust MP based on the flag
           MP = case_when(
             flag > 0 & MP == 1 ~ 2, # If flag is active and MP is 1, change it to 2
             TRUE ~ MP           # Otherwise, keep MP as it is
           )
    ) %>%
    ungroup() %>%
    dplyr::select(-flag) 
  
  misc <- calculate_m_matrix(pop,1)
  apply_flexsurv_base(pop_ms,misc, LCA_obj, result_folder,nsim)
  apply_msm(pop_ms,misc, LCA_obj, result_folder,nsim)
  gc()
  apply_nhm(pop_ms,misc, LCA_obj,result_folder,nsim)
  gc()
  
}
run_multistate_wrapper <- function(dataset,result_folder,LCA_obj,nsim) {
  # assign LCA class based on LCA model
  dataset<- apply_LCA(dataset,LCA_obj)
  run_multistate(dataset, function(x)result_folder(x,result_folder,LCA_obj,nsim))
}

run_analysis <- function(study_type,nsim,LCA_obj){
  #load data
  load(paste0("Data Simulation/schema_",study_type,"_",nsim,"_all.RData"))
  
  ########### split by dataset_id #########
  # run in parallel, group pop by dataset_id
  cores <- min(detectCores() -5, 20) #modify number of cores
  plan(multisession, workers=cores)  # Or `plan(multiprocess)` for cross-platform compatibility
  
  datasets <- split(data_mm, data_mm$dataset_id)
  
  # create result folder using timestamp
  result_folder <-  file.path("results",create_unique_folder(paste0("results"),study_type))
  dir.create(result_folder, recursive = TRUE)
  
  # call function to apply nhm, msm, and flexsurv in parallel
  tic("Multistate in parallel")
  options(future.globals.maxSize=1*1e9)
  lapply(datasets, function(x)run_multistate_wrapper(x,result_folder,LCA_obj,nsim))
  #future_lapply(datasets, FUN = function(x)run_multistate_wrapper(x,result_folder,LCA_obj,nsim))
  toc()
  
  
}