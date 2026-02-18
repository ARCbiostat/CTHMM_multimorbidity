run_sim_baseline2 <- function(db,threshold=0.05){
  
  
  require(magrittr)
  require(tidyverse)
  require(msm)
  require(MMLCA)
  death_rows <- db %>%
    filter(dht == 1) %>%
    group_by(subject_id) %>%
    slice(1) %>%
    mutate(
      visit_number = 99,
      age = age_exit
    ) %>%
    dplyr:: select(colnames(db))
  
  db_clean <- db %>% #set death at 0 in visit's rows
    mutate(dht = 0)
  
  db <- bind_rows(db_clean, death_rows) %>% # row-bind datasets with visit rows and death rows
    arrange(subject_id, visit_number)
  
  # -----> Keep only those with at least 2 visits or 1 visit + dead 
  db <- db %>%
    group_by(subject_id) %>%
    filter(n() >= 2) %>%
    ungroup() 
  
  # -------> Add dis_ prefix to easily create the disease matrix
  
  names(db)[13:71] <- paste0("dis_", names(db)[13:71]) # Add dis_ prefix
  
  # -----> prepare matrix: one with only baseline diseases (Xb) and the other with all diseases (X)
  
  Xb <- prepare_data(db %>% filter(visit_number == 1), dis_string = "dis_", keepmm = T) #only baseline
  X <- prepare_data(db, dis_string = "dis_", keepmm = T)
  
  # ------> Set the threshold to include only those diseases with at least 2% prevalence
  

  disease_names <- select_conditions(Xb,
                                     threshold = 0.02)
  
  # -------> Run the LCA with different number of classes and compare the metrics:
  
  set.seed(1234)
  res <- select_number_LCA(
    nclasses = 2:5,
    X = Xb,
    conditions = disease_names,
    nrep = 15
  )
  
  
  n_class <-as.numeric(res$metrics[which.min(as.numeric(res$metrics[,8])),2])
  print(res$metrics)
  if(n_class!=2) return(c(0,NA))
  
  mm_pattern <- assign_LCA(res$obj[[1]],X)
  table(mm_pattern)
  db$mm_pattern <- mm_pattern
  
  # -------> Add death status
  db %<>% mutate(mm_pattern_death = case_when(dht==1~3,
                                              TRUE~mm_pattern)) 
  
  q_matrix <- rbind(c(0,1,1), # from state 1 to 2 and 3
                    c(1,0,1), # from state 2 to 1 and 3
                    c(0,0,0)) # from state 3 to nothing (absorbing state!)
 
  
  # Step 2: compute empirical transition frequencies
  emp_f <- statetable.msm(mm_pattern_death, subject=subject_id, data= db)/apply(statetable.msm(mm_pattern_death, subject = subject_id, data = db),1,sum)

  
  # Step 3: copy the original qmatrix

  
  # Step 4: Apply threshold 
  mask <- emp_f < threshold
  
  # Step 5: Apply threshold mask to q0.0 
  # NOTE: emp_f does not contain the absorbing state (death), therefore mask needs to be applied only to rows 1 and 2
  q_matrix[1:2, 1:3][mask] <- 0
  q_matrix
  ntrans <- which(as.numeric(q_matrix)!=0)
  gc()
  if(ntrans!=3)return(c(1,0))
  else return(c(1,1))
}