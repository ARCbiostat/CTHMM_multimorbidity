

apply_LCA <- function(pop, scenario_obj) {
  pop %<>% ungroup() 
  X <- pop %>%  dplyr::select(any_of(tolower(colnames(scenario_obj$pattern_obj$obj$y)))) %>%
    dplyr::mutate_all(function(x)
    x + 1) %>% as.data.frame()
  
  if(ncol(X) !=ncol(scenario_obj$pattern_obj$obj$y))stop("Some diseases in the LCA model not found in the data!")
  
  post <- poLCA::poLCA.posterior(scenario_obj$pattern_obj$obj,
                          y = X,
                          x = as.matrix(cbind(pop$age), ncol = 1))
  
  pop$MP <- as.numeric(apply(post, 1, which.max))# MODE
  
  
  # renaming MP with right order:
  #scenario_obj$right_order # for renaming MP
  n <- dim(scenario_obj$tmat)[1]
  recode_map <- setNames(1:length(scenario_obj$right_order),
                         scenario_obj$right_order)
  pop <- pop %>%
    mutate(MP = recode(MP, !!!recode_map))
  # move MP col next to MP_sim
  pop <- pop %>%
    dplyr::select(1:which(names(.) == "MP_sim"), "MP", setdiff(names(.), c("MP", "MP_sim")))
  
  return(pop)
  
}
