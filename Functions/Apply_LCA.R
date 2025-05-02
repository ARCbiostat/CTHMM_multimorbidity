
apply_LCA <- function(pop,scenario_obj, scenario){
  
  colnames(scenario_obj$pattern_obj$obj$y)<- tolower(colnames(scenario_obj$pattern_obj$obj$y))
  X <- pop %>%  dplyr::select(any_of(colnames(scenario_obj$pattern_obj$obj$y))) %>% dplyr::mutate_all(function(x)x+1)
  if (scenario %in% c("A", "Av2")){
    post <- poLCA.posterior(scenario_obj$pattern_obj$obj,y=X)
  } else if (scenario %in% c("B")){
    post <- poLCA.posterior(scenario_obj$pattern_obj$obj,y=X,x=as.matrix(cbind(pop$age), ncol=1))
  }
  pop$MP <-as.numeric(apply(post,1,which.max))# MODE
  
  
  # renaming MP with right order:
  #scenario_obj$right_order # for renaming MP
  n <- dim(scenario_obj$tmat)[1]
  recode_map <- setNames(1:length(scenario_obj$right_order), scenario_obj$right_order)
  pop <- pop %>%
    mutate(MP = recode(MP, !!!recode_map))
  # move MP col next to MP_sim
  pop <- pop %>%
    dplyr::select(
      1:which(names(.) == "MP_sim"), "MP", setdiff(names(.), c("MP", "MP_sim"))   
    )
  
  return(pop)
  
}



