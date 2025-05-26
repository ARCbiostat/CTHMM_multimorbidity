fix_table <- function(pm_data){
  table <- split(pm_data, pm_data$trans)
  table <- lapply(table, function(x)
    x %>% dplyr::select(parameter, model, estimate, bias,se_bias, coverage))
  
  table <- lapply(table, function(x)
    x %>% mutate_at(3:5, round, 3) %>% mutate_at(6, round, 2)  %>% 
      mutate_all(as.character) %>% 
      filter(parameter%in%c("beta_cov1","beta_cov2","beta_cov3")))
  
  table <- do.call("rbind",lapply(1:3,function(y){
    res <- split(table[[y]],table[[y]]$parameter)
    tab <- do.call("rbind",lapply(1:3, function(x){
      res[[x]] %<>% dplyr::select(-parameter)
      res2 <- rbind(rep(paste0("**beta_",x,"**"),5),res[[x]] %>% slice(c(5,2,1,4,3)))
      return(res2)
    }))
    tab <- cbind(tab,trans=rep(y,nrow(tab)))
  }
  
  ))
  
}
