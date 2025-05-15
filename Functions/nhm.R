library(nhm)
library(msm)
library(tictoc)



apply_nhm <- function(pop_ms,misc, sim_obj, result_folder,nsim){
    model_0 <- NULL
    model_misc <- NULL
    pop_ms$state <- pop_ms$MP
    pop_ms <- as.data.frame(pop_ms)
    q_nhm <- rbind(c(0,1,2),
                   c(0,0,3),
                   c(0,0,0))
    colnames(q_nhm)<-c(1:dim(q_nhm)[1])
    rownames(q_nhm)<-c(1:dim(q_nhm)[1])
    # same structure for q_mat and nonh
    nonh <- q_nhm
    covm <- list(
      cov1 = rbind(c(0,1,2), c(0,0,3), c(0,0,0)),
      cov2 = rbind(c(0,4,5), c(0,0,6), c(0,0,0)),
      cov3 = rbind(c(0,7,8), c(0,0,9), c(0,0,0)),
      )
    tic("nhm gompertz with cov")
    model_obj_gomp<- model.nhm(state ~ age, subject=subject_id, type='gompertz', data=pop_ms, trans= q_nhm,
                               nonh= nonh,
                               covariates= c("cov1", "cov2", "cov3"),
                               covm = covm,
                               death=T, death.states = 3)
    tryCatch({
    model_0 <- nhm(model_obj_gomp,
                   #initial= nhm_init,
                   gen_inits = TRUE, #not converging
                   control=nhm.control(splits = c(60,61,63,65,75,76,85,95,99,102,103,105,109),
                                       ncores = 4
                                       #, verbose=TRUE
                   )
    )
    }, error = function(e) {
      print(paste("Error during model fitting:", e$message))
    })
    
    t<-toc()
    
    if (!is.null(model_0)){
      model_obj_0 <- list(
        model = model_0,
        time = t$toc - t$tic
      )
      save(model_obj_0, file =paste0(result_folder,"/TIMM_scenario_",scenario,"_", nsim,"_",unique(pop_ms$dataset_id),".RData"))
    } else {
      cat("Model nhm_base is NULL for dataset_id ",unique(pop_ms$dataset_id),". Model was not saved", "\n")
     }
    
    non_diagonal <- misc
    diag(non_diagonal) <- 0
    misc_param <- non_diagonal[non_diagonal != 0]
    
    #nhm_init <- c(model_0$par, misc_param)
    nhm_init <- c(model_0$par, log(misc_param))
    
    misc_shape <- rbind(c(0,2,0),
                        c(1,0,0),
                        c(0,0,0))
    tic("nhm misc gompertz with cov")
    model_obj_gomp_misc<- model.nhm(state ~ age, subject=subject_id, type='gompertz', data=pop_ms, trans= q_nhm,
                                    nonh= nonh,
                                    emat = misc_shape,
                                    covariates= c("cov1", "cov2", "cov3"),
                                    covm = covm,
                                    #firstobs="exact",
                                    #initp = c(1,2,2,2,2,0),
                                    #initp_value = c(0.6,0.4,0),
                                    #centre_time=72,   
                                    death=T, death.states = 3)
    
    
    # misc parameters need to be passed to initial and then fixed in fixedpar
    n1<- length(model_0$par) + 1
    n2<- length(nhm_init)
    split_points <- c(60,61,63,65,75,76,85,95,99,102,103,105,106,109)
    
    
    tryCatch({
      model_misc <- nhm(
        model_obj_gomp_misc,
        initial = nhm_init,
        fixedpar = c(n1:n2),
        control = nhm.control(
          splits = split_points,
          ncores = 4, 
          rtol = 1e-8, 
          atol = 1e-8
        )
      )
    }, error = function(e) {
      print(paste("Error during model fitting:", e$message))
    })
    t<-toc()
    if (!is.null(model_misc)) {
      model_obj_m <- list(
        model = model_misc,
        time = t$toc - t$tic
      )
      save(model_obj_m, file =paste0(result_folder,"/TIHMM_scenario_", nsim,"_",unique(pop_ms$dataset_id),".RData"))
      
    } else {
      cat("Model nhm_misc is NULL for dataset_id ",unique(pop_ms$dataset_id),". Model was not saved")
    
    }
    
  print("nhm executed")
}

