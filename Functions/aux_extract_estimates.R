library(msm)
library(flexsurv)
library(nhm)

##
### calculate rate and shape from msm:

hazards_mine <- function(x, b.covariates, no.years, trans = NULL,SE=FALSE, CI = FALSE,
                         age.shift = 0) {

  model <- x  # Fitted MSM model

  # Handle transition matrix (convert vector to matrix if needed)
  if (is.vector(trans)) {
    trans <- matrix(trans, 1, 2)
  }

  # Create the transition matrix if not provided
  if (!is.null(trans)) {
    ntrans <- nrow(trans)  # Number of transitions
  } else {
    Q <- model$qmodel$imatrix  # Transition rate matrix
    ntrans <- sum(Q)  # Count the number of transitions
    trans <- matrix(NA, ntrans, 2)  # Initialize matrix for transitions
    index <- 1
    for (i in 1:nrow(Q)) {
      for (j in 1:nrow(Q)) {
        if (Q[i, j]) {
          trans[index, ] <- c(i, j)  # Add transitions to matrix
          index <- index + 1
        }
      }
    }
  }

  # Number of states in the model
  nstates <- nrow(model$Qmatrices$baseline)

  # Create a null transition matrix with zeros on the diagonal
  Q.null <- matrix(as.numeric(model$Qmatrices$baseline != 0), nstates, nstates)
  diag(Q.null) <- 0
  ntrans <- sum(Q.null)  # Count the number of transitions

  # Number of covariates (including age and other covariates)
  ncovs <- 1 + length(b.covariates)
  nbeta <- max(which(names(model$estimates) == "qcov"))

  if (nbeta != (ntrans * ncovs)) {
    stop("\nAll covariates in msm model should be specified.\n\n")
  }

  # Extract the age variable from the covariates
  b.age <- b.covariates$age

  # Create a sequence of ages for computation
  age.grid <- seq(b.age, b.age + no.years, by = 1/12)

  # Initialize an empty list to store the hazards
  hazards <- list()

  # Loop over each transition to calculate the hazards
  for (i in 1:ntrans) {
    haz <- rep(NA, length(age.grid))
    hazLB <- rep(NA, length(age.grid))
    hazUB <- rep(NA, length(age.grid))
    hazSE <- rep(NA, length(age.grid))

    # For each age, calculate the hazard for the current transition
    for (j in 1:length(age.grid)) {
      if (ncovs == 2) {
        covariates <- list(age = age.grid[j])
      } else {
        covariates <- c(age = age.grid[j], b.covariates[2:length(b.covariates)])
      }

      # Compute the hazard for the current transition (without confidence intervals)
      haz[j] <- qmatrix.msm(model, covariates = covariates, ci = "none")[trans[i, 1], trans[i, 2]]

      if (CI || SE) {
        # Compute confidence intervals if requested
        Q.CI <- qmatrix.msm(model, covariates = covariates, ci = "delta", cores=4)
        hazLB[j] <- Q.CI$L[trans[i, 1], trans[i, 2]]
        hazUB[j] <- Q.CI$U[trans[i, 1], trans[i, 2]]
        hazSE[j] <- Q.CI$SE[trans[i, 1], trans[i, 2]]
      }
    }

    # Store the hazards (and optionally the confidence intervals)
    hazard_data <- list(hazard = haz)
    if (CI) {
      hazard_data$hazardLB <- hazLB
      hazard_data$hazardUB <- hazUB
    }
    if (SE){
      hazard_data$hazardSE <-hazSE
    }

    # Name the transition for easier identification later
    hazard_name <- paste("Transition", trans[i, 1], "to", trans[i, 2])
    hazards[[hazard_name]] <- hazard_data
  }

  # Return the computed hazards (and confidence intervals if requested)
  return(hazards)
}
# 

get_estimates <- function(model, model_name, i,j=0){
  # 
  rate<- numeric()
  rate_se <-numeric()
  shape <- numeric()
  shape_se <- numeric()
  trans <-numeric()
  beta_cov1 <- numeric()
  beta_cov2 <-numeric()
  beta_cov3 <-numeric()
  beta_cov1_se <- numeric()
  beta_cov2_se <- numeric()
  beta_cov3_se <- numeric()

  if(model_name %in% c("benchmark_model"))
  {   
    endings <- names(model$transmodel)
    for (ending in endings) {
      current_model <- model$transmodel[[ending]]
      trans<- c(trans, ending)
      rate <-c(rate, current_model$res.t["rate", "est"])
      rate_se <- c(rate_se, current_model$res.t["rate", "se"])
      shape <- c(shape, current_model$res.t["shape", "est"])
      shape_se <- c(shape_se, current_model$res.t["shape", "se"])
      beta_cov1 <- c(beta_cov1, current_model$res.t["cov1", "est"])
      beta_cov1_se <- c(beta_cov1_se, current_model$res.t["cov1", "se"])
      beta_cov2 <- c(beta_cov2, current_model$res.t["cov2", "est"])
      beta_cov2_se <- c(beta_cov2_se, current_model$res.t["cov2", "se"])
      beta_cov3 <- c(beta_cov3, current_model$res.t["cov3", "est"])
      beta_cov3_se <- c(beta_cov3_se, current_model$res.t["cov3", "se"])
    }

  } else if (model_name %in% c("ApproxTIMM", "ApproxTIHMM")){
    # simil gompertz
    
    haz <- hazards_mine(model$model,SE=TRUE, b.covariates = list(age = 60, cov1 = 0, cov2 = 0, cov3=0), no.years = 40)
    #haz <- hazards_mine(model$model, b.covariates = list(age = 0, educ_el = 0, dm_sex = 0), no.years = 40)
    
    # Assuming this hazards come from fitting a gompertz model I wanna retrieve for each transition
    # shape and rate value, parameters of the distribution
    # h(t)=rate*exp(shape*t)
    # log(h(t))= log(rate) + shape*t
    # can be seen as y(t)= a+b*t
    
    # >0 element in Qmatrices
    ntrans <-  sum(model$model$qmodel$imatrix)
    trans <- rename_trans(names(haz))
    min_age <- min(model$model$data[[1]]$age)
    max_age <- max(model$model$data[[1]]$age)
    
    # SE for rate and shape are computed with bootstrapping
    n_boot <- 1000
    for (i in 1:ntrans){
      age_grid <- seq(min_age, max_age , length.out = length(unlist(haz[[1]]$hazard)))
      y <- log(as.numeric(unlist(haz[[i]]$hazard)))
      reg_model <- lm(y ~ age_grid)
      rate <- c(rate,reg_model$coefficients[1])
      shape<- c(shape, reg_model$coefficients[2])
      haz_se <- unlist(haz[[i]]$hazardSE)
      
      bootstrap_results <- replicate(n_boot, {
        y_noisy <- log(as.numeric(unlist(haz[[i]]$hazard))) + rnorm(length(y), 0, haz_se)
        indices <- sample(1:length(y), replace = TRUE)
        y_boot <- y_noisy[indices]
        age_boot <- age_grid[indices]
      
        boot_model <- lm(y_boot ~ age_boot)
        boot_model$coefficients
      })
      
      rate_se <- c(rate_se, sd(bootstrap_results[1, ]))
      shape_se <- c(shape_se, sd(bootstrap_results[2, ]))

      #rate_se <- c(rate_se, )
      #shape_se <- rep(0,length(trans))
    }
    std_errors <- sqrt(diag(model$model$covmat))
    beta_cov1 <- model$model$estimates[((ntrans*2)+1):(3*ntrans)]
    beta_cov1_se <- std_errors[((ntrans*2)+1):(3*ntrans)]
    beta_cov2 <- model$model$estimates[((ntrans*3)+1):(4*ntrans)]
    beta_cov2_se <- std_errors[((ntrans*3)+1):(4*ntrans)]
    beta_cov3 <- model$model$estimates[((ntrans*4)+1):(5*ntrans)]
    beta_cov3_se <- std_errors[((ntrans*4)+1):(5*ntrans)]
   

  } else if (model_name %in% c("TIMM", "TIHMM")){
    np <- model$model$nstate
    trans <- model$model$parnames[1:np]
    cov_matrix <- solve(model$model$hess)
    std_errors <- sqrt(diag(cov_matrix))
    trans <- gsub("^Base: ", "", trans)
    rate <- model$model$par[1:np]
    rate_se <- std_errors[1:np]
    shape <- model$model$par[(np+1):(2*np)]
    shape_se <- std_errors[(np+1):(2*np)]
    beta_cov1 <- model$model$par[((np*2)+1):(3*np)]
    beta_cov1_se <- std_errors[((np*2)+1):(3*np)]
    beta_cov2 <- model$model$par[((np*3)+1):(4*np)]
    beta_cov2_se <- std_errors[((np*3)+1):(4*np)]
    beta_cov3 <- model$model$par[((np*4)+1):(5*np)]
    beta_cov3_se <- std_errors[((np*4)+1):(5*np)]
  }
  if (model_name=="flexsurv_imp"){
    est_obj <- data.frame(
      model = rep(model_name, length(trans)),
      dataset_id = rep(i,length(trans)),
      imp_id = rep(j,length(trans)),
      trans = trans,
      rate = rate,
      rate_se = rate_se,
      shape = shape,
      shape_se = shape_se,
      beta_cov1 = beta_cov1,
      beta_cov1_se = beta_cov1_se,
      beta_cov2 = beta_cov2,
      beta_cov2_se = beta_cov2_se,
      beta_cov3 = beta_cov3,
      beta_cov3_se = beta_cov3_se
    )
  }else {
    est_obj <- data.frame(
    model = rep(model_name, length(trans)),
    dataset_id = rep(i,length(trans)),
    trans = trans,
    rate = rate,
    rate_se = rate_se,
    shape = shape,
    shape_se = shape_se,
    beta_cov1 = beta_cov1,
    beta_cov1_se = beta_cov1_se,
    beta_cov2 = beta_cov2,
    beta_cov2_se = beta_cov2_se,
    beta_cov3 = beta_cov3,
    beta_cov3_se = beta_cov3_se
  )
  }
  
  return(est_obj)
}


rename_trans<- function(names){
  
  new_names <- gsub("Transition (\\d+) to (\\d+)", "\\1->\\2", names)
  return(new_names)
}


load_and_extract_est <- function(results_path, nsim, N, st, models, load_model=FALSE){
  if (load_model){
    pt <-load(paste0("results/results_estimates/model_est_",st,"_",nsim,".RData"))
    param <- get(pt)
    param_df <- param$param_df
    conv_time <- param$time_df
  }else{
    param_df <-data.frame(dataset_id = integer(), 
                          trans=character(),
                          rate=numeric(),
                          shape=numeric(),
                          rate_se=numeric(),
                          shape_se=numeric(),
                          beta_cov1 = numeric(),
                          beta_cov1_se = numeric(),
                          beta_cov2 = numeric(),
                          beta_cov2_se = numeric(),
                          beta_cov3 = numeric(),
                          beta_cov3_se = numeric()
    )
    
    # extract rate, shape, their se and the transition name and store it in param for each dataset
    conv_time <- data.frame(dataset_id = integer(), 
                            model_name = character(),
                            time = numeric()
    )
    
    for(model_name in models){
      print(model_name)
      for(i in 1:N){
        model_file <- paste0(results_path, model_name, "_", nsim, "_", i, ".RData")
        model_n <- tryCatch({
          load(model_file)
        }, error = function(e) {
          message(paste("File not found or cannot load:", model_file))
          return(NULL) 
        })
        
        if (is.null(model_n)) {
          next
        }
        model_n <-load(paste0(results_path,model_name,"_",nsim,"_",i,".RData"))
        model <- get(model_n)
        conv_time <- rbind(conv_time, data.frame(dataset_id = i, model_name = model_name, time=model$time))
        #print(param_df)
        param_df <- rbind(param_df,get_estimates(model, model_name,i))
        
        model_est <-list(
          param_df = param_df,
          time_df = conv_time
        )
      }
    }
  }
  
  if (!load_model){
    save(model_est, file =paste0(result_folder,"/model_est_",st,"_", nsim,".RData"))
  }
}
