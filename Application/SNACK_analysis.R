library(magrittr)
library(tidyverse)
library(gtsummary)
library(ggplot2)
library(ggprism)
library(ggpubr)
library(poLCA)
library(sBIC)
library(dplyr)
library(future.apply)
library(tidyr)
library(msm)
library(flexsurv)
library(nhm)
library(tictoc)
#install.packages("psych")  # If not installed
#install.packages("polycor")
library(psych)
library(polycor)

source("Functions//misc_matrix.R")
source("Functions//helper_f_snack.R")
########################################
#### import and preprocess the data ####
########################################
scenario <- "B"
load(paste0("Data SNAC-K/Data MM pattern LCA.RData"))
load("Data SNAC-K/coviates SNac K.RData") #sex and educ_level
load("Data SNAC-K/Socioecon_behaviors_SNACK.RData") 
sim_obj_name <- load(paste0("Data SNAC-K/scenario_B.RData"))
sim_obj <- get(sim_obj_name)
load("Data SNAC-K/01_LCA.RData")
pattern_obj<- res2_covs
covs$lopnr <- as.numeric(as.character(covs$lopnr))
covs$dm_sex <-ifelse(as.character(covs$dm_sex)=="Females",1,0)
cov_tot <- covs %>% left_join(covs_sel)
cov_tot <- cov_tot %>%
  mutate(sei_long_cat=haven::as_factor(sei_long_cat), finstrain = haven::as_factor(finstrain))

dat <- dat %>% left_join(cov_tot)
length(unique(dat$lopnr))

snack_2 <- dat %>% drop_na(Age) %>% ungroup() # 8 subjects removed


# modify from 5 to 2 states
snack_death <- snack_2 %>% filter(MP=="Death")
tapply(snack_death$Age, snack_death$dm_sex, summary) # summarize death statistics
snack_death$MP <- 3
snack_death$MP_prob <-NA
snack_2 <- snack_2 %>%
  filter(!if_any(any_of(colnames(pattern_obj$obj$y)), is.na)) # ci sono molte ultime osservazioni con NA (dropout?)


X <- snack_2 %>% ungroup() %>% dplyr::select(any_of(colnames(pattern_obj$obj$y))) %>% mutate_all(function(x)x+1)


post <- poLCA.posterior(pattern_obj$obj,y=X, x=as.matrix(cbind(snack_2$Age), ncol=1))

#length(which(is.na(post[,1])))
snack_2$MP <-as.numeric(apply(post,1,which.max))
snack_2$MP_prob <- as.numeric(apply(post, 1, max))
length(unique(snack_2$lopnr))

snack_2 <- bind_rows(snack_2,snack_death)
# eliminate those with no follow ups
obs_count <- table(snack_2$lopnr)
single_obs_subjects <- names(obs_count[obs_count == 1])

snack_2_filtered <- snack_2[!snack_2$lopnr %in% single_obs_subjects,]
# modify transitions from 2->1 
snack_2_filtered <- snack_2_filtered%>% arrange(lopnr) 

snack_2_filtered <- snack_2_filtered%>% group_by(lopnr) %>%
  mutate(flag = cumsum(MP == 2),
       # Adjust MP based on the flag
       MP = case_when(
         flag > 0 & MP == 1 ~ 2, # If flag is active and MP is 1, change it to 2
         TRUE ~ MP           # Otherwise, keep MP as it is
       )
  ) %>%
  ungroup() %>%
  dplyr::select(-flag) 
length(unique(snack_2_filtered$lopnr))


# folder where estimates will be saved
result_folder <-"Application/results"

if (!dir.exists(result_folder)) {
  dir.create(result_folder)
  message("Folder 'results_estimates' created.")
} else {
  message("Folder 'result_estimate' already exists.")
}

#######################
### Baseline Stat.#####
#######################
# 
snack_base <- snack_2_filtered %>%
  group_by(lopnr) %>%
  slice_head()

tapply(snack_base$Age, snack_base$dm_sex, summary) 
table(snack_base$MP)

snack_nhm <- snack_2_filtered
snack_nhm$state <- snack_nhm$MP

# factor covs are not supported by nhm -> convert to binary dummy variables
snack_nhm <- snack_nhm%>% 
  mutate(finstrain_dummy = ifelse(finstrain == levels(finstrain)[2], 1, 0),
         sei_long_cat_dummy = ifelse(sei_long_cat == levels(sei_long_cat)[2], 1, 0),
         if_ever_smoke= ifelse(R_SMOKE %in% c("former", "current"), 1,0 ),
         heavy_alcool= ifelse(ALCO_CONSUMP %in% c("heavy", "light/mod"), 1,0 )
  )

snack_nhm <- as.data.frame(snack_nhm)
statetable.msm(state, subject=lopnr, data= snack_nhm)

# mean number of visits
snack_3 <- snack_2_filtered %>%
  group_by(lopnr) %>%
  summarise(n_visits = n()) %>%
  ungroup()

summary(snack_3$n_visits)
###################
### Covariates ###
##################
hetcor(cov_tot[2:11])


####################################
######## FITTING THE MODELS ########


statetable.msm(MP, subject=lopnr, data= snack_2_filtered)
q_matrix <- rbind(c(0,1,1),
                  c(0,0,1),
                  c(0,0,0))
#inizialization
q0 <- crudeinits.msm(MP ~ Age, subject=lopnr, data=snack_2_filtered, qmatrix= q_matrix)
q0

# compute misc matrix:
source("Functions//misc_matrix.R")

X2 <- snack_base %>%ungroup() %>% dplyr::select(any_of(colnames(sim_obj$pattern_obj$obj$y))) %>% mutate_all(function(x)x+1)
misc <- get_internal_validation_matrix(sim_obj$pattern_obj$obj, X2, covs =as.matrix(cbind(snack_base$Age), ncol=1) )
n <- dim(sim_obj$tmat)[1]
misc <- add_death(misc, n-1)
rownames(misc) <- rownames(sim_obj$tmat)
colnames(misc) <- rownames(sim_obj$tmat)
misc

#################################
# covariates only on trans 1->2 #
#################################
#### ApproxTIMM: 9 covariates ####
model_age_cov5 <- msm(MP ~ Age, subject=lopnr, data=snack_nhm, qmatrix= q0, deathexact = 3,
                      covariates = list("1-2" = ~ Age + educ_el + dm_sex + no_pa + life_alone + heavy_alcool+ if_ever_smoke+ fin_strain_early+finstrain_dummy+sei_long_cat_dummy, "1-3" = ~ Age, "2-3" = ~ Age), 
                      method = "BFGS", control = list(fnscale = 14000, maxit= 10000))
hazard.msm(model_age_cov5)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save(model_age_cov5, file =paste0(result_folder,"/ApproxTIMM_9cov_",timestamp,".RData"))

#### ApproxTIHMM: 9 covariates ####
qmatrix_init <- model_age_cov5$Qmatrices$baseline 
covinits_init <-model_age_cov5$estimates
n1<- sum(q0 > 0)
n2 <- length(model_age_cov5$estimates)
covinits_init <-as.list(covinits_init[(n1+1):n2])
#n1<- sum(q0 > 0) + sum(covinits_init!=0)
covinits_init <- covinits_init[covinits_init != 0]

n_fixed <- sum(misc > 0 & row(misc) != col(misc))
n1 <- sum(q0 > 0) +length(covinits_init)
n2<- n1+ n_fixed


# filter out NAs
snack_nhm_filtered <- snack_nhm %>% drop_na(fin_strain_early,finstrain_dummy,sei_long_cat_dummy)
model_5.2_misc <- msm(MP ~ Age, subject=lopnr, data=snack_nhm_filtered, qmatrix= qmatrix_init, 
                      covinits = covinits_init, 
                      deathexact = 3,
                      initprobs = c(0.8,0.2,0),
                      covariates = list("1-2" = ~ Age + educ_el + dm_sex + no_pa + life_alone + heavy_alcool+ if_ever_smoke+ fin_strain_early+finstrain_dummy+sei_long_cat_dummy, "1-3" = ~ Age, "2-3" = ~ Age), 
                      ematrix=misc, 
                      #fixedpars= c((n1+1):n2),
                      fixedpars= c(11:12),
                      control = list(fnscale = 20000, maxit= 10000))
hazard.msm(model_5.2_misc) 
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save(model_5.2_misc, file =paste0(result_folder,"/ApproxTIHMM_9cov_",timestamp,".RData"))


##################################
###  ApproxTIMM and ApproxTIHMM: 6 covariates #################
# removing covariates with too many NAs
q_matrix <- rbind(c(0,1,1),
                  c(0,0,1),
                  c(0,0,0))
#inizialization
q0 <- crudeinits.msm(MP ~ Age, subject=lopnr, data=snack_nhm, qmatrix= q_matrix)
q0
model_age_cov6 <- msm(MP ~ Age, subject=lopnr, data=snack_nhm, qmatrix= q0, deathexact = 3,
                      covariates = list("1-2" = ~ Age + educ_el + dm_sex + no_pa + life_alone + heavy_alcool+ if_ever_smoke, "1-3" = ~ Age, "2-3" = ~ Age), 
                      method = "BFGS", control = list(fnscale = 14000, maxit= 10000))
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save(model_age_cov6, file =paste0(result_folder,"/ApproxTIMM_6cov_",timestamp,".RData"))
qmatrix_init <- model_age_cov6$Qmatrices$baseline 
covinits_init <-model_age_cov6$estimates
n1<- sum(q0 > 0)
n2 <- length(model_age_cov6$estimates)
covinits_init <-as.list(covinits_init[(n1+1):n2])
covinits_init <- covinits_init[covinits_init != 0]
n1 <- sum(q0 > 0) +length(covinits_init)
n2<- n1+ n_fixed

model_6_misc <- msm(MP ~ Age, subject=lopnr, data=snack_nhm, qmatrix= qmatrix_init, 
                      covinits = covinits_init, # passare lista valori age a convergenza
                      deathexact = 3,
                      initprobs = c(0.8,0.2,0),
                      covariates = list("1-2" = ~ Age + educ_el + dm_sex + no_pa + life_alone + heavy_alcool+ if_ever_smoke, "1-3" = ~ Age, "2-3" = ~ Age), 
                      ematrix=misc, 
                      fixedpars= c((n1+1):n2), 
                      control = list(fnscale = 20000, maxit= 10000))
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save(model_6_misc, file =paste0(result_folder,"/ApproxTIHMM_6cov_",timestamp,".RData"))


#####################################################
###########TIMM and TIHMM: 9 covariates################

q_nhm <- rbind(c(0,1,2),
               c(0,0,3),
               c(0,0,0))
colnames(q_nhm)<-c(1:dim(q_nhm)[1])
rownames(q_nhm)<-c(1:dim(q_nhm)[1])
# same structure for q_mat and nonh
nonh <- q_nhm

covm1 <- list(
  educ_el = rbind(c(0,1,0), c(0,0,0), c(0,0,0)),
  dm_sex = rbind(c(0,2,0), c(0,0,0), c(0,0,0)),
  no_pa = rbind(c(0,3,0), c(0,0,0), c(0,0,0)),
  life_alone = rbind(c(0,4,0), c(0,0,0), c(0,0,0)),
  heavy_alcool = rbind(c(0,5,0), c(0,0,0), c(0,0,0)),
  if_ever_smoke = rbind(c(0,6,0), c(0,0,0), c(0,0,0)),
  fin_strain_early = rbind(c(0,7,0), c(0,0,0), c(0,0,0)),
  finstrain_dummy = rbind(c(0,8,0), c(0,0,0), c(0,0,0)),
  sei_long_cat_dummy = rbind(c(0,9,0), c(0,0,0), c(0,0,0))
)

tic("nhm gompertz with 9 cov")
model_obj_gomp2<- model.nhm(state ~ Age, subject=lopnr, type='gompertz', data=snack_nhm, trans= q_nhm,
                            nonh= nonh,
                            covariates= c("educ_el", "dm_sex","no_pa" , "life_alone",  "heavy_alcool","if_ever_smoke", "fin_strain_early","finstrain_dummy","sei_long_cat_dummy"),
                            covm = covm1,
                            death=T, death.states = 3)

model_2 <- nhm(model_obj_gomp2,
               #initial= nhm_init,
               gen_inits = TRUE, #not converging
               control=nhm.control(splits = c(60,61,63,65,75,76,85,95,99,102,103,105,109),
                                   ncores = detectCores()-2
                                   #, verbose=TRUE
               )
)
toc() # nhm gompertz with 9 cov: 12586.949 sec elapsed
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save(model_2, file =paste0(result_folder,"/TIMM_9cov_",timestamp,".RData"))

misc_shape <- rbind(c(0,2,0),
                    c(1,0,0),
                    c(0,0,0))
tic()
model_obj_gomp_misc<- model.nhm(state ~ Age, subject=lopnr, type='gompertz', data=snack_nhm, trans= q_nhm,
                                nonh= nonh,
                                emat = misc_shape,
                                covariates= c("educ_el", "dm_sex","no_pa" , "life_alone",  "heavy_alcool","if_ever_smoke", "fin_strain_early","finstrain_dummy","sei_long_cat_dummy"),
                                covm = covm1,
                                death=T, death.states = 3)


# misc parameters need to be passed to initial and then fixed in fixedpar
non_diagonal <- misc
diag(non_diagonal) <- 0
misc_param <- non_diagonal[non_diagonal != 0]
nhm_init <- c(model_2$par, log(misc_param))
n1<- length(model_2$par) + 1
n2<- length(nhm_init)
split_points <- c(60,61,63,65,75,76,85,95,99,102,103,105,106,109)

model_misc <- nhm(
  model_obj_gomp_misc,
  initial = nhm_init,
  fixedpar = c(n1:n2),
  control = nhm.control(
    splits = split_points,
    ncores = detectCores()-2, 
    rtol = 1e-6, 
    atol = 1e-6
  )
)
toc() # 39420.092

tic("nhm misc preinizializzato")
nhm_init <- c(model_misc$par, log(misc_param))
model_misc2 <- nhm(
  model_obj_gomp_misc,
  initial = nhm_init,
  fixedpar = c(n1:n2),
  control = nhm.control(
    splits = split_points,
    ncores = detectCores()-2, 
    rtol = 1e-6, 
    atol = 1e-6
  )
)
toc()

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save(model_misc2, file =paste0(result_folder,"/THIMM_9cov_",timestamp,".RData"))


####################################################
########## TIMM and TIHMM: 6 covariates ############
covm1 <- list(
  educ_el = rbind(c(0,1,0), c(0,0,0), c(0,0,0)),
  dm_sex = rbind(c(0,2,0), c(0,0,0), c(0,0,0)),
  no_pa = rbind(c(0,3,0), c(0,0,0), c(0,0,0)),
  life_alone = rbind(c(0,4,0), c(0,0,0), c(0,0,0)),
  if_ever_smoke = rbind(c(0,5,0), c(0,0,0), c(0,0,0)),
  heavy_alcool = rbind(c(0,6,0), c(0,0,0), c(0,0,0))
)

tic("nhm gompertz with cov")
model_obj_gomp1<- model.nhm(state ~ Age, subject=lopnr, type='gompertz', data=snack_nhm, trans= q_nhm,
                            nonh= nonh,
                            covariates= c("educ_el", "dm_sex","no_pa" , "life_alone", "if_ever_smoke", "heavy_alcool"),
                            covm = covm1,
                            death=T, death.states = 3)

model_1 <- nhm(model_obj_gomp1,
               #initial= nhm_init,
               gen_inits = TRUE, #not converging
               control=nhm.control(splits = c(60,61,63,65,75,76,85,95,99,102,103,105,109),
                                   ncores = 4
                                   #, verbose=TRUE
               )
)
toc()

nhm_init <- c(model_1$par)
model_1.1 <- nhm(model_obj_gomp1,
                 initial= nhm_init,
                 #gen_inits = TRUE, #not converging
                 control=nhm.control(splits = c(60,61,63,65,75,76,85,95,99,102,103,105,109),
                                     ncores = 4
                                     #, verbose=TRUE
                 )
)
#

# misc model
non_diagonal <- misc
diag(non_diagonal) <- 0
misc_param <- non_diagonal[non_diagonal != 0]
nhm_init <- c(model_1.1$par, log(misc_param))

misc_shape <- rbind(c(0,2,0),
                    c(1,0,0),
                    c(0,0,0))
covm2
tic("nhm misc gompertz with cov")
model_obj_gomp_misc<- model.nhm(tate ~ Age, subject=lopnr, type='gompertz', data=snack_nhm, trans= q_nhm,
                                nonh= nonh,
                                emat = misc_shape,
                                covariates= c("educ_el", "dm_sex","no_pa" , "life_alone", "if_ever_smoke", "heavy_alcool"),
                                covm = covm1,
                                death=T, death.states = 3)


# misc parameters need to be passed to initial and then fixed in fixedpar
n1<- length(model_0$par) + 1
n2<- length(nhm_init)
split_points <- c(60,61,63,65,75,76,85,95,99,102,103,105,106,109)

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


########################################
############ baseline hazard plots ##########
# obtain approx Gompertz param:


bp_reduced <- data.frame(model= character(),
                         trans=character(),
                         rate=numeric(),
                         shape=numeric(),
                         rate_se=numeric(),
                         shape_se=numeric()
)
tic("Approx gompertz TIMM")
model <- model_age_cov6
bp_reduced <- rbind(bp_reduced,baseline_param2(model, "ApproxTIMM"))
toc()
baseline_param2(model, "ApproxTIMM")
tic("Approx gompertz TIHMM")
model <- model_6_misc
bp_reduced <- rbind(bp_reduced,baseline_param2(model, "ApproxTIHMM"))
toc()

plot_hazard_snack(bp_reduced,t_vals)
plot_hazard_snack_ci(bp_reduced,t_vals)


haz <- hazards_snack(model,SE=TRUE, b.covariates = list(Age = 60, dm_sex = 0,no_pa=0,life_alone=0), no.years = 30)
ntrans <-  sum(model$qmodel$imatrix)
trans <- rename_trans(names(haz))
min_age <- min(model$data[[1]]$Age)
max_age <- max(model$data[[1]]$Age)
age_grid <- seq(min_age, max_age , length.out = length(unlist(haz[[1]]$hazard)))
plot(age_grid, haz[[2]]$hazard, type = "l", main = paste("Transition 1 ->3"), xlab = "Age", ylab = "Hazard")


base_param <- data.frame(model= character(),
                         trans=character(),
                         rate=numeric(),
                         shape=numeric(),
                         rate_se=numeric(),
                         shape_se=numeric()
)
hazards_snack <- function(x, b.covariates, no.years, trans = NULL,SE=FALSE, CI = FALSE,
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
  b.age <- b.covariates$Age
  
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
        covariates <- list(Age = age.grid[j])
      } else {
        covariates <- c(Age = age.grid[j], b.covariates[2:length(b.covariates)])
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

tic("Approx gompertz TIMM")
model <- model_age_cov5
base_param <- rbind(base_param,baseline_param(model, "ApproxTIMM"))
toc()

tic("Approx gompertz TIHMM")
model <- model_5.2_misc
base_param <- rbind(base_param,baseline_param(model, "ApproxTIHMM"))
toc()

baseline_nhm <- function(model,model_name){
  np <- model$nstate
  trans <- model$parnames[1:np]
  cov_matrix <- solve(model$hess)
  std_errors <- sqrt(diag(cov_matrix))
  trans <- gsub("^Base: ", "", trans)
  rate <- model$par[1:np]
  rate_se <- std_errors[1:np]
  shape <- model$par[(np+1):(2*np)]
  shape_se <- std_errors[(np+1):(2*np)]
  base_est <- data.frame(
    model = rep(model_name, length(trans)),
    trans = trans,
    rate = rate,
    rate_se = rate_se,
    shape = shape,
    shape_se = shape_se
  )
  return(base_est)
}
model <- model_2
base_param <- rbind(base_param,baseline_nhm(model, "TIMM"))

# singular hessian -> some se are NaN
model <- model_misc2
base_param <- rbind(base_param,baseline_nhm(model, "TIHMM"))
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save(base_param, file =paste0("results/snack_baseline_",timestamp,".RData"))


# plots


t_vals <- seq(60, 90, length.out = 100)
plot_cum_hazard_snack(base_param,t_vals)
plot_hazard_snack(base_param,t_vals)
#plot_hazard_snack_ci(base_param,t_vals)



########################
### Goodness of fit ####
########################
# all covariates
AIC(model_age_cov5)
AIC(model_5.2_misc)
AIC(model_2)
+2*model_2$value + 2*model_2$npar
AIC(model_misc2) # -2 logLik + 2*n_param
AIC_TIHMM <- 2*model_misc2$value + 2*(model_misc2$npar)
AIC_TIHMM

#summary(model_age_cov5)
# loglikelihood
logLik(model_age_cov5) 
logLik(model_5.2_misc) 
-model_2$value
-model_misc2$value

# BIC :  -2 logLik + log(n)*n_param. n = number of observations
n <- dim(snack_nhm_filtered)[1]
-2*logLik(model_age_cov5) + log(n)*length(model_age_cov5$estimates[model_age_cov5$estimates!=0])
-2*logLik(model_5.2_misc) + log(n)*length(model_5.2_misc$estimates[model_5.2_misc$estimates!=0])
+2*model_2$value + log(n)*model_2$npar
+2*model_misc2$value + log(n)*model_misc2$npar

#
library(msm)
pearson.msm(model_age_cov5)
pearson.msm(model_5.2_misc)

# only covariates without NAs
AIC(model_age_cov6)
AIC(model_6_misc)
AIC(model_3)
+2*model_3$value + 2*model_3$npar
AIC(model_misc3) # -2 logLik + 2*n_param
AIC_TIHMM <- 2*model_misc3$value + 2*(model_misc3$npar)
AIC_TIHMM


# loglikelihood
logLik(model_age_cov6) 
logLik(model_6_misc) 
-model_3$value
-model_misc3$value

# BIC :  -2 logLik + log(n)*n_param. n = number of observations
n <- dim(snack_nhm)[1]
-2*logLik(model_age_cov6) + log(n)*length(model_age_cov6$estimates[model_age_cov6$estimates!=0])
-2*logLik(model_6_misc) + log(n)*length(model_6_misc$estimates[model_6_misc$estimates!=0])
+2*model_3$value + log(n)*model_3$npar
+2*model_misc3$value + log(n)*model_misc3$npar


######### Plots ###########
#### TIMM plots ###
plot(model_2, time0=60)
plot(model_2, time0=60, times= t_vals)
plot(model_2,what="intensities" , time0=60)

# female vs male
plot(model_2, time0=60, covvalue = c(0,1,0,0,0,0,1,0,0))
plot(model_2, time0=60, covvalue = c(0,0,0,0,0,0,1,0,0))
#predict(model_2, time0=60, times= t_vals)

#### TIHMM plots ###
plot.nhm.mine(model_misc2, time0=60, xlab="Age")
plot.nhm.mine(model_misc3, time0=60)
plot.nhm.mine(model_misc2, what= "intensities",time0=60, times= t_vals, main_arg= "Hazard function", xlab="Age")
plot.nhm.mine(model_misc3, what= "intensities",time0=60, times= t_vals, main_arg= "Hazard function")
plot.nhm.mine(model_2, time0=60, covvalue = c(0,1,0,0,0,0,1,0,0))
plot.nhm.mine(model_2, time0=60, covvalue = c(0,0,0,0,0,0,1,0,0))


# non-param estimates for hazard 1->3
data <- msm2Surv(data = snack_nhm, subject = "lopnr", time= "Age", state="MP", Q = q_matrix)
data_2 <- data %>% filter(trans==2)
model_flex <- flexsurvreg(formula = Surv(Tstart,Tstop, status) ~ educ_el + dm_sex+no_pa + life_alone + heavy_alcool + if_ever_smoke +fin_strain_early +finstrain_dummy+ sei_long_cat_dummy, subset=(trans==2),
                       data = data_2, dist = "gompertz")
plot(model_flex,type= "cumhaz", est= TRUE, ci=FALSE, t=t_vals, xlim= c(60,90))
plot_cumhazard_with_ci <- function(x, mp, col.obs = "black", xlim, lty.obs = 1, lwd.obs = 1, ylim = NULL, t = NULL, start = 0, est = TRUE, ci = NULL, B = 1000, cl = 0.95, col = "red", 
                                   lty = 1, lwd = 2, col.ci = "black", lty.ci = 2, lwd.ci = 0.5, 
                                   add = FALSE) {
  dat <- x$data
  mf <- model.frame(x)
  mm <- as.data.frame(model.matrix(x))
  if (isTRUE(attr(dat$Y, "type") == "interval")) {
    form <- "Surv(dat$Y[,\"time1\"],dat$Y[,\"time2\"],type=\"interval2\") ~ "
  } else {
    form <- "Surv(dat$Y[,\"start\"],dat$Y[,\"stop\"],dat$Y[,\"status\"]) ~ "
  }
  form <- paste(form, if (x$ncovs > 0 && all(x$covdata$isfac)) 
    paste("mm[,", 1:x$ncoveffs, "]", collapse = " + ")
    else 1)
  form <- as.formula(form)
  
  #plot(survfit(form, data = mm, weights = dat$m$`(weights)`), fun = "cumhaz", col = col.obs, lty = lty.obs, lwd = lwd.obs, ylim = ylim, xlim=xlim)
  surv_fit <- survfit(form, data = mm, weights = dat$m$`(weights)`)
  cumhaz_data <- data.frame(time = surv_fit$time, cumhaz = -log(surv_fit$surv),
                            lcl = -log(surv_fit$upper),
                            ucl = -log(surv_fit$lower))
  model_ids <- unique(mp$model)
  model_hazards <- list()
  
  for (j in seq_along(model_ids)) {
    model_data <- mp %>% filter(model == model_ids[j])
    shape_param <- model_data$shape
    rate_param <- model_data$rate
    
    if (length(shape_param) == 1 && length(rate_param) == 1) {
      model_hazards[[j]] <- tibble(
        t = t_vals,
        hazard = Hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE),
        model = paste0("Model ", model_ids[j])
      )
    }
  }
  
  p <- ggplot(cumhaz_data, aes(x = time, y = cumhaz)) +
    geom_line(color = col.obs, linetype = lty.obs, size = lwd.obs) +
    geom_line(aes(y = lcl), color = col.ci, linetype = lty.ci, size = lwd.ci) +
    geom_line(aes(y = ucl), color = col.ci, linetype = lty.ci, size = lwd.ci) +
    labs(x = "Age", y = "Cumulative Hazard") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(color = "lightgray", fill = NA, linewidth = 0.5),
                  plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                  legend.position = "bottom" )
  
  if (!is.null(xlim)) {
    p <- p + xlim(xlim)
  }
  
  print(p)
}
plot_cumhazard_with_ci <- function(x, mp, col.obs = "black", xlim = NULL, lty.obs = 1, lwd.obs = 1, ylim = NULL, t = NULL, start = 0, est = TRUE, ci = NULL, B = 1000, cl = 0.95, col = "red", 
                                   lty = 1, lwd = 2, col.ci = "black", lty.ci = 2, lwd.ci = 0.5, 
                                   add = FALSE) {
  dat <- x$data
  mf <- model.frame(x)
  mm <- as.data.frame(model.matrix(x))
  
  form <- as.formula(paste0("Surv(dat$Y[,\"start\"], dat$Y[,\"stop\"], dat$Y[,\"status\"]) ~ ",
                            if (x$ncovs > 0 && all(x$covdata$isfac)) paste("mm[,", 1:x$ncoveffs, "]", collapse = " + ") else 1))
  
  surv_fit <- survfit(form, data = mm, weights = dat$m$`(weights)`)
  cumhaz_data <- data.frame(time = surv_fit$time, cumhaz = -log(surv_fit$surv),
                            lcl = -log(surv_fit$upper),
                            ucl = -log(surv_fit$lower))
  if (!is.null(xlim)) {
    cumhaz_data <- cumhaz_data %>% filter(time >= xlim[1])
  }
  
  # Build cumulative hazard curves for mp models
  model_ids <- unique(mp$model)
  model_hazards <- list()
  t_vals <- seq(min(cumhaz_data$time), max(cumhaz_data$time), length.out = 500)
  
  for (j in seq_along(model_ids)) {
    model_data <- mp %>% filter(model == model_ids[j])
    shape_param <- model_data$shape
    rate_param <- model_data$rate
    
    if (length(shape_param) == 1 && length(rate_param) == 1) {
      model_hazards[[j]] <- tibble(
        time = t_vals,
        cumhaz = Hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE),
        model = paste0("Model ", model_ids[j])
      )
    }
  }
  
  model_hazards_df <- bind_rows(model_hazards)
  
  p <- ggplot(cumhaz_data, aes(x = time, y = cumhaz)) +
    geom_line(aes(color = "Nonparametric"), linetype = lty.obs, linewidth = lwd.obs) +
    geom_line(aes(y = lcl, color = "LCL"), linetype = lty.ci, linewidth = lwd.ci) +
    geom_line(aes(y = ucl, color = "UCL"), linetype = lty.ci, linewidth = lwd.ci) +
    geom_line(data = model_hazards_df, aes(x = time, y = cumhaz, color = model), linetype = "dashed", linewidth = lwd) +
    scale_color_manual(
      values = c("Nonparametric" = col.obs, 
                 "LCL" = col.ci, 
                 "UCL" = col.ci, 
                 setNames(rainbow(length(model_ids)), paste0("Model ", model_ids)))
    ) +
    labs(x = "Time", y = "Cumulative Hazard", color = "Models", title = "Cumulative Hazard transition 1->3") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "lightgray", fill = NA, linewidth = 0.5),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = "bottom")
  
  if (!is.null(xlim)) {
    p <- p + xlim(xlim)
  }
  
  if (!is.null(ylim)) {
    p <- p + ylim(ylim)
  }
  
  print(p)
}

base_param_2 <- base_param %>% filter(trans=="1->3")
plot_cumhazard_with_ci(model_flex, base_param_2, t=t_vals, xlim= c(60,90))


library(muhaz)

data_3 <-  snack_nhm %>%
  arrange(lopnr, Age) %>%
  group_by(lopnr) %>%
  mutate(
    next_state = lead(state),
    event_1_to_3 = ifelse(state == 1 & next_state == 3, 1, 0),
    remove = ifelse(state==2 | (state==2 & next_state == 3),1,0)
  ) %>%
  ungroup()

data_3 <- data_3 %>%
  group_by(lopnr) %>%
  filter(!any(remove == 1)) %>%
  ungroup()

survival_data <- data_3 %>%
  group_by(lopnr) %>%
  summarize(
    start_time = min(Age),
    event_time_idx = which(event_1_to_3 == 1)[1],  # Index of the event if exists
    stop_time = if (!is.na(event_time_idx)) Age[event_time_idx + 1] else max(Age),
    event = ifelse(!is.na(event_time_idx), 1, 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-event_time_idx)


# Estimate non-parametric hazard using muhaz
haz_est <- muhaz(
  times = survival_data$stop_time,
  delta = survival_data$event,
  min.time = min(survival_data$start_time, na.rm = TRUE),
  #max.time = max(survival_data$stop_time, na.rm=TRUE)
  #min.time=60,
  max.time=90
)

# Plot the estimated hazard function
plot(haz_est, main = "Non parametric hazard function: 1->3", xlab = "Age", ylab = "Hazard")

########
plot_hazard_snack_np<-function(pm_df,  t_vals, np_est, nsim=NULL, ST=NULL) {
  
  plot_list <- list()  # Store individual plots
  
  for (i in unique(pm_df$trans)) {
    mp <- pm_df %>% filter(trans == i)
    
    model_ids <- unique(mp$model)
    
    if (!is.null(ST)){
      mp <- mp %>%
        filter(model=="ETES" |
                 (ST == "PS" & grepl("^PS_", model)) |
                 (ST == "IV" & grepl("^IV_", model)) |
                 !(ST %in% c("PS", "IV")) # Keep all models for other study types
        )
    }
    
    model_ids <- unique(mp$model)
    model_hazards <- list()
    
    for (j in seq_along(model_ids)) {
      model_data <- mp %>% filter(model == model_ids[j])
      shape_param <- model_data$shape
      rate_param <- model_data$rate
      
      if (length(shape_param) == 1 && length(rate_param) == 1) {
        model_hazards[[j]] <- tibble(
          t = t_vals,
          hazard = hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE),
          model = paste0("Model ", model_ids[j])
        )
      }
    }
    
    if (i=="1->3"){
      np_hazard <- tibble(
        t= t_vals,
        hazard = np_est$haz.est[1:length(t_vals)],
        model = "Nonparametric Estimate"
      )
      
      hazard_data <- bind_rows(model_hazards,np_hazard)
      
    }
    else{
      hazard_data <- bind_rows(model_hazards)
    }
    model_colors <- c("Nonparametric Estimate" = "black", setNames(rainbow(length(model_ids)), paste0("Model ", model_ids)))
    model_linetypes <- c("Nonparametric Estimate" = "solid", setNames(rep("dashed", length(model_ids)), paste0("Model ", model_ids)))
    
    p <- ggplot(hazard_data, aes(x = t, y = hazard, color = model, linetype = model)) +
      geom_line(data = hazard_data %>% filter(model == "Nonparametric Estimate"), 
                aes(color = model, linetype = model), linewidth = 1.0) +
      # Plot other models on top
      geom_line(data = hazard_data %>% filter(model != "Nonparametric Estimate"), 
                aes(color = model, linetype = model), linewidth = 1) +
      scale_color_manual(values = model_colors) +
      scale_linetype_manual(values = model_linetypes) +
      labs(
        title = paste0("Gompertz Hazard Function - Transition ", i),
        x = "Age",
        y = "Hazard Function",
        color = "Model",
        linetype = "Model"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 12,  hjust = 0.5)) # Set legend to bottom for easy extraction
    
    plot_list[[as.character(i)]] <- p  # Store plot in the list
  }
  
  # Extract legend from one plot
  legend <- get_legend(plot_list[[1]] + theme(legend.position = "right"))
  
  # Remove legends from individual plots
  plots_no_legend <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
  
  combined_plot <- plot_grid(
    plots_no_legend[[1]], plots_no_legend[[2]], 
    plots_no_legend[[3]], legend, 
    nrow = 2, ncol = 2, 
    rel_widths = c(1, 1), rel_heights = c(1, 1)
  )
  final_plot <- ggdraw() +
    draw_label(paste("Gompertz Hazard Function Comparison"), fontface = "bold", size = 16, x = 0.5, y = 0.95) + 
    draw_plot(combined_plot, x = 0, y = 0, width = 1, height = 0.9)
  
  print(final_plot)
  
}

plot_hazard_snack_np<-function(pm_df,  t_vals, np_est, nsim=NULL, ST=NULL) {
  
  plot_list <- list()  # Store individual plots
  
  for (i in unique(pm_df$trans)) {
    mp <- pm_df %>% filter(trans == i)
    
    model_ids <- unique(mp$model)
    
    if (!is.null(ST)){
      mp <- mp %>%
        filter(model=="ETES" |
                 (ST == "PS" & grepl("^PS_", model)) |
                 (ST == "IV" & grepl("^IV_", model)) |
                 !(ST %in% c("PS", "IV")) # Keep all models for other study types
        )
    }
    
    model_ids <- unique(mp$model)
    model_hazards <- list()
    
    for (j in seq_along(model_ids)) {
      model_data <- mp %>% filter(model == model_ids[j])
      shape_param <- model_data$shape
      rate_param <- model_data$rate
      
      if (length(shape_param) == 1 && length(rate_param) == 1) {
        model_hazards[[j]] <- tibble(
          t = t_vals,
          hazard = hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE),
          model = paste0("Model ", model_ids[j])
        )
      }
    }
    
    if (i=="1->3"){
      np_hazard <- tibble(
      t= t_vals,
      hazard = np_est$haz.est[1:length(t_vals)],
      model = "Nonparametric Estimate"
    )
    
    hazard_data <- bind_rows(model_hazards,np_hazard)
      
    }
    else{
      hazard_data <- bind_rows(model_hazards)
    }
    model_colors <- c("Nonparametric Estimate" = "black", setNames(rainbow(length(model_ids)), paste0("Model ", model_ids)))
    model_linetypes <- c("Nonparametric Estimate" = "solid", setNames(rep("dashed", length(model_ids)), paste0("Model ", model_ids)))
    
    p <- ggplot(hazard_data, aes(x = t, y = hazard, color = model, linetype = model)) +
      geom_line(data = hazard_data %>% filter(model == "Nonparametric Estimate"), 
                aes(color = model, linetype = model), linewidth = 1.0) +
      # Plot other models on top
      geom_line(data = hazard_data %>% filter(model != "Nonparametric Estimate"), 
                aes(color = model, linetype = model), linewidth = 1) +
      scale_color_manual(values = model_colors) +
      scale_linetype_manual(values = model_linetypes) +
      labs(
        title = paste0("Gompertz Hazard Function - Transition ", i),
        x = "Age",
        y = "Hazard Function",
        color = "Model",
        linetype = "Model"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 12,  hjust = 0.5)) # Set legend to bottom for easy extraction
    
    plot_list[[as.character(i)]] <- p  # Store plot in the list
  }
  
  # Extract legend from one plot
  legend <- get_legend(plot_list[[1]] + theme(legend.position = "right"))
  
  # Remove legends from individual plots
  plots_no_legend <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
  
  combined_plot <- plot_grid(
    plots_no_legend[[1]], plots_no_legend[[2]], 
    plots_no_legend[[3]], legend, 
    nrow = 2, ncol = 2, 
    rel_widths = c(1, 1), rel_heights = c(1, 1)
  )
  final_plot <- ggdraw() +
    draw_label(paste("Gompertz Hazard Function Comparison"), fontface = "bold", size = 16, x = 0.5, y = 0.95) + 
    draw_plot(combined_plot, x = 0, y = 0, width = 1, height = 0.9)
  
  print(final_plot)
  
}
plot_hazard_snack_np(base_param,t_vals, haz_est)

###### HR plots ########
library(ggplot2)
library(ggprism)

# helper functions in helper_f_snack.R
HR_est_9cov <- data.frame()
#HR_est_9cov <- extract_hr_estimates(model_2, "TIMM", "nhm")
HR_est_9cov <- rbind(HR_est_9cov, 
                     extract_hr_estimates(model_2, "TIMM", "nhm"),
                     extract_hr_estimates(model_misc2, "TIHMM", "nhm"),
                     extract_hr_estimates(model_age_cov5, "ApproxTIMM", "msm"),
                     extract_hr_estimates(model_5.2_misc, "ApproxTIMM", "msm")
                     )

HR_est_9cov <- HR_est_9cov %>%
  mutate(Variable = recode(Variable,
                           "if_ever_smoke"        = "smoke (Y/N)",
                           "dm_sex"               = "sex (F/M)",
                           "no_pa"                = "physical activity (N/Y)",
                           "life_alone"           = "living alone (Y/N)",
                           "fin_strain_early"     = "financial strain early (Y/N)",
                           "finstrain_dummy"      = "financial strain (Y/N)",
                           "educ_el"              = "education (L/H)",
                           "heavy_alcool"         = "alcohol (Y/N)",
                           "sei_long_cat_dummy"   = "manual occupation (Y/N)" 
  ))
HR_est_6cov <- HR_est_6cov %>%
  mutate(Variable = recode(Variable,
                           "if_ever_smoke"        = "smoke (Y/N)",
                           "dm_sex"               = "sex (F/M)",
                           "no_pa"                = "physical activity (N/Y)",
                           "life_alone"           = "living alone (Y/N)",
                           "educ_el"              = "education (L/H)",
                           "heavy_alcool"         = "alcohol (Y/N)"
  ))

plot_HR(HR_est_9cov)
plot_HR(HR_est_6cov)


