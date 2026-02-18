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

#####################################################
###########TIMM and TIHMM: covariates################

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
  if_ever_smoke = rbind(c(0,5,0), c(0,0,0), c(0,0,0)),
  heavy_alcool = rbind(c(0,6,0), c(0,0,0), c(0,0,0)),
  sei_long_cat_dummy = rbind(c(0,7,0), c(0,0,0), c(0,0,0))
)

tic("nhm gompertz with 6 cov")
model_obj_gomp1<- model.nhm(state ~ Age, subject=lopnr, type='gompertz', data=snack_nhm, trans= q_nhm,
                            nonh= nonh,
                            covariates= c("educ_el", "dm_sex","no_pa" , "life_alone", "if_ever_smoke", "heavy_alcool","sei_long_cat_dummy"),
                            covm = covm1,
                            death=T, death.states = 3)

model_1 <- nhm(model_obj_gomp1,
               #initial= nhm_init,
               gen_inits = TRUE, #not converging
               control=nhm.control(splits = c(60,61,63,65,75,76,85,95,99,102,103,105,109),
                                   ncores = 1
                                   #, verbose=TRUE
               )
)
toc()

nhm_init <- c(model_1$par)
model_1.1 <- nhm(model_obj_gomp1,
                 initial= nhm_init,
                 #gen_inits = TRUE, #not converging
                 control=nhm.control(splits = c(60,61,63,65,75,76,85,95,99,102,103,105,109),
                                     ncores = 1
                                     #, verbose=TRUE
                 )
)
#
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save(model_1.1, file =paste0(result_folder,"/TIMM_6cov_",timestamp,".RData"))


# misc model
non_diagonal <- misc
diag(non_diagonal) <- 0
misc_param <- non_diagonal[non_diagonal != 0]
nhm_init <- c(model_1.1$par, log(misc_param))



tic("nhm misc gompertz with 6 cov")
model_obj_gomp_misc<- model.nhm(state ~ Age, subject=lopnr, type='gompertz', data=snack_nhm, trans= q_nhm,
                                nonh= nonh,
                                emat = misc_shape,
                                covariates= c("educ_el", "dm_sex","no_pa" , "life_alone", "if_ever_smoke", "heavy_alcool","sei_long_cat_dummy"),
                                covm = covm1,
                                death=T, death.states = 3)


# misc parameters need to be passed to initial and then fixed in fixedpar
n1<- length(model_1.1$par) + 1
n2<- length(nhm_init)
split_points <- c(60,61,63,65,75,76,85,95,99,102,103,104,105,106,109)

model_misc3 <- nhm(
  model_obj_gomp_misc,
  initial = nhm_init,
  fixedpar = c(n1:n2),
  control = nhm.control(
    splits = split_points,
    ncores = 1, 
    rtol = 1e-8, 
    atol = 1e-8
  )
)


save.image("fits_snack.RData")

######### Plots ###########

t_vals <- seq(60, 90, length.out = 100)

# female (blu) vs male (red)
plot.nhm.mine(model_misc, what= "intensity",trans = 1, covvalue = c(0,1,0,0,0,0), colours = c("blue","red"), colours_fill = c("lightblue","pink"),
               labels = c("Female","Male"), time0=60, times= t_vals, main_arg= "Hazard function", xlab="Age")
# no physical activity (orange) vs physical activity (green)
plot.nhm.mine(model_misc, what= "intensity",trans = 1, covvalue = c(0,0,1,0,0,0), colours = c("darkorange","forestgreen"), colours_fill = c("moccasin","#ccffcc"),
               labels = c("No physical activity","Physical activity"), time0=60, times= t_vals, main_arg= "Hazard function", xlab="Age")

plot.nhm.mine(model_misc, time0=60, xlab="Age")
plot.nhm.mine(model_misc, what= "intensities",time0=60, times= t_vals, main_arg= "Hazard function", xlab="Age")
plot.nhm.mine(model_misc, what= "probabilities",time0=60, times= t_vals, main_arg= "Hazard function", xlab="Age")


###### HR plots ########
library(ggplot2)
library(ggprism)

#helper functions in helper_f_snack.R
HR_est_cov <- data.frame()
HR_est_cov <- rbind(HR_est_cov,
                     extract_hr_estimates(model_misc3, "TIMM", "nhm")
                     )

HR_est_cov <- HR_est_cov %>%
  mutate(Variable = recode(Variable,
                           "if_ever_smoke"        = "Smoking (Y/N)",
                           "dm_sex"               = "Sex (F/M)",
                           "no_pa"                = "Sedentarism (N/Y)",
                           "life_alone"           = "Living alone (Y/N)",
                           "educ_el"              = "Elementary Education (Y/N)",
                           "heavy_alcool"         = "Alcohol (Y/N)",
                           "sei_long_cat_dummy"   = "Manual occupation (Y/N)"
  ))


plot_HR(HR_est_cov)

ggsave(plot_HR(HR_est_cov),file="Application/Figures/HR_plot.jpeg",dpi=300)

