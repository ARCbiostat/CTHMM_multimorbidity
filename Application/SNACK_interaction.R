model_int_sex <- msm(MP ~ Age, subject=lopnr, data=snack_nhm, qmatrix= q0, deathexact = 3,
                       covariates = list("1-2" = ~ Age + educ_el:dm_sex + dm_sex + no_pa + no_pa:dm_sex+ life_alone + life_alone:dm_sex+ heavy_alcool:dm_sex+ if_ever_smoke:dm_sex, "1-3" = ~ Age, "2-3" = ~ Age), 
                       method = "BFGS", control = list(fnscale = 14000, maxit= 10000))
model_int_sex


model_int_sex_6 <- msm(MP ~ Age, subject=lopnr, data=snack_nhm, qmatrix= q0, deathexact = 3,
                      covariates = list("1-2" = ~ Age + educ_el + educ_el:dm_sex + dm_sex + no_pa + no_pa:dm_sex+ life_alone + life_alone:dm_sex+ heavy_alcool+ heavy_alcool:dm_sex+ if_ever_smoke:dm_sex, "1-3" = ~ Age, "2-3" = ~ Age), 
                      method = "BFGS", control = list(fnscale = 14000, maxit= 10000))
model_int_sex_6

model_int_sex_9 <- msm(MP ~ Age, subject=lopnr, data=snack_nhm_filtered, qmatrix= q0, 
                                         deathexact = 3,
                                         covariates = list("1-2" = ~ Age + educ_el*dm_sex + no_pa*dm_sex + life_alone*dm_sex + heavy_alcool*dm_sex+ if_ever_smoke*dm_sex+ fin_strain_early*dm_sex+finstrain_dummy*dm_sex+sei_long_cat_dummy*dm_sex, "1-3" = ~ Age, "2-3" = ~ Age), 
                                         control = list(fnscale = 20000, maxit= 10000))

model_int_sex_9

model_int_sex_f1 <- msm(MP ~ Age, subject=lopnr, data=snack_nhm, qmatrix= q0, 
                       deathexact = 3,
                       covariates = list("1-2" = ~ Age + no_pa*dm_sex + life_alone*dm_sex + heavy_alcool*dm_sex+ if_ever_smoke*dm_sex+sei_long_cat_dummy*dm_sex, "1-3" = ~ Age, "2-3" = ~ Age), 
                       control = list(fnscale = 20000, maxit= 10000))

model_int_sex_f1
# no_pa significant, interaction with sex not sign
# heavy_alcool significant interaction
model_int_sex_f2 <- msm(MP ~ Age, subject=lopnr, data=snack_nhm, qmatrix= q0, 
                        deathexact = 3,
                        covariates = list("1-2" = ~ Age + no_pa+dm_sex+  life_alone*dm_sex + sei_long_cat_dummy*dm_sex, "1-3" = ~ Age, "2-3" = ~ Age), 
                        control = list(fnscale = 20000, maxit= 10000))

model_int_sex_f2

########### altre prove ###########
model_p1 <- msm(MP ~ Age, subject=lopnr, data=snack_nhm, qmatrix= q0, 
                        deathexact = 3,
                        covariates = list("1-2" = ~ Age + no_pa+dm_sex+ + educ_el + life_alone*dm_sex + sei_long_cat_dummy*dm_sex+ heavy_alcool+ if_ever_smoke, "1-3" = ~ Age, "2-3" = ~ Age), 
                        control = list(fnscale = 20000, maxit= 10000))
model_p1

model_p2 <- msm(MP ~ Age, subject=lopnr, data=snack_nhm, qmatrix= q0, 
                deathexact = 3,
                covariates = list("1-2" = ~ Age + no_pa+dm_sex+ + educ_el + life_alone*dm_sex + sei_long_cat_dummy*dm_sex+ heavy_alcool+ if_ever_smoke + finstrain_dummy, "1-3" = ~ Age, "2-3" = ~ Age), 
                control = list(fnscale = 20000, maxit= 10000))
model_p2


###### NHM #######
q_nhm <- rbind(c(0,1,2),
               c(0,0,3),
               c(0,0,0))
colnames(q_nhm)<-c(1:dim(q_nhm)[1])
rownames(q_nhm)<-c(1:dim(q_nhm)[1])
# same structure for q_mat and nonh
nonh <- q_nhm

snack_nhm$sei_long_cat_dummy_sex <- snack_nhm$dm_sex*snack_nhm$sei_long_cat_dummy
covm_int <- list(
  educ_el = rbind(c(0,1,0), c(0,0,0), c(0,0,0)),
  dm_sex = rbind(c(0,2,0), c(0,0,0), c(0,0,0)),
  no_pa = rbind(c(0,3,0), c(0,0,0), c(0,0,0)),
  life_alone = rbind(c(0,4,0), c(0,0,0), c(0,0,0)),
  heavy_alcool = rbind(c(0,5,0), c(0,0,0), c(0,0,0)),
  if_ever_smoke = rbind(c(0,6,0), c(0,0,0), c(0,0,0)),
  sei_long_cat_dummy = rbind(c(0,7,0), c(0,0,0), c(0,0,0)),
  sei_long_cat_dummy_sex = rbind(c(0,8,0), c(0,0,0), c(0,0,0))
)

tic("nhm 7 cov interaction")
model_obj_gomp2<- model.nhm(state ~ Age, subject=lopnr, type='gompertz', data=snack_nhm, trans= q_nhm,
                            nonh= nonh,
                            covariates= c("educ_el", "dm_sex","no_pa" , "life_alone",  "heavy_alcool","if_ever_smoke", "sei_long_cat_dummy", "sei_long_cat_dummy_sex"),
                            covm = covm_int,
                            death=T, death.states = 3)

model_int_sex_7 <- nhm(model_obj_gomp2,
               #initial= nhm_init,
               gen_inits = TRUE, #not converging
               control=nhm.control(splits = c(60,61,63,65,75,76,85,95,99,102,103,105,109),
                                   ncores = 16
                                   #, verbose=TRUE
               )
)
toc()

misc_shape <- rbind(c(0,2,0),
                    c(1,0,0),
                    c(0,0,0))
tic()
model_int_sex_7_misc<- model.nhm(state ~ Age, subject=lopnr, type='gompertz', data=snack_nhm, trans= q_nhm,
                                nonh= nonh,
                                emat = misc_shape,
                                covariates= c("educ_el", "dm_sex","no_pa" , "life_alone",  "heavy_alcool","if_ever_smoke", "sei_long_cat_dummy"),
                                covm = covm1,
                                death=T, death.states = 3)
toc()