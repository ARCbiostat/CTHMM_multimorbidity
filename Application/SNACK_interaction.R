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


