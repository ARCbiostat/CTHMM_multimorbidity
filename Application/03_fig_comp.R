library(mvtnorm)
library(MMLCA)
library(ggplot2)
library(ggprism)


### alluvial plot

snack_nhm_7$Age_disc<- snack_nhm_7$Age/3
snack_nhm_7$MP_label <- factor(as.factor(snack_nhm_7$MP), levels = as.character(c(1,2,3)), labels = c("Mild MM", "Complex MM","Death"))

all_snack <- ggalluvial(
  snack_nhm_7,
  time_var = "Age_disc",
  "lopnr",
  "MP_label",
  colors= c("lightblue", "dodgerblue", "grey"),
  space = 20,
  "Death"
)

all_snack$plot+scale_x_continuous("Age",breaks = seq(20,35,by=2),labels = c("60-63","66-69", "72-75", "78-81", "84-87", "90-93", "96-99", "100+"))



### transition to Death
t_vals <- seq(60, 90, length.out = 100)

preds <- predict.nhm.mine(model_misc3, times=t_vals, flag="Transition Probability", time0=60)
nstate <- model_misc3$nstate
prevalence <- preds$probabilities
times <- preds$times
prevL <- preds$lower
prevU <- preds$upper
# main_arg <- "trans prob"
# ci <- TRUE
# par(mfrow = c(2, ceiling(nstate/2)))
# for (i in 1:nstate) {
#   plot(times, prevalence[, i], type = "l", xlab = "Age", 
#        ylab = "Probability", main = paste(main_arg, 
#                                           i, sep = ""), ylim = c(0, 1))
#   if (ci) {
#     lines(times, prevL[, i], lty = 2)
#     lines(times, prevU[, i], lty = 2)
#   }
#   else {
#     prevL <- prevU <- NULL
#   }
# }

trans_death_data=data.frame(age=rep(times,2),
                            prob=c(prevalence[, 2],prevalence[,3]),
                            lower=c(prevL[, 2],prevL[,3]),
                            upper=c(prevU[, 2],prevU[,3]),
                            from=c(rep("From Mild MM",length(times)),rep("From Complex MM",length(times))))


trans_death_data$from=factor(as.factor(trans_death_data$from),levels = c("From Mild MM", "From Complex MM"))
Fig_appl_4b=ggplot(trans_death_data)+
  geom_line(aes(age,prob,group=from,color=from),linewidth=1)+
  geom_ribbon(aes(age,ymin=lower,ymax=upper,group=from,fill=from),alpha=0.3)+
  scale_fill_manual(values = c("lightblue", "dodgerblue"))+
  scale_color_manual(values = c("lightblue", "dodgerblue"))+
  scale_x_continuous("Age (years)")+
  scale_y_continuous("Probability")+
  ggtitle("Transition to Death")+
  theme_prism()+
  theme(legend.position = "bottom")

Fig_appl_4b

############ HR ratio
# helper functions in helper_f_snack.R
HR_est_7cov <- data.frame()
HR_est_7cov <- rbind(HR_est_7cov, 
                     extract_hr_estimates(model_misc3, "TIHMM", "nhm"))


HR_est_7cov["dm_sex",2] <- 1/HR_est_7cov["dm_sex",2]
HR_est_7cov["dm_sex", c(3, 4)] <- 1 / rev(HR_est_7cov["dm_sex", c(3, 4)])

HR_est_7cov <- HR_est_7cov %>%
  mutate(Variable = recode(Variable,
                           "if_ever_smoke"        = "Smoking (Y/N)",
                           "dm_sex"               = "Sex (M/F)",
                           "no_pa"                = "Sedentarism (Y/N)",
                           "life_alone"           = "Living alone (Y/N)",
                           "educ_el"              = "Elementary Education (Y/N)",
                           "heavy_alcool"         = "Alcohol (Y/N)",
                           "sei_long_cat_dummy"   = "Manual occupation (Y/N)" 
  ))



HR_est_7cov$Variable <- factor(HR_est_7cov$Variable, levels = rev(unique(HR_est_7cov$Variable)))

Fig_appl_4a <- ggplot(HR_est_7cov, aes(x = HR, y = Variable)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = L, xmax = U), height = 0.2, size=1.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  scale_x_continuous(limits = c(0.5, 2)) +
  labs(
    x = "Hazard Ratio", 
    y = "", 
    title = "Transition from Mild to Complex MM"
  ) +
  theme_prism() +
  theme(
    axis.text.y = element_text(size = 10)
  )
Fig_appl_4a

########## Transition mild to complex MM by covariates

fix_pred_complex=function(model,preds,lab){
  nstate <- model$nstate
  prevalence <- preds$probabilities
  times <- preds$times
  prevL <- preds$lower
  prevU <- preds$upper
  
  pred_data=data.frame(age=times,
                              prob=prevalence[, 1],
                              lower=prevL[, 1],
                              upper=prevU[, 1],
                              group=lab)
  return(pred_data)
  
}

covvalue=model_misc3$covmeans
covvalue[2]=0
pred_sex_male <- predict.nhm.mine(model_misc3, 
                             times=t_vals, 
                             flag="Transition Probability", time0=60,
                             covvalue =covvalue)

data_pred_male=fix_pred_complex(model_misc3,pred_sex_male,lab = "Males")

covvalue=model_misc3$covmeans
covvalue[2]=1
pred_sex_female <- predict.nhm.mine(model_misc3, 
                                  times=t_vals, 
                                  flag="Transition Probability", time0=60,
                                  covvalue =covvalue)


data_pred_female=fix_pred_complex(model_misc3,pred_sex_female,lab = "Females")
data_pred_sex=rbind(data_pred_female,data_pred_male)



Fig_appl_res2a=ggplot(data_pred_sex)+
  geom_line(aes(age,prob,group=group,color=group),linewidth=1)+
  geom_ribbon(aes(age,ymin=lower,ymax=upper,group=group,fill=group),alpha=0.3)+
  scale_fill_manual(values = c("purple", "red"))+
  scale_color_manual(values = c("purple", "red"))+
  scale_x_continuous("Age (years)")+
  scale_y_continuous("Probability",limits = c(0,0.7))+
  ggtitle("Sex")+
  theme_prism()+
  theme(legend.position = "bottom")

Fig_appl_res2a


covvalue=model_misc3$covmeans
covvalue[3]=1

pred_nopa<- predict.nhm.mine(model_misc3, 
                                  times=t_vals, 
                                  flag="Transition Probability", time0=60,
                                  covvalue =covvalue)

data_pred_nopa=fix_pred_complex(model_misc3,pred_nopa,lab = "Yes")

covvalue[3]=0
pred_pa <- predict.nhm.mine(model_misc3, 
                                    times=t_vals, 
                                    flag="Transition Probability", time0=60,
                                    covvalue =covvalue)


data_pred_pa=fix_pred_complex(model_misc3,pred_pa,lab = "No")
data_pred_paall=rbind(data_pred_nopa,data_pred_pa)



Fig_appl_res2b=ggplot(data_pred_paall)+
  geom_line(aes(age,prob,group=group,color=group),linewidth=1)+
  geom_ribbon(aes(age,ymin=lower,ymax=upper,group=group,fill=group),alpha=0.3)+
  scale_fill_manual(values = c("forestgreen", "orange"))+
  scale_color_manual(values = c("forestgreen", "orange"))+
  scale_x_continuous("Age (years)")+
  scale_y_continuous("Probability",limits = c(0,0.7))+
  ggtitle("Sedentarism")+
  theme_prism()+
  theme(legend.position = "bottom")

Fig_appl_res2b


ggpubr::ggarrange(Fig_appl_res2a,Fig_appl_res2b)





# pred_mean <- predict.nhm.mine(model_misc3, 
#                                     times=t_vals, 
#                                     flag="Transition Probability", time0=60)
# 
# 
# data_pred_mean=fix_pred_complex(model_misc3,pred_mean,lab = "Mean")
# 
# 
# Fig_appl_res2=ggplot(data_pred_mean)+
#   geom_line(aes(age,prob),linewidth=1)+
#   geom_ribbon(aes(age,ymin=lower,ymax=upper),alpha=0.3)+
#   scale_x_continuous("Age (years)")+
#   scale_y_continuous("Probability",limits = c(0,0.7))+
#   ggtitle("From Mild to Complex MM")+
#   theme_prism()+
#   theme(legend.position = "bottom")
# 
# Fig_appl_res2



