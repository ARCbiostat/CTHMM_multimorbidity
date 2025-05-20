library(mvtnorm)

t_vals <- seq(60, 90, length.out = 100)

preds <- predict.nhm.mine(model_misc3, times=t_vals, flag="Transition Probability", time0=60)
#plot.nhm.mine2(model_misc3, times=seq(60,95,1), flag="Transition Probability")
#plot.nhm.mine(model_misc3, flag="Transition Probability", what= "probabilities",time0=60, times= t_vals, main_arg= "Hazard function", xlab="Age")

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
Fig_appl_res3=ggplot(trans_death_data)+
  geom_line(aes(age,prob,group=from,color=from),linewidth=1)+
  geom_ribbon(aes(age,ymin=lower,ymax=upper,group=from,fill=from),alpha=0.3)+
  scale_fill_manual(values = c("lightblue", "dodgerblue"))+
  scale_color_manual(values = c("lightblue", "dodgerblue"))+
  scale_x_continuous("Age (years)")+
  scale_y_continuous("Probability")+
  ggtitle("Transition to Death")+
  theme_prism()+
  theme(legend.position = "bottom")

Fig_appl_res3
#########################################################################################################


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

data_pred_nopa=fix_pred_complex(model_misc3,pred_nopa,lab = "No")

covvalue[3]=0
pred_pa <- predict.nhm.mine(model_misc3, 
                                    times=t_vals, 
                                    flag="Transition Probability", time0=60,
                                    covvalue =covvalue)


data_pred_pa=fix_pred_complex(model_misc3,pred_pa,lab = "Yes")
data_pred_paall=rbind(data_pred_nopa,data_pred_pa)



Fig_appl_res2b=ggplot(data_pred_paall)+
  geom_line(aes(age,prob,group=group,color=group),linewidth=1)+
  geom_ribbon(aes(age,ymin=lower,ymax=upper,group=group,fill=group),alpha=0.3)+
  scale_fill_manual(values = c("forestgreen", "orange"))+
  scale_color_manual(values = c("forestgreen", "orange"))+
  scale_x_continuous("Age (years)")+
  scale_y_continuous("Probability",limits = c(0,0.7))+
  ggtitle("Physical Activity")+
  theme_prism()+
  theme(legend.position = "bottom")

Fig_appl_res2b


ggpubr::ggarrange(Fig_appl_res2a,Fig_appl_res2b)







pred_mean <- predict.nhm.mine(model_misc3, 
                                    times=t_vals, 
                                    flag="Transition Probability", time0=60)


data_pred_mean=fix_pred_complex(model_misc3,pred_mean,lab = "Mean")


Fig_appl_res2=ggplot(data_pred_mean)+
  geom_line(aes(age,prob),linewidth=1)+
  geom_ribbon(aes(age,ymin=lower,ymax=upper),alpha=0.3)+
  scale_x_continuous("Age (years)")+
  scale_y_continuous("Probability",limits = c(0,0.7))+
  ggtitle("From Mild to Complex MM")+
  theme_prism()+
  theme(legend.position = "bottom")

Fig_appl_res2


ggsave(Fig_appl_res2,filename="")