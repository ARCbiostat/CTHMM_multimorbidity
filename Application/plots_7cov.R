library(ggplot2)
library(ggprism)
library(dplyr)
library(msm)
source("Functions/helper_f_snack.R")

# helper functions in helper_f_snack.R
HR_est_9cov <- data.frame()
HR_est_9cov <- rbind(HR_est_9cov, 
                     extract_hr_estimates(model_2, "TIMM", "nhm"),
                     extract_hr_estimates(model_misc2, "TIHMM", "nhm"),
                     extract_hr_estimates(model_age_cov5, "ApproxTIMM", "msm"),
                     extract_hr_estimates(model_5.2_misc, "ApproxTIHMM", "msm")
)

HR_est_9cov <- HR_est_9cov %>%
  mutate(Variable = recode(Variable,
                           "if_ever_smoke"        = "Smoke (Y/N)",
                           "dm_sex"               = "Sex (F/M)",
                           "no_pa"                = "Physical activity (N/Y)",
                           "life_alone"           = "Living alone (Y/N)",
                           "fin_strain_early"     = "Financial strain early (Y/N)",
                           "finstrain_dummy"      = "Financial strain (Y/N)",
                           "educ_el"              = "Education (L/H)",
                           "heavy_alcool"         = "Alcohol (Y/N)",
                           "sei_long_cat_dummy"   = "Manual occupation (Y/N)" 
  ))
HR_est_6cov <- data.frame()
HR_est_6cov <- rbind(HR_est_6cov, 
                     extract_hr_estimates(model_1, "TIMM", "nhm"),
                     extract_hr_estimates(model_misc, "TIHMM", "nhm"),
                     extract_hr_estimates(model_age_cov6, "ApproxTIMM", "msm"),
                     extract_hr_estimates(model_6_misc, "ApproxTIHMM", "msm")
)
HR_est_6cov <- HR_est_6cov %>%
  mutate(Variable = recode(Variable,
                           "if_ever_smoke"        = "Smoke (Y/N)",
                           "dm_sex"               = "Sex (F/M)",
                           "no_pa"                = "Physical activity (N/Y)",
                           "life_alone"           = "Living alone (Y/N)",
                           "educ_el"              = "Elementary education (N/Y)",
                           "heavy_alcool"         = "Alcohol (Y/N)"
  ))

plot_HR(HR_est_9cov)
plot_HR(HR_est_6cov)




HR_est_7cov <- data.frame()
HR_est_7cov <- rbind(HR_est_7cov, 
                     extract_hr_estimates(model_misc3, "TIHMM", "nhm"))
HR_est_7cov <- HR_est_7cov %>%
  mutate(Variable = recode(Variable,
                           "if_ever_smoke"        = "Smoking (Y/N)",
                           "dm_sex"               = "Sex (F/M)",
                           "no_pa"                = "Physical activity (N/Y)",
                           "life_alone"           = "Living alone (Y/N)",
                           "educ_el"              = "Education (L/H)",
                           "heavy_alcool"         = "Alcohol (Y/N)",
                           "sei_long_cat_dummy"   = "Manual occupation (Y/N)" 
  ))

plot_HR(HR_est_7cov)
ggsave(plot_HR(HR_est_7cov),file="Application/Figures/HR_plot.jpeg",dpi=300)


plot_fm <-plot.nhm.mine2(model_misc3, what= "intensity",trans = 1, covvalue = c(0,1,0,0,0,0,0), colours = c("blue","red"),
                         labels = c("Female","Male"), time0=60, times= t_vals, main_arg= "Transition Mild to Complex MM", xlab="Age")
plot_fm
# no physical activity (orange) vs physical activity (green)
plot_pa <- plot.nhm.mine2(model_misc3, what= "intensity",trans = 1, covvalue = c(0,0,1,0,0,0,0), colours = c("darkorange","forestgreen"),
              labels = c("No physical activity","Physical activity"), time0=60, times= t_vals, main_arg= "Transition Mild to Complex MM", xlab="Age")


