


compute_pm <- function(param_df, true_param) {
  
  names(true_param) <- c("trans", "rate", "shape", "beta_cov1", "beta_cov2", "beta_cov3")
  
  param_long <- param_df %>%
    pivot_longer(cols = c(rate, shape, beta_cov1, beta_cov2,beta_cov3),
                 names_to = "parameter", values_to = "estimate_value") %>%
    pivot_longer(cols = c(rate_se, shape_se, beta_cov1_se, beta_cov2_se,beta_cov3_se),
                 names_to = "se_param", values_to = "se_value") %>%
    mutate(se_param = sub("_se$", "", se_param)) %>%
    filter(parameter == se_param) %>%
    dplyr::select(-se_param)
  param_long

  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_cov1, beta_cov2, beta_cov3),
                 names_to = "parameter", values_to = "true_value")
  
  
  results <- param_long %>%
    left_join(true_param_long, by = c("trans", "parameter")) %>%
    group_by(model, trans, parameter) %>%
    summarise(
      theta= mean(true_value),     
      estimate = mean(estimate_value),
      se_estimate = sqrt(mean(se_value^2, na.rm = TRUE) + var(estimate_value, na.rm = TRUE)/n()),
      bias = mean(estimate_value - true_value, na.rm = TRUE),
      se_bias = sd(estimate_value - true_value, na.rm = TRUE) / sqrt(n()),
      rel_bias = mean((estimate_value - true_value) / abs(true_value), na.rm = TRUE),
      se_rel_bias = sd((estimate_value - true_value) / abs(true_value), na.rm = TRUE) / sqrt(n()),
      mse = mean((estimate_value - true_value)^2, na.rm = TRUE),
      coverage = calculate_coverage(estimate_value, se_value, true_value),
      .groups = "drop"
    )
  
  return(results)
}

compute_pm2 <- function(param_df,true_param){
  
  names(true_param) <- c("trans","t_rate","t_shape","t_beta_cov1","t_beta_cov2","t_beta_cov3")
  results <- param_df %>%
    left_join(true_param, by = "trans") %>%
    dplyr::select(-dataset_id) %>%
    group_by(model, trans) %>%
    summarise(
      bias_rate = mean(rate - t_rate, na.rm = TRUE),
      se_bias_rate = sqrt(var(rate - t_rate, na.rm = TRUE) / n() + rate_se^2),
      rel_bias_rate = mean((rate - t_rate) / abs(t_rate), na.rm = TRUE),
      se_rel_bias_rate = sqrt((var(rate - t_rate, na.rm = TRUE) / n() + rate_se^2)/t_rate, nar.rm=TRUE),
      mse_rate = mean((rate - t_rate)^2, na.rm = TRUE),
      coverage_rate = calculate_coverage(rate, rate_se, t_rate),
      
      bias_shape = mean(shape - t_shape, na.rm = TRUE),
      se_bias_shape = sd(shape - t_shape, na.rm = TRUE) / sqrt(n()),
      rel_bias_shape = mean((shape - t_shape) / abs(t_shape), na.rm = TRUE),
      se_rel_bias_shape = sd((shape - t_shape) / abs(t_shape), na.rm = TRUE) / sqrt(n()),
      mse_shape = mean((shape - t_shape)^2, na.rm = TRUE),
      coverage_shape = calculate_coverage(shape, shape_se, t_shape),
      
      bias_beta_cov1 = mean(beta_cov1 - t_beta_cov1, na.rm = TRUE),
      se_bias_beta_cov1 = sd(beta_cov1 - t_beta_cov1, na.rm = TRUE) / sqrt(n()),
      rel_bias_beta_cov1 = mean((beta_cov1 - t_beta_cov1) / abs(t_beta_cov1), na.rm = TRUE),
      se_rel_bias_beta_cov1 = sd((beta_cov1 - t_beta_cov1) / abs(t_beta_cov1), na.rm = TRUE) / sqrt(n()),
      mse_beta_cov1 = mean((beta_cov1 - t_beta_cov1)^2, na.rm = TRUE),
      coverage_beta_cov1 = calculate_coverage(beta_cov1, beta_cov1_se, t_beta_cov1),
      
      bias_beta_cov2 = mean(beta_cov2 - t_beta_cov2, na.rm = TRUE),
      se_bias_beta_cov2 = sd(beta_cov2 - t_beta_cov2, na.rm = TRUE) / sqrt(n()),
      rel_bias_beta_cov2 = mean((beta_cov2 - t_beta_cov2) / abs(t_beta_cov2), na.rm = TRUE),
      se_rel_bias_beta_cov2 = sd((beta_cov2 - t_beta_cov2) / abs(t_beta_cov2), na.rm = TRUE) / sqrt(n()),
      mse_beta_cov2 = mean((beta_cov2 - t_beta_cov2)^2, na.rm = TRUE),
      coverage_beta_cov2 = calculate_coverage(beta_cov2, beta_cov2_se, t_beta_cov2),
      
      bias_beta_cov3 = mean(beta_cov3 - t_beta_cov3, na.rm = TRUE),
      se_bias_beta_cov3 = sd(beta_cov3 - t_beta_cov3, na.rm = TRUE) / sqrt(n()),
      rel_bias_beta_cov3 = mean((beta_cov3 - t_beta_cov3) / abs(t_beta_cov3), na.rm = TRUE),
      se_rel_bias_beta_cov3 = sd((beta_cov3 - t_beta_cov3) / abs(t_beta_cov3), na.rm = TRUE) / sqrt(n()),
      mse_beta_cov3 = mean((beta_cov3 - t_beta_cov3)^2, na.rm = TRUE),
      coverage_beta_cov3 = calculate_coverage(beta_cov3, beta_cov3_se, t_beta_cov3),
      
      
      
      .groups = "drop" 
    )
  return(results)
  
  
}

calculate_coverage <- function(estimates, se, true_param, confidence_level = 0.95) {
  alpha <- 1 - confidence_level
  z <- qnorm(1 - alpha / 2)  
  
  lower_ci <- estimates - z * se
  upper_ci <- estimates + z * se
  within_ci <- (true_param >= lower_ci) & (true_param <= upper_ci)
  
  coverage <- mean(within_ci, na.rm = TRUE)
  
  
  return(coverage)
}


