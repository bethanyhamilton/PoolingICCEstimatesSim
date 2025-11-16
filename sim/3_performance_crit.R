# performance criteria related to the ICC estimate

calc_performance <- function(results) {
  
  performance_measures <- 
    results |> 
    group_by(method, 
             var_icc_name, 
             icc_est_n, 
             nj_size, 
             n_bar_size, 
             n_bar_prop, 
             var_combo,
             tau) |> 
    mutate(
           l1_es = (residvar_mean - l1_var) / residvar_sd,
           l2_es = (clustervar_mean - l2_var) / clustervar_sd,
           se_ICC_true = ifelse(var_icc_name != "Fisher TF", sd(pooled_icc_est), NA),
           pooled_icc_est2 = ifelse(var_icc_name == "Fisher TF", .5 * log((1 + (n_weighted - 1) * pooled_icc_est) / (1 - pooled_icc_est)), NA),
           se_ICC_true= ifelse(var_icc_name == "Fisher TF", sd(pooled_icc_est2), se_ICC_true)#,
         #  CI_lower = pooled_icc_est - (qt(p = .975, df = (icc_est_n - 1)) * se_pooled_icc_est),
          # CI_upper = pooled_icc_est + (qt(p = .975, df = (icc_est_n - 1)) * se_pooled_icc_est),
           ) |>
    summarize(
              # number of iterations
              K = n(),
      
              # getting one SE value for each condition. 
              se_ICC_true = mean(se_ICC_true),
              
              # ICC estimate
              variance = var(pooled_icc_est),
              ICC_true = (mean(l2_var) / (mean(l2_var) + mean(l1_var))),
              
              # PB -- MCSE
              bias_icc = mean(pooled_icc_est) - ICC_true,
              bias_icc_mcse = sqrt(variance / K),
              
              # MSE
              MSE = sum((pooled_icc_est - ICC_true)^2) / K,
              
              # RMSE
              RMSE = sqrt(MSE),
              ## alt RMSE that is slightly rounded
              RMSE_rel_bias_icc = sqrt(var(pooled_icc_est) + bias_icc^2),
              
              # RPB
              rel_bias_icc = (mean(pooled_icc_est) - ICC_true) / (ICC_true),
              rel_bias_icc_mcse = sqrt(variance / K)*(1/ICC_true),
              
              # RPB -- ratio
              rel_bias_icc_alt = mean(pooled_icc_est) / (ICC_true),
              
              # SE of ICC
              rel_bias_se_ICC = (mean(se_pooled_icc_est) - se_ICC_true) / (se_ICC_true),
              
              # skewness
              skewness = (1/(K*sqrt(variance)^3)*sum((pooled_icc_est-mean(pooled_icc_est))^3)),
              
              # kurtosis
              kurtosis = (1/(K*variance^2))*sum((pooled_icc_est-mean(pooled_icc_est))^4),
              
              # coverage
             # coverage = mean(CI_lower <= ICC_true & ICC_true <= CI_upper),
              
              .groups = 'drop')
  
  return(performance_measures)
  
}

# performance criteria related to the variance components

calc_performance_var <- function(results) {
  
  performance_measures <- 
    results |>  
    select(-method, 
           -var_icc_name, 
           -pooled_icc_est, 
           -se_pooled_icc_est, 
           -tau_est, 
           -n_weighted, 
           -L2_new) |>  
    distinct() |> 
    group_by(icc_est_n, 
             nj_size, 
             n_bar_size, 
             n_bar_prop, 
             var_combo,
             tau) |> 
  summarize(K = n(),
            l2_var = mean(l2_var),
            l1_var = mean(l1_var),
            
            rel_bias_clus = (mean(clustervar_median) - l2_var) / (l2_var),
            
            rel_bias_clus_ratio = mean(clustervar_median) / l2_var,
            
            bias_clus = (mean(clustervar_median) - l2_var),

            RMSE_rel_bias_clus = sqrt((var(clustervar_median)) + bias_clus^2),

            rel_bias_res = (mean(residvar_median) - l1_var) / (l1_var),
            
            rel_bias_res_ratio = mean(residvar_median) / (l1_var),
            
            bias_res = (mean(residvar_median) - l1_var),

            RMSE_rel_bias_res = sqrt((var(residvar_median)) + bias_res^2),
            
            .groups = 'drop') |>  distinct()
  
  return(performance_measures)
  
}





