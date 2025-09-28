calc_performance <- function(results) {
  
  # performance criteria related to the ICC estimate

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
    mutate(K = n(),
           l1_es = (residvar_mean - l1_var) / residvar_sd,
           l2_es = (clustervar_mean - l2_var) / clustervar_sd,
           se_ICC_true = ifelse(var_icc_name != "Fisher TF", sd(pooled_icc_est), NA),
           pooled_icc_est2 = ifelse(var_icc_name == "Fisher TF", .5 * log((1 + (n_weighted - 1) * pooled_icc_est) / (1 - pooled_icc_est)), NA),
           se_ICC_true= ifelse(var_icc_name == "Fisher TF", sd(pooled_icc_est2), se_ICC_true)) |>
    summarize(se_ICC_true = mean(se_ICC_true),
              
              # ICC estimate
              variance = var(pooled_icc_est),
              ICC_true = (mean(l2_var) / (mean(l2_var) + mean(l1_var))),
              
              # PB -- MCSE
              bias_icc = mean(pooled_icc_est) - ICC_true,
              bias_icc_mcse = sqrt(variance / K),
              
              # RMSE
              RMSE = sqrt(sum(pooled_icc_est - ICC_true)^2 / K),
              
              # RPB
              rel_bias_icc = (mean(pooled_icc_est) - ICC_true) / (ICC_true),
              
              # Relative RMSE
              RMSE_rel_bias_icc = sqrt((var(pooled_icc_est)) + rel_bias_icc^2),
              
              # SE of ICC
              rel_bias_se_ICC = (mean(se_pooled_icc_est) - se_ICC_true) / (se_ICC_true),
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

            RMSE_rel_bias_clus = sqrt((var(clustervar_median)) + rel_bias_clus^2),

            rel_bias_res = (mean(residvar_median) - l1_var) / (l1_var),

            RMSE_rel_bias_res = sqrt((var(residvar_median)) + rel_bias_res^2),
            
            .groups = 'drop') |>  distinct()
  
  return(performance_measures)
  
}





