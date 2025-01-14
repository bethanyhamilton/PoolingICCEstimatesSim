calc_performance <- function(results) {
  
  # require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  # require(tidyr, quietly = TRUE, warn.conflicts = FALSE)
  
  performance_measures <- 
    results %>%
    group_by(method, var_icc_name, icc_est_n, nj_size, n_bar_size, n_bar_prop, var_combo,
             tau
    ) %>%
    mutate(K = n(),
           # test = mean(pooled_icc_est),
           
           # For jackknife RMSE
           T_j = (1/(K-1)) *((K*mean(pooled_icc_est))-pooled_icc_est),
           var_T_j = (1/(K-2))*(((K-1)*var(pooled_icc_est))-((K/(K-1))*(pooled_icc_est - mean(pooled_icc_est))^2)),
           RMSE_rel_bias_icc_j = sqrt( ((T_j-(l2_var/(l2_var+l1_var)))^2+var_T_j)/ (l2_var/(l2_var+l1_var))^2),
           
           # SMD per replication of L1 and L2 variance. Want to see how many SDs off.
           l1_es = (residvar_mean - l1_var)/residvar_sd,
           l2_es = (clustervar_mean - l2_var)/clustervar_sd,
           
           
    ) %>%
    summarize(K = n(),
              # ICC estimate
              variance = var(pooled_icc_est),
              ICC_true = (mean(l2_var)/(mean(l2_var)+mean(l1_var))),
              rel_bias_icc2 = (mean(pooled_icc_est))/(ICC_true),
              rel_bias_icc = (mean(pooled_icc_est) - ICC_true)/(ICC_true),
              rel_bias_icc_median = (median(pooled_icc_est) - ICC_true)/(ICC_true),
              
              
              rel_bias_icc_mcse = sqrt(variance/(K*(ICC_true^2))),
              bias_icc_mcse = sqrt(variance/(K)),
              
              RMSE_rel_bias_icc2 = sqrt( (var(pooled_icc_est)/ (ICC_true^2) ) + rel_bias_icc^2),
              RMSE_rel_bias_icc_mcse = sqrt(((K-1)/K)*sum((RMSE_rel_bias_icc_j-RMSE_rel_bias_icc2)^2)),
              
              
              RMSE_rel_bias_icc = sqrt( (var(pooled_icc_est)) + rel_bias_icc^2),
              
              RMSE_rel_bias_icc_median = sqrt( (var(pooled_icc_est)) + rel_bias_icc_median^2),
              
              # SE of ICC
              
              
              se_ICC_true= ifelse(var_icc_name!= "Fisher TF", sd(pooled_icc_est), NA),
              
              pooled_icc_est2 = ifelse(var_icc_name == "Fisher TF", .5*log((1+(n_weighted-1)*pooled_icc_est)/(1-pooled_icc_est)), NA),
              
              se_ICC_true= ifelse(var_icc_name == "Fisher TF", sd(pooled_icc_est2), se_ICC_true),
              
              rel_bias_se_ICC = (mean(se_pooled_icc_est) - se_ICC_true)/(se_ICC_true),
              
              RMSE_rel_bias_se_ICC = sqrt( (var(se_pooled_icc_est) ) + rel_bias_se_ICC^2),
              
              rel_bias_se_ICC_median = (median(se_pooled_icc_est) - se_ICC_true)/(se_ICC_true),
              
              RMSE_rel_bias_se_ICC_median = sqrt( (var(se_pooled_icc_est) ) + rel_bias_se_ICC_median^2),
              
              
              
              .groups = 'drop')
  
  return(performance_measures)
}


calc_performance_var <- function(results) {
  
  # require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  # require(tidyr, quietly = TRUE, warn.conflicts = FALSE)
  
  performance_measures <- 
    results %>% select(-method, -var_icc_name, -pooled_icc_est, -se_pooled_icc_est) %>% distinct() %>%
    group_by( icc_est_n, nj_size, n_bar_size, n_bar_prop, var_combo,
              tau
    ) %>%
    # mutate(K = n(),
    # test = mean(pooled_icc_est),
    
    # For jackknife RMSE
    #   T_j = (1/(K-1)) *((K*mean(pooled_icc_est))-pooled_icc_est),
    #   var_T_j = (1/(K-2))*(((K-1)*var(pooled_icc_est))-((K/(K-1))*(pooled_icc_est - mean(pooled_icc_est))^2)),
    #   RMSE_rel_bias_icc_j = sqrt( ((T_j-(l2_var/(l2_var+l1_var)))^2+var_T_j)/ (l2_var/(l2_var+l1_var))^2),
    
    # SMD per replication of L1 and L2 variance. Want to see how many SDs off.
    #   l1_es = (residvar_mean - l1_var)/residvar_sd,
    #   l2_es = (clustervar_mean - l2_var)/clustervar_sd,
  
  
  # ) %>%
  summarize(K = n(),
            
            ### Maybe change this later -- Should talk to Tasha
            
            variance_res = var(residvar_median),
            variance_clus = var(clustervar_median),
            
            rel_bias_clus2 = (mean(clustervar_median))/(l2_var),
            rel_bias_clus = (mean(clustervar_median) - l2_var)/(l2_var),
            rel_bias_clus3 = (median(clustervar_median) - l2_var)/(l2_var),
            
            rel_bias_clus_mcse = sqrt(variance_clus/(K*(l2_var^2))),
            RMSE_rel_bias_clus2 = sqrt( (var(clustervar_median)/ (l2_var^2) ) + rel_bias_clus^2),
            RMSE_rel_bias_clus = sqrt( (var(clustervar_median) ) + rel_bias_clus^2),
            
            RMSE_rel_bias_clus_median = sqrt( (var(clustervar_median) ) + rel_bias_clus3^2),
            
            
            rel_bias_res2 = (mean(residvar_median))/(l1_var),
            rel_bias_res = (mean(residvar_median) - l1_var)/(l1_var),
            rel_bias_res3 = (median(residvar_median) - l1_var)/(l1_var),
            
            rel_bias_res_mcse = sqrt(variance_res/(K*(l1_var^2))),
            RMSE_rel_bias_res2 = sqrt( (var(residvar_median)/ (l1_var^2) ) + rel_bias_res^2),
            RMSE_rel_bias_res = sqrt( (var(residvar_median) ) + rel_bias_res^2),
            
            RMSE_rel_bias_res_median = sqrt( (var(residvar_median) ) + rel_bias_res3^2),
            
            
            
            
            
            ## mean and sd L2 SMD 
            
            
            .groups = 'drop') %>% distinct()
  
  return(performance_measures)
}





