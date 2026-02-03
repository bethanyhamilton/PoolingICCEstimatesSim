#-------------------------------------------------------------------------------
# Estimation Functions
#-------------------------------------------------------------------------------
# library(metafor)
# library(clubSandwich)
# library(robumeta)
# library(dplyr)
# library(nlme)
# library(lmeInfo)
# library(purrr)
# library(tidyr)
# source("sim/1_gen_sample_ICC.R")



rve_estimation <- function(icc_est_sample, icc_value, var_icc_est, var_icc_name = NULL){
  tryCatch( { 

    if (!is.null(icc_est_sample)) {
      var_icc_est_call <- substitute(var_icc_est)
      icc_value_call <- substitute(icc_value)


      env <- list2env(icc_est_sample, parent = parent.frame())


      var_icc_est <- eval(var_icc_est_call, env)
      icc_value <- eval(icc_value_call, env)


    }
    
    
    # Number missing iccs
    test2 <- with(icc_est_sample,!is.na(icc_value))
    all_converge <- all(test2)
    num_no_cong_icc = sum(test2==FALSE)
    
    rve_fit_robu <- robumeta::robu(formula = icc_value  ~ 1, studynum = study_id, var.eff.size = var_icc_est, data = icc_est_sample, small = TRUE)
    n <-  mean(icc_est_sample$n_0, na.rm = TRUE)
    
    if(var_icc_name == "Fisher TF"){
      
      pooled_icc_est = (exp(2 *as.numeric(rve_fit_robu$reg_table[2])) - 1)/(n - 1 + exp(2*as.numeric(rve_fit_robu$reg_table[2])))
      
    } else {
      
      pooled_icc_est =  as.numeric(rve_fit_robu$reg_table[2])
    }
    
    tau_est  = sqrt(rve_fit_robu$mod_info$tau.sq)
    
    return(results = data.frame(method = "RVE",
                                pooled_icc_est = pooled_icc_est,
                                se_pooled_icc_est = as.numeric(rve_fit_robu$reg_table[3]),
                                var_icc_name = var_icc_name,
                                tau_est  = tau_est,
                                n_weighted = n, 
                                all_converge = all_converge,
                                num_no_cong_icc = num_no_cong_icc
    ))
    
    
  }
  , error = function(w) { return(result = data.frame(method = "RVE",
                                                     pooled_icc_est = NA,
                                                     se_pooled_icc_est = NA,
                                                     var_icc_name = var_icc_name,
                                                     tau_est  = NA,
                                                     n_weighted = NA, 
                                                     all_converge = all_converge,
                                                     num_no_cong_icc = num_no_cong_icc
  ))})
  
  
}



rma_estimation <- function(icc_est_sample, icc_value, var_icc_est, var_icc_name = NULL){
  
  tryCatch( { 
    
    
    if (!is.null(icc_est_sample)) {
      var_icc_est_call <- substitute(var_icc_est)
      icc_value_call <- substitute(icc_value)

      env <- list2env(icc_est_sample, parent = parent.frame())


      var_icc_est <- eval(var_icc_est_call, env)
      icc_value <- eval(icc_value_call, env)


    }
    

    # Number missing iccs
    test2 <- with(icc_est_sample,!is.na(icc_value))
    all_converge <- all(test2)
    num_no_cong_icc = sum(test2==FALSE)
    
    reml_fit <- metafor::rma.uni(yi = icc_value, vi = var_icc_est, data = icc_est_sample, method = "REML", test="knha") # test
    

    # assumption -- used the mean of the n_weighted across all of the primary studies. 
    n <-  mean(icc_est_sample$n_0, na.rm = TRUE)
    
    if(var_icc_name == "Fisher TF"){
      
      pooled_icc_est = (exp(2*as.numeric(reml_fit$beta)) - 1)/(n - 1 + exp(2*as.numeric(reml_fit$beta)))
      
    } else {
      
      pooled_icc_est = as.numeric(reml_fit$beta)
    }
    
    
    tau_est  <- sqrt(reml_fit$tau2)
    
    
    # extra that I may add in at some point:
    
    # weights <- 1 / (var_icc_est + tau_est^2)
    # se_alt <- sqrt(1 / sum(weights))
    # var_KH <- (1 / (length(weights) - 1)) * sum(weights * (icc_value - pooled_icc_est)^2)
    # se_khna <- sqrt(var_KH / sum(weights))
    # CI_k_lower <- pooled_icc_est - (qt(p = .975, df = (k - 1)) * se_khna)
    # CI_k_upper <- pooled_icc_est + (qt(p = .975, df = (k - 1)) * se_khna)
    # reml_fit$ci.lb
    # reml_fit$ci.ub
    # reml_fit$QE
    # reml_fit$QEp  # to look at power of the test
    # tau_sq_ci <- confint(reml_fit, type="PL")
    # tau_sq_ci_lower <- tau_sq_ci$random[1,2]
    # tau_sq_ci_upper <- tau_sq_ci$random[1,3]
    
    return(results = data.frame(method = "REML",
                                pooled_icc_est = pooled_icc_est,
                                se_pooled_icc_est = as.numeric(reml_fit$se),
                                var_icc_name = var_icc_name,
                                tau_est  = tau_est,
                                n_weighted = n, 
                                all_converge = all_converge,
                                num_no_cong_icc = num_no_cong_icc
                                
    ))
    
  }, error = function(w) { return(result = data.frame(method = "REML",
                                                     pooled_icc_est = NA,
                                                     se_pooled_icc_est = NA,
                                                     var_icc_name = var_icc_name,
                                                     tau_est = NA,
                                                     n_weighted = NA, 
                                                     all_converge = all_converge,
                                                     num_no_cong_icc = num_no_cong_icc
  ))})
  
  
}



analysis <- function(icc_est_sample){

  #possibly_rve_estimation <- possibly(.f = rve_estimation, otherwise = NULL)
  #possibly_rma_estimation <- possibly(.f = rma_estimation, otherwise = NULL)
  
  analysis_results <-   dplyr::bind_rows( 
    rve_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est, var_icc_est = var_hedges, var_icc_name = "Hedges"),
    rve_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est, var_icc_est = var_donner, var_icc_name = "Donner"),
    rve_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est, var_icc_est = var_fisher, var_icc_name = "Fisher"),
    rve_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est, var_icc_est = var_swiger, var_icc_name = "Swiger"),
    rve_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est, var_icc_est = var_smith, var_icc_name = "Smith"),
    rve_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est_fisher_tf,  var_icc_est = var_fisher_tf, var_icc_name = "Fisher TF"),
    
    
    rma_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est, var_icc_est = var_hedges, var_icc_name = "Hedges"),
    rma_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est, var_icc_est = var_donner, var_icc_name = "Donner"),
    rma_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est, var_icc_est = var_fisher, var_icc_name = "Fisher"),
    rma_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est, var_icc_est = var_swiger, var_icc_name = "Swiger"),
    rma_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est, var_icc_est = var_smith, var_icc_name = "Smith"),
    rma_estimation(icc_est_sample = icc_est_sample, icc_value = icc_est_fisher_tf,  var_icc_est = var_fisher_tf, var_icc_name = "Fisher TF"))
  
  # estimated l1 and l2 variances -- mean and sd
  analysis_results$L2_new_mean <- mean( icc_est_sample$L2_new)
  analysis_results$L2_new_sd <- sd( icc_est_sample$L2_new)
  analysis_results$L2_new_min <- min( icc_est_sample$L2_new)
  analysis_results$L2_new_max <- max( icc_est_sample$L2_new)
  analysis_results$L2_new_median <- median( icc_est_sample$L2_new)
  
  
  
  analysis_results$clustervar_mean <- mean( icc_est_sample$clustervar)
  analysis_results$clustervar_sd <- sd( icc_est_sample$clustervar)
  analysis_results$clustervar_min <- min( icc_est_sample$clustervar)
  analysis_results$clustervar_max <- max( icc_est_sample$clustervar)
  analysis_results$clustervar_median <- median( icc_est_sample$clustervar)
  
  
  analysis_results$residvar_mean <-  mean( icc_est_sample$residvar)
  analysis_results$residvar_sd <-  sd( icc_est_sample$residvar)
  analysis_results$residvar_min <- min( icc_est_sample$residvar)
  analysis_results$residvar_max <- max( icc_est_sample$residvar)
  analysis_results$residvar_median <- median( icc_est_sample$residvar)
 
  analysis_results$clustervar <- list( icc_est_sample$clustervar)
  analysis_results$residvar <- list( icc_est_sample$residvar)
  analysis_results$L2_new <- list( icc_est_sample$L2_new)
  
  analysis_results$l1_var <- unique(icc_est_sample$l1_var)
  analysis_results$l2_var <- unique(icc_est_sample$l2_var)
  
  if(sum(is.na(analysis_results$pooled_icc_est)) > 0) {
    
    analysis_results$sample_conv <- "no"
    
  }else{
    
    analysis_results$sample_conv <- "yes"
    
  }
  
  
  return(analysis_results)
  
}

# test
# 
# ICC_test_dist <- gen_icc_unbalanced(
#                    icc_est_n= 30,
#                    nj_size = "small",
#                    n_bar_size = "small",
#                    n_bar_prop= .5,
#                    var_combo= "small_large",
#                    tau= .01)
# # 
#  reml_est <- rma_estimation(ICC_test_dist, icc_value= icc_est, var_icc_est = var_hedges, var_icc_name = "Hedges")
#  reml_fit <- metafor::rma.uni(yi =  icc_est, vi = var_hedges, data = ICC_test_dist, method = "REML", test="knha") # test
# # 
# # all.equal(as.numeric(reml_fit$beta), unlist(reml_est[[2]]))
# # all.equal(as.numeric(reml_fit$se), unlist(reml_est[[3]]))
# 
# 
# 
# weights <- 1 / (ICC_test_dist$var_hedges + reml_fit$tau2)
# se_alt <- sqrt(1 / sum(weights))
# var_KH <- (1 / (length(weights) - 1)) * sum(weights * (ICC_test_dist$icc_est - reml_fit$beta)^2)
# se_khna <- sqrt(var_KH / sum(weights))
# CI_k_lower <- reml_fit$beta - (qt(p = .975, df = (length(weights) - 1)) * se_khna)
# CI_k_upper <- reml_fit$beta + (qt(p = .975, df = (length(weights) - 1)) * se_khna)
#  
# all.equal(as.numeric(CI_k_lower), reml_fit$ci.lb)
# all.equal(as.numeric(CI_k_upper), reml_fit$ci.b)

# 
# #
  # rve_est <- rve_estimation(icc_est_sample = ICC_test_dist, icc_value = icc_est, var_icc_est = var_hedges, var_icc_name = "Hedges")
  # rve_fit <- robumeta::robu(formula = icc_est  ~ 1, studynum = study_id, var.eff.size = var_hedges, data = ICC_test_dist, small = TRUE)
  # res.CA  <- rma(icc_est, var_hedges, method="HE", data=ICC_test_dist)
  # res.CA2 <- rma(icc_est, var_hedges, method="GENQ", weights=1/(var_hedges + res.CA$tau2), data=ICC_test_dist)
  # 
  # rve_rma_ce <- rma.mv(icc_est, var_hedges, random = ~ icc_est | study_id,  method="MoM", struct="UN", data=ICC_test_dist)
  # 
  # coef_RVE1 <-  metafor::robust(reml_fit, cluster = study_id, clubSandwich = TRUE)
  # 
  # coef_RVE2 <-  metafor::robust(res.CA, cluster = study_id, clubSandwich = TRUE)
#
#
#  all.equal(as.numeric(rve_fit$reg_table[2]), unlist(rve_est[[2]]))
#  all.equal(as.numeric(rve_fit$reg_table[3]), unlist(rve_est[[3]]))
  # all.equal(as.numeric(rve_fit$reg_table[2]), as.numeric(res.CA2$beta))
  # all.equal(as.numeric(rve_fit$reg_table[2]), as.numeric(coef_RVE1$beta))
  # all.equal(as.numeric(rve_fit$reg_table[2]), as.numeric(coef_RVE2$beta))
#  
 #as.numeric(rve_fit$reg_table[2]) - (qt(p = .975, df = (30 - 1)) * as.numeric(rve_fit$reg_table[3]))
 #as.numeric(rve_fit$reg_table[2]) - (qt(p = .975, df = 23.5) * as.numeric(rve_fit$reg_table[3]))


### Test

#tm5 <- system.time(test <- analysis(ICC_test_dist) )

 
# test <- analysis(ICC_test_dist) 
 
 # ICC_test_dist_neg <- ICC_test_dist |>  mutate(rownum = row_number(),
 #                                           icc_est = ifelse(rownum == 20, icc_est*-1, icc_est),
 #                                           icc_est_fisher_tf = ifelse(rownum == 20,icc_est, icc_est_fisher_tf),
 #                                           var_hedges = ifelse(rownum == 20,NA, var_hedges),
 #                                           var_donner = ifelse(rownum == 20,NA, var_donner),
 #                                           var_fisher = ifelse(rownum == 20,NA, var_fisher),
 #                                           var_swiger = ifelse(rownum == 20,NA, var_swiger),
 #                                           var_smith = ifelse(rownum == 20,NA, var_smith),
 #                                           var_fisher_tf = ifelse(rownum == 20,NA, var_fisher_tf),
 # 
 #                                           icc_est = ifelse(rownum == 22, icc_est*-1, icc_est),
 #                                           icc_est_fisher_tf = ifelse(rownum == 22,icc_est, icc_est_fisher_tf),
 #                                           var_hedges = ifelse(rownum == 22,NA, var_hedges),
 #                                           var_donner = ifelse(rownum == 22,NA, var_donner),
 #                                           var_fisher = ifelse(rownum == 22,NA, var_fisher),
 #                                           var_swiger = ifelse(rownum == 22,NA, var_swiger),
 #                                           var_smith = ifelse(rownum == 22,NA, var_smith),
 #                                           var_fisher_tf = ifelse(rownum == 22,NA, var_fisher_tf)
 # 
 #                                           ) |>  select(-rownum)
  # tm6 <- system.time(test2 <- analysis(ICC_test_dist_neg) )

