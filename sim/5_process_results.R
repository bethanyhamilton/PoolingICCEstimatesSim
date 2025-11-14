rm(list = ls())
library(tidyverse)
library(kableExtra)
library(ggplot2)

load("sim/data/ICC_sim_results_1.Rdata")
results1 <- results
load("sim/data/ICC_sim_results_2.Rdata")
results2 <- results
load("sim/data/ICC_sim_results_3.Rdata")
results3 <- results
load("sim/data/ICC_sim_results_4.Rdata")
results4 <- results
load("sim/data/ICC_sim_results_5.Rdata")
results5 <- results
load("sim/data/ICC_sim_results_6.Rdata")
results6 <- results
load("sim/data/ICC_sim_results_7.Rdata")
results7 <- results
load("sim/data/ICC_sim_results_8.Rdata")
results8 <- results
load("sim/data/ICC_sim_results_9.Rdata")
results9 <- results
load("sim/data/ICC_sim_results_10.Rdata")
results10 <- results


load("sim/data/ICC_sim_results_update1.Rdata")
results11 <- results
load("sim/data/ICC_sim_results_update2.Rdata")
results12 <- results
load("sim/data/ICC_sim_results_update3.Rdata")
results13 <- results
load("sim/data/ICC_sim_results_update4.Rdata")
results14 <- results
load("sim/data/ICC_sim_results_update5.Rdata")
results15 <- results
load("sim/data/ICC_sim_results_update6.Rdata")
results16 <- results
load("sim/data/ICC_sim_results_update7.Rdata")
results17 <- results

results <- bind_rows(results1, results2, results3, results4, results5, results6, results7, results8, results9, results10, results11, results12, results13, results14, results15, results16, results17)

rm(results1, results2, results3, results4, results5, results6, results7, results8, results9, results10, results11, results12, results13, results14, results15, results16, results17)

source("sim/3_performance_crit.R")

# convergence
test <- results |>  
  filter(sample_conv == "no") |>  
  filter(all_converge == TRUE)
#sum(is.na(results$pooled_icc_est))


results |>  filter(is.na(pooled_icc_est)) |>  
  group_by(nj_size) |>  
  tally() |> 
  knitr::kable() |>
  kableExtra::kable_styling() 


converge_sum <- results |>  
  group_by(method, var_icc_name, icc_est_n, nj_size, 
           n_bar_size, n_bar_prop, var_combo, tau)  |> 
  summarise(k = n(),
            num_no_converge = sum(is.na(pooled_icc_est)),
            num_converge = k - num_no_converge,
            num_converge2 = sum(!is.na(pooled_icc_est)),
            converge_rate = (num_converge/k)*100,
            nonconverge_rate = (num_no_converge/k)*100, 
            .groups = "drop"
  ) |>  
  ungroup()



#min(converge_sum$converge_rate)
#which(converge_sum$converge_rate < 100)


# number of samples that did not converge per method.
converge_sum |>
  group_by(method, var_icc_name) |>
  summarise(total_non_convergence = sum(num_no_converge, na.rm = TRUE),
            total_converged = sum(num_converge, na.rm= TRUE),
            converge_rate = (total_converged/(total_non_convergence + total_converged))*100, .groups = "drop" ) |>
  ungroup() |> 
  knitr::kable(digits = 3) |>
  kableExtra::kable_styling() 


results <- results |>  
  group_by(method, var_icc_name, icc_est_n, nj_size, n_bar_prop, n_bar_size, tau, var_combo)  |> 
  mutate(combination_id = cur_group_id()) |> 
  ungroup()





results |>  group_by(combination_id) |>  tally()

results |> filter(sample_conv != "no") |> filter(all_converge == TRUE) |>  group_by(combination_id) |>  tally() |> filter(n < 1000)


results_all_converged <- results |>  
  filter(sample_conv != "no") |>  
  filter(all_converge == TRUE) |>  
  group_by(combination_id) |>  
  slice(1:1000) #|> ungroup()





# calculate performance criteria
perf_crit <- calc_performance(results_all_converged)

perf_crit2 <- calc_performance_var(results_all_converged)



# add in labels for graphs

## raw converged results

# results_all_converged <- results_all_converged |> 
#   tidyr::unite("method2",  c(method, var_icc_name), remove = FALSE) |> 
#   mutate(method2 = str_replace(method2, "_", "\n"))
# 
# results_all_converged <- results_all_converged  |>  mutate(
#   var_combo_graph = case_when(
#     var_combo == 'small_large' ~ 0.05,
#     var_combo == 'icc_0.1' ~ 0.1,
#     var_combo == 'medium_large' ~ 0.15,
#     var_combo == 'large_large' ~ 0.25,
#     var_combo == 'icc_0.5' ~ 0.5,
#     var_combo == 'icc_0.9' ~ 0.9
#   ),
#   nj_size = case_when(
#     nj_size == "small" ~'U[30, 50]',
#     
#     nj_size == "large" ~'U[50, 100]',
#   ),
#   n_bar_size = case_when(
#     n_bar_size == "small" ~'U[10, 30]',
#     
#     n_bar_size == "large" ~'U[30, 50]',
#   ),
#   n_bar_graph = paste("n_bar_size = ", n_bar_size),
#   nj_graph = paste("nj_size = ", nj_size),
#   icc_est_n_graph = paste("k = ", icc_est_n)
# ) |>  
#   separate(method2, into = c("model", "variance"), sep = "\n", remove = FALSE 
#   )  |>  
#   mutate(variance = ifelse( variance == "Fisher TF", "Fisher\nTF", variance))
# 
# 
# results_all_converged$icc_est_n_graph <- factor(results_all_converged$icc_est_n_graph,levels=c(
#   "k =  20",
#   "k =  50",
#   "k =  100"))

## performance criteria of the pooled ICC

perf_crit <- perf_crit |> 
  tidyr::unite("method2",  c(method, var_icc_name), remove = FALSE) |> 
  mutate(method2 = str_replace(method2, "_", "\n"))

perf_crit$var_combo_graph <- factor(perf_crit$var_combo,labels=c(
  'icc_0.9' = 0.9,
  'icc_0.5' = 0.5,
  'large_large' = 0.25,
  'medium_large' = 0.15,
  'icc_0.1' = 0.1,
  'small_large' = 0.05)) 

perf_crit <- perf_crit  |>  mutate(
  var_combo_graph = case_when(
    var_combo == 'small_large' ~ 0.05,
    var_combo == 'icc_0.1' ~ 0.1,
    var_combo == 'medium_large' ~ 0.15,
    var_combo == 'large_large' ~ 0.25,
    var_combo == 'icc_0.5' ~ 0.5,
    var_combo == 'icc_0.9' ~ 0.9
  ),
  nj_size = case_when(
    nj_size == "small" ~'U[30, 50]',
    
    nj_size == "large" ~'U[50, 100]',
  ),
  n_bar_size = case_when(
    n_bar_size == "small" ~'U[10, 30]',
    
    n_bar_size == "large" ~'U[30, 50]',
  ),
  n_bar_graph = paste("n_bar_size = ", n_bar_size), 
  nj_graph = paste("nj_size = ", nj_size), 
  icc_est_n_graph = paste("k = ", icc_est_n)
) |>  
  separate(method2, into = c("model", "variance"), sep = "\n", remove = FALSE 
  )  |>  
  mutate(variance = ifelse( variance == "Fisher TF", "Fisher\nTF", variance))


perf_crit$icc_est_n_graph <- factor(perf_crit$icc_est_n_graph,levels=c(
  "k =  20",
  "k =  50",
  "k =  100"))

perf_crit$var_combo_graph_f <- factor(as.character(perf_crit$var_combo_graph), levels = c("0.9", "0.5", "0.25", "0.15", "0.1", "0.05"))



## performance criteria of the variance components

perf_crit2$var_combo_graph <- factor(perf_crit2$var_combo,labels=c(
  'icc_0.9' = 0.9,
  'icc_0.5' = 0.5,
  'large_large' = 0.25,
  'medium_large' = 0.15,
  'icc_0.1' = 0.1,
  'small_large' = 0.05)) 

perf_crit2 <- perf_crit2  |>  mutate(
  var_combo_graph = case_when(
    var_combo == 'small_large' ~ 0.05,
    var_combo == 'icc_0.1' ~ 0.1,
    var_combo == 'medium_large' ~ 0.15,
    var_combo == 'large_large' ~ 0.25,
    var_combo == 'icc_0.5' ~ 0.5,
    var_combo == 'icc_0.9' ~ 0.9
  ),
  nj_size = case_when(
    nj_size == "small" ~'U[30, 50]',
    
    nj_size == "large" ~'U[50, 100]',
  ),
  n_bar_size = case_when(
    n_bar_size == "small" ~'U[10, 30]',
    
    n_bar_size == "large" ~'U[30, 50]',
  ),
  n_bar_graph = paste("n_bar_size = ", n_bar_size), 
  nj_graph = paste("nj_size = ", nj_size), 
  icc_est_n_graph = paste("k = ", icc_est_n)
) 

perf_crit2$icc_est_n_graph <- factor(perf_crit2$icc_est_n_graph,levels=c(
  "k =  20",
  "k =  50",
  "k =  100"))


perf_crit2$var_combo_graph_f <- factor(as.character(perf_crit2$var_combo_graph), levels = c("0.9", "0.5", "0.25", "0.15", "0.1", "0.05"))



# save results

#save(results_all_converged, file = "sim/data/results_conv.RData")

save(perf_crit, file = "sim/data/results_pooled_est.RData")

save(perf_crit2, file = "sim/data/results_var_comp.RData")


