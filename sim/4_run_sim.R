rm(list=ls())
library(dplyr)
library(tidyr)
library(lmeInfo)
library(nlme)
library(purrr)
library(metafor)
library(robumeta)


source("sim/1_gen_sample_ICC.R")
source("sim/2_estimate_pooledICC.R")
source("sim/3_performance_crit.R")


# ----------------------------------------------------------------------



run_sim <- function(iterations,
                    icc_est_n, 
                    nj_size, 
                    n_bar_size, 
                    n_bar_prop, 
                    var_combo, 
                    tau,
                    seed = NULL,
                    summarize_results = FALSE){
  
  if (!is.null(seed)) set.seed(seed)
  

  
  
  results <-
    
    rerun(iterations,{
      
      # generate data --------------------------------------------------------- 
       dat <- gen_icc_unbalanced(icc_est_n = icc_est_n,
                                 nj_size = nj_size,
                                 n_bar_size = n_bar_size,
                                 n_bar_prop = n_bar_prop,
                                 var_combo = var_combo,
                                 tau = tau) 
       
       
       est_res <-  analysis(dat)
      
      
    }) %>% dplyr::bind_rows()
  
  
 if (summarize_results) {
   performance <- calc_performance(results = results)
   return(performance)
 } else {
   return(results)
 }
  
}



# ----------------------------------------------------------------------

set.seed(20251005) 



design_factors <- list(
  icc_est_n = c(20, 50, 100),
  nj_size = c("small", "large"),
  n_bar_size = c("small", "large"),
  n_bar_prop = c(.1, .5),
  tau= c(.01,.02),
#  var_combo = c("small_large", "medium_large", "large_large", "icc_0.9", "icc_0.5", "icc_0.1")
var_combo = c("icc_0.9", "icc_0.5", "icc_0.1")
)


batches <- 10
total_reps <- 1500


lengths(design_factors)

params <- expand.grid(c(design_factors, list(batch = 1:batches)))


params$iterations <- total_reps / batches
params$seed <- round(runif(1) * 2^30) + 1:nrow(params)

params |> group_by(seed) |> tally()


source_obj <- ls()

# ----------------------------------------------------------------------
batch_file <-  10

#batch_file <-  1
params2 <- params %>% filter(batch == batch_file)
params2$batch <- NULL

# params2 <- params2 |> slice(1:1)

library(future)
library(furrr)

run_date1 <- date()
run_date1


options(error=recover)

no_cores <- 4 - 1


plan("multisession", workers = no_cores) 



tm <- system.time(
  results <-
    params2 %>%
    mutate(res = future_pmap(., .f = run_sim, .options = furrr_options(seed=NULL))) %>%
    unnest(cols = res)
) 



tm 


# ----------------------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

FileName <- paste("sim/data/ICC_sim_results_update", batch_file,".Rdata",sep="")

save(tm, params2, results, session_info, run_date, file =  FileName)

