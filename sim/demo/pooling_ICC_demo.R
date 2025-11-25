rm(list=ls())

library(tidyverse)
library(nlme)
library(purrr)
library(lmeInfo)


# Source estimation function
source("sim/2_estimate_pooledICC.R")

# Read in data 
data1 <- read.table("sim/demo/data/g0706nd.txt",
                    sep="\t", header=TRUE, fill = TRUE, colClasses=c(rep("character", 114)))

data2 <- read.table("sim/demo/data/g0806nd.txt",
                    sep="\t", header=TRUE, fill = TRUE, colClasses=c(rep("character", 281)))

data3 <- read.table("sim/demo/data/g1006nd.txt",
                    sep="\t", header=TRUE, fill = TRUE, colClasses=c(rep("character", 281)))

data4 <- read.table("sim/demo/data/g0306nd.txt",
                    sep="\t", header=TRUE, fill = TRUE, colClasses=c(rep("character", 114)))

data5 <- read.table("sim/demo/data/g0406nd.txt",
                    sep="\t", header=TRUE, fill = TRUE, colClasses=c(rep("character", 114)))

data6 <- read.table("sim/demo/data/g0506nd.txt",
                    sep="\t", header=TRUE, fill = TRUE, colClasses=c(rep("character", 161)))

data7 <- read.table("sim/demo/data/g0606nd.txt",
                    sep="\t", header=TRUE, fill = TRUE, colClasses=c(rep("character", 114)))


# set empty cells as NA


all <- data2 |> 
  bind_rows(data1, data3, data4, data5,data6,data7) |>  
  mutate(cluster = schname) |> 
  filter(grade != "10") |> # removed grade 10
  mutate(sprp_sch = as.factor(sprp_sch),
         sprp_dis = as.factor(sprp_dis),
         stu_id = row_number())

all[all == " "] <- NA

# take a look at the data -- breakdown by grade
all |>  group_by(grade) |>  
  summarise(n_districts = length(unique(sprp_dis)), 
            n_districts2 = length(unique(district)), 
            n_schools = length(unique(sprp_sch)),
            n_schools2 = length(unique(schname)),
            students = sum(!is.na(erawsc)),
            students2 = sum(!is.na(mrawsc)),
            .groups = 'drop')




# function to estimate one-way random effects ICCs  
est_icc <- function(dat, outcome) {
  
  fixed_formula <- as.formula(paste(outcome, "~ 1"))
    
  
    mod <- lme(fixed = fixed_formula, random =  ~ 1 | sprp_sch, data = dat)  
    
    # obtaining variances
    info <- extract_varcomp(mod)
    clustervar <- as.numeric(info$Tau$sprp_sch) # cluster variance
    residvar <- as.numeric(info$sigma_sq) # residual variance
    
    # ICC estimate
    icc_est <- clustervar / (clustervar + residvar)
    
    # need var of clustervar for hedges ICC var formula. use fisher information matrix.
    I <- Fisher_info(mod, type = "expected")
    v2 <- diag(solve(I))[1]
    v1 <- diag(solve(I))[2]
    
    # Hedges
    var_hedges <- as.numeric(((1 - icc_est)^2 * v2) / (clustervar + residvar)^2)
    
    # mean sample size
    n <- as.numeric(dat |> 
                      group_by(sprp_sch) |> 
                      summarise(ni = length(stu_id), .groups = 'drop') |> 
                      summarise(mean = mean(ni), .groups = 'drop')
    )
    
    # number of clusters
    m <- length(unique(dat$sprp_sch))
    
    # Fisher
    var_fisher <- (2 * (1 - icc_est)^2 * (1 + (n - 1) * icc_est)^2) / (n * (n - 1) * (m - 1))
    
    # list of sample sizes 
    list <- mod$groups |> 
      group_by(sprp_sch) |> 
      tally() |> 
      select(n)
    
    W <- 1 + (list$n - 1) * icc_est
    V <- 1 + (list$n - 1) * icc_est^2
    
    # Total Sample Size
    N <- mod$dims$N 
    
    # Donner
    var_donner <- (2 * N * (1 - icc_est)^2) / (N * sum(list$n * (list$n - 1) * V * (W^(-2))) - (icc_est^2 * (sum(list$n * (list$n - 1) * (W^(-1))))^2))
    
    #N = length(uncond_data$y)
    
    # weighted mean cluster size
    n_0 <- (1/(m - 1)) * (N - sum(list$n^2 / N))
    
    # Swiger
    # uses weighted mean cluster size
    var_swiger <- (2 * (N - 1) * (1 - icc_est)^2 * (1 + (n_0 - 1) * icc_est)^2)/(n_0^2 * (N - m) * (m - 1))
    
    
    ## Smith
    pt1 <- (2 * ((1 - icc_est)^2)/(n_0^2))
    pt2a <- (((1 + (icc_est * (n_0 - 1)))^2)/(N - m))
    pt2b <- (m - 1) * (1 - icc_est) * (1 + (icc_est * (2 * n_0 - 1)))
    pt2c <- (icc_est^2) * (sum(list$n^2) - ((2 * (N^(-1))) * sum(list$n^3)) + ((N^(-2)) * ((sum(list$n^2))^2)))
    pt2d <- (m - 1)^2
    var_smith <- pt1 * ((pt2a + ((pt2b + pt2c) / pt2d)))
    
    
    ## Fisher transformed 
    icc_est_fisher_tf <- .5 * log((1 + (n_0 - 1) * icc_est) / (1 - icc_est))
    var_fisher_tf <- .5 * (((m - 1)^(-1)) + ((N - m)^(-1)))
    
    return(
      data.frame(
        icc_est = icc_est,
        icc_est_fisher_tf = icc_est_fisher_tf,
        var_hedges = var_hedges,
        var_donner = var_donner,
        var_fisher = var_fisher,
        var_swiger = var_swiger,
        var_smith = var_smith,
        var_fisher_tf = var_fisher_tf,
        clustervar = clustervar,
        residvar = residvar,
        clustervar_var = as.numeric(v2),
        residvar_var = as.numeric(v1), 
        n_0 =  n_0,
        n_arith = n
        )
      )
    
}  



#rve est
rve_estimation <- function(icc_est_sample, icc_value, var_icc_est, var_icc_name = NULL){
  tryCatch( { 
    
    if (!is.null(icc_est_sample)) {
      var_icc_est_call <- substitute(var_icc_est)
      icc_value_call <- substitute(icc_value)
      
      
      env <- list2env(icc_est_sample, parent = parent.frame())
      
      
      var_icc_est <- eval(var_icc_est_call, env)
      icc_value <- eval(icc_value_call, env)
      
      
    }
    
    
    
    rve_fit_robu <- robumeta::robu(formula = icc_value  ~ 1, studynum = study_id, var.eff.size = var_icc_est, data = icc_est_sample, small = TRUE)
    n_0 <-  mean(icc_est_sample$n_0, na.rm = TRUE)
    
    if(var_icc_name == "Fisher TF"){
      
      pooled_icc_est = (exp(2*as.numeric(rve_fit_robu$reg_table[2])) - 1)/(n_0 - 1 + exp(2*as.numeric(rve_fit_robu$reg_table[2])))
      ci_lower = (exp(2*as.numeric(rve_fit_robu$reg_table$CI.L)) - 1)/(n_0 - 1 + exp(2*as.numeric(rve_fit_robu$reg_table$CI.L)))
      ci_upper = (exp(2*as.numeric(rve_fit_robu$reg_table$CI.U)) - 1)/(n_0 - 1 + exp(2*as.numeric(rve_fit_robu$reg_table$CI.U)))
      
    } else {
      
      pooled_icc_est =  as.numeric(rve_fit_robu$reg_table[2])
      ci_lower = rve_fit_robu$reg_table$CI.L
      ci_upper = rve_fit_robu$reg_table$CI.U  
      
    }
    
    tau_est  = sqrt(rve_fit_robu$mod_info$tau.sq)
    
    return(results = tibble(method = "RVE",
                            pooled_icc_est = pooled_icc_est,
                            se_pooled_icc_est = as.numeric(rve_fit_robu$reg_table[3]),
                            var_icc_name = var_icc_name,
                            tau_est  = tau_est,
                            n_weighted = n_0
    ))
    
    
  }
  , error = function(w) { return(result = tibble(method = "RVE",
                                                     pooled_icc_est = NA,
                                                     se_pooled_icc_est = NA,
                                                     var_icc_name = var_icc_name,
                                                     tau_est  = NA,
                                                     n_weighted = NA
  ))})
  
  
}


#rma est
rma_estimation <- function(icc_est_sample, icc_value, var_icc_est, var_icc_name = NULL){
  
  tryCatch( { 
    
    
    if (!is.null(icc_est_sample)) {
      var_icc_est_call <- substitute(var_icc_est)
      icc_value_call <- substitute(icc_value)
      
      env <- list2env(icc_est_sample, parent = parent.frame())
      
      
      var_icc_est <- eval(var_icc_est_call, env)
      icc_value <- eval(icc_value_call, env)
      
      
    }
    
  
    
    reml_fit <- metafor::rma.uni(yi = icc_value, vi = var_icc_est, data = icc_est_sample, method = "REML", test="knha") # test
    
    
    # assumption -- used the mean of the n_weighted across all of the primary studies. 
    n_0 <-  mean(icc_est_sample$n_0, na.rm = TRUE)
    
    if(var_icc_name == "Fisher TF"){
      
      pooled_icc_est = (exp(2*as.numeric(reml_fit$beta)) - 1)/(n_0 - 1 + exp(2*as.numeric(reml_fit$beta)))
      
    } else {
      
      pooled_icc_est = as.numeric(reml_fit$beta)
    }
    
    
    tau_est  <- sqrt(reml_fit$tau2)
    
    
    return(results = tibble(method = "REML",
                                pooled_icc_est = pooled_icc_est,
                                se_pooled_icc_est = as.numeric(reml_fit$se),
                                var_icc_name = var_icc_name,
                                tau_est  = tau_est,
                                n_weighted = n_0
                                
    ))
    
  }, error = function(w) { return(result = tibble(method = "REML",
                                                      pooled_icc_est = NA,
                                                      se_pooled_icc_est = NA,
                                                      var_icc_name = var_icc_name,
                                                      tau_est = NA,
                                                      n_weighted = NA 

  ))})
  
  
}


#analysis function
analysis <- function(icc_est_sample){
  
  
 # possibly_rve <- possibly(.f = rve_estimation, otherwise = NULL)
  
 # possibly_reml <- possibly(.f = rma_estimation, otherwise = NULL)
  
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
  

  
  return(analysis_results)
  
}

# Graphs

library(patchwork) 

  nested_districts_eng <- all |> 
    mutate(erawsc = as.numeric(erawsc)) |> 
    group_by(district, grade) |> 
    mutate(num_school = length(unique(sprp_sch))) |> 
    filter(!is.na(erawsc)) |> 
    filter (num_school > 4) |>  
    nest() |> 
    dplyr::mutate(icc_dat = purrr::map(data, est_icc, outcome = "erawsc")) |>  
    select(-data) |> 
    unnest(icc_dat) |>
    mutate(study_id = district) |> # call district "study_id" to use the analysis() function
    nest(data = -grade) |>
    mutate(grade= paste("g_", grade, sep = ""),
           data = set_names(data, grade)) |>     
    dplyr::mutate(pooled = purrr::map(data, analysis))
  
  
  nested_districts_math <- all |> 
    mutate(mrawsc = as.numeric(mrawsc)) |> 
    group_by(district, grade) |> 
    mutate(num_school = length(unique(sprp_sch))) |> 
    filter(!is.na(mrawsc)) |> 
    filter (num_school > 4) |>  
    nest() |> 
    dplyr::mutate(icc_dat = purrr::map(data, est_icc, outcome = "mrawsc")) |>  
    select(-data) |> 
    unnest(icc_dat) |>
    mutate(study_id = district) |> # call district "study_id" to use the analysis() function
    nest(data = -grade) |>
    mutate(grade= paste("g_", grade, sep = ""),
           data = set_names(data, grade)) |>     
    dplyr::mutate(pooled = purrr::map(data, analysis))
  
  
  

  
  eng_dis <- nested_districts_eng |>
    select(-pooled) |>
    unnest(data) |>
    ggplot(aes(x = icc_est, colour = grade)) +
    geom_density() + theme_bw() + scale_colour_discrete(name = "Grade", labels = c("Grade 3", "Grade 4", "Grade 5", "Grade 6", "Grade 7", "Grade 8")) +
    labs(y = "Density", x = "ICC Estimate")
  
  math_dis <- nested_districts_math |>
    select(-pooled) |>
    unnest(data) |>
    ggplot(aes(x = icc_est, colour = grade)) +
    geom_density() + theme_bw() + scale_colour_discrete(name = "Grade", labels = c("Grade 3", "Grade 4", "Grade 5", "Grade 6", "Grade 7", "Grade 8")) +
    labs(y = "Density", x = "ICC Estimate")
  
  
  all_dis <- eng_dis + math_dis + plot_layout(axis_titles = "collect", guides = "collect")  & theme(legend.position = 'bottom')
  
  ggsave("sim/demo/demo_density.png", width = 6.8, height = 5, units = "in")
  
    
  nested_districts_eng |>
    select(-pooled) |>
    unnest(data) |>
    ungroup()|>
    select(grade, icc_est) |>
    group_by(grade) |>
    summarise(mean_icc = mean(icc_est),
           sd_icc = sd(icc_est),
           n_icc = n(), .groups = 'drop')
  
  nested_districts_math |>
    select(-pooled) |>
    unnest(data) |>
    ungroup()|>
    select(grade, icc_est) |>
    group_by(grade) |>
    summarise(mean_icc = mean(icc_est),
              sd_icc = sd(icc_est),
              n_icc = n(), .groups = 'drop')
  
  
 means_eng <- nested_districts_eng |>
    select(-pooled) |>
    unnest(data) |>
    ungroup()|>
    select(grade, icc_est) |>
    group_by(grade) |>
    summarise(mean_icc = mean(icc_est), 
              .groups = 'drop')
 
 means_math <- nested_districts_eng |>
   select(-pooled) |>
   unnest(data) |>
   ungroup()|>
   select(grade, icc_est) |>
   group_by(grade) |>
   summarise(mean_icc = mean(icc_est),
             .groups = 'drop')
  
  
  #nested_districts$pooled[[1]] |> select(method, pooled_icc_est, se_pooled_icc_est, var_icc_name)

  #dummy single dataset without intercept grade for each dummy coded. 
  
  grade.labs <- c("Grade 3", "Grade 4", "Grade 5", "Grade 6", "Grade 7", "Grade 8")
  names(grade.labs) <- c("g_03", "g_04", "g_05", "g_06", "g_07", "g_08")
  
  

  
  mean_eng_graph <- nested_districts_eng |>
    select(-data) |>
    unnest(pooled) |>
    ggplot(aes(y = pooled_icc_est, x = var_icc_name, group = var_icc_name, 
               shape = method, color=var_icc_name)) +
    geom_point(aes(shape = method, color = var_icc_name)) + 
    geom_hline(data = means_eng, aes(yintercept = mean_icc), color="red", linetype="dashed" ) +
    facet_wrap(~grade,
               labeller = labeller(grade = grade.labs))  +
    labs(shape = "Meta-Analytic Pooling Method", color = "Variance ICC", y = "Pooled ICC Estimate") + theme_bw()+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) 
  
  mean_math_graph <- nested_districts_math |>
    select(-data) |>
    unnest(pooled) |>
    ggplot(aes(y = pooled_icc_est, x = var_icc_name, group = var_icc_name)) +
    geom_point(aes(shape = method, color = var_icc_name)) + 
    geom_hline(data = means_math, aes(yintercept = mean_icc), color="red", linetype="dashed" ) +
    facet_wrap(~grade,
               labeller = labeller(grade = grade.labs))  +
    labs(shape = "Meta-Analytic Pooling Method", color = "Variance ICC", y = "Pooled ICC Estimate") + theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) 
  
  mean_all_graph <- mean_eng_graph + mean_math_graph + plot_layout(axis_titles = "collect", guides = "collect")  & theme(legend.position = 'bottom', legend.box="vertical")
  
  
  
  ggsave("sim/demo/demo_results_pooledICC.png", width = 6.8, height = 5, units = "in")
  
  
  