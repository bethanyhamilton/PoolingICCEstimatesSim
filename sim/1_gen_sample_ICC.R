# library(nlme)
# library(lmeInfo)
# library(purrr)
# library(dplyr)
# library(tidyr)

# Generates raw data
## data from unconditonal model 
gen_uncond_data <- function(nj, l1_var, l2_var, n_bar = NULL, n_bar_prop = NULL, tau){
  
  icc_true = l2_var/(l2_var+l1_var)
  
  # between-study heterogeneity
  

  
  lower = (-2*tau) #
  upper = (2*tau)
  
  u <-  runif(1, min = lower, max = upper)
  

  # add in check to make sure it is not less 0. throw away sample and make note.. this sample was below 0.
  icc_est <- icc_true + u
  
  # find new L2 by using fixed L1. So, when looking at RPB for L2 will this be factored in?
  l2_var_new= l1_var*icc_est/(1-icc_est)
  

  stopifnot( n_bar_prop >= 0 && n_bar_prop < 1 )
  
  # generate site sizes 
  n_min = ceiling( n_bar * (1 - n_bar_prop) )
  n_max = floor( n_bar * (1 + n_bar_prop) )
  # num_ind <- sample( n_min:n_max, nj, replace=TRUE )
  # num_ind <- round(runif(nj, min = n_min, max = n_max))
  
  if (n_min < n_max) {
    num_ind <- sample( n_min:n_max, nj, replace=TRUE )
  } else {
    num_ind <- rep(n_bar,nj)
  }
  
  # cluster and individual ID values
  uncond_data <- data.frame(cluster = rep(seq(1,nj), times = num_ind))
  
  # Obtaining individual and cluster residuals
  uncond_data <- uncond_data %>%
    group_by(cluster) %>%
    mutate(individual = row_number(),
           cluster_resid = rnorm(1, 0, sqrt(l2_var_new)),
           ind_resid = rnorm(length(individual), 0, sqrt(l1_var)),
           y = cluster_resid + ind_resid)
          
  uncond_data$u <- u
  uncond_data$L2_new <- l2_var_new
  
  
  return(uncond_data)
}


# Estimate ICC
est_icc <- function(uncond_data) {
  
  
  u = unique(uncond_data$u)
  L2_new = unique(uncond_data$L2_new)
  
  
  possibly_model <- possibly(.f = lme, otherwise = NULL)
  
  
  model <- possibly_model(fixed = y ~ 1, random =  ~ 1 | cluster, data = uncond_data) #uncond
  
  if(is.null(model)){
    converged = "no"
    
  }else{
    converged = "yes"
    model <- lme(fixed = y ~ 1, random =  ~ 1 | cluster, data = uncond_data) #uncond
    
    
  }

  if(converged == "yes"){
    
  
  ## obtaining variances
  info <-  extract_varcomp(model)
  clustervar <- as.numeric(info$Tau$cluster) # cluster variance
  residvar <- as.numeric(info$sigma_sq) # residual variance
  
  ## ICC estimate
  icc_est <- clustervar / (clustervar + residvar)
  
  ## need se of clustervar for hedges ICC var formulae. use fisher information matrix.
  I <- Fisher_info(model, type = "expected")
  v2 <- diag(solve(I))[1]
  v1 <- diag(solve(I))[2]
  
  ## Hedges
  var_hedges <- as.numeric(((1 - icc_est) ^ 2 * v2) / (clustervar + residvar) ^ 2)
  
  ## Fishers
  n <- as.numeric(
    uncond_data %>% group_by(cluster) %>%
      summarise(ni = length(individual), .groups = 'drop') %>%
      summarise(mean = mean(ni), .groups = 'drop')
  )
  
  m <- length(unique(uncond_data$cluster))
  var_fisher <- (2 * (1 - icc_est) ^ 2 * (1 + (n - 1) * icc_est) ^ 2) / (n * (n - 1) * (m - 1))
  
  ## Donner
  list <- model$groups %>% group_by(cluster) %>% tally() %>% select(n)
  
  W <- 1 + (list$n -1)*icc_est
  V <- 1+ (list$n -1)*icc_est^2
  
  var_donner <- (2 * model$dims$N * (1 - icc_est) ^ 2) /
     (model$dims$N * sum(list$n*(list$n-1)*V*(W^(-2)))-(icc_est^2*(sum(list$n*(list$n-1)*(W^(-1))))^2))
  

  ## Smith
  N = length(uncond_data$y)
  #weighted mean cluster size
  n_0 = (1/(m -1))*(N - sum(list$n^2/N))
  
  
  ## Swiger
  # weighted mean cluster size not arithmetic mean cluster size 
  var_swiger = (2*(N-1)*(1-icc_est)^2*(1+(n_0-1)*icc_est)^2)/(n_0^2*(N-m)*(m-1))
  
  
  ## Smith
  
  pt1 = (2*((1-icc_est)^2)/(n_0^2))
  
  pt2a = (((1+(icc_est*(n_0-1)))^2)/(N-m))
  pt2b =(m-1)*(1-icc_est)*(1+(icc_est*(2*n_0-1)))
  pt2c =(icc_est^2)*(sum(list$n^2)-((2*(N^(-1)))*sum(list$n^3))+((N^(-2))*((sum(list$n^2))^2)))
  pt2d =(m-1)^2
  
  var_smith = pt1 * ((pt2a + ((pt2b + pt2c)/pt2d)))
  
  
  ## fisher transformed 
  
  icc_est_fisher_tf = .5*log((1+(n_0-1)*icc_est)/(1-icc_est))
  
  
  var_fisher_tf = .5*( ((m-1)^(-1)) + ((N-m)^(-1)) )
  
  return(
    data.frame(
      icc_est = icc_est,
      u =u,
      L2_new = L2_new,
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
      n_arith = n,
      converged = "yes"
    )
  )
  
  }else{
    
    return(
      data.frame(
        icc_est = NA,
        u=NA,
        L2_new = NA,
        icc_est_fisher_tf = NA,
        var_hedges = NA,
        var_donner = NA,
        var_fisher = NA,
        var_swiger = NA,
        var_smith = NA,
        var_fisher_tf = NA,
        clustervar = NA,
        residvar = NA,
        clustervar_var = NA,
        residvar_var = NA, 
        n_0 =  NA,
        n_arith = NA,
        converged = "no"
      ))
    
  }
  
}



### Generates data sets repeatedly and estimated ICC for each data set. returns only ICC estimates. how many do we create for a population??
gen_icc <- function(icc_est_n, nj, n_bar, n_bar_prop, var_combo, tau){
  # different conditions for level 1 and level 2 variances
  if(var_combo == "small_large"){
    l1_var <- 95
    l2_var <- 5
  }else if(var_combo == "medium_large"){
    l1_var <- 85
    l2_var <- 15
  }else{
    l1_var <- 75
    l2_var <- 25
  }

  icc_est <- replicate(icc_est_n,
                       {

                         possibly_gen_uncond_data <- possibly(.f = gen_uncond_data, otherwise = NULL)
                         possibly_est_icc <- possibly(.f = est_icc, otherwise = NULL)

                         dat <- possibly_gen_uncond_data(nj = nj, n_bar = n_bar, n_bar_prop = n_bar_prop, l1_var = l1_var, l2_var = l2_var)

                      #   dat <- gen_uncond_data(nj, n_bar, n_bar_prop, l1_var, l2_var) # or gen_cond_data
                         cbind(possibly_est_icc(uncond_data = dat, tau= tau),nj,  n_bar)
                       }, simplify = FALSE) %>%
    dplyr::bind_rows() %>% mutate(study_id = row_number(), l1_var = l1_var, l2_var = l2_var)



  return(icc_est)
}





gen_icc_unbalanced <- function(icc_est_n, nj_size, n_bar_size, n_bar_prop, var_combo, tau){
  if(nj_size == "small"){
    n_min <- 30
    n_max <- 50
  }else{
    n_min <- 50
    n_max <- 100
  }
  
  
  if(n_bar_size == "small"){
    n_bar_min <- 10
    n_bar_max <- 30
  }else{
    n_bar_min <- 30
    n_bar_max <- 50
  }
  

  
  # different conditions for level 1 and level 2 variances
  if(var_combo == "small_large"){
    l1_var <- 95
    l2_var <- 5
  }else if(var_combo == "medium_large"){
    l1_var <- 85
    l2_var <- 15
  }else{
    l1_var <- 75
    l2_var <- 25
  }
  
  
  icc_est <- replicate(icc_est_n,
                       { nj <- sample( n_min:n_max, 1, replace=TRUE )
                         
                         n_bar <- sample( n_bar_min:n_bar_max, 1, replace=TRUE )

                         dat <- gen_uncond_data(nj = nj, n_bar = n_bar, n_bar_prop = n_bar_prop, l1_var = l1_var, l2_var = l2_var, tau = tau)
                         
                         cbind(est_icc(uncond_data = dat), nj,  n_bar,  tau)
                         
                       }, simplify = FALSE) %>%
 bind_rows() %>% mutate(study_id = row_number(), l1_var = l1_var, l2_var = l2_var)
  
  
  
  return(icc_est)
  
  
  
}

### Test

# tm1 <- system.time( test_uncond_data <- gen_uncond_data(nj = 34, n_bar = 15 , n_bar_prop = .5, l1_var = 95, l2_var = 5, tau=.01) )
# tm1


### Test Estimating an ICC function
# tm2  <- system.time(test <-  est_icc(test_uncond_data))
# tm2




### Test Generating ICC functions




# tm3 <-
#   system.time(
#     ICC_test_dist2 <-
#       gen_icc(
#         icc_est_n = 100,
#         nj = 100,
#         n_bar = 50,
#         n_bar_prop = .5,
#         var_combo = "medium_large",
#         tau = sqrt(.01)
#       )
#   )
# 
# hist(ICC_test_dist2$icc_est)
# 
# 
# 
#  tm4 <-
#    system.time(
# ICC_test_dist <- gen_icc_unbalanced(
#                    icc_est_n= 150,
#                    nj_size = "small",
#                    n_bar_size = "small",
#                    n_bar_prop= .5,
#                    var_combo= "small_large",
#                    tau= .01)      )
#  tm4 <-
#    system.time(
# ICC_test_dist_neg <- gen_icc_unbalanced(
#                    icc_est_n= 150,
#                    nj_size = "small",
#                    n_bar_size = "small",
#                    n_bar_prop= .5,
#                    var_combo= "medium_large",
#                    tau= .1)      )

#hist(ICC_test_dist$icc_est)
# 
#  tm3
#  tm4
