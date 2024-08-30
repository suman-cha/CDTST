rm(list=ls())
setwd("~/WORKSPACE/PROJECTS/CDT/simulation/Ours")  
suppressPackageStartupMessages({
  library(MASS)
  library(glmnet)
  library(ANN2)
  library(CVST)
  library(kernlab)
  library(pbapply)
  library(ggplot2)
  library(tmvtnorm)
  library(RCIT)
})
source("all_tests.R")
set.seed(1203)

# Function to generate covariates using truncated multivariate normal distributions
generate_data <- function(n, p, group) {
  mu <- c(1, 1, -1, -1, rep(0, p - 4))
  sigma <- diag(1, p)

  lb <- rep(-.5, p)
  ub <- rep(.5, p)

  if (group == 1) {
    x <- rtmvnorm(n, mean = mu, sigma = sigma, lower = lb, upper = ub, algorithm="gibbs")
  } else if (group == 2) {
    x <- rtmvnorm(n, mean = rep(0, p), sigma = sigma, lower = lb, upper = ub, algorithm="gibbs")
  }
  return(x)
}

# generate_data <- function(n, p, group, a = 0.1, b = 0.9) {
#   if (group == 1) {
#     x <- matrix(rbeta(n * p, shape1 = 4, shape2 = 2), n, p)
#   } else if (group == 2) {
#     x <- matrix(rbeta(n * p, shape1 = 4, shape2 = 6), n, p)
#   }
# 
#   # Apply truncation
#   x <- pmax(a, pmin(b, x))
#   x <- matrix(x, n, p)
#   return(x)
# }

generate_y <- function(x, is_null = TRUE, sigma = 2) {
  n <- nrow(x)
  epsilon <- rt(n, df = sigma)
  f0 <- x %*% c(1, -1, -1, 1, rep(0, dim(x)[2] - 4))
  mean_shift <- if (is_null) 0 else .5
  y <- f0 + epsilon + mean_shift
  return(y)
}

# List of test functions to apply
c2st_test_functions <- list(
  # LinearMMD_test = LinearMMD_test,
  # CLF_test = CLF_test,
  # CP_test = CP_test,
  debiased_test = debiased_test,
  CV_LinearMMD_test = CV_LinearMMD_test,
  CV_CLF_test = CV_CLF_test
)

cit_test_functions <- list(
  RCIT_test = RCIT_test,
  RCoT_test = RCoT_test,
  PCM_test = PCM_test,
  GCM_test = GCM_test
  # WGSC_test = WGSC_test
  # WGCM_test = WGCM_test
)

# Define parameters
n_values <- c(2000)
n_sims <- 500
alpha <- 0.05
d_values <- c(10)
estimators <- c("LL", "QL", "KLR")

results <- list()

print("Model A' is running...")  
for (n in n_values) {
  for (d in d_values) {
    for (is_null in c(TRUE, FALSE)) {
      h_label <- if (is_null) "Null" else "Alternative"
      
      for (test_type in c("C2ST", "CIT")) {
        test_functions <- if (test_type == "C2ST") c2st_test_functions else cit_test_functions
        for (test_name in names(test_functions)) {
          for (extra_param in if (test_type == "C2ST") estimators else c(TRUE, FALSE)) {
            cat("[Test] ", test_name, "\n")
            cat("[Settings] ", "sample size: ", n, " | dimension: ", d, " | ", if (test_type == "C2ST") "estimator: " else "algorithm 1: ", extra_param, " | under ", h_label, "\n")
            
            result <- pbapply::pbsapply(1:n_sims, function(sim) {
              seed <- 1203 + sim
              set.seed(seed)
              
              # Generate data for Group 1
              x1 <- generate_data(n, d, group=1)
              y1 <- generate_y(x1, is_null = TRUE)
              
              # Generate data for Group 2
              seed <- seed + n_sims
              set.seed(seed)
              x2 <- generate_data(n, d, group=2)
              y2 <- generate_y(x2, is_null = is_null)
              
              test_args <- list(x1, x2, y1, y2, seed = seed)
              if (test_type == "C2ST") {
                test_args$est.method <- extra_param
              } else {
                
                test_args$regr.method <- lm_reg_method
                test_args$binary.regr.method <- lm_reg_method_binary
                test_args$alg1 <- extra_param
              }
              
              result <- do.call(test_functions[[test_name]], test_args)
              return(result)
            }, simplify = "array")
            
            mean_result <- mean(result)
            results[[paste(test_type, "_", test_name, extra_param, n, d, h_label, sep = "_")]] <- mean_result
            
            cat("[Result]: ", mean_result, "\n")
            cat(rep("-", 50), '\n')
          }
        }
      }
    }
  }
}

print(results)