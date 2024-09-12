rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(MASS)
  library(glmnet)
  library(ANN2)
  library(CVST)
  library(kernlab)
  library(pbapply)
  library(data.table)
  library(ggplot2)
  library(ranger)
  library(xgboost)
})
tag <- "S1U"
source("all_tests.R")

# Data generation functions
generate_data <- function(n, p, group) {
  mu <- if (group == 1) c(1, 1, -1, -1, rep(0, p - 4)) else rep(0, p)
  sigma <- diag(1, p)
  x <- mvrnorm(n, mu = mu, Sigma = sigma)
  return(x)
}

generate_y <- function(x, is_null = TRUE, sigma = 2) {
  n <- nrow(x)
  epsilon <- rt(n, df = sigma)
  f0 <- x %*% c(1, -1, -1, 1, rep(0, dim(x)[2] - 4))
  mean_shift <- if (is_null) 0 else .5
  y <- f0 + epsilon + mean_shift
  return(y)
}

# Test functions
c2st_test_functions <- list(
  # LinearMMD_test = LinearMMD_test,
  # CLF_test = CLF_test,
  # CP_test = CP_test,
  # debiased_test = debiased_test,
  # CVLinearMMD_test = CV_LinearMMD_test,
  # CVCLF_test = CV_CLF_test
)

cit_test_functions <- list(
  # RCIT_test = RCIT_test,
  # RCoT_test = RCoT_test,
  PCM_test = PCM_test,
  GCM_test = GCM_test,
  WGSC_test = WGSC_test
)

# Parameters
n_values <- c(200, 500, 1000, 2000)
n_sims <- 500
alpha <- 0.05
d_values <- c(10)
estimators <- c("LL", "QL", "KLR")
regressors <- c("rf")

results_list <- list()

print("Scenario1(U) is running...")

for (n in n_values) {
  for (d in d_values) {
    for (is_null in c(TRUE, FALSE)) {
      h_label <- if (is_null) "Null" else "Alternative"
      
      for (test_type in c("C2ST", "CIT")) {
        test_functions <- if (test_type == "C2ST") c2st_test_functions else cit_test_functions
        for (test_name in names(test_functions)) {
          for (extra_param in if (test_type == "C2ST") estimators else c(TRUE, FALSE)) {
            if (test_name %in% c("PCM_test", "GCM_test", "WGSC_test")) {
              for (regressor in regressors) {
                cat("[Test] ", test_name, "\n")
                cat("[Settings] ", "sample size: ", n, " | dimension: ", d, " | ", 
                    if (test_type == "C2ST") "estimator: " else "algorithm 1: ", extra_param, 
                    " | regressor: ", regressor, " | under ", h_label, "\n")
                
                result <- pbapply::pbsapply(1:n_sims, function(sim) {
                  seed <- 1203 + sim
                  set.seed(seed)
                  
                  # Generate data for Group 1
                  x1 <- generate_data(n, d, group = 1)
                  y1 <- generate_y(x1, is_null = TRUE)
                  
                  # Generate data for Group 2
                  seed <- seed + n_sims
                  set.seed(seed)
                  x2 <- generate_data(n, d, group = 2)
                  y2 <- generate_y(x2, is_null = is_null)
                  
                  test_args <- list(x1, x2, y1, y2, seed = seed)
                  
                  if (regressor == "lm") {
                    test_args$regr.method <- lm_reg_method
                    test_args$binary.regr.method <- lm_reg_method_binary
                  } else if (regressor == "rf") {
                    test_args$regr.method <- ranger_reg_method
                    test_args$binary.regr.method <- ranger_reg_method_binary
                  } else if (regressor == "xgboost") {
                    test_args$regr.method <- xgboost_reg_method
                    test_args$binary.regr.method <- xgboost_reg_method_binary
                  }
                  test_args$alg1 <- extra_param
                  
                  result <- do.call(test_functions[[test_name]], test_args)
                  return(result)
                }, simplify = "array")
                
                mean_result <- mean(result, na.rm = TRUE)
                results_list[[length(results_list) + 1]] <- data.table(
                  test_type = test_type,
                  test_name = test_name,
                  extra_param = extra_param,
                  regressor = regressor,
                  n = n,
                  d = d,
                  h_label = h_label,
                  mean_result = mean_result
                )
                
                cat("[Result]: ", mean_result, "\n")
                cat(rep("-", 50), '\n')
              }
            } else {
              cat("[Test] ", test_name, "\n")
              cat("[Settings] ", "sample size: ", n, " | dimension: ", d, " | ", 
                  if (test_type == "C2ST") "estimator: " else "algorithm 1: ", extra_param, 
                  " | under ", h_label, "\n")
              
              result <- pbapply::pbsapply(1:n_sims, function(sim) {
                seed <- 1203 + sim
                set.seed(seed)
                
                # Generate data for Group 1
                x1 <- generate_data(n, d, group = 1)
                y1 <- generate_y(x1, is_null = TRUE)
                
                # Generate data for Group 2
                seed <- seed + n_sims
                set.seed(seed)
                x2 <- generate_data(n, d, group = 2)
                y2 <- generate_y(x2, is_null = is_null)
                
                test_args <- list(x1, x2, y1, y2, seed = seed)
                if (test_type == "CIT") {
                  test_args$alg1 <- extra_param
                } else {
                  test_args$est.method <- extra_param
                }
                
                result <- do.call(test_functions[[test_name]], test_args)
                return(result)
              }, simplify = "array")
              
              mean_result <- mean(result, na.rm = TRUE)
              results_list[[length(results_list) + 1]] <- data.table(
                test_type = test_type,
                test_name = test_name,
                extra_param = extra_param,
                regressor = NA,
                n = n,
                d = d,
                h_label = h_label,
                mean_result = mean_result
              )
              
              cat("[Result]: ", mean_result, "\n")
              cat(rep("-", 50), '\n')
            }
          }
        }
      }
    }
  }
}

# Combine all results into a single data.table
results_dt <- rbindlist(results_list)

# Define filename and save to CSV
filename <- paste0("results/cit_simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")