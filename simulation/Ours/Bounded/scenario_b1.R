# Clear workspace
rm(list = ls())

# Suppress startup messages and load packages
suppressPackageStartupMessages({
  library(MASS)
  library(glmnet)
  library(ANN2)
  library(CVST)
  library(kernlab)
  library(foreach)
  library(doParallel)
  library(data.table)  
  library(ggplot2)
  library(tmvtnorm)
  library(RCIT)
})

# Source additional scripts
source("all_tests.R")

set.seed(1203)

# Function to generate covariates using truncated multivariate normal distributions
generate_data <- function(n, p, group) {
  mu <- c(1, 1, -1, -1, rep(0, p - 4))
  sigma <- diag(1, p)
  lb <- rep(-.5, p)
  ub <- rep(.5, p)
  
  if (group == 1) {
    x <- rtmvnorm(n, mean = mu, sigma = sigma, lower = lb, upper = ub, algorithm = "gibbs")
  } else {
    x <- rtmvnorm(n, mean = rep(0, p), sigma = sigma, lower = lb, upper = ub, algorithm = "gibbs")
  }
  return(x)
}

# Function to generate response variable
generate_y <- function(x, is_null = TRUE, sigma = 2) {
  n <- nrow(x)
  epsilon <- rt(n, df = sigma)
  f0 <- x %*% c(1, -1, -1, 1, rep(0, dim(x)[2] - 4))
  mean_shift <- if (is_null) 0 else .5
  y <- f0 + epsilon + mean_shift
  return(y)
}

# Define test functions
c2st_test_functions <- list(
  LinearMMD_test = LinearMMD_test,
  CLF_test = CLF_test,
  CP_test = CP_test,
  debiased_test = debiased_test,
  CV_LinearMMD_test = CV_LinearMMD_test,
  CV_CLF_test = CV_CLF_test
)

cit_test_functions <- list(
  RCIT_test = RCIT_test,
  RCoT_test = RCoT_test,
  PCM_test = PCM_test,
  GCM_test = GCM_test,
  WGSC_test = WGSC_test
)

# Define parameters
n_values <- c(200, 500, 1000, 2000)
n_sims <- 500
alpha <- 0.05
d_values <- c(10)
estimators <- c("LL", "QL", "KLR")

# Set up parallel backend
num_cores <- detectCores() - 1  # Use all but one core to avoid system overload
cat("Number of available cores:", num_cores, "\n")
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Explicitly source the necessary R scripts on each worker node
clusterEvalQ(cl, {
  source("all_tests.R")
  library(MASS)
  library(glmnet)
  library(ANN2)
  library(CVST)
  library(kernlab)
  library(tmvtnorm)
  library(RCIT)
  library(data.table)
})

results_dt <- data.table()

print("Model A' is running...")

# Parallelize the computations
results_dt <- foreach(n = n_values, .combine = 'rbind', .packages = c("MASS", "glmnet", "ANN2", "CVST", "kernlab", "tmvtnorm", "RCIT", "data.table", "foreach", "doParallel")) %:%
  foreach(d = d_values, .combine = 'rbind') %:%
  foreach(is_null = c(TRUE, FALSE), .combine = 'rbind') %:%
  foreach(test_type = c("C2ST", "CIT"), .combine = 'rbind') %:%
  foreach(test_name = if (test_type == "C2ST") names(c2st_test_functions) else names(cit_test_functions), .combine = 'rbind') %:%
  foreach(extra_param = if (test_type == "C2ST") estimators else c(TRUE, FALSE), .combine = 'rbind') %dopar% {
    
    # Choose the appropriate test function list
    test_functions <- if (test_type == "C2ST") c2st_test_functions else cit_test_functions
    h_label <- if (is_null) "Null" else "Alternative"
    
    result <- foreach(sim = 1:(n_sims / 50), .combine = 'c', .packages = c("MASS", "glmnet", "ANN2", "CVST", "kernlab", "tmvtnorm", "RCIT")) %dopar% {
      batch_result <- numeric(50)
      for (i in 1:50) {
        seed <- 1203 + (sim - 1) * 50 + i
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
        if (test_type == "C2ST") {
          test_args$est.method <- extra_param
        } else {
          test_args$regr.method <- ranger_reg_method
          test_args$binary.regr.method <- ranger_reg_method_binary
          test_args$alg1 <- extra_param
        }
        
        batch_result[i] <- do.call(test_functions[[test_name]], test_args)
      }
      return(batch_result)
    }
    
    mean_result <- mean(result)
    return(data.table(
      test_type = test_type,
      test_name = test_name,
      extra_param = extra_param,
      n = n,
      d = d,
      h_label = h_label,
      mean_result = mean_result
    ))
  }

# Convert results to data frame and print
results_df <- as.data.frame(results_dt)
print(results_df)

# Save the data frame to a CSV file
fwrite(results_df, "simulation_results.csv", row.names = FALSE)

# Stop the parallel cluster
stopCluster(cl)
