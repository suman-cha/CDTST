rm(list=ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(MASS)
  library(data.table)     
  library(tmvtnorm)
})
source("./experiments/all_tests.R")
tag <- "S3U_regression_comparison"

# Data generation functions
generate_data <- function(n, p, group) {
  mu <- if (group == 1) c(1, 1, -1, -1, rep(0, p - 4)) else rep(0, p)
  sigma <- diag(1, p)
  x <- mvrnorm(n, mu = mu, Sigma = sigma)
  return(x)
}


generate_y <- function(x, is_null) {
  n <- nrow(x)
  # List of possible transformations
  transformations <- list(
    function(z) z,
    function(z) z^2,
    function(z) z^3,
    function(z) sin(z),
    function(z) tanh(z)
  )
  # Randomly choose a transformation
  random_transformation <- sample(transformations, 1)[[1]]
  
  if (is_null) {
    y <- cos(rowSums(x) + 2 * rnorm(n))
    
  } else {
    y<- random_transformation(rowSums(x) + 2 * rnorm(n))
    
  }
  
  return(y)
}

# Test functions
cit_test_functions <- list(
  PCM_test = PCM_test,
  GCM_test = GCM_test,
  WGSC_test = WGSC_test
)

# Parameters
n_values <- c(200, 500, 1000, 2000)
n_sims <- 500
alpha <- 0.05
d <- 10
alg1_list <- c(TRUE, FALSE)
regression_methods <- c("lm", "rf", "xgboost")
results_list <- list()

# Simulation loop
for (n in n_values) {
  for (is_null in c(TRUE, FALSE)) {
    h_label <- if (is_null) "Null" else "Alternative"
    
    for (test_name in names(cit_test_functions)) {
      for (regressor in regression_methods) {
        for (alg1 in alg1_list) {
          cat("[Test] ", test_name, "\n")
          cat("[Settings] n:", n, "| regression method:", regressor, "| alg1:", alg1, "| under", h_label, "\n")
          
          result <- pbapply::pbsapply(1:n_sims, function(sim) {
            seed <- 1203 + sim
            set.seed(seed)
            
            # Generate data for Group 1 and 2
            x1 <- generate_data(n, d, group = 1)
            y1 <- generate_y(x1, is_null = TRUE)
            set.seed(seed + n_sims)
            x2 <- generate_data(n, d, group = 2)
            y2 <- generate_y(x2, is_null = is_null)
            
            test_args <- list(x1, x2, y1, y2, seed = seed, alg1 = alg1)
            
            # Assign regression methods
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
            
            do.call(cit_test_functions[[test_name]], test_args)
          }, simplify = "array")
          
          mean_result <- mean(result)
          results_list[[length(results_list) + 1]] <- data.table(
            test_name = test_name,
            regression_method = regressor,
            n = n,
            h_label = h_label,
            alg1 = alg1,
            rejection_rate = mean_result
          )
          
          cat("[Result]:", mean_result, "\n", strrep("-", 50), "\n")
        }
      }
    }
  }
}

results_dt <- rbindlist(results_list)

# Save the results
filename <- paste0("results/ablation_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")
