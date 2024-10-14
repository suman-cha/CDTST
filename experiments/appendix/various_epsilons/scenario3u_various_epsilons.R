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
tag <- "S3U_various_epsilon"
source("all_tests.R")

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
test_functions <- list(
  RCIT_test = RCIT_test,
  PCM_test = PCM_test,
  GCM_test = GCM_test,
  WGSC_test = WGSC_test
)

# Parameters
n_values <- c(200, 500, 1000, 2000)
n_sims <- 500
alpha <- 0.05
d <- 10
epsilon_types <- c("1/n", "1/sqrt(log(n))", "1/log(n)", "1/sqrt(n)")

results_list <- list()

for (n in n_values) {
  for (is_null in c(TRUE, FALSE)) {
    h_label <- if (is_null) "Null" else "Alternative"
    
    for (epsilon_type in epsilon_types) {
      for (test_name in names(test_functions)) {
        cat("[Test] ", test_name, "\n")
        cat("[Settings] n = ", n, ", epsilon = ", epsilon_type, ", ", h_label, "\n")
        
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          seed <- 1203 + sim
          set.seed(seed)
          
          # Generate data for Group 1 (Null)
          x1 <- generate_data(n, d, group = 1)
          y1 <- generate_y(x1, is_null = TRUE)
          
          # Generate data for Group 2 (Null or Alternative)
          seed <- seed + n_sims
          set.seed(seed)
          x2 <- generate_data(n, d, group = 2)
          y2 <- generate_y(x2, is_null = is_null)
          
          # Create test argument list
          test_args <- list(x1, x2, y1, y2, seed = seed)
          
          # Apply epsilon based on the type
          test_args$epsilon <- switch(epsilon_type,
                                      "1/n" = 1/n,
                                      "1/sqrt(log(n))" = 1/sqrt(log(n)),
                                      "1/log(n)" = 1/log(n),
                                      "1/sqrt(n)" = 1/sqrt(n))
          
          # Run the test
          result <- do.call(test_functions[[test_name]], test_args)
          return(result)
          
        }, simplify = "array")
        
        # Collect results
        mean_result <- mean(result)
        results_list[[length(results_list) + 1]] <- data.table(
          test_name = test_name,
          extra_param = TRUE,
          regressor = NA,  # Regressor column is no longer relevant
          n = n,
          d = d,
          epsilon_type = epsilon_type,
          h_label = h_label,
          mean_result = mean_result
        )
        
        cat("[Result]: ", mean_result, "\n")
        cat(rep("-", 50), '\n')
      }
    }
  }
}

results_dt <- rbindlist(results_list)

# Define filename and save to CSV
filename <- paste0("results/ablation_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")