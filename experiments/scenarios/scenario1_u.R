rm(list = ls())
set.seed(1203)
setwd("/cloud/project/CDTST/simulation/Ours")
suppressPackageStartupMessages({
  library(MASS)
  library(glmnet)
  library(ANN2)
  library(CVST)
  library(kernlab)
  library(doParallel)     
  library(foreach)        
  library(data.table)     
  library(ggplot2)
})

all_tests_path <- normalizePath("all_tests.R", mustWork = TRUE)
utils_path <- normalizePath("utils.R", mustWork = TRUE)
source(all_tests_path)
source(utils_path)

cores <- detectCores() - 1  
cat("Number of available cores:", cores, "\n")
cl <- makeCluster(cores)    
registerDoParallel(cl)      

clusterEvalQ(cl, {
  source("all_tests.R")
  source("utils.R")
})

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
  LinearMMD_test = LinearMMD_test,
  CLF_test = CLF_test,
  CP_test = CP_test,
  debiased_test = debiased_test,
  CVLinearMMD_test = CV_LinearMMD_test,
  CVCLF_test = CV_CLF_test
)

cit_test_functions <- list(
  RCIT_test = RCIT_test,
  RCoT_test = RCoT_test,
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

results <- list()
results_dt <- data.table()

print("Scenario1(U) is running...")

for (n in n_values) {
  for (d in d_values) {
    for (is_null in c(FALSE, TRUE)) {
      h_label <- if (is_null) "Null" else "Alternative"
      
      for (test_type in c("C2ST", "CIT")) {
        test_functions <- if (test_type == "C2ST") c2st_test_functions else cit_test_functions
        for (test_name in names(test_functions)) {
          for (extra_param in if (test_type == "C2ST") estimators else c(TRUE, FALSE)) {
            cat("[Test] ", test_name, "\n")
            cat("[Settings] ", "sample size: ", n, " | dimension: ", d, " | ", if (test_type == "C2ST") "estimator: " else "algorithm 1: ", extra_param, " | under ", h_label, "\n")
            
            chunk_size <- 50  
            result <- foreach(sim = 1:(n_sims / chunk_size), .combine = 'c', 
                              .packages = c("MASS", "glmnet", "ANN2", "CVST", "kernlab"),
                              .export = c("generate_data", "generate_y")) %dopar% {
                                batch_result <- numeric(chunk_size)
                                for (i in 1:chunk_size) {
                                  seed <- 1203 + (sim - 1) * chunk_size + i
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
            results_dt <- rbindlist(list(results_dt, data.table(
              test_type = test_type,
              test_name = test_name,
              extra_param = extra_param,
              n = n,
              d = d,
              h_label = h_label,
              mean_result = mean_result
            )))
            
            cat("[Result]: ", mean_result, "\n")
            cat(rep("-", 50), '\n')
          }
        }
      }
    }
  }
}

results_df <- as.data.frame(results_dt)
print(results_df)

# Save the data frame to a CSV file
fwrite(results_df, "simulation_results.csv", row.names = FALSE)
stopCluster(cl)
