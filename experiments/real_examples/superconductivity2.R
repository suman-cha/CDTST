rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(pbapply)
  library(MASS)
  library(dplyr)
  library(glmnet)
  library(ANN2)
  library(CVST)
  library(kernlab)
  library(sn)
  library(RCIT)
})
source("all_tests.R")

cur_wd <- getwd()
file_path <- "/real_examples/data/superconductivity.csv"
full_path <- file.path(cur_wd, file_path)
data <- read.csv(full_path)

X <- as.matrix(data[, -ncol(data)])
Y <- data[, "critical_temp"]

# Remove missing data
X <- na.omit(X)
Y <- Y[complete.cases(Y)]

# Normalize 
X_norm <- scale(X)
Y_norm <- scale(Y)

# Define test functions
c2st_tests <- list(
  LinearMMD_test = LinearMMD_test,
  CLF_test = CLF_test,
  CP_test = CP_test,
  CV_LinearMMD_test = CV_LinearMMD_test,
  CV_CLF_test = CV_CLF_test,
  debiased_test = debiased_test
)

cit_tests <- list(
  RCIT_test = RCIT_test,
  GCM_test = GCM_test,
  WGSC_test = WGSC_test,
  PCM_test = PCM_test,
  RCoT_test = RCoT_test
)

# Parameters
n_vals <- c(200, 400, 800, 1200, 1600, 2000)
n_sims <- 500
estimators <- c("LL")

results_df <- data.frame(
  TestType = character(),
  Test = character(),
  Extraparam = character(),
  SampleSize = numeric(),
  Hypothesis = character(),
  Result = numeric(),
  stringsAsFactors = FALSE
)

# Run C2ST tests
for (test_name in names(c2st_tests)) {
  for (estimator in estimators) {
    for (is_null in c(FALSE, TRUE)) {
      h_label <- if (is_null) "Null" else "Alt"
      for (n in n_vals) {
        cat("[C2ST Test] ", test_name, "\n")
        cat("[Settings] Size: ", n, " | Estimator: ", estimator, " | Hypothesis: ", h_label, "\n")
        
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          seed <- 1203 + sim
          set.seed(seed)
          
          # Uniform sampling for x1
          idx1 <- sample(1:nrow(X_norm), n, replace = FALSE)
          x1 <- X_norm[idx1, , drop = FALSE]
          
          # Sampling for x2 to create covariate shift
          feature_to_bias <- X_norm[, ncol(X)-13]  
          center <- mean(feature_to_bias)
          bandwidth <- sd(feature_to_bias)
          prob <- dt(feature_to_bias, df = 2)
          prob <- prob / sum(prob)
          idx2 <- sample(1:nrow(X_norm), n, replace = FALSE, prob = prob)
          x2 <- X_norm[idx2, , drop = FALSE]
          
          if (is_null) {
            y1 <- Y_norm[idx1]
            y2 <- Y_norm[idx2]
          } else {
            # For y1: higher probability for values in the middle
            mid_50 <- Y_norm >= quantile(Y_norm, 0.25) & Y_norm <= quantile(Y_norm, 0.75)
            prob_y1 <- ifelse(mid_50, 0.5, 0.5)
            prob_y1 <- prob_y1 / sum(prob_y1)
            y1 <- sample(Y_norm, size = n, replace = FALSE, prob = prob_y1)
            
            # For y2: higher probability for values in the tails
            tails <- Y_norm < quantile(Y_norm, 0.25) | Y_norm > quantile(Y_norm, 0.75)
            prob_y2 <- ifelse(tails, 0.7, 0.3)
            prob_y2 <- prob_y2 / sum(prob_y2)
            y2 <- sample(Y_norm, size = n, replace = FALSE, prob = prob_y2)
          }
          
          test_args <- list(x1, x2, y1, y2, est.method = estimator, seed = seed)
          result <- do.call(c2st_tests[[test_name]], test_args)
          return(result)
        }, simplify = "array")
        
        mean_result <- mean(result)
        new_row <- data.frame(
          TestType = "C2ST",
          Test = test_name,
          Extraparam = estimator,
          SampleSize = n,
          Hypothesis = h_label,
          Result = mean_result,
          stringsAsFactors = FALSE
        )
        
        results_df <- rbind(results_df, new_row)
        
        cat("[Result]: ", mean_result, "\n")
        cat(rep("-", 40), '\n')
      }
    }
  }
}

# Run CIT tests
for (test_name in names(cit_tests)) {
  for (alg1 in c(TRUE, FALSE)) {
    for (is_null in c(FALSE, TRUE)) {
      h_label <- if (is_null) "Null" else "Alt"
      for (n in n_vals) {
        cat("[CIT Test] ", test_name, "\n")
        cat("[Settings] Size: ", n, " | Algorithm1: ", alg1, " | Hypothesis: ", h_label, "\n")
        
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          # Uniform sampling for x1
          seed <- 1203 + sim
          set.seed(seed)
          idx1 <- sample(1:nrow(X_norm), n, replace = FALSE)
          x1 <- X_norm[idx1, , drop = FALSE]
          
          # Sampling for x2 to create covariate shift
          feature_to_bias <- X_norm[, ncol(X)-13]  
          center <- mean(feature_to_bias)
          bandwidth <- sd(feature_to_bias)
          prob <- dt(feature_to_bias, df = 2)
          prob <- prob / sum(prob)
          idx2 <- sample(1:nrow(X_norm), n, replace = FALSE, prob = prob)
          x2 <- X_norm[idx2, , drop = FALSE]
          
          if (is_null) {
            y1 <- Y_norm[idx1]
            y2 <- Y_norm[idx2]
          } else {
            # For y1: Higher probability for values in the middle 50%
            mid_50 <- Y_norm >= quantile(Y_norm, 0.25) & Y_norm <= quantile(Y_norm, 0.75)
            prob_y1 <- ifelse(mid_50, 0.5, 0.5)
            prob_y1 <- prob_y1 / sum(prob_y1)
            y1 <- sample(Y_norm, size = n, replace = FALSE, prob = prob_y1)
            
            # For y2: Higher probability for values in the tails
            tails <- Y_norm < quantile(Y_norm, 0.25) | Y_norm > quantile(Y_norm, 0.75)
            prob_y2 <- ifelse(tails, 0.7, 0.3)
            prob_y2 <- prob_y2 / sum(prob_y2)
            y2 <- sample(Y_norm, size = n, replace = FALSE, prob = prob_y2)
          }
          
          test_args <- list(x1, x2, y1, y2, regr.method = ranger_reg_method, binary.regr.method = ranger_reg_method_binary, alg1 = alg1, seed = seed)
          result <- do.call(cit_tests[[test_name]], test_args)
          return(result)
        }, simplify = "array")
        
        mean_result <- mean(result)
        new_row <- data.frame(
          TestType = "CIT",
          Test = test_name,
          Extraparam = alg1,
          SampleSize = n,
          Hypothesis = h_label,
          Result = mean_result,
          stringsAsFactors = FALSE
        )
        
        results_df <- rbind(results_df, new_row)
        
        cat("[Result]: ", mean_result, "\n")
        cat(rep("-", 40), '\n')
      }
    }
  }
}

# Save the results to CSV
write.csv(results_df, file = "results/simulation_results_high_dim.csv", row.names = FALSE)

print(results_df)