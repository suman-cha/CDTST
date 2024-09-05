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

# Sampling indices
X1_idx <- sample(1:nrow(X), size = 10000, replace = FALSE)
X1 <- X[X1_idx, , drop = FALSE]
re_idx <- setdiff(1:nrow(X), X1_idx)
X_re <- X[re_idx, , drop = FALSE]
X_re_mean <- rowMeans(X_re)
X_re_var <- rowSums((X_re - X_re_mean)^2)
X_re_prob <- exp(- X_re_mean / sqrt(X_re_var))
X_re_prob <- X_re_prob / sum(X_re_prob)

# Sampling for X2
set.seed(1204)
X2_idx <- sample(1:nrow(X_re), size = 10000, replace = FALSE, prob = X_re_prob)

# Normalize using scale()
X_norm <- scale(X)
X1 <- X_norm[X1_idx, , drop = FALSE]
X2 <- X_norm[X2_idx, , drop = FALSE]
Y_norm <- scale(Y)
Y1 <- Y_norm[X1_idx]
Y2 <- Y_norm[X2_idx]

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
          idx1 <- sample(1:nrow(X1), n, replace = FALSE)
          seed <- 1203 + sim + n_sims
          set.seed(seed)
          idx2 <- sample(1:nrow(X2), n, replace = FALSE)
          x1 <- X1[idx1, , drop = FALSE]
          x2 <- X2[idx2, , drop = FALSE]
          y1 <- if (is_null) Y1[idx1] else {
            prob_y1 <- ifelse(Y1 < quantile(Y1, 0.25) | Y1 > quantile(Y1, 0.75), 0.5, 0.5)
            prob_y1 <- prob_y1 / sum(prob_y1)
            Y1[sample(1:length(Y1), size = length(idx1), replace = FALSE, prob = prob_y1)]
          }
          y2 <- if (is_null) Y2[idx2] else {
            prob_y2 <- ifelse(Y2 < quantile(Y2, 0.25) | Y2 > quantile(Y2, 0.75), 0.6, 0.4)
            prob_y2 <- prob_y2 / sum(prob_y2)
            Y2[sample(1:length(Y2), size = length(idx2), replace = FALSE, prob = prob_y2)]
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
          seed <- 1203 + sim
          set.seed(seed)
          idx1 <- sample(1:nrow(X1), n, replace = FALSE)
          seed <- 1203 + sim + n_sims
          set.seed(seed)
          idx2 <- sample(1:nrow(X2), n, replace = FALSE)
          x1 <- X1[idx1, , drop = FALSE]
          x2 <- X2[idx2, , drop = FALSE]
          y1 <- if (is_null) Y1[idx1] else {
            prob_y1 <- ifelse(Y1 < quantile(Y1, 0.25) | Y1 > quantile(Y1, 0.75), 0.5, 0.5)
            prob_y1 <- prob_y1 / sum(prob_y1)
            Y1[sample(1:length(Y1), size = length(idx1), replace = FALSE, prob = prob_y1)]
          }
          y2 <- if (is_null) Y2[idx2] else {
            prob_y2 <- ifelse(Y2 < quantile(Y2, 0.25) | Y2 > quantile(Y2, 0.75), 0.6, 0.4)
            prob_y2 <- prob_y2 / sum(prob_y2)
            Y2[sample(1:length(Y2), size = length(idx2), replace = FALSE, prob = prob_y2)]
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