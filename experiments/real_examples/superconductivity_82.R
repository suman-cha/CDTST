# Clean up the environment
rm(list = ls())
setwd("/cloud/project/CDTST/simulation/Ours")

# Load necessary libraries
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

# Load the tests source file
source("all_tests.R")

# Set seed for reproducibility
set.seed(1203)

# Load the Superconductivity Dataset
file_path <- "/cloud/project/CDTST/simulation/Ours/Real_examples/data/superconductivty/superconductivty/train.csv"
superconductivity_data <- read.csv(file_path)

# Prepare data: Select features and the target variable
X <- as.matrix(superconductivity_data[, -ncol(superconductivity_data)])
Y <- superconductivity_data[, "critical_temp"]

# Remove missing data
X <- na.omit(X)
Y <- Y[complete.cases(Y)]

# Step 1: Uniformly sample X1
sample_indices_X1 <- sample(1:nrow(X), size = nrow(X) / 2, replace = FALSE)
X1 <- X[sample_indices_X1, , drop = FALSE]
re_idx <- setdiff(1:nrow(X), sample_indices_X1)
X_remaining <- X[re_idx, , drop = FALSE]
mean_X1 <- colMeans(X1)
linear_projection <- X_remaining %*% mean_X1
quadratic_potential <- rowSums((X_remaining - mean_X1)^2)
prob_X_remaining <- exp(-quadratic_potential / sd(quadratic_potential))
prob_X_remaining <- prob_X_remaining / sum(prob_X_remaining)  # Normalize probabilities

set.seed(1204)
sample_indices_X2 <- sample(re_idx, size = nrow(X) / 2, replace = FALSE, prob = prob_X_remaining)


# Normalize function: z-score normalization
normalize <- function(data) {
  return((data - min(data)) / (max(data) - min(data)))
}

# Normalize the entire X matrix
X_normalized <- apply(X, 2, normalize)

# After normalization, split the data again into X1 and X2
X1 <- X_normalized[sample_indices_X1, , drop=FALSE]
X2 <- X_normalized[sample_indices_X2, , drop=FALSE]

# Normalize the response variable Y
Y_normalized <- normalize(Y)

# Split the normalized Y into Y1 and Y2
Y1 <- Y_normalized[sample_indices_X1]
Y2 <- Y_normalized[sample_indices_X2]

# List of test functions to apply
c2st_test_functions <- list(
  # LinearMMD_test = LinearMMD_test,
  # CLF_test = CLF_test,
  # CP_test = CP_test,
  # CV_LinearMMD_test = CV_LinearMMD_test,
  # CV_CLF_test = CV_CLF_test,
  # debiased_test = debiased_test
)

cit_test_functions <- list(
  RCIT_test = RCIT_test,
  GCM_test = GCM_test,
  WGSC_test = WGSC_test,
  PCM_test = PCM_test,
  RCoT_test = RCoT_test
)

# Define parameters
n_values <- c(200, 400, 800, 1200, 1600, 2000)
n_sims <- 500
estimators <- c("LL")

results <- list()
c2st_results_df <- data.frame()
cit_results_df <- data.frame()

# Run C2ST tests
for (c2st_test in names(c2st_test_functions)) {
  for (estimator in estimators) {
    for (is_null in c(TRUE, FALSE)) {
      h_label <- if (is_null) "Null" else "Alternative"
      for (n in n_values) {
        cat("[Test] ", c2st_test, "\n")
        cat("[Settings] ", "sample size: ", n, " | estimator: ", estimator, " | under ", h_label, "\n")
        
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          seed <- 1203 + sim
          set.seed(seed)
          idx1 <- sample(1:nrow(X1), n, replace = FALSE)
          seed <- 1203 + sim + n_sims
          set.seed(seed)
          idx2 <- sample(1:nrow(X2), n, replace = FALSE)
          x1 <- X1[idx1,,drop=FALSE]
          x2 <- X2[idx2,,drop=FALSE]
          if (is_null) {
            y1 <- Y1[idx1]
            y2 <- Y2[idx2]
          } else {
            seed <- 1203 + sim
            set.seed(seed)
            # Uniform sampling for y1 (unbiased sampling from Y1)
            q1 <- quantile(Y1, 0.25)
            q3 <- quantile(Y1, 0.75)
            prob_y1 <- ifelse(Y1 < q1 | Y1 > q3, 0.5, 0.5)  # Assign higher weights to tails
            prob_y1 <- prob_y1 / sum(prob_y1)  # Normalize
            y1 <- Y1[sample(1:length(Y1), size = length(idx1), replace = FALSE, prob = prob_y1)]
            
            # y2 sampling: Emphasize tails based on quantiles
            set.seed(1203 + sim + n_sims)
            q1 <- quantile(Y2, 0.25)
            q3 <- quantile(Y2, 0.75)
            prob_y2 <- ifelse(Y2 < q1 | Y2 > q3, 0.6, 0.4)  # Assign higher weights to tails
            prob_y2 <- prob_y2 / sum(prob_y2)  # Normalize
            y2 <- Y2[sample(1:length(Y2), size = length(idx2), replace = FALSE, prob = prob_y2)]
          }
          
          test_args <- list(x1, x2, y1, y2, est.method = estimator, seed = 1203 + sim)
          result <- do.call(c2st_test_functions[[c2st_test]], test_args)
          return(result)
        }, simplify = "array")
        
        mean_result <- mean(result)
        results[[paste("C2ST_", c2st_test, estimator, n, h_label, sep = "_")]] <- mean_result
        
        c2st_results_df <- rbind(c2st_results_df, data.frame(
          TestType = "C2ST",
          Test = c2st_test,
          Estimator = estimator,
          SampleSize = n,
          Hypothesis = h_label,
          Result = mean_result
        ))
        
        cat("[Result]: ", mean_result, "\n")
        cat(rep("-", 40), '\n')
      }
    }
  }
}

# Run CIT tests
for (cit_test in names(cit_test_functions)) {
  for (alg1 in c(TRUE, FALSE)) {
    for (is_null in c(FALSE, TRUE)) {
      h_label <- if (is_null) "Null" else "Alternative"
      for (n in n_values) {
        cat("[Test] ", cit_test, "\n")
        cat("[Settings] ", "sample size: ", n, " | algorithm 1: ", alg1, " | under ", h_label, "\n")
        
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          seed <- 1203 + sim
          set.seed(seed)
          idx1 <- sample(1:nrow(X1), n, replace = FALSE)
          seed <- 1203 + sim + n_sims
          set.seed(seed)
          idx2 <- sample(1:nrow(X2), n, replace = FALSE)
          x1 <- X1[idx1,,drop=FALSE]
          x2 <- X2[idx2,,drop=FALSE]
          if (is_null) {
            y1 <- Y1[idx1]
            y2 <- Y2[idx2]
          } else {
            seed <- 1203 + sim
            set.seed(seed)
            # Uniform sampling for y1 (unbiased sampling from Y1)
            q1 <- quantile(Y1, 0.25)
            q3 <- quantile(Y1, 0.75)
            prob_y1 <- ifelse(Y1 < q1 | Y1 > q3, 0.5, 0.5)  # Assign higher weights to tails
            prob_y1 <- prob_y1 / sum(prob_y1)  # Normalize
            y1 <- Y1[sample(1:length(Y1), size = length(idx1), replace = FALSE, prob = prob_y1)]
            
            # y2 sampling: Emphasize tails based on quantiles
            set.seed(1203 + sim + n_sims)
            q1 <- quantile(Y2, 0.25)
            q3 <- quantile(Y2, 0.75)
            prob_y2 <- ifelse(Y2 < q1 | Y2 > q3, 0.6, 0.4)  # Assign higher weights to tails
            prob_y2 <- prob_y2 / sum(prob_y2)  # Normalize
            y2 <- Y2[sample(1:length(Y2), size = length(idx2), replace = FALSE, prob = prob_y2)]
          }
          
          test_args <- list(x1, x2, y1, y2, regr.method = ranger_reg_method, binary.regr.method = ranger_reg_method_binary, alg1 = alg1, seed = 1203 + sim)
          result <- do.call(cit_test_functions[[cit_test]], test_args)
          return(result)
        }, simplify = "array")
        
        mean_result <- mean(result)
        results[[paste("CIT_", cit_test, n, alg1, h_label, sep = "_")]] <- mean_result
        
        cit_results_df <- rbind(cit_results_df, data.frame(
          TestType = "CIT",
          Test = cit_test,
          Algorithm1 = alg1,
          SampleSize = n,
          Hypothesis = h_label,
          Result = mean_result
        ))
        
        cat("[Result]: ", mean_result, "\n")
        cat(rep("-", 40), '\n')
      }
    }
  }
}

# Save results to CSV
write.csv(c2st_results_df, file = "c2st_test_results_high_dim.csv", row.names = FALSE)
write.csv(cit_results_df, file = "cit_test_results_high_dim.csv", row.names = FALSE)

print(results)
