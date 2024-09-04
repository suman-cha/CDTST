# Clean up the environment and set working directory
rm(list = ls())
# setwd("/cloud/project/CDTST/simulation/Ours")

# Load necessary libraries
suppressPackageStartupMessages({
  library(MASS)
  library(glmnet)
  library(ANN2)
  library(CVST)
  library(kernlab)
  library(pbapply)
  library(dplyr)
  library(ggplot2)
})
source("all_tests.R")

# Set seed for reproducibility
set.seed(1203)
data("diamonds")
mydata <- diamonds

s <- 6
X <- as.matrix(mydata[, c("carat", "depth", "table", "x", "y", "z")], nrow=nrow(mydata), ncol=s)
colnames(X) <- c("V1", "V2", "V3", "V4", "V5", "V6")

Y <- mydata$price


# Normalize function: min-max normalization
normalize <- function(data) {
  return((data - min(data)) / (max(data) - min(data)))
}

# Normalize the entire X matrix and Y vector
X_normalized <- apply(X, 2, normalize)
Y_normalized <- normalize(Y)

# Function to sample data
sample_data <- function(X, Y, n, is_null = TRUE, is_x1 = TRUE) {
  if (is_x1) {
    # Sampling method for x1: Uniform sampling
    idx <- sample(1:nrow(X), n, replace = FALSE)
    x <- X[idx, , drop = FALSE]
  } else {
    # Sampling method for x2: Mixture of two normals
    means <- apply(X, 2, mean)
    sds <- apply(X, 2, sd)
    mixture <- rbinom(n, 1, 0.5)
    x <- matrix(0, nrow = n, ncol = ncol(X))
    for (i in 1:ncol(X)) {
      x[,i] <- ifelse(mixture == 1, 
                      rnorm(n, means[i] - sds[i], sds[i]/2), 
                      rnorm(n, means[i] + sds[i], sds[i]/2))
    }
    x <- pmin(pmax(x, 0), 1)  # Ensure values are between 0 and 1
  }
  
  if (is_null) {
    # For null hypothesis, use inverse CDF sampling for both y1 and y2
    u <- runif(n)
    y <- quantile(Y, u)
  } else {
    if (is_x1) {
      # For y1 under alternative hypothesis: emphasize tails
      u <- rbeta(n, 0.5, 0.5)  # U-shaped distribution
      y <- quantile(Y, u)
    } else {
      # For y2 under alternative hypothesis: emphasize center
      u <- rbeta(n, 2, 2)  # Bell-shaped distribution
      y <- quantile(Y, u)
    }
  }
  
  return(list(x = x, y = y))
}

# List of test functions to apply
c2st_test_functions <- list(
  # LinearMMD_test = LinearMMD_test,
  # CLF_test = CLF_test,
  # CP_test = CP_test,
  # CV_LinearMMD_test = CV_LinearMMD_test,
  # CV_CLF_test = CV_CLF_test,
  debiased_test = debiased_test
)

cit_test_functions <- list(
  # RCIT_test = RCIT_test,
  # GCM_test = GCM_test,
  # WGSC_test = WGSC_test,
  # PCM_test = PCM_test,
  # RCoT_test = RCoT_test
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
    for (is_null in c(FALSE, TRUE)) {
      h_label <- if (is_null) "Null" else "Alternative"
      for (n in n_values) {
        cat("[Test] ", c2st_test, "\n")
        cat("[Settings] ", "sample size: ", n, " | estimator: ", estimator, " | under ", h_label, "\n")
        
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          set.seed(1203 + sim)
          data_1 <- sample_data(X_normalized, Y_normalized, n, is_null, is_x1 = TRUE)
          set.seed(1203 + sim + n_sims)
          data_2 <- sample_data(X_normalized, Y_normalized, n, is_null, is_x1 = FALSE)
          
          test_args <- list(data_1$x, data_2$x, data_1$y, data_2$y, est.method = estimator, seed = 1203 + sim)
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
          set.seed(1203 + sim)
          data_1 <- sample_data(X_normalized, Y_normalized, n, is_null, is_x1 = TRUE)
          set.seed(1203 + sim + n_sims)
          data_2 <- sample_data(X_normalized, Y_normalized, n, is_null, is_x1 = FALSE)
          
          test_args <- list(data_1$x, data_2$x, data_1$y, data_2$y, regr.method = ranger_reg_method, binary.regr.method = ranger_reg_method_binary, alg1 = alg1, seed = 1203 + sim)
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
write.csv(c2st_results_df, file = "c2st_test_results_low_dim.csv", row.names = FALSE)
write.csv(cit_results_df, file = "cit_test_results_low_dim.csv", row.names = FALSE)

print(results)