rm(list = ls())
set.seed(1203)
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

data("diamonds")
data <- diamonds

s <- 6
X <- as.matrix(data[, c("carat", "depth", "table", "x", "y", "z")], nrow=nrow(data), ncol=s)
colnames(X) <- c("V1", "V2", "V3", "V4", "V5", "V6")
Y <- data$price

normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

X_norm <- apply(X, 2, normalize)
Y_norm <- normalize(Y)

sample_data <- function(X, Y, n, is_null = TRUE, is_x1 = TRUE) {
  if (is_x1) {
    idx <- sample(1:nrow(X), n, replace = FALSE)
    x <- X[idx, , drop = FALSE]
  } else {
    means <- apply(X, 2, mean)
    sds <- apply(X, 2, sd)
    mixture <- rbinom(n, 1, 0.5)
    x <- matrix(0, nrow = n, ncol = ncol(X))
    for (i in 1:ncol(X)) {
      x[,i] <- ifelse(mixture == 1, 
                      rnorm(n, means[i] - sds[i], sds[i]), 
                      rnorm(n, means[i] + sds[i], sds[i]))
    }
  }
  if (is_null) {
    y <- sample(Y, n, replace = FALSE)
  } else {
    if (is_x1) {
      weights <- dbeta(Y, 4, 1.5)
      weights <- weights / sum(weights)
      y <- sample(Y, n, replace = FALSE, prob = weights)
    } else {
      weights <- dbeta(Y, 3.5, 2)
      weights <- weights / sum(weights)
      y <- sample(Y, n, replace = FALSE, prob = weights)
    }
  }
  
  return(list(x = x, y = y))
}

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

for (test_name in names(c2st_tests)) {
  for (estimator in estimators) {
    for (is_null in c(TRUE, FALSE)) {
      h_label <- if (is_null) "Null" else "Alt"
      for (n in n_vals) {
        cat("[Test] ", test_name, "\n")
        cat("[Settings] ", "Sample size: ", n, " | Estimator: ", estimator, " | Hypothesis: ", h_label, "\n")
        
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          set.seed(1203 + sim)
          d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
          set.seed(1203 + sim + n_sims)
          d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
          
          test_args <- list(d1$x, d2$x, d1$y, d2$y, est.method = estimator, seed = 1203 + sim)
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

for (test_name in names(cit_tests)) {
  for (alg1 in c(TRUE, FALSE)) {
    for (is_null in c(TRUE, FALSE)) {
      h_label <- if (is_null) "Null" else "Alt"
      for (n in n_vals) {
        cat("[CIT Test] ", test_name, "\n")
        cat("[Settings] Size: ", n, " | Algorithm1: ", alg1, " | Hypothesis: ", h_label, "\n")
        
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          set.seed(1203 + sim)
          d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
          set.seed(1203 + sim + n_sims)
          d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
          
          test_args <- list(d1$x, d2$x, d1$y, d2$y, regr.method = ranger_reg_method, binary.regr.method = ranger_reg_method_binary, alg1 = alg1, seed = 1203 + sim)
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

write.csv(results_df, file = "results/simulation_results_low_dim.csv", row.names = FALSE)
print(results_df)
