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
      y <- sample(Y, n, replace = FALSE)
    } else {
      Y_sorted <- sort(Y)
      weights <- dbeta(Y, 2, 3)
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

n_vals <- c(200, 400)
n_sims <- 500
estimators <- c("LL")

results_df <- data.frame()

for (test in names(c2st_tests)) {
  for (est in estimators) {
    for (null in c(FALSE, TRUE)) {
      label <- if (null) "Null" else "Alt"
      for (n in n_vals) {
        cat("[Test] ", test, "\n")
        cat("[Settings] ", "Sample size: ", n, " | Estimator: ", est, " | Hypothesis: ", label, "\n")
        
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          set.seed(1203 + sim)
          d1 <- sample_data(X_norm, Y_norm, n, null, TRUE)
          set.seed(1203 + sim + n_sims)
          d2 <- sample_data(X_norm, Y_norm, n, null, FALSE)
          
          test_args <- list(d1$x, d2$x, d1$y, d2$y, est.method = est, seed = 1203 + sim)
          result <- do.call(c2st_tests[[test]], test_args)
          return(result)
        }, simplify = "array")
        
        mean_result <- mean(result)
        results_df <- rbind(results_df, data.frame(
          TestType = "C2ST",
          Test = test,
          Estimator = est,
          SampleSize = n,
          Hypothesis = label,
          Result = mean_result
        ))
        
        cat("[Result]: ", mean_result, "\n")
        cat(rep("-", 40), '\n')
      }
    }
  }
}

for (test in names(cit_tests)) {
  for (alg in c(TRUE, FALSE)) {
    for (null in c(TRUE, FALSE)) {
      label <- if (null) "Null" else "Alt"
      for (n in n_vals) {
        cat("[Test] ", test, "\n")
        cat("[Settings] ", "Sample size: ", n, " | Algorithm 1: ", alg, " | Hypothesis: ", label, "\n")
        
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          set.seed(1203 + sim)
          d1 <- sample_data(X_norm, Y_norm, n, null, TRUE)
          set.seed(1203 + sim + n_sims)
          d2 <- sample_data(X_norm, Y_norm, n, null, FALSE)
          
          test_args <- list(d1$x, d2$x, d1$y, d2$y, regr.method = ranger_reg_method, binary.regr.method = ranger_reg_method_binary, alg1 = alg, seed = 1203 + sim)
          result <- do.call(cit_tests[[test]], test_args)
          return(result)
        }, simplify = "array")
        
        mean_result <- mean(result)
        results_df <- rbind(results_df, data.frame(
          TestType = "CIT",
          Test = test,
          Algorithm1 = alg,
          SampleSize = n,
          Hypothesis = label,
          Result = mean_result
        ))
        
        cat("[Result]: ", mean_result, "\n")
        cat(rep("-", 40), '\n')
      }
    }
  }
}

write.csv(results_df, file = "results/simulation_results_low_dim.csv", row.names = FALSE)
print(results_df)
