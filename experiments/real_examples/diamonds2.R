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

# Normalize 
X_norm <- scale(X)
Y_norm <- scale(Y)

sample_data <- function(X, Y, n, is_null = TRUE, is_x1 = TRUE) {
  if (is_x1) {
    # 샘플링할 인덱스 선택
    idx1 <- sample(1:nrow(X), n, replace = FALSE)
    x <- X[idx1, , drop = FALSE]
    y <- Y[idx1]
    
    if (!is_null) {
      prob_y1_norm1 <- dnorm(Y, mean(x)/2, 8*sd(x))
      prob_y1_norm2 <- dnorm(Y, -mean(x)/4, 4*sd(x))
      prob_y1 <- 0.5 * prob_y1_norm1 + 0.5 * prob_y1_norm2
      prob_y1 <- prob_y1 / sum(prob_y1)  
      y <- sample(Y, n, replace = FALSE, prob = prob_y1)
    }
    
  } else {
    feature_to_bias <- X[, "V3"]
    prob <- dt(feature_to_bias, df = 2)
    prob <- prob / sum(prob)
    idx2 <- sample(1:nrow(X), n, replace = FALSE, prob = prob)
    x <- X[idx2, , drop = FALSE]
    y <- Y[idx2]
    
    if (!is_null) {
      mean_x2 <- mean(x)
      sd_x2 <- sd(x)
      
      prob_y2_norm1 <- dnorm(Y, mean_x2/4, sd_x2/4)
      prob_y2_norm2 <- dnorm(Y, -mean_x2/2, sd_x2/8)
      prob_y2 <- 0.5 * prob_y2_norm1 + 0.5 * prob_y2_norm2
      prob_y2 <- prob_y2 / sum(prob_y2) 
      y <- sample(Y, n, replace = FALSE, prob = prob_y2)
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

all_tests <- c(c2st_tests, cit_tests)

results_df <- data.frame(
  TestType = character(),
  Test = character(),
  Extraparam = character(),
  SampleSize = numeric(),
  Hypothesis = character(),
  Result = numeric(),
  stringsAsFactors = FALSE
)

for (test_name in names(all_tests)) {
  test_type <- if (test_name %in% names(c2st_tests)) "C2ST" else "CIT"
  
  for (param in if(test_type == "C2ST") estimators else c(TRUE, FALSE)) {
    param_name <- if(test_type == "C2ST") "Estimator" else "Algorithm1"
    
    for (is_null in c(TRUE, FALSE)) {
      h_label <- if (is_null) "Null" else "Alt"
      for (n in n_vals) {
        cat("[", test_type, " Test] ", test_name, "\n")
        cat("[Settings] Size: ", n, " | ", param_name, ": ", param, " | Hypothesis: ", h_label, "\n")
        
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          seed <- 1203 + sim
          set.seed(seed)
          
          d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
          d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
          
          if (test_type == "C2ST") {
            test_args <- list(d1$x, d2$x, d1$y, d2$y, est.method = param, seed = seed)
          } else {
            test_args <- list(d1$x, d2$x, d1$y, d2$y, regr.method = ranger_reg_method, 
                              binary.regr.method = ranger_reg_method_binary, alg1 = param, seed = seed)
          }
          
          result <- do.call(all_tests[[test_name]], test_args)
          return(result)
        }, simplify = "array")
        
        mean_result <- mean(result)
        new_row <- data.frame(
          TestType = test_type,
          Test = test_name,
          Extraparam = as.character(param),
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