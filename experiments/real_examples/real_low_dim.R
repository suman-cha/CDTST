rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(MASS)
  library(pbapply)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(parallel)
})
source("all_tests.R")

tag <- "real_low_dim_ver2"
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
    # uniform sampling for X1
    X_idx <- sample(1:nrow(X), nrow(X)%/%2, replace = FALSE)
    X1 <- X[X_idx, , drop = FALSE]
    Y_subset <- Y[X_idx]
    
    # subsample from X1 to construct x1
    x_idx <- sample(1:nrow(X1), n, replace=FALSE)
    x <- X1[x_idx,,drop=FALSE]
  } else {
    # biased sampling for X2 based on normal distribution
    feature_to_bias <- X[, 1]  
    prob <- dnorm(feature_to_bias, 0, 1)
    prob <- prob / sum(prob)  
    X_idx <- sample(1:nrow(X), nrow(X)%/%2, replace = FALSE, prob = prob)
    X2 <- X[X_idx, , drop = FALSE]
    Y_subset <- Y[X_idx]  
    
    # subsample from X2 to construct x2
    x_idx <- sample(1:nrow(X2), n, replace=FALSE)
    x <- X2[x_idx,,drop=FALSE]
  }
  
  if (is_null) {
    # Null hypothesis: uniform sampling from Y values
    y <- sample(Y_subset, size=n, replace=FALSE)
  } else {
    # Alternative hypothesis: introduce bias in Y1 and Y2
    if (is_x1) {
      u <- dunif(Y_subset, 0, 1)
    } else {
      u <- exp(-Y_subset)
    }
    u <- u / sum(u)
    y <- sample(Y_subset, size = n, prob = u, replace = FALSE)
  }
  
  return(list(x = x, y = y))
}

# Define test functions
drt_test_functions <- list(
  LinearMMD_test = LinearMMD_test,
  CLF_test = CLF_test,
  CP_test = CP_test,
  CV_LinearMMD_test = CV_LinearMMD_test,
  CV_CLF_test = CV_CLF_test,
  debiased_test = debiased_test
)

cit_test_functions <- list(
  RCIT_test = RCIT_test,
  GCM_test = GCM_test,
  WGSC_test = WGSC_test,
  PCM_test = PCM_test
)

n_values <- c(200, 400, 800, 1200, 1600, 2000)
n_sims <- 500
estimators <- c("LL", "KLR")
results_list <- list()

cl <- makeCluster(detectCores() - 6)
pbapply::pboptions(cl = cl)

clusterExport(cl, c("sample_data", "X_norm", "Y_norm", "normalize", "drt_test_functions", "cit_test_functions", "n_sims"))

clusterEvalQ(cl, {
  source("utils.R")  # 클러스터에 utils.R 파일을 로드하여 apply_alg1 함수 사용 가능
  source("all_tests.R")  # all_tests.R 파일도 필요하면 로드
  library(data.table)
})

for (n in n_values) {
  for (is_null in c(FALSE, TRUE)) {
    h_label <- if (is_null) "Null" else "Alternative"
    
    for (test_type in c("DRT", "CIT")) {
      test_functions <- if (test_type == "DRT") drt_test_functions else cit_test_functions
      
      for (test_name in names(test_functions)) {
        if (test_type == "DRT") {
          for (est in estimators) {
            result <- parSapply(cl, 1:n_sims, function(sim, n, is_null, est, test_functions, test_name) {
              seed <- 1203 + sim
              set.seed(seed)
              
              # 데이터 생성
              d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
              set.seed(seed + n_sims)
              d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
              
              test_args <- list(d1$x, d2$x, d1$y, d2$y, est.method = est, seed = seed)
              do.call(test_functions[[test_name]], test_args)
            }, n = n, is_null = is_null, est = est, test_functions = test_functions, test_name = test_name, simplify = "array")
            
            mean_result <- mean(result)
            results_list[[length(results_list) + 1]] <- data.table(
              test_type = test_type,
              test_name = test_name,
              n = n,
              h_label = h_label,
              estimator = est,
              rejection_rate = mean_result
            )
            
            cat("[Test]", test_name, "| n:", n, "| Estimator:", est, "|", h_label, "| Rejection Rate:", mean_result, "\n", strrep("-", 80), "\n")
          }
        } else {
          result <- parSapply(cl, 1:n_sims, function(sim, n, is_null, test_functions, test_name) {
            seed <- 1203 + sim
            set.seed(seed)
            
            # 데이터 생성
            d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
            set.seed(seed + n_sims)
            d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
            
            epsilon <- 1 / sqrt(log(n))
            test_args <- list(d1$x, d2$x, d1$y, d2$y, alg1 = TRUE, epsilon = epsilon, seed = seed)
            
            do.call(test_functions[[test_name]], test_args)
          }, n = n, is_null = is_null, test_functions = test_functions, test_name = test_name, simplify = "array")
          
          mean_result <- mean(result)
          results_list[[length(results_list) + 1]] <- data.table(
            test_type = test_type,
            test_name = test_name,
            n = n,
            h_label = h_label,
            estimator = NA,
            rejection_rate = mean_result
          )
          
          cat("[Test]", test_name, "| n:", n, "|", h_label, "| Rejection Rate:", mean_result, "\n", strrep("-", 80), "\n")
        }
      }
    }
  }
}

# 결과 병합
results_dt <- rbindlist(results_list)

# 클러스터 종료
stopCluster(cl)


# Save the results
filename <- paste0("results/simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")