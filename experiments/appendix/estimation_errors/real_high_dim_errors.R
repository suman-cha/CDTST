rm(list=ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(pbapply)
  library(MASS)
  library(dplyr)
  library(data.table)
  library(CVST)
  library(parallel)
})
source("utils.R")
tag <- "density_ratio_errors_real_high_dim"

cur_wd <- getwd()
file_path <- file.path(cur_wd, "real_examples", "data", "superconductivity.csv")
data <- read.csv(file_path)

X <- as.matrix(data[, !names(data) %in% c("critical_temp")])
Y <- data[, "critical_temp"]

normalize <- function(x){
  (x - min(x)) / (max(x) - min(x))
}

X_norm <- apply(X, 2, normalize)
Y_norm <- normalize(Y)

mse <- function(true, est) {
  mean((true - est)^2)
}

# True Marginal Density Ratio
true_marginal_density_ratio <- function(X1, X2, x_subset, is_x1 = TRUE) {
  prob_X1 <- rep(1 / nrow(X), nrow(X))  
  
  feature_to_bias <- X2[, 1]
  prob_X2 <- dnorm(feature_to_bias, mean = 0, sd = 1)
  prob_X2 <- prob_X2 / sum(prob_X2)
  
  if (is_x1) {
    indices <- match(x_subset[,1], X1[,1])
    ratio <- prob_X1[indices] / prob_X2[indices]
  } else {
    indices <- match(x_subset[,1], X2[,1])
    ratio <- prob_X1[indices] / prob_X2[indices]
  }
  
  return(ratio)
}


# True Conditional Density Ratio 
true_conditional_density_ratio <- function(Y_subset, is_null = TRUE, is_x1 = TRUE) {
  if (is_null) {
    prob_Y_given_X1 <- rep(1 / length(Y_subset), length(Y_subset))
    prob_Y_given_X2 <- rep(1 / length(Y_subset), length(Y_subset))
  } else {
    prob_Y_given_X1 <- dunif(Y_subset, min = 0, max = 1)
    prob_Y_given_X1 <- prob_Y_given_X1 / sum(prob_Y_given_X1)  
    
    prob_Y_given_X2 <- exp(-Y_subset)  
    prob_Y_given_X2 <- prob_Y_given_X2 / sum(prob_Y_given_X2)  
  }
  
  conditional_ratio <- prob_Y_given_X1 / prob_Y_given_X2
  
  return(conditional_ratio)
}

sample_data <- function(X, Y, n, is_null = TRUE, is_x1 = TRUE) {
  
  X_idx <- sample(1:nrow(X), nrow(X) %/% 2, replace = FALSE)
  X1 <- X[X_idx, , drop = FALSE]
  Y1 <- Y[X_idx]
  
  feature_to_bias <- X[, 1]
  prob <- dnorm(feature_to_bias, mean = 0, sd = 1)
  prob <- prob / sum(prob)
  
  X2_idx <- sample(1:nrow(X), nrow(X) %/% 2, replace = FALSE, prob = prob)
  X2 <- X[X2_idx, , drop = FALSE]
  Y2 <- Y[X2_idx]
  
  if (is_x1) {
    x_idx <- sample(1:nrow(X1), n, replace = FALSE)
    x <- X1[x_idx, , drop = FALSE]
    Y_subset <- Y1[x_idx]
  } else {
    x_idx <- sample(1:nrow(X2), n, replace = FALSE)
    x <- X2[x_idx, , drop = FALSE]
    Y_subset <- Y2[x_idx]
  }
  
  if (is_null) {
    y <- sample(Y_subset, size = n, replace = FALSE)  # Uniform sampling under null hypothesis
  } else {
    u <- if (is_x1) dunif(Y_subset, min = 0, max = 1) else exp(-Y_subset)
    u <- u / sum(u)
    y <- sample(Y_subset, size = n, prob = u, replace = FALSE)
  }
  
  return(list(x = x, y = y, X1 = X1, X2 = X2, y_subset = Y_subset))
}


estimate_r <- function(x11, x12, x21, x22, y11, y12, y21, y22, est.method = "LL", seed = NULL){
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  n11 <- length(y11); n12 <- length(y12)
  n21 <- length(y21); n22 <- length(y22)
  label.fit <- factor(c(rep(0, n11), rep(1, n21)))
  
  if (est.method == "LL"){
    xy.fit <- cbind(rbind(x11, x21), c(y11, y21))
    fit.joint <- glm(label.fit~., data=as.data.frame(xy.fit), family = binomial())
    x.fit <- rbind(x11, x21)
    fit.marginal <- glm(label.fit~., data=as.data.frame(x.fit), family = binomial())
    
    marginal_new_data <- rbind(x12, x22)
    prob.marginal <- predict(fit.marginal, newdata=as.data.frame(marginal_new_data), type="response")
    
    g12.est <- prob.marginal[1:n12]/(1-prob.marginal[1:n12])*n11/n21
    g22.est <- prob.marginal[(n12+1):(n12+n22)]/(1-prob.marginal[(n12+1):(n12+n22)])*n11/n21
    
    joint_new_data <- cbind(marginal_new_data, c(y12, y22))
    prob.joint <- predict(fit.joint, newdata=as.data.frame(joint_new_data), type="response")
    
    v12.est <- (1-prob.joint[1:n12])/prob.joint[1:n12]*g12.est
    v22.est <- (1-prob.joint[(n12+1):(n12+n22)])/prob.joint[(n12+1):(n12+n22)]*g22.est
  } else if (est.method == "KLR"){
    xy.fit <- cbind(rbind(x11, x21), c(y11, y21))
    data.fit <- constructData(xy.fit, label.fit)
    klrlearner <- constructKlogRegLearner()
    params <- list(kernel='rbfdot', sigma=0.005, lambda=0.0005, tol=1e-6, maxiter=500)
    fit.joint <- klrlearner$learn(data.fit, params)
    
    x.fit <- rbind(x11, x21)
    data.fit <- constructData(x.fit, label.fit)
    fit.marginal <- klrlearner$learn(data.fit, params)
    
    newdata <- rbind(x12, x22)
    K <- kernelMult(fit.marginal$kernel, newdata, fit.marginal$data, fit.marginal$alpha)
    pi <- 1 / (1 + exp(-as.vector(K)))
    
    g12.est <- pi[1:n12]/(1-pi[1:n12])*n11/n21
    g22.est <- pi[(n12+1):(n12+n22)]/(1-pi[(n12+1):(n12+n22)])*n11/n21
    
    newdata <- cbind(rbind(x12, x22), c(y12, y22))
    K <- kernelMult(fit.joint$kernel, newdata, fit.joint$data, fit.joint$alpha)
    pi <- 1 / (1 + exp(-as.vector(K)))
    
    v12.est <- (1-pi[1:n12])/pi[1:n12]*g12.est
    v22.est <- (1-pi[(n12+1):(n12+n22)])/pi[(n12+1):(n12+n22)]*g22.est
  }
  
  list(g12.est = g12.est, g22.est = g22.est, v12.est = v12.est, v22.est = v22.est)
}

run_simulation <- function(X, Y, n, is_null, estimator, seed) {
  set.seed(seed)
  
  d1 <- sample_data(X, Y, n, is_null, TRUE) 
  set.seed(seed + 500)
  d2 <- sample_data(X, Y, n, is_null, FALSE)
  
  split_idx <- sample(1:n, n/2)
  x11 <- d1$x[split_idx, ]      
  y11 <- d1$y[split_idx]        
  x12 <- d1$x[-split_idx, ]     
  y12 <- d1$y[-split_idx]       
  x21 <- d2$x[split_idx, ]      
  y21 <- d2$y[split_idx]        
  x22 <- d2$x[-split_idx, ]     
  y22 <- d2$y[-split_idx]       
  
  true_g12 <- true_marginal_density_ratio(d1$X1, d2$X2, x12)
  true_g22 <- true_marginal_density_ratio(d1$X1, d2$X2, x22)
  
  true_v12 <- true_conditional_density_ratio(d1$y_subset, is_null = is_null, is_x1 = TRUE)
  true_v22 <- true_conditional_density_ratio(d2$y_subset, is_null = is_null, is_x1 = FALSE)
  
  estimated_r <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est.method = estimator)
  
  # MSE
  mse_g12 <- mse(true_g12, estimated_r$g12.est)
  mse_g22 <- mse(true_g22, estimated_r$g22.est)
  mse_v12 <- mse(true_v12, estimated_r$v12.est)
  mse_v22 <- mse(true_v22, estimated_r$v22.est)
  
  return(c(mse_g12 = mse_g12, mse_g22 = mse_g22, mse_v12 = mse_v12, mse_v22 = mse_v22))
}

# Set up parallel processing
cl <- makeCluster(detectCores() - 1)
pbapply::pboptions(cl = cl)

# Parameter settings
n_values <- c(200, 400, 800, 1200, 1600, 2000)
n_sims <- 500
estimators <- c("LL", "KLR")
results_list <- list() 

for (n in n_values){
  for (is_null in c(TRUE, FALSE)){
    h_label <- if(is_null) "Null" else "Alternative"
    
    for (est in estimators){
      result <- pbapply::pbsapply(1:n_sims, function(sim) {
        seed <- 1203 + sim 
        set.seed(seed)
        
        run_simulation(X_norm, Y_norm, n, is_null, est, seed)
      }, simplify = "array")
      
      median_result <- apply(result, MARGIN = 1, median)
      
      results_list[[length(results_list) + 1]] <- data.table(
        n = n,
        h_label = h_label,
        estimator = est,
        mse_g12 = median_result["mse_g12"],
        mse_g22 = median_result["mse_g22"],
        mse_v12 = median_result["mse_v12"],
        mse_v22 = median_result["mse_v22"]
      )
      
      # Print results
      cat("[Estimator]", est, "| n:", n, "|", h_label, "| MSE (g12, g22, v12, v22):", 
          median_result["mse_g12"], median_result["mse_g22"], median_result["mse_v12"], median_result["mse_v22"], "\n", strrep("-", 80), "\n")
    }
  }
}


results_dt <- rbindlist(results_list)
stopCluster(cl)

# Save the results
filename <- paste0("results/ablation_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")