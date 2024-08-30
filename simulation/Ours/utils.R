# Install the required packages
required_pkgs <- c("glmnet",
                   "caret",
                   "ranger",
                   "mgcv",
                   "nnet",
                   "xgboost",
                   "brulee",
                   "recipes",
                   "kernlab",
                   "SuperLearner",
                   "CVST",
                   "densratio")

install_pkgs <- function(pkgs){
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs)){
    install.packages(new_pkgs, repos="http://cran.us.r-project.org")
  }
}

install_pkgs(required_pkgs)
source("FunFiles.R")

# Gaussian kernel 
gaussian.kernel <- function(x, y = NULL, h=1) {
  if (is.null(y)) {
    res <- (exp(-0.5 * (x/h)^2) / (h * sqrt(2 * pi)))
  } else {
    dist <- sum((x - y)^2)
    res <- exp(-dist / (2 * h^2))
  }
  return(res)
}

# Laplacian kernel
laplacian.kernel <- function(x, y, h=1){
  dist <- sum(abs(x-y))
  res <- exp(-dist / h)
  return(res)
}

# Linear kernel function
linear.kernel <- function(x, y) {
  return(sum(x * y))
}

# Polynomial kernel function
polynomial.kernel <- function(x, y, degree=3, alpha=1, c=0) {
  return((alpha * sum(x * y) + c)^degree)
}

delta.kernel <- function(x, y = NULL) {
  if (is.null(y)) {
    stop("Delta kernel requires two inputs for comparison.")
  } else {
    res <- as.numeric(x != y)
  }
  return(res)
}

# Median heuristic
median_heuristic <- function(x12, x22){
  x <- rbind(x12, x22)
  dists <- as.vector(dist(x))
  return(median(dists))
}

cv_bandwidth <- function(x12, x22) {
  # Combine x12 and x22 into a single matrix
  X <- rbind(x12, x22)
  
  # Perform cross-validation to select the optimal bandwidth
  h_values <- seq(0.1, 2, by=0.1)
  best_h <- h_values[1]
  best_score <- Inf
  
  for (h in h_values) {
    score <- 0
    for (i in 1:nrow(X)) {
      x_train <- X[-i, , drop=FALSE]
      x_test <- X[i, , drop=FALSE]
      dist <- rowSums((x_train - matrix(x_test, nrow(x_train), ncol(x_train), byrow=TRUE))^2)
      score <- score + mean(exp(-dist / (2 * h^2)))
    }
    if (score < best_score) {
      best_score <- score
      best_h <- h
    }
  }
  return(best_h)
}

# Conditional Energy Distance(CED) 
CED <- function(x1, x2, y1, y2, h=1){
  n1 <- length(y1); n2 <- length(y2)
  x1 <- matrix(x1, nrow=n1)
  x2 <- matrix(x2, nrow=n2)
  
  # kernel density estimation 
  Kx1 <- sapply(x1, function(x) rowSums(gaussian.kernel(as.matrix(outer(x, x1, "-"), h))))
  Kx2 <- sapply(x2, function(x) rowSums(gaussian.kernel(as.matrix(outer(x, x2, "-"), h))))
  
  # U-stats for CED 
  U_stat <- function(Kx1, Kx2, y1, y2){
    n1 <- length(y1); n2 <- length(y2)
    
    term1 <- 2/(n1*n2) * sum(outer(y1, y2, function(a,b) abs(a-b)) * Kx1 * Kx2)
    term2 <- 1/n1^2 * sum(outer(y1, y1, function(a,b) abs(a-b)) * Kx1 * Kx1)
    term3 <- 1/n2^2 * sum(outer(y2, y2, function(a,b) abs(a-b)) * Kx2 * Kx2)
    
    return (term1  - term2 - term3)
  }
}

# Conformity score function
conformity_score <- function(X, U){
  as.numeric(rowSums(X) < rowSums(U))
}

# Linear MMD
MMDl <- function(x12, x22, y12, y22, h_x=1, h_y=1, r_hat, kernel.type="gaussian", seed=NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # x12, x22, y12, y22 are splitted data 2n -> n
  stopifnot(length(y12) == length(y22))
  n <- length(y12)
  m <- floor(n/2)
  S_hat_values <- numeric(m)

  for(i in 1:m){
    if (kernel.type == "discrete") { 
      k_zz <- gaussian.kernel(x12[i, ], x12[i+m, ], h_x) * delta.kernel(y12[i], y12[i+m])
      k_ww <- gaussian.kernel(x22[i, ], x22[i+m, ], h_x) * delta.kernel(y22[i], y22[i+m])
      k_wz <- gaussian.kernel(x22[i, ], x12[i+m, ], h_x) * delta.kernel(y22[i], y12[i+m])
      k_zw <- gaussian.kernel(x12[i, ], x22[i+m, ], h_x) * delta.kernel(y12[i], y22[i+m])
    } else if (kernel.type == "gaussian") {
      k_zz <- gaussian.kernel(x12[i, ], x12[i+m, ], h_x) * gaussian.kernel(y12[i], y12[i+m], h_y)
      k_ww <- gaussian.kernel(x22[i, ], x22[i+m, ], h_x) * gaussian.kernel(y22[i], y22[i+m], h_y)
      k_wz <- gaussian.kernel(x22[i, ], x12[i+m, ], h_x) * gaussian.kernel(y22[i], y12[i+m], h_y)
      k_zw <- gaussian.kernel(x12[i, ], x22[i+m, ], h_x) * gaussian.kernel(y12[i], y22[i+m], h_y)
    } else if (kernel.type == "linear") { 
      k_zz <- linear.kernel(x12[i, ], x12[i+m, ]) * linear.kernel(y12[i], y12[i+m])
      k_ww <- linear.kernel(x22[i, ], x22[i+m, ]) * linear.kernel(y22[i], y22[i+m])
      k_wz <- linear.kernel(x22[i, ], x12[i+m, ]) * linear.kernel(y22[i], y12[i+m])
      k_zw <- linear.kernel(x12[i, ], x22[i+m, ]) * linear.kernel(y12[i], y22[i+m])
    } else if (kernel.type == "polynomial") {
      k_zz <- polynomial.kernel(x12[i, ], x12[i+m, ]) * polynomial.kernel(y12[i], y12[i+m])
      k_ww <- polynomial.kernel(x22[i, ], x22[i+m, ]) * polynomial.kernel(y22[i], y22[i+m])
      k_wz <- polynomial.kernel(x22[i, ], x12[i+m, ]) * polynomial.kernel(y22[i], y12[i+m])
      k_zw <- polynomial.kernel(x12[i, ], x22[i+m, ]) * polynomial.kernel(y12[i], y22[i+m])
    } else {
      stop("Unsupported kernel type")
    }
    
    S_hat_values[i] <- k_zz + r_hat[i]*r_hat[i+m]*k_ww - r_hat[i]*k_wz - r_hat[i+m]*k_zw
  }
  S_bar <- mean(S_hat_values)
  sigma_hat <- sum((S_hat_values - S_bar)^2) / (m - 1)
  MMDl_hat2 <- sqrt(m) * S_bar / sqrt(sigma_hat)
  return(MMDl_hat2)
}

estimate_r <- function(x11, x12, x21, x22, y11, y12, y21, y22, est.method="LL", seed=NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n11 <- length(y11); n12 <- length(y12)
  n21 <- length(y21); n22 <- length(y22)
  label.fit <- as.factor(c(rep(0,n11), rep(1,n21)))
  
  if (est.method == "LL"){
    xy.fit <- data.frame(x.fit=rbind(x11, x21), y.fit=c(y11, y21))
    fit.joint <- glm(label.fit~., data=xy.fit, family="binomial")
    x.fit <- data.frame(x.fit=rbind(x11,x21))
    fit.marginal <- glm(label.fit~., data=x.fit, family="binomial")
    prob.marginal <- predict(fit.marginal, newdata=data.frame(x.fit=rbind(x12,x22)), type="response")
    prob.marginal[prob.marginal<0.01] <- 0.01; prob.marginal[prob.marginal>0.99] <- 0.99
    g12.est <- prob.marginal[1:n12]/(1-prob.marginal[1:n12])*n11/n21
    g22.est <- prob.marginal[(n12+1):(n12+n22)]/(1-prob.marginal[(n12+1):(n12+n22)])*n11/n21
    prob.joint <- predict(fit.joint, newdata=data.frame(x.fit=rbind(x12,x22), y.fit=c(y12,y22)), type="response")
    prob.joint[prob.joint<0.01] <- 0.01; prob.joint[prob.joint>0.99] <- 0.99
    v12.est <- (1-prob.joint[1:n12])/prob.joint[1:n12]*g12.est
    v22.est <- (1-prob.joint[(n12+1):(n12+n22)])/prob.joint[(n12+1):(n12+n22)]*g22.est
    
  } else if (est.method == "QL"){ 
    xy.fit = cbind(rbind(x11,x21), c(y11,y21))
    xy.fit = data.frame(poly(as.matrix(xy.fit), degree = 2, raw = T))
    fit.joint = glm(label.fit~., data=xy.fit, family="binomial")
    x.fit <- rbind(x11,x21)
    x.fit = data.frame(poly(as.matrix(x.fit), degree = 2, raw = T))
    fit.marginal <- glm(label.fit~., data=x.fit, family="binomial")
    newdata = data.frame(poly(as.matrix(rbind(x12,x22)), degree = 2, raw = T))
    prob.marginal <- predict(fit.marginal, newdata=newdata, type="response")
    prob.marginal[prob.marginal < 0.01] <- 0.01
    prob.marginal[prob.marginal > 0.99] <- 0.99
    
    g12.est <- prob.marginal[1:n12] / (1 - prob.marginal[1:n12]) * n11 / n21
    g22.est <- prob.marginal[(n12 + 1):(n12 + n22)] / (1 - prob.marginal[(n12 + 1):(n12 + n22)]) * n11 / n21
    
    newdata = data.frame(poly(as.matrix(cbind(rbind(x12, x22), c(y12, y22))), degree = 2, raw = T))
    prob.joint <- predict(fit.joint, newdata = newdata, type = "response")
    prob.joint[prob.joint < 0.01] <- 0.01
    prob.joint[prob.joint > 0.99] <- 0.99
    
    v12.est <- (1 - prob.joint[1:n12]) / prob.joint[1:n12] * g12.est
    v22.est <- (1 - prob.joint[(n12 + 1):(n12 + n22)]) / prob.joint[(n12 + 1):(n12 + n22)] * g22.est
    
  } else if (est.method == "KLR"){
    xy.fit <- cbind(rbind(x11, x21), c(y11, y21))
    data.fit <- constructData(xy.fit, label.fit)
    klrlearner <- constructKlogRegLearner()
    params <- list(kernel='rbfdot', sigma=0.005, lambda=0.05/getN(data.fit), tol=10e-6, maxiter=500)
    fit.joint <- klrlearner$learn(data.fit, params)
    x.fit <- rbind(x11, x21)
    data.fit <- constructData(x.fit, label.fit)
    klrlearner <- constructKlogRegLearner()
    params <- list(kernel='rbfdot', sigma=0.005, lambda=0.05/getN(data.fit), tol=10e-6, maxiter=500)
    fit.marginal <- klrlearner$learn(data.fit, params)
    
    newdata <- rbind(x12, x22)
    K = kernelMult(fit.marginal$kernel, newdata, fit.marginal$data, fit.marginal$alpha)
    pi = 1 / (1 + exp(-as.vector(K))) # predicted probabilities
    pi[pi<0.01] <- 0.01; pi[pi>0.99] <- 0.99
    g12.est <- pi[1:n12]/(1-pi[1:n12])*n11/n21
    g22.est <- pi[(n12+1):(n12+n22)]/(1-pi[(n12+1):(n12+n22)])*n11/n21
    newdata <- cbind(rbind(x12, x22), c(y12, y22))
    K = kernelMult(fit.joint$kernel, newdata, fit.joint$data, fit.joint$alpha)
    pi = 1 / (1 + exp(-as.vector(K))) # predicted probabilities
    pi[pi<0.01] <- 0.01; pi[pi>0.99] <- 0.99
    v12.est <- (1-pi[1:n12])/pi[1:n12]*g12.est
    v22.est <- (1-pi[(n12+1):(n12+n22)])/pi[(n12+1):(n12+n22)]*g22.est
    
  } else if (est.method == "NN"){
    hidden.layers <- c(10,10)
    learn.rates <- 0.001
    n.epochs <- 500
    x.fit <- data.frame(x=rbind(x11, x21))
    newdata1 <- data.frame(x=x12)
    newdata2 <- data.frame(x=x22)
    temp <- NNfun(x.fit, label.fit, newdata1, newdata2, nnrep = 5, hidden.layers = hidden.layers,
                  n.epochs = n.epochs, learn.rates = learn.rates)
    g12.est <- temp$prob1.fit/(1-temp$prob1.fit)*n11/n21
    g22.est <- temp$prob2.fit/(1-temp$prob2.fit)*n11/n21
    
    xy.fit <- data.frame(x=rbind(x11, x21), y=c(y11, y21))
    newdata1 <- data.frame(x=x12, y=y12)
    newdata2 <- data.frame(x=x22, y=y22)
    temp <- NNfun(xy.fit, label.fit, newdata1, newdata2, nnrep = 5, hidden.layers = hidden.layers,
                  n.epochs = n.epochs, learn.rates = learn.rates)
    v12.est <- (1-temp$prob1.fit)/temp$prob1.fit*g12.est
    v22.est <- (1-temp$prob2.fit)/temp$prob2.fit*g22.est
  } else if (est.method == "KLIEP"){
    xy1 <- data.frame(x=x11, y=y11)
    xy2 <- data.frame(x=x21, y=y21)
    fit.joint <- densratio(xy2, xy1, method="KLIEP", verbose=FALSE) # q(x,y)/p(x,y)
    fit.marginal <- densratio(x21, x11, method="KLIEP", verbose=FALSE) # q(x)/p(x)
    g12.est <- fit.marginal$compute_density_ratio(x12)
    g22.est <- fit.marginal$compute_density_ratio(x22)
    v12.est <- fit.joint$compute_density_ratio(data.frame(x=x12, y=y12))*g12.est
    v22.est <- fit.joint$compute_density_ratio(data.frame(x=x22, y=y22))*g22.est
  }
  return(list(g12.est = g12.est, g22.est = g22.est, v12.est=v12.est, v22.est=v22.est))
  
}


# From Chen and Lei (2024)

# Estimates marginal density ratio given probability
# @param eta estimated probability of being in class 1
# @param n0 sample size of class 0
# @param n1 sample size of class 1
marg <- function(eta,n0,n1){
  return((n0/n1)*(eta/(1-eta)))
}

# Function to compute joint density ratio given probability
joint <- function(eta) {
  return((1 - eta) / eta)
}

# Returns a function that computes estimated marginal density ratio at a point
# @param data data to estimate density ratio
# @param n0 number of points from class 0
# @param n1 number of points from class 1
# @param type "ld" for low dimensional "hd" for high dimensional
estimate_marginal_ratio <- function(data, n0, n1, type) {
  if(type == "LL") {
    # Linear Logistic Regression
    model <- glm(class ~. -y, data = data, family = binomial())
    marg_ratio <- function(x) {
      eta <- predict(model, newdata = x, type = "response")
      eta[eta < 0.01] <- 0.01
      eta[eta > 0.99] <- 0.99
      return(sapply(eta, function(x) { return(marg(x, n0, n1)) }))
    }
  } else if(type == "QL") {
    # Quadratic Logistic Regression
    data_poly <- data.frame(poly(as.matrix(data[, !names(data) %in% "class"]), degree = 2, raw = TRUE))
    data_poly$class <- data$class
    model <- glm(class ~ ., data = data_poly, family = "binomial")
    
    marg_ratio <- function(x) {
      x_poly <- data.frame(poly(as.matrix(x[, !names(x) %in% "class"]), degree = 2, raw = TRUE))
      eta <- predict(model, newdata = x_poly, type = "response")
      eta[eta < 0.01] <- 0.01
      eta[eta > 0.99] <- 0.99
      return(sapply(eta, function(x) { return(marg(x, n0, n1)) }))
    }
  } else if(type == "KLR") {
    # Kernel Logistic Regression using kernelMult and kernellearner
    xy.fit <- cbind(as.matrix(data[, !names(data) %in% "class"]))
    label.fit <- as.factor(data$class)
    
    # Construct the KLR learner
    klrlearner <- constructKlogRegLearner()
    params <- list(kernel = 'rbfdot', sigma = 0.005, lambda = 0.05 / nrow(data), 
                   tol = 1e-6, maxiter = 500)
    
    # Train the model
    fit.marginal <- klrlearner$learn(constructData(xy.fit, label.fit), params)
    
    marg_ratio <- function(x) {
      # Predict using the learned model
      new_x <- as.matrix(x[, !names(x) %in% "class"])
      K <- kernelMult(fit.marginal$kernel, new_x, fit.marginal$data, fit.marginal$alpha)
      pi <- 1 / (1 + exp(-as.vector(K)))  # predicted probabilities
      pi[pi < 0.01] <- 0.01
      pi[pi > 0.99] <- 0.99
      return(sapply(pi, function(x) { return(marg(x, n0, n1)) }))
    }
  } else if(type == "superlearner") {
    # SuperLearner approach
    xgb_grid <- create.SL.xgboost(tune = list(ntrees = c(100, 200, 500, 1000), 
                                              max_depth = c(2, 4, 6), 
                                              shrinkage = 0.3, 
                                              minobspernode = 1), 
                                  detailed_names = FALSE, env = .GlobalEnv,
                                  name_prefix = "SL.xgb")
    model <- SuperLearner(Y = data$class, 
                          X = data[, !names(data) %in% c("class", "y")],
                          SL.library = c("SL.ranger", "SL.lm", "SL.ksvm", xgb_grid$names), 
                          family = binomial())
    
    marg_ratio <- function(x) {
      eta <- predict(model, x[, !names(x) %in% c("class", "y")], onlySL = TRUE)$pred
      eta[eta < 0.01] <- 0.01
      eta[eta > 0.99] <- 0.99
      return(sapply(eta, function(x) { return(marg(x, n0, n1)) }))
    }
  } else if(type == "nn") {
    # Neural Network approach
    rec <- recipe(class ~., data = data[,-5]) %>%
      step_normalize(all_numeric_predictors())
    model <- brulee_mlp(rec, data = data, epochs = 15000, hidden_units = c(10, 10), learn_rate = 0.001, verbose = TRUE)
    
    marg_ratio <- function(x) {
      eta <- predict(model, x[, 1:4], type = "prob")$.pred_1
      eta[eta < 0.01] <- 0.01
      eta[eta > 0.99] <- 0.99
      return(sapply(eta, function(x) { return(marg(x, n0, n1)) }))
    }
  }
  return(marg_ratio)
}

# Function to estimate joint density ratio
estimate_joint_ratio <- function(data, type) {
  if(type == "LL") {
    # Linear Logistic Regression
    model <- glm(class ~ ., data = data, family = binomial())
    joint_ratio <- function(point) {
      eta <- predict(model, newdata = point, type = "response")
      eta[eta < 0.01] <- 0.01
      eta[eta > 0.99] <- 0.99
      return(sapply(eta, joint))
    }
  } else if(type == "QL") {
    # Quadratic Logistic Regression
    data_poly <- data.frame(poly(as.matrix(data[, !names(data) %in% "class"]), degree = 2, raw = TRUE))
    data_poly$class <- data$class
    model <- glm(class ~ ., data = data_poly, family = "binomial")
    
    joint_ratio <- function(point) {
      point_poly <- data.frame(poly(as.matrix(point[, !names(point) %in% "class"]), degree = 2, raw = TRUE))
      eta <- predict(model, newdata = point_poly, type = "response")
      eta[eta < 0.01] <- 0.01
      eta[eta > 0.99] <- 0.99
      return(sapply(eta, joint))
    }
  } else if(type == "KLR") {
    # Kernel Logistic Regression using kernelMult and kernellearner
    xy.fit <- cbind(as.matrix(data[, !names(data) %in% "class"]))
    label.fit <- as.factor(data$class)
    
    # Construct the KLR learner
    klrlearner <- constructKlogRegLearner()
    params <- list(kernel = 'rbfdot', sigma = 0.005, lambda = 0.05 / nrow(data), 
                   tol = 1e-6, maxiter = 500)
    
    # Train the joint model
    fit.joint <- klrlearner$learn(constructData(xy.fit, label.fit), params)
    
    joint_ratio <- function(point) {
      # Predict using the learned joint model
      new_point <- as.matrix(point[, !names(point) %in% "class"])
      K <- kernelMult(fit.joint$kernel, new_point, fit.joint$data, fit.joint$alpha)
      pi <- 1 / (1 + exp(-as.vector(K)))  # predicted probabilities
      pi[pi < 0.01] <- 0.01
      pi[pi > 0.99] <- 0.99
      return(sapply(pi, joint))
    }
  } else if(type == "superlearner") {
    # SuperLearner approach
    xgb_grid <- create.SL.xgboost(tune = list(ntrees = c(100, 200, 500), 
                                              max_depth = c(2, 6), 
                                              shrinkage = 0.3, 
                                              minobspernode = 1), 
                                  detailed_names = FALSE, env = .GlobalEnv,
                                  name_prefix = "SL.xgb")
    model <- SuperLearner(Y = data$class, 
                          X = data[, !names(data) %in% "class"],
                          SL.library = c("SL.ranger", "SL.lm", xgb_grid$names),
                          family = binomial())
    
    joint_ratio <- function(point) {
      eta <- predict(model, point[, !names(point) %in% "class"], onlySL = TRUE)$pred
      eta[eta < 0.01] <- 0.01
      eta[eta > 0.99] <- 0.99
      return(sapply(eta, joint))
    }
  }
  return(joint_ratio)
}


