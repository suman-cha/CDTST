# Load necessary libraries
library(ggplot2)
library(dplyr)

# Prepare data
data <- data.frame(
  test = rep(c("LinearMMD_LL", "LinearMMD_QL", "CV_LinearMMD_LL", "CV_LinearMMD_QL", 
               "CP_LL", "CP_QL", "PCM_alg1_TRUE", "PCM_alg1_FALSE", 
               "GCM_alg1_TRUE", "GCM_alg1_FALSE", "WGCM_alg1_TRUE", "WGCM_alg1_FALSE"), times = 5),
  sample_size = rep(c(400, 800, 1000, 1600, 2000), each = 12),
  value = c(0.058, 0.06, 0.046, 0.056, 0.05, 0.398, 0.53, 0.622, 0.144, 0.172, 0.452, 0.512,
            0.046, 0.052, 0.05, 0.058, 0.038, 0.674, 0.912, 0.924, 0.306, 0.34, 0.818, 0.858,
            0.04, 0.038, 0.076, 0.066, 0.052, 0.78, 0.946, 0.984, 0.36, 0.38, 0.912, 0.938,
            0.062, 0.06, 0.048, 0.054, 0.038, 0.932, 1, 0.998, 0.624, 0.646, 0.996, 0.998,
            0.052, 0.05, 0.052, 0.06, 0.048, 0.986, 1, 1, 0.708, 0.746, 0.998, 1)
)

# Plot all results in one graph
ggplot(data, aes(x = sample_size, y = value, color = test)) +
  geom_line() +
  geom_point() +
  labs(title = "Alternative", x = "Sample Size", y = "Average Rejection Rate") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ylim(0,1)

# Prepare data for null hypothesis
data_null <- data.frame(
  test = rep(c("LinearMMD_LL", "LinearMMD_QL", "CV_LinearMMD_LL", "CV_LinearMMD_QL", 
               "CP_LL", "CP_QL", "PCM_alg1_TRUE", "PCM_alg1_FALSE", 
               "GCM_alg1_TRUE", "GCM_alg1_FALSE", "WGCM_alg1_TRUE", "WGCM_alg1_FALSE"), times = 5),
  sample_size = rep(c(400, 800, 1000, 1600, 2000), each = 12),
  value = c(0.06, 0.062, 0.052, 0.05, 0.086, 0.084, 0.096, 0.066, 0.066, 0.074, 0.084, 0.096,
            0.05, 0.054, 0.05, 0.058, 0.086, 0.14, 0.068, 0.072, 0.048, 0.052, 0.054, 0.054,
            0.04, 0.036, 0.076, 0.056, 0.114, 0.166, 0.066, 0.054, 0.048, 0.05, 0.064, 0.064,
            0.058, 0.06, 0.058, 0.052, 0.142, 0.182, 0.072, 0.074, 0.066, 0.056, 0.06, 0.07,
            0.058, 0.054, 0.058, 0.056, 0.2, 0.264, 0.058, 0.076, 0.05, 0.046, 0.044, 0.048),
  hypothesis = "Null"
)

# Plot all results in one graph
ggplot(data_null, aes(x = sample_size, y = value, color = test)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +
  labs(title = "Null", x = "Sample Size", y = "Average Rejection Rate") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ylim(0, 1)
