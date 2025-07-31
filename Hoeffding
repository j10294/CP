########################################################################
# Tolerance Interval Simulation with Risk Control and Localization
# 
# This script includes:
# - Data generation for various noise models
# - Localized quantile computation
# - Risk-controlled lambda estimation
# - Content-based interval evaluation
# - Monte Carlo simulation for pointwise coverage
# - Visualization of coverage and NA proportion
#
# Author: Jisu Kim
# Date: 20250713
########################################################################


########## Library Dependencies ##########
library(gam)
library(dplyr)
library(ggplot2)

########## Function Definitions ##########

# Function to generate synthetic data according to 8 different noise models
generate_data <- function(model, n) {
  x <- seq(0, 10, length.out=n)
  base <- 3 * cos(x) - 5 * (x / 15)^2

  if (model == 1) y <- base + rnorm(n, 0, sqrt(2))
  if (model == 2) y <- base + rnorm(n, 0, sqrt(2 * x))
  if (model == 3) y <- base + 4 - rgamma(n, shape = 2, scale = 2)
  if (model == 4) y <- base + rt(n, df = 3)
  if (model == 5) {
    sd_x <- 0.3 + 3 * exp(-(x - 5)^2 / (2 * 1.5^2))
    y <- base + rnorm(n, 0, sd_x)
  }
  if (model == 6) y <- base + rnorm(n, 0, 1) + abs(rnorm(n, 0, 1)) - sqrt(2/pi)
  if (model == 7) {
    noise <- ifelse(runif(n) < 0.5, rnorm(n, -3, 1), rnorm(n, 3, 1))
    y <- base + noise
  }
  if (model == 8) {
    noise <- ifelse(x < 5, rnorm(n, 0, 1 + 0.5 * x), rnorm(n, 0, 2 - 0.3 * (x - 5)))
    y <- base + noise
  }
  return(data.frame(x = x, y = y))
}

# Function to compute theoretical content (true coverage probability)
content_function <- function(model, lower, upper, x) {
  y <- 3 * cos(x) - 5 * (x / 15)^2
  if (model == 1) return(pnorm(upper - y, 0, sqrt(2)) - pnorm(lower - y, 0, sqrt(2)))
  if (model == 2) return(pnorm(upper - y, 0, sqrt(2 * x)) - pnorm(lower - y, 0, sqrt(2 * x)))
  if (model == 3) return(pgamma(4 - (lower - y), 2, 2) - pgamma(4 - (upper - y), 2, 2))
  if (model == 4) return(pt(upper - y, df = 3) - pt(lower - y, df = 3))
  if (model == 5) {
    sd_x <- 0.3 + 3 * exp(-(x - 5)^2 / (2 * 1.5^2))
    return(pnorm(upper - y, 0, sd_x) - pnorm(lower - y, 0, sd_x))
  }
  if (model == 6) {
    N <- 1e5
    z1 <- rnorm(N)
    z2 <- abs(rnorm(N)) - sqrt(2 / pi)
    eps <- z1 + z2
    return(mean((eps >= lower - y) & (eps <= upper - y)))
  }
  if (model == 7) {
    N <- 1e5
    u <- runif(N)
    eps <- ifelse(u < 0.5, rnorm(N, -3, 1), rnorm(N, 3, 1))
    return(mean((eps >= lower - y) & (eps <= upper - y)))
  }
  if (model == 8) {
    sd_x <- ifelse(x < 5, 1 + 0.5 * x, 2 - 0.3 * (x - 5))
    return(pnorm(upper - y, 0, sd_x) - pnorm(lower - y, 0, sd_x))
  }
}

# Lambda estimation using empirical risk and truncation
find_lambda_hat <- function(alpha, delta, y, pred, variance, n,
                            lower = 0, upper = 30, step = 1e-3) {
  threshold <- max(0, alpha - sqrt(log(1 / delta) / (2 * n)))
  lambda_seq <- seq(lower, upper, by = step)
  for (lambda in lambda_seq) {
    risk <- mean(y <= pred - sqrt(variance) * lambda |
                 y >= pred + sqrt(variance) * lambda, na.rm = TRUE)
    if (risk <= threshold) return(lambda)
  }
  return(NA)
}

# Lambda estimation using localization
find_lambda_local <- function(x0, x_cal, y_cal, pred_fun, var_fun,
                              alpha, delta, kernel = "gaussian", h = 0.5,
                              fallback_lambda = NA) {
  dist <- abs(x_cal - x0)
  if (kernel == "gaussian") {
    w <- exp(-dist^2 / (2 * h^2))
  } else {
    w <- as.numeric(dist <= h)
  }
  w <- w / sum(w)
  
  eff_w2 <- sum(w^2)
  if (eff_w2 < 1e-6) return(fallback_lambda)
  
  threshold <- alpha - sqrt(log(1 / delta) / (2 * eff_w2))
  threshold <- max(threshold, 0)
  
  lambda_seq <- seq(0, 30, by = 0.01)
  pred <- pred_fun(x_cal)
  var_hat <- var_fun(x_cal)
  
  for (lambda in lambda_seq) {
    lower <- pred - lambda * sqrt(var_hat)
    upper <- pred + lambda * sqrt(var_hat)
    risk <- sum(w * ((y_cal < lower) | (y_cal > upper)))
    if (risk <= threshold) return(lambda)
  }
  return(fallback_lambda)
}


########## Simulation ##########

models <- 1:8
sample_sizes <- c(50, 100, 200)
alpha <- 0.10
delta <- 0.05
M <- 100

pointwise_coverage_result <- list()
lambda_na_count_result <- list()

for (model in models) {
  for (m in sample_sizes) {
    start_time <- Sys.time()
    
    content_matrix <- matrix(NA, nrow = m, ncol = M)
    lambda_na_mat <- matrix(0, nrow = m, ncol = M)
    
    for (i in 1:M) {
      set.seed(i)
      data <- generate_data(model, m)
      x <- data$x
      y <- data$y
      
      n_train <- m / 2
      n_cal <- m / 2
      
      idx <- sample(m)
      train_index <- idx[1:n_train]
      cal_index <- idx[(n_train + 1):m]
      
      train_x <- x[train_index]
      train_y <- y[train_index]
      cal_x <- x[cal_index]
      cal_y <- y[cal_index]
      
      fit <- smooth.spline(train_x, train_y)
      train_res <- train_y - predict(fit, train_x)$y
      log_res2 <- log(train_res^2 + 1e-6)
      gam_model <- gam(log_res2 ~ s(train_x))
      
      for (k in 1:n_cal) {
        x0 <- cal_x[k]
        
        lambda_hat <- find_lambda_local(
          x0 = x0,
          x_cal = cal_x,
          y_cal = cal_y,
          pred_fun = function(xx) predict(fit, xx)$y,
          var_fun = function(xx) exp(predict(gam_model, newdata = data.frame(train_x = xx))),
          alpha = alpha,
          delta = delta,
          h = 0.5,
          fallback_lambda = NA
        )
        
        
        mu_hat <- predict(fit, x0)$y
        var_hat <- exp(predict(gam_model, newdata = data.frame(train_x = x0)))
        TI_lower <- mu_hat - sqrt(var_hat) * lambda_hat
        TI_upper <- mu_hat + sqrt(var_hat) * lambda_hat
        j <- which(x == x0)
        
        content_matrix[j, i] <- content_function(model, TI_lower, TI_upper, x0)
      }
    }
    
    na_proportion <- rowMeans(lambda_na_mat)
    coverage_proportion <- rowMeans(content_matrix >= 0.90, na.rm = TRUE)
    
    pointwise_coverage_result[[paste0("Model_", model, "_m_", m)]] <- data.frame(
      x = x,
      coverage = coverage_proportion,
      na_proportion = na_proportion
    )
    
    duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat("[END] Model:", model, "Sample size:", m, "—", round(duration, 2), "seconds\n")
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)


coverage_df <- bind_rows(lapply(names(pointwise_coverage_result), function(name) {
  df <- pointwise_coverage_result[[name]]
  df$model <- str_extract(name, "Model_\\d+") %>% str_remove("Model_") %>% as.integer()
  df$m <- str_extract(name, "m_\\d+") %>% str_remove("m_") %>% as.integer()
  return(df)
}), .id = "key")

coverage_summary <- coverage_df %>%
  group_by(model, m) %>%
  summarize(mean_coverage = mean(coverage, na.rm = TRUE), .groups = "drop")


# marginal 시각화
p <- ggplot(coverage_summary, aes(x = m, y = mean_coverage, group = model, color = as.factor(model))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = unique(coverage_summary$m)) +
  labs(
    title = "Coverage by Model and Sample Size (m)",
    x = "Sample Size (m)",
    y = "Mean Coverage",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  ) + ylim(0,1); p

# pointwise coverage 시각화
p2 <- ggplot(coverage_df, aes(x = x, y = coverage)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  facet_grid(rows = vars(model), cols = vars(m), labeller = label_both) +
  labs(
    title = "Pointwise Coverage Proportion for Hoeffding (content ≥ 0.90)",
    x = "x", y = "Coverage Proportion (for 100 Simulations)"
  ) + ylim(0, 1)+
  theme_minimal(base_size = 13); p2


# NA 비율 시각화 (lambda_hat 실패율)
ggplot(coverage_df, aes(x = x, y = na_proportion)) +
  geom_line(color = "gray40", size = 1) +
  facet_grid(rows = vars(model), cols = vars(m), labeller = label_both) +
  labs(
    title = "Proportion of NA (λ_hat estimation failure)",
    x = "x", y = "NA Proportion"
  ) +
  theme_minimal(base_size = 13)
