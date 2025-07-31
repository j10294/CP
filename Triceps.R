########################################################################
# Triceps Dataset - Risk-based Tolerance Interval Estimation
# Including 'find_lambda' function for Hoeffding-style risk control
########################################################################

#-------------------------------------------------------------#
# Lambda Estimation Function
# - Uses Hoeffding-style inequality to find smallest lambda 
#   that satisfies the desired risk control level.
#-------------------------------------------------------------#
find_lambda <- function(x_cal, y_cal, pred, var_fit,
                        alpha, delta, fallback_lambda = NA,
                        lambda_seq = seq(0, 30, by = 0.01)) {
  n <- length(x_cal)
  w <- rep(1 / n, n)          # uniform weights (global average)
  
  # Hoeffding risk threshold
  threshold <- alpha - sqrt(log(1 / delta) / (2 * n))
  threshold <- max(threshold, 0)
  
  # Search for lambda satisfying the risk constraint
  for (lambda in lambda_seq) {
    lower <- pred - lambda * sqrt(var_fit)
    upper <- pred + lambda * sqrt(var_fit)
    risk <- sum(w * ((y_cal < lower) | (y_cal > upper)))
    if (risk <= threshold) return(lambda)
  }
  
  return(fallback_lambda)  # fallback if no lambda found
}

#-------------------------------------------------------------#
# Main Procedure: Load data, split, model, and evaluate TIs
#-------------------------------------------------------------#

# Load necessary libraries
library(MultiKink)
library(zoo)
library(gam)
library(ggplot2)

# Load Triceps dataset
data(triceps)
df <- triceps
x <- df$age
y <- df$lntriceps

# Split data into train, calibration, and test sets
set.seed(2)
n <- length(x)
idx <- sample(n)
n_train <- floor(n * 0.3)
n_cal <- floor(n * 0.4)
train_idx <- idx[1:n_train]
cal_idx <- idx[(n_train + 1):(n_train + n_cal)]
test_idx <- idx[(n_train + n_cal + 1):n]

x_train <- x[train_idx]; y_train <- y[train_idx]
x_cal   <- x[cal_idx];   y_cal <- y[cal_idx]
x_test  <- x[test_idx];  y_test <- y[test_idx]

# Fit smoothing spline and estimate residual variance
fit <- smooth.spline(x_train, y_train)
train_residual <- y_train - predict(fit, x_train)$y
var_fit <- gam(train_residual^2 ~ s(x_train))

# Define tolerance interval parameters
alpha <- 0.10
delta <- 0.05

# Predict mean and variance on calibration set
pred <- predict(fit, x_cal)$y
var_hat <- predict(var_fit, newdata = data.frame(x_train = x_cal))

# Estimate lambda
lambda_hat <- find_lambda(x_cal = data.frame(x_cal), y_cal = data.frame(y_cal),
                          pred = pred, var_fit = var_hat,
                          alpha = alpha, delta = delta)

# Construct tolerance intervals
upper <- pred + sqrt(var_hat) * lambda_hat
lower <- pred - sqrt(var_hat) * lambda_hat

# Empirical risk on calibration set
emp_risk <- mean((y_cal >= upper) | (y_cal <= lower))

# Plot on calibration set
plot_df <- data.frame(x = x_cal, y = y_cal, lower = lower, upper = upper, pred = pred)
ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(color = 'black', alpha = 0.5) +
  geom_line(aes(y = pred), color = 'blue') +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = 'skyblue') +
  ggtitle(paste0('Tolerance Interval with Empirical Risk = ', round(emp_risk, 4))) +
  theme_bw()

# Interpolate TI bounds on test set
upper_interp <- approx(x = x_cal, y = upper, xout = x_test, rule = 2)$y
lower_interp <- approx(x = x_cal, y = lower, xout = x_test, rule = 2)$y
ours_in_TI <- (y_test >= lower_interp) & (y_test <= upper_interp)
ours_risk <- mean(!ours_in_TI)

# Plot on test set
test_plot_df <- data.frame(x = x_test, y = y_test,
                           lower = lower_interp, upper = upper_interp)
ggplot(test_plot_df, aes(x = x, y = y)) +
  geom_point(color = 'black') +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = 'blue') +
  ggtitle(paste0('Tolerance Interval (Ours) with Risk = ', round(ours_risk, 4))) +
  theme_bw()
