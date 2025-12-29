find_lambda <- function(alpha, delta, y, pred, variance, n, lower=0, upper=30, step=1e-3){
  threshold <- max(0, alpha-sqrt(log(1/delta)/(2*n)))
  lambda_seq <- seq(lower, upper, by=step)
  for (lambda in lambda_seq){
    risk <- mean(y <= pred - sqrt(variance)*lambda | y >= pred + sqrt(variance)*lambda, na.rm=TRUE)
    if (risk <= threshold) return(lambda)
  } 
  return(NA)
}

find_k_factor <- function(nu, norm_lx_h, P, gamma){
  sqrt(nu*qchisq(p=P, df=1, ncp=norm_lx_h^2) / qchisq(p=1-gamma, df=nu))
}

compute_norm <- function(vector){
  return(sqrt(sum(vector^2)))
}

compute_norm_fast <- function(B, A, BtB){
  # returns vector: ||S_{:,j}|| for j = 1,...,n
  # S_{:,j} = B %*% A %*% t(B[j, ])
  
  M <- A %*% BtB %*% A          # small matrix (df x df)
  
  # 각 j에 대해: B[j, ] %*% M %*% t(B[j, ])
  sqrt(rowSums((B %*% M) * B))
}

tf  <- function(y) log1p(y)
itf <- function(z) expm1(z)



library(readr)
library(gam)
library(dplyr)
library(ggplot2)
library(MASS)
library(tidyr)

local_path <- 'photoz_data'

# data 
data_path <- file.path(local_path, 'Happy', 'happy_A')
data <- read_delim(data_path, delim=' ', comment='#', col_names=FALSE, trim_ws=TRUE)
colnames(data) <- c(
  "id", "mag_r", "u_g", "g_r", "r_i", "i_z",
  "z_spec", "feat1", "feat2", "feat3", "feat4", "feat5"
)

# data visualization
ggplot(data, aes(x=mag_r, y=z_spec))+
  geom_point(alpha=0.2)

data <- data.frame(x=data$mag_r, y=data$z_spec)

# data sampling
set.seed(123)
data_sample <- data %>% sample_n(5000)

x <- data_sample$x
y <- data_sample$y
n <- length(x)

n_train <- floor(n*0.5)
n_cal <- floor(n*0.5)

idx <- sample(n)
train_idx <- idx[1:n_train]
cal_idx <- idx[(n_train + 1):(n_train + n_cal)]
train_x <- x[train_idx]; train_y <- y[train_idx]
cal_x <- x[cal_idx]; cal_y <- y[cal_idx]


#test data
test_data <- read_delim(file='photoz_data/Happy/happy_B', delim=' ', comment='#', col_names=FALSE)
colnames(test_data) <- c(
  "id", "mag_r", "u_g", "g_r", "r_i", "i_z",
  "z_spec", "feat1", "feat2", "feat3", "feat4", "feat5"
)
test_x <- test_data$mag_r
test_y <- test_data$z_spec

# (1) Hoeffding (log1p scale)
alpha <- 0.10
delta <- 0.05

# transform y
train_z <- tf(train_y)
cal_z   <- tf(cal_y)

# fit mean on transformed scale
fit <- smooth.spline(train_x, train_z)

# residuals on transformed scale
train_res <- train_z - predict(fit, train_x)$y

# variance model on transformed scale
var_fit <- gam(train_res^2 ~ s(train_x))

# predictions on cal set (transformed scale)
mean_pred_z <- predict(fit, cal_x)$y
var_pred    <- predict(var_fit, newdata = data.frame(train_x = cal_x))
var_pred[var_pred < 0] <- 1e-6

lambda_hat <- find_lambda(
  alpha = alpha, delta = delta,
  y = cal_z, pred = mean_pred_z, variance = var_pred,
  n = n_cal
)

upper_z <- mean_pred_z + sqrt(var_pred) * lambda_hat
lower_z <- mean_pred_z - sqrt(var_pred) * lambda_hat

# back-transform to original y scale
upper <- itf(upper_z)
lower <- itf(lower_z)

# optional: clip lower at 0 (recommended if you want nonnegative lower bound)
lower <- pmax(lower, 0)

cal_df <- data.frame(x = cal_x, y = cal_y, upper = upper, lower = lower,
                     pred = itf(mean_pred_z))

ggplot(cal_df, aes(x = x, y = y)) +
  geom_point(alpha = 0.3) +
  geom_line(aes(y = pred), color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

# interpolate to test_x
upper_interp <- approx(x = cal_x, y = upper, xout = test_x, rule = 1)$y
lower_interp <- approx(x = cal_x, y = lower, xout = test_x, rule = 1)$y

plot_df <- data.frame(x = test_x, y = test_y, upper = upper_interp, lower = lower_interp)

ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(alpha = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "skyblue")

ours_content <- mean((test_y >= lower_interp) & (test_y <= upper_interp), na.rm = TRUE)
ours_content


# (2) GY (log1p scale)
proper_train_idx <- idx[1:(n_train + n_cal)]
gy_x <- x[proper_train_idx]
gy_y <- y[proper_train_idx]

# transform response
gy_z <- tf(gy_y)

# mean fit on transformed scale (for residuals / variance estimation)
gy_fit <- smooth.spline(x = gy_x, y = gy_z)
gy_train_residual <- gy_z - predict(gy_fit, gy_x)$y

# variance fit on transformed scale
gy_var_fit <- gam(gy_train_residual^2 ~ s(gy_x))
gy_variance_hat <- predict(gy_var_fit, newdata = data.frame(gy_x = gy_x))
gy_variance_hat[gy_variance_hat < 0] <- 1e-6

# standardized response on transformed scale
transform_y <- gy_z / sqrt(gy_variance_hat)

# spline on standardized response
gy_fit_trans <- smooth.spline(x = gy_x, y = transform_y)

# basis computations
B <- bs(gy_x, df = gy_fit_trans$df)
D <- diff(diag(ncol(B)), differences = 2)
BtB <- crossprod(B)
A <- ginv(BtB + gy_fit_trans$lambda * crossprod(D))

# residual variance on standardized scale
residuals2 <- transform_y - predict(gy_fit_trans, gy_x)$y
est_var <- as.numeric(crossprod(residuals2) / (nrow(B) - 1))

nu <- length(gy_x) - 1

norm_lx_h_values <- compute_norm_fast(B, A, BtB)
k_factors <- find_k_factor(nu, norm_lx_h_values, P = 1 - alpha, gamma = delta)

gy_mean_pred <- predict(gy_fit_trans, gy_x)$y

# IMPORTANT: elementwise multiplication (your original line used %*% and will inflate)
gy_upper_t <- gy_mean_pred + sqrt(est_var) * k_factors
gy_lower_t <- gy_mean_pred - sqrt(est_var) * k_factors

# back to transformed-y scale (Z = log1p(Y))
gy_upper_z <- gy_upper_t * sqrt(gy_variance_hat)
gy_lower_z <- gy_lower_t * sqrt(gy_variance_hat)

# back to original y scale
gy_upper <- itf(gy_upper_z)
gy_lower <- itf(gy_lower_z)

# optional lower clipping
gy_lower <- pmax(gy_lower, 0)

gy_train_df <- data.frame(
  x = gy_x, y = gy_y,
  upper = gy_upper, lower = gy_lower,
  pred = itf(predict(gy_fit, gy_x)$y)  # mean curve on original scale
)

ggplot(data = gy_train_df, aes(x = x, y = y)) +
  geom_point(alpha = 0.3) +
  geom_line(aes(y = pred), color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

gy_upper_interp <- approx(x = gy_x, y = gy_upper, xout = test_x, rule = 1)$y
gy_lower_interp <- approx(x = gy_x, y = gy_lower, xout = test_x, rule = 1)$y

gy_test_df <- data.frame(x = test_x, y = test_y, upper = gy_upper_interp, lower = gy_lower_interp)

ggplot(data = gy_test_df, aes(x = x, y = y)) +
  geom_point(alpha = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "skyblue")

gy_content <- mean((test_y >= gy_lower_interp) & (test_y <= gy_upper_interp), na.rm = TRUE)
gy_content

# (3) 합쳐서 그리기
ours_df <- data.frame(x = test_x, y = test_y,
                      lower = lower_interp, upper = upper_interp,
                      method = "Ours")

gy_df <- data.frame(x = test_x, y = test_y,
                    lower = gy_lower_interp, upper = gy_upper_interp,
                    method = "GY")
combined_df <- bind_rows(ours_df, gy_df)

paste0("Ours content : ", round(ours_content,3),
       "\nGY content : ", round(gy_content,3))

p <- ggplot(combined_df, aes(x = x, y = y)) +
  geom_point(alpha = 0.05, size = 0.6, color = "grey30") +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = method),
    alpha = 0.35
  ) +
  scale_fill_manual(values = c("Ours" = "tomato", "GY" = "skyblue")) +
  labs(title = "Tolerance Intervals", x = "x", y = "y") +
  theme_bw()

p
