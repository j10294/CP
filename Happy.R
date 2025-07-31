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




library(readr)
library(gam)
library(dplyr)
library(ggplot2)
library(MASS)
library(tidyr)

local_path <- 'photoz_data'

# data 
data_path <- file.path(local_path, 'Happy', 'Happy_A')
data <- read_delim(data_path, delim=' ', comment='#', col_names=FALSE, trim_ws=TRUE)
colnames(data) <- c(
  "id", "mag_r", "u_g", "g_r", "r_i", "i_z",
  "z_spec", "feat1", "feat2", "feat3", "feat4", "feat5"
)

# data visualization
ggplot(data, aes(x=mag_r, y=z_spec))+
  geom_point(alpha=0.2)

data <- data.frame(x=data$mag_r, y=data$z_spec)

# oulier 제거
ggplot(data, aes(x=x))+
  geom_density()
q_x <- quantile(data$x, probs = c(0.01, 0.99))
data <- data %>% filter(x >= q_x[1], x <= q_x[2])

# data sampling
set.seed(123)
data_sample <- data %>% sample_n(5000)

ggplot(data_sample, aes(x = x)) +
  geom_density() +
  ggtitle("Sampled Data (5000개)")

x <- data_sample$x
y <- data_sample$y
n <- length(x)

n_train <- floor(n*0.3)
n_cal <- floor(n*0.4)
n_test <- n - n_train - n_cal

idx <- sample(n)
train_idx <- idx[1:n_train]
cal_idx <- idx[(n_train + 1):(n_train + n_cal)]
test_idx <- idx[(n_train + n_cal + 1):n]

train_x <- x[train_idx]; train_y <- y[train_idx]
cal_x <- x[cal_idx]; cal_y <- y[cal_idx]
test_x <- x[test_idx]; test_y <- y[test_idx]



# (1) Hoeffding
fit <- smooth.spline(train_x, train_y)
train_res <- train_y - predict(fit, train_x)$y
var_fit <- gam(train_res^2 ~ s(train_x))

alpha <- 0.10
delta <- 0.05

mean_pred <- predict(fit, cal_x)$y
var_pred <- predict(var_fit, newdata = data.frame(train_x = cal_x))
var_pred[var_pred < 0] <- 1e-6
plot(var_pred)


lambda_hat <- find_lambda(alpha, delta, cal_y, mean_pred, var_pred, 
                          n = n_cal)

upper <- mean_pred + sqrt(var_pred)*lambda_hat
lower <- mean_pred - sqrt(var_pred)*lambda_hat

cal_df <- data.frame(x=cal_x, y=cal_y, upper=upper, lower=lower, pred = mean_pred)

ggplot(cal_df, aes(x=x, y=y))+
  geom_point(alpha=0.3)+
  geom_line(aes(y=pred), color='blue')+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3)

upper_interp <- approx(x=cal_x, y=upper, xout=test_x, rule=1)$y
lower_interp <- approx(x=cal_x, y=lower, xout=test_x, rule=1)$y

plot_df <- data.frame(x=test_x, y=test_y, upper=upper_interp, lower=lower_interp)
ggplot(plot_df, aes(x=x, y=y))+
  geom_point(alpha=0.2)+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill='skyblue')

ours_content <- mean((test_y >= lower_interp)&(test_y <= upper_interp), na.rm=TRUE)
ours_content



# (2) GY
proper_train_idx <- idx[1:(n_train+n_cal)]
gy_x <- x[proper_train_idx]; gy_y <- y[proper_train_idx]
gy_fit <- smooth.spline(x=gy_x, y=gy_y)
gy_train_residual <- gy_y - predict(gy_fit, gy_x)$y

gy_var_fit <- gam(gy_train_residual^2 ~ s(gy_x))
gy_variance_hat <- predict(gy_var_fit, newdata=data.frame(gy_x = gy_x))
gy_variance_hat[gy_variance_hat<0] <- 1e-6
transform_y <- gy_y / sqrt(gy_variance_hat)
plot(transform_y)

gy_fit_trans <- smooth.spline(x=gy_x, y=transform_y)

B <- bs(gy_x, df = gy_fit_trans$df)
D <- diff(diag(ncol(B)), differences = 2)
S_inv <- ginv(t(B) %*% B + gy_fit_trans$lambda * t(D) %*% D)
S <- B %*% S_inv %*% t(B)
I_n <- diag(length(gy_x))
R <- I_n - S

residuals2 <- transform_y - predict(gy_fit_trans, gy_x)$y
est_var <- (t(residuals2) %*% residuals2) / sum(diag(t(R) %*% R))

nu <- length(gy_x)-1
norm_lx_h_values <- sapply(1:ncol(S), function(j) compute_norm(S[, j]))
k_factors <- find_k_factor(nu, norm_lx_h_values, P=1-alpha, gamma=delta)

gy_mean_pred <- predict(gy_fit_trans, gy_x)$y
gy_upper <- as.vector(gy_mean_pred + c(sqrt(est_var)) %*% k_factors)
gy_lower <- as.vector(gy_mean_pred - c(sqrt(est_var)) %*% k_factors)

gy_upper <- gy_upper * sqrt(gy_variance_hat)
gy_lower <- gy_lower * sqrt(gy_variance_hat)

gy_train_df <- data.frame(x=gy_x, y=gy_y, upper=gy_upper, lower=gy_lower, pred = predict(gy_fit,gy_x)$y)
ggplot(data=gy_train_df, aes(x=x, y=y))+
  geom_point(alpha=0.3)+
  geom_line(aes(y=pred), color='blue')+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

gy_upper_interp <- approx(x=gy_x, y=gy_upper, xout=test_x, rule=1)$y
gy_lower_interp <- approx(x=gy_x, y=gy_lower, xout=test_x, rule=1)$y

gy_test_df <- data.frame(x=test_x, y=test_y, upper=gy_upper_interp, lower=gy_lower_interp)
ggplot(data=gy_test_df, aes(x=x, y=y))+
  geom_point(alpha = 0.2)+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill='skyblue')

gy_content <- mean((test_y >= gy_lower_interp) & (test_y <= gy_upper_interp), na.rm=TRUE); gy_content

# (3) 합쳐서 그리기
ours_df <- data.frame(x = test_x, y = test_y,
                      lower = lower_interp, upper = upper_interp,
                      method = "Ours")

gy_df <- data.frame(x = test_x, y = test_y,
                    lower = gy_lower_interp, upper = gy_upper_interp,
                    method = "GY")
combined_df <- bind_rows(ours_df, gy_df)

p <- ggplot(combined_df, aes(x=x, y=y))+
  geom_point(alpha=0.3)+
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=method), alpha=0.3)+
  scale_fill_manual(values=c("Ours" = 'skyblue', 'GY' = 'tomato'))+
  labs(title = 'Tolerance Intervals', x='x', y='y')+
  annotate('text', x=Inf, y=Inf, hjust=1.1, vjust=1.5,
           label = paste0("Ours content : ", round(ours_content,3),
                          "\nGY content : ", round(gy_content,3)),
           size=4.5, fontface='italic'); p

ggsave(filename='0731_result.png', plot=p, width=10, height=8)
