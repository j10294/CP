library(gam)
library(dplyr)
library(MASS)
library(parallel)
# 데이터 생성 함수 정의
generate_data <- function(model, n) {
  #x <- sort(runif(n, 0, 10)) # x도 랜덤하게 생성
  x <- seq(0, 10, length.out=n)
  base <- 3 * cos(x) - 5 * (x / 15)^2
  
  if (model == 1) {
    y <- base + rnorm(n, mean = 0, sd = sqrt(2))
  }
  if (model == 2) {
    y <- base + rnorm(n, mean = 0, sd = sqrt(2 * x))
  }
  if (model == 3) {
    y <- base + 4 - rgamma(n, shape = 2, scale = 2)
  }
  if (model == 4) {
    y <- base + rt(n, df = 3)
  }
  if (model == 5) {
    sd_x <- 0.3 + 3 * exp(-(x - 5)^2 / (2 * 1.5^2))
    y <- base + rnorm(n, mean = 0, sd = sd_x)
  }
  if (model == 6) {
    y <- base + rnorm(n, mean = 0, sd = 1) + abs(rnorm(n, 0, 1)) - sqrt(2/pi)
  }
  if (model == 7) {
    noise <- ifelse(runif(n) < 0.5, rnorm(n, mean = -3, sd = 1), rnorm(n, mean = 3, sd = 1))
    y <- base + noise
  }
  if (model == 8) {
    noise <- ifelse(x < 5,
                    rnorm(n, mean = 0, sd = 1 + 0.5 * x),
                    rnorm(n, mean = 0, sd = 2 - 0.3 * (x - 5)))
    y <- base + noise
  }
  
  return(data.frame(x = x, y = y))
}

content_function <- function(model,lower,upper,x){
  y = 3*cos(x) - 5*(x/15)^2
  if (model == 1){
    content <- pnorm(upper-y, mean=0, sd=sqrt(2)) - pnorm(lower-y, mean=0, sd=sqrt(2)) #분산 2로 수정. guo and young model 1과 동일 
  }
  if (model == 2){
    content <- pnorm(upper-y, mean=0, sd=sqrt(2 * x)) - pnorm(lower-y, mean=0, sd=sqrt(2 * x))
  }
  if (model == 3){
    content <- pgamma(4 - (lower - y), shape = 2, scale = 2) -
      pgamma(4 - (upper - y), shape = 2, scale = 2)
  }
  if (model == 4){
    content <- pt(upper-y, df=3) - pt(lower-y, df=3)
  }
  if (model == 5){
    sd_x <- 0.3 + 3 * exp(-(x - 5)^2 / (2 * 1.5^2))
    content <- pnorm(upper-y, mean=0, sd=sd_x)- pnorm(lower-y, mean=0, sd=sd_x)
  }
  if (model == 6){
    #empirical content 계산 
    N <- 1e+5
    z1 <- rnorm(N, 0, 1)
    z2 <- abs(rnorm(N, 0, 1)) - sqrt(2 / pi)
    eps <- z1 + z2 #에러항 계산 후 
    content <- sapply(seq_along(x), function(i) {
      y_i <- 3 * cos(x[i]) - 5 * (x[i]/15)^2
      mean((eps >= lower[i] - y_i) & (eps <= upper[i] - y_i))
    })
  }
  
  if (model == 7) {
    N <- 1e+5
    u <- runif(N)
    eps <- ifelse(u < 0.5, 
                  rnorm(N, mean = -3, sd = 1), 
                  rnorm(N, mean = 3, sd = 1))
    content <- sapply(seq_along(x), function(i) {
      y_i <- 3 * cos(x[i]) - 5 * (x[i]/15)^2
      mean((eps >= lower[i] - y_i) & (eps <= upper[i] - y_i))
    })
  } 
  
  if (model == 8){
    sd_x <- ifelse(x < 5, 1 + 0.5 * x, 2 - 0.3 * (x - 5))
    content <- pnorm(upper - y, mean = 0, sd = sd_x) - pnorm(lower - y, mean = 0, sd = sd_x)
  }
  return(content)
}
# --- GY 방법 관련 함수들 ---
compute_norm <- function(vector) sqrt(sum(vector^2))

compute_probability <- function(nu, t, P, k, norm_lx_h) {
  tryCatch({
    q_val <- qchisq(P, df = 1, ncp = pmin(1e4,t^2))
    out <- numeric(length(q_val))
    for (i in seq_along(q_val)) {
      q <- q_val[i]
      if (is.nan(q) || is.na(q)) { out[i] <- 0; next }
      threshold <- (nu * q) / (k^2)
      prob <- pchisq(threshold, df = nu, lower.tail = FALSE)
      out[i] <- ifelse(is.nan(prob) || is.na(prob), 1e-6, prob)
    }
    return(out)
  }, error = function(e) rep(0, length(t)))
} #prop 3.1 의 Pr부분 계산하는 함수. t값을 받아 확률을 계산함 

integrand <- function(t, k, nu, P, norm_lx_h) {
  exp_term <- exp(-t^2 / (2 * norm_lx_h^2))
  prob_term <- compute_probability(nu, t, P, k, norm_lx_h)
  if (log(exp_term) < -20) return(1e-3) #너무 작으면 -inf이 나오므로 이부분 방지
  return(exp_term * prob_term)
}

find_k_factor <- function(nu, norm_lx_h, P=0.90, gamma=0.95){
  c <- norm_lx_h^2
  q1 <- qchisq(P, df=1, ncp=c)
  q2 <- qchisq(1-gamma, df=nu)
  
  if (is.nan(q1) || is.nan(q2) || q2 == 0) return(NA)
  
  k_sq <- (nu * q1) / q2
  return(sqrt(k_sq))
}


sample_sizes <- c(50, 100, 200)
models <- 1:8
n_sim <- 100
results <- expand.grid(model = models, m = sample_sizes)
results$mean_coverage <- NA
results$mean_content <- NA

coverage_list <- list() 

for (idx in seq_len(nrow(results))) {
  model <- results$model[idx]
  m <- results$m[idx]
  coverages <- numeric(n_sim)
  contents <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    set.seed(i)
    data <- generate_data_2(model, m)
    x <- data$x
    y <- data$y
    fit <- smooth.spline(x = x, y = y, cv = FALSE)
    residuals <- y - predict(fit, x)$y
    log_variance_fit <- smooth.spline(x, log(residuals^2 + 1e-6), cv = FALSE)
    variance_hat <- exp(predict(log_variance_fit, x)$y)
    
    transform_y <- y / sqrt(variance_hat)
    fit_trans <- smooth.spline(x, transform_y, cv = FALSE)
    
    B <- bs(x, df = fit_trans$df)
    D <- diff(diag(ncol(B)), differences = 2)
    S_inv <- ginv(t(B) %*% B + fit_trans$lambda * t(D) %*% D)
    S <- B %*% S_inv %*% t(B)
    I_n <- diag(length(x))
    R <- I_n - S
    
    residuals2 <- transform_y - predict(fit_trans, x)$y
    est_var <- (t(residuals2) %*% residuals2) / sum(diag(t(R) %*% R))
    num <- sum(diag(t(R) %*% R))^2
    den <- sum(diag((t(R) %*% R)^2))
    nu <- num / den
    
    compute_norm <- function(vector) sqrt(sum(vector^2))
    norm_lx_h_values <- sapply(1:ncol(S), function(j) compute_norm(S[, j]))
    k_factors <- unlist(mclapply(norm_lx_h_values, function(nlh) {
      find_k_factor(nu, nlh)
    }, mc.cores = detectCores(logical = FALSE)))
    
    pred_seq <- predict(fit_trans, x)$y
    gy_upper <- ( pred_seq + sqrt(as.vector(est_var)) * k_factors ) * sqrt(variance_hat)
    gy_lower <- ( pred_seq - sqrt(as.vector(est_var)) * k_factors ) * sqrt(variance_hat)
    
    content <- content_function_2(model=model, lower=gy_lower, upper=gy_upper, x=x)
    pointwise_covered <- as.numeric(content >= 0.90)
    interval_widths <- gy_upper - gy_lower
    
    coverages[i] <- mean(pointwise_covered)
    contents[i] <- mean(content)
    
    # ✅ 각 simulation의 pointwise coverage 저장
    coverage_list[[length(coverage_list) + 1]] <- data.frame(
      x = x,
      coverage = pointwise_covered,
      interval_width = interval_widths,
      model = model,
      m = m,
      sim = i
    )
  }
  
  results$mean_coverage[idx] <- mean(coverages)
  results$mean_content[idx] <- mean(contents)
  print(results[idx, ])
}

library(ggplot2)

#평균 커버리지 시각화 (0.90 이상인 비율 의미!)
p <- ggplot(results, aes(x = m, y = mean_coverage, color = factor(model))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  labs(
    title = "Coverage by Model and Sample Size (m)",
    x = "Sample Size (m)",
    y = "Mean Coverage",
    color = "Model"
  ) + ylim(0,1)+
  theme_minimal(base_size = 14) +
  theme(legend.position = "right"); p


ggsave(
  filename = "0729_GY_marginal2.png",
  plot = p,
  width = 8,       # 너비 (inch 단위)
  height = 5,      # 높이
  dpi = 300        # 해상도 (논문용이면 300 이상 권장)
)


coverage_df_all <- do.call(rbind, coverage_list)

library(dplyr)

coverage_df <- coverage_df_all |>
  group_by(model, m, x) |>
  summarise(mean_coverage = mean(coverage), .groups = "drop")

p2 <- ggplot(coverage_df, aes(x = x, y = mean_coverage)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
  facet_grid(model ~ m, labeller = label_both) +
  labs(
    title = "Pointwise Coverage Proportion for GY (content ≥ 0.90)",
    x = "x",
    y = "Coverage Proportion (100 Simulations)"
  ) + ylim(0,1)+
  theme_minimal(base_size = 13); p2
ggsave(filename = '0729_GY_pointwise2.png', 
       plot=p2, width=14, height=10, dpi=300)

# Pointwise 평균 NMPIW 계산 및 시각화
width_df <- do.call(rbind, coverage_list) |>
  group_by(model, m, x) |>
  summarise(mean_width = mean(interval_width, na.rm=TRUE), .groups = "drop")

p3 <- ggplot(width_df, aes(x = x, y = mean_width)) +
  geom_line(color = "darkgreen") +
  facet_grid(model ~ m, labeller = label_both) +
  labs(
    title = "Pointwise Mean Interval Width (NMPIW)",
    x = "x",
    y = "Mean Interval Width"
  ) +
  theme_minimal(base_size = 13)
p3

ggsave(filename = '0729_GY_pointwise_NMPIW.png', 
       plot = p3, width = 14, height = 10, dpi = 300)
