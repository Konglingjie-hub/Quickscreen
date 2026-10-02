# 估计均值（Q1, median, Q3）
est_mean_qmq <- function(q1, median, q3) {
  (q1 + median + q3) / 3
}

# 估计SD（Wan 2014，基于IQR并带样本量修正）
est_sd_wan <- function(q1, q3, n) {
  pL <- (0.25 * n - 0.125) / (n + 0.25)
  pU <- (0.75 * n - 0.125) / (n + 0.25)
  zL <- qnorm(pL)
  zU <- qnorm(pU)
  (q3 - q1) / (zU - zL)
}

# 简化近似（不带样本量修正）
est_sd_simple <- function(q1, q3) {
  (q3 - q1) / 1.349
}

# STN 组
q1 <- 33; med <- 40.44; q3 <- 51; n <- 9
mean_stn <- est_mean_qmq(q1, med, q3)
sd_stn_wan <- est_sd_wan(q1, q3, n)
sd_stn_simple <- est_sd_simple(q1, q3)

mean_stn
sd_stn_wan
sd_stn_simple


# 假设你已经知道每组的信息
mean_vec <- c(18.6, 23.2)
sd_vec   <- c(11.6, 13.5)
n_vec    <- c(13,13)
# 合并后的均值
M <- sum(mean_vec * n_vec) / sum(n_vec)
# 合并后的标准差
S <- sqrt(sum(n_vec * ((sd_vec)^2 + (mean_vec - M)^2)) / sum(n_vec))
M
S

