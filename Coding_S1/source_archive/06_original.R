# ---------------------------------------------
# Crossover / paired differences calculator
# 输入：每位受试者在Co-stimulation与STN-only下的差值（Co − STN）
# 输出：Mean difference、SD、SE、95%CI、t值、p值
# ---------------------------------------------

# Step 1. 输入你的数据（这里是示例，你替换成自己6个差值）
data <- read.csv("<ORIGINAL_LOCAL_PATH_REDACTED>")
d <-  data$Co_Base # 例如每位患者的差值
colnames(data)
# Step 2. 计算基本统计量
n  <- length(d)
md <- mean(d)         # 平均差值
sd <- sd(d)           # 差值的标准差
se <- sd / sqrt(n)    # 标准误
df <- n - 1           

# Step 3. 计算置信区间（95%）
t_crit <- qt(0.975, df)
ci_low <- md - t_crit * se
ci_high<- md + t_crit * se

# Step 4. 计算t值及p值（配对t检验）
t_val <- md / se
p_val <- 2 * (1 - pt(abs(t_val), df))

# Step 5. 打印结果
cat("Within-subject mean difference:\n")
cat("n =", n, "\nMean =", round(md,2), "\nSD =", round(sd,2),
    "\nSE =", round(se,2), "\nt =", round(t_val,2),
    "\np =", round(p_val,4), "\n95% CI = [",
    round(ci_low,2), ", ", round(ci_high,2), "]\n")


md_se_indep <- function(mean_S, sd_S, n_S, mean_C, sd_C, n_C) {
  md <- mean_C - mean_S
  se <- sqrt(sd_S^2 / n_S + sd_C^2 / n_C)
  list(md = md, se = se)
}
res1 <- md_se_indep(mean_S = 41.48, sd_S = 13.34, n_S = 9,
                    mean_C = 43.37, sd_C = 6.86, n_C = 9)
print(res1)
md_se_paired_corr <- function(mean_S, sd_S, mean_C, sd_C, n, rho = 0) {
  # ----
  # mean_S: 刺激条件（S）的均值
  # sd_S:   S条件的SD
  # mean_C: 对照条件（C1 或 C2）的均值
  # sd_C:   对照条件SD
  # n:      配对样本数
  # rho:    两条件的相关系数（若无实测数据，可取0.5或文献经验值）
  # ----
  
  # 配对均值差（within-subject difference）
  MD_paired <- mean_C - mean_S  # 或反过来取符号，看你的方向定义
  
  # 配对均值差的方差（近似）
  Var_MD_paired <- (sd_S^2 + sd_C^2 - 2 * rho * sd_S * sd_C) / n
  
  SE_MD_paired <- sqrt(Var_MD_paired)
  
  CI95 <- c(MD_paired - 1.96 * SE_MD_paired,
            MD_paired + 1.96 * SE_MD_paired)
  
  list(
    MD_paired   = MD_paired,
    Var_paired  = Var_MD_paired,
    SE_paired   = SE_MD_paired,
    CI95        = CI95
  )
}
res1 <- md_se_paired_corr(mean_S = 16.17, sd_S = 3.83,
                          mean_C = 14.5, sd_C = 4.89,
                          n = 12, rho = 0.5)

print(res1)
combine_md_two <- function(md1, se1, md2, se2) {
  v1 <- se1^2; v2 <- se2^2
  md_comb <- (md1 / v1 + md2 / v2) / (1 / v1 + 1 / v2)
  var_comb <- 1 / (1 / v1 + 1 / v2)
  se_comb <- sqrt(var_comb)
  ci95 <- c(md_comb - 1.96 * se_comb, md_comb + 1.96 * se_comb)
  list(MD_combined = md_comb, Var_combined = var_comb, SE_combined = se_comb, CI95 = ci95)
}

n <- 13         # 用你的实际配对样本量
df <- n - 1

res1 <- md_se_paired(md = -0.99, n = n, p = 0.642, df = df)  # S vs C1
res2 <- md_se_paired(md = -2.73, n = n, p = 0.191, df = df)  # S vs C2
combine_md_two(-4, 3.448813806,-13,2.555634298)

2.01	3.346537987
0.17	2.553429589

-4	3.448813806
-13	2.555634298


#################################################################
# --------------------------------------------
# 输入：每位受试者在不同频率下的差值（co-stimulation − STN）
# 每行一个受试者，每列一个频率差值
# 示例：10例 × 3频率
# --------------------------------------------
diffs <- matrix(c(
  # 每行替换为你的实际差值（例：C1, C2, C3）
  -2, -1, -3,
  -5, -3, -4,
  -4, -2, -3,
  -1,  0, -2,
  -3, -4, -5,
  -2, -1, -2,
  0, -1, -2,
  -3, -4, -3,
  -2, -1, -1,
  -1,  0, -2
), nrow=10, byrow=TRUE)

# Step 1. 对每个受试者取平均差（合并三个频率）
d_avg <- rowMeans(diffs, na.rm=TRUE)
data <- read.csv("<ORIGINAL_LOCAL_PATH_REDACTED>")
d_avg <- data$Mean1
# Step 2. 计算统计指标
n  <- length(d_avg)
MD <- mean(d_avg)
SD <- sd(d_avg)
SE <- SD / sqrt(n)
df <- n - 1

# Step 3. 95% CI & p-value
t_val <- MD / SE
p_val <- 2 * (1 - pt(abs(t_val), df))
CI <- MD + qt(c(0.025, 0.975), df) * SE

# Step 4. 输出结果
cat("Combined within-subject estimate (across 3 frequencies):\n")
cat("n =", n, "\nMean Difference =", round(MD,2),
    "\nSD =", round(SD,2), "\nSE =", round(SE,3),
    "\nt =", round(t_val,2), "\np =", round(p_val,4),
    "\n95% CI = [", round(CI[1],2), ",", round(CI[2],2), "]\n")


#####################################
#############################
n <- 10
mean_STN <- 63.2
sd_STN   <- 12.52

mean_co  <- c(39.9, 39.4, 42.5)     # 30, 70, 110 Hz
sd_co    <- c( 9.86, 13.55, 15.86)

combine_crossover <- function(mean_STN, sd_STN, mean_co, sd_co, n, rho=0.5) {
  MDs <- mean_co - mean_STN
  Var <- (sd_STN^2 + sd_co^2 - 2*rho*sd_STN*sd_co) / n
  SEs <- sqrt(Var)
  w   <- 1 / (SEs^2)
  MD_star <- sum(w * MDs) / sum(w)
  SE_star <- sqrt(1 / sum(w))
  CI <- MD_star + c(-1,1) * 1.96 * SE_star
  list(MD_each=MDs, SE_each=SEs, w_each=w,
       MD=MD_star, SE=SE_star, CI=CI)
}
res_r05 <- combine_crossover(mean_STN, sd_STN, mean_co, sd_co, n, rho=0.5)
res_r05 

##############合并#####
MD_on  <- 0.17
SE_on  <- 2.32
MD_off <- -2
SE_off <- 2.79

# 计算权重
w_on  <- 1 / SE_on^2
w_off <- 1 / SE_off^2

# 合并效应
MD_comb <- (MD_on * w_on + MD_off * w_off) / (w_on + w_off)

# 合并SE
SE_comb <- sqrt(1 / (w_on + w_off))

# 95% CI
CI_low  <- MD_comb - 1.96 * SE_comb
CI_high <- MD_comb + 1.96 * SE_comb

cat("Combined MD =", round(MD_comb,3),
    "\nCombined SE =", round(SE_comb,3),
    "\n95% CI =", round(CI_low,3), "to", round(CI_high,3))


