library(dplyr)
library(readr)
library(metafor)
library(forestplot)

# 读入你的文件
dat <- read_csv("<ORIGINAL_LOCAL_PATH_REDACTED>")
colnames(dat)
dat_overall <- read_csv("<ORIGINAL_LOCAL_PATH_REDACTED>")
make_effect <- function(row){
  N  <- as.numeric(row["N"])
  MeanCo <- as.numeric(row["Meanco"])
  MeanCtrl <- as.numeric(row["MeanCtrl"])
  SDco <- as.numeric(row["SDco"])
  SDctrl <- as.numeric(row["SDctrl"])
  md_within <- as.numeric(row["Within_mean_co_STN"])
  se_within <- as.numeric(row["Within_SE_co_STN"])
  md_ind <- as.numeric(row["Mean_Difference_CO_STN"])
  se_ind <- as.numeric(row["SE_CO_STN"])
  
  # 优先使用配对差
  if(!is.na(md_within) & !is.na(se_within)){
    yi <- md_within
    sei <- se_within
  } else if(!is.na(md_ind) & !is.na(se_ind)){
    yi <- md_ind
    sei <- se_ind
  } else if(!is.na(MeanCo) & !is.na(MeanCtrl) & !is.na(SDco) & !is.na(SDctrl) & !is.na(N)){
    yi <- MeanCo - MeanCtrl
    sei <- sqrt(SDco^2/N + SDctrl^2/N)
  } else {
    yi <- NA; sei <- NA
  }
  return(c(yi, sei))
}
dat_overall_2 <- dat_overall %>%
  rowwise() %>%
  mutate(tmp = list(make_effect(cur_data()))) %>%
  mutate(yi = tmp[1], sei = tmp[2]) %>%
  ungroup()

dat2 <- dat %>%
  rowwise() %>%
  mutate(tmp = list(make_effect(cur_data()))) %>%
  mutate(yi = tmp[1], sei = tmp[2]) %>%
  ungroup()

updrs <- dat_overall_2 %>%
  filter(Assessment == "UPDRS_3",
         !is.na(yi), !is.na(sei))

# UPDRS-III medication ON
updrs_on <- dat2 %>%
  filter(Assessment == "UPDRS_3",
         Medication == "ON",
         !is.na(yi), !is.na(sei))

# UPDRS-III medication OFF
updrs_off <- dat2 %>%
  filter(Assessment == "UPDRS_3",
         Medication == "OFF",
         !is.na(yi), !is.na(sei))

library(metafor)
res_on  <- rma.uni(yi, sei, data = updrs_on, method = "REML")
res_off <- rma.uni(yi, sei, data = updrs_off, method = "REML")
res_updrs <- rma.uni(yi = updrs$yi,
                     sei = updrs$sei,
                     method = "REML",
                     test = "knha")

# ============================================================# ============================================================
# ★★★ 新增 1: Medication-state interaction test
# ============================================================

# --- 方法A: Z-test for subgroup interaction ---
# (两亚组估计值之差 / 差值的SE)
diff_b      <- as.numeric(res_on$b) - as.numeric(res_off$b)
se_diff     <- sqrt(res_on$se^2 + res_off$se^2)
z_interact  <- diff_b / se_diff
p_interact_z <- 2 * pnorm(-abs(z_interact))

cat("====== Medication-State Interaction (Z-test) ======\n")
cat(sprintf("  Difference (ON - OFF) = %.2f\n", diff_b))
cat(sprintf("  SE of difference      = %.2f\n", se_diff))
cat(sprintf("  z = %.3f, p = %.4f\n\n", z_interact, p_interact_z))

# --- 方法B: Random-effects meta-regression with Knapp-Hartung adjustment ---
updrs_combined <- bind_rows(
  updrs_on %>% mutate(Medication = "ON"),
  updrs_off %>% mutate(Medication = "OFF")
)

res_mod <- rma.uni(yi  = updrs_combined$yi,
                   sei = updrs_combined$sei,
                   mods = ~ factor(Medication),
                   data = updrs_combined,
                   method = "REML",
                   test = "knha")   # Knapp-Hartung adjustment

cat("====== Medication-State Interaction (Meta-regression, Knapp-Hartung) ======\n")
cat(sprintf("  Moderator coefficient (OFF vs ON): %.2f (95%% CI %.2f to %.2f)\n",
            res_mod$b[2], res_mod$ci.lb[2], res_mod$ci.ub[2]))
cat(sprintf("  t = %.3f, p = %.4f\n", res_mod$zval[2], res_mod$pval[2]))
cat(sprintf("  Omnibus test of moderator: QM = %.3f, p = %.4f\n\n",
            res_mod$QM, res_mod$QMp))

# ============================================================
# ★★★ 新增 2: Overall pooled MD 的 P 值
# ============================================================
p_overall <- res_updrs$pval   # Knapp-Hartung adjusted p-value

cat("====== Overall Pooled Effect ======\n")
cat(sprintf("  Pooled MD = %.2f (95%% CI %.2f to %.2f)\n",
            as.numeric(res_updrs$b), res_updrs$ci.lb, res_updrs$ci.ub))
cat(sprintf("  t = %.3f, p = %.4f\n\n", res_updrs$zval, p_overall))

# 1. 提取 rma 结果
# ============================================================
b_on     <- as.numeric(res_on$b)
ci_lb_on <- res_on$ci.lb
ci_ub_on <- res_on$ci.ub

b_off     <- as.numeric(res_off$b)
ci_lb_off <- res_off$ci.lb
ci_ub_off <- res_off$ci.ub

b_all <- as.numeric(res_updrs$b)
ci_lb_all<- res_updrs$ci.lb
ci_ub_all <- res_updrs$ci.ub




# ============================================================
# 2. 提取异质性统计量
# ============================================================
# --- ON ---
I2_on   <- round(res_on$I2, 1)
tau2_on <- round(res_on$tau2, 2)
Q_on    <- round(res_on$QE, 2)
Qp_on   <- res_on$QEp

# --- OFF ---
I2_off   <- round(res_off$I2, 1)
tau2_off <- round(res_off$tau2, 2)
Q_off    <- round(res_off$QE, 2)
Qp_off   <- res_off$QEp

# --- Overall ---
I2_all   <- round(res_updrs$I2, 1)
tau2_all <- round(res_updrs$tau2, 2)
Q_all    <- round(res_updrs$QE, 2)
Qp_all   <- res_updrs$QEp

# ★ 预测区间计算
# ON
pred_on  <- predict(res_on)
pred_off <- predict(res_off)
pred_all <- predict(res_updrs)

pi_lb_on  <- pred_on$pi.lb;  pi_ub_on  <- pred_on$pi.ub
pi_lb_off <- pred_off$pi.lb; pi_ub_off <- pred_off$pi.ub
pi_lb_all <- pred_all$pi.lb; pi_ub_all <- pred_all$pi.ub


# ============================================================
# 3. 格式化 p 值的辅助函数
# ============================================================
fmt_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01)  return(sprintf("p = %.3f", p))
  return(sprintf("p = %.2f", p))
}

# ============================================================
# 4. 构建异质性文本
# ============================================================
het_on_text <- sprintf(
  "    (I\u00B2 = %.1f%%, \u03C4\u00B2 = %.2f, Q = %.2f, %s)",
  I2_on, tau2_on, Q_on, fmt_p(Qp_on)
)

het_off_text <- sprintf(
  "    (I\u00B2 = %.1f%%, \u03C4\u00B2 = %.2f, Q = %.2f, %s)",
  I2_off, tau2_off, Q_off, fmt_p(Qp_off)
)

het_all_text <- sprintf(
  "    (I\u00B2 = %.1f%%, \u03C4\u00B2 = %.2f, Q = %.2f, %s)",
  I2_all, tau2_all, Q_all, fmt_p(Qp_all)
)

# ============================================================
# 5. 构建数据行（含异质性信息行）
# ============================================================

# --- ON 亚组标题 ---
on_label <- tibble(
  Research = "Medication ON", Medication = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = TRUE
)

# --- ON 个体研究 ---
on_data <- updrs_on %>%
  transmute(
    Research, Medication,
    mean  = yi,
    lower = yi - 1.96 * sei,
    upper = yi + 1.96 * sei,
    is_summary = FALSE
  )

# --- ON subtotal（菱形） ---
on_subtotal <- tibble(
  Research = "  Subtotal (Medication ON)", Medication = "",
  mean = b_on, lower = ci_lb_on, upper = ci_ub_on,
  is_summary = TRUE
)

# ★ ON 预测区间行
on_pi_row <- tibble(
  Research = sprintf("  Prediction Interval: "),
  Medication = "",
  mean = b_on, lower = pi_lb_on, upper = pi_ub_on,
  is_summary = FALSE
)

# --- ON 异质性信息行 ---
on_het_row <- tibble(
  Research = het_on_text, Medication = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = FALSE
)

# --- 空行 ---
blank_row <- tibble(
  Research = "", Medication = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = FALSE
)

# --- OFF 亚组标题 ---
off_label <- tibble(
  Research = "Medication OFF", Medication = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = TRUE
)

# --- OFF 个体研究 ---
off_data <- updrs_off %>%
  transmute(
    Research, Medication,
    mean  = yi,
    lower = yi - 1.96 * sei,
    upper = yi + 1.96 * sei,
    is_summary = FALSE
  )

# --- OFF subtotal（菱形） ---
off_subtotal <- tibble(
  Research = "  Subtotal (Medication OFF)", Medication = "",
  mean = b_off, lower = ci_lb_off, upper = ci_ub_off,
  is_summary = TRUE
)

# ★ OFF 预测区间行
off_pi_row <- tibble(
  Research = sprintf("  Prediction Interval: "),
  Medication = "",
  mean = b_off, lower = pi_lb_off, upper = pi_ub_off,
  is_summary = FALSE
)

# --- OFF 异质性信息行 ---
off_het_row <- tibble(
  Research = het_off_text, Medication = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = FALSE
)

# --- Overall（菱形） ---
overall_row <- tibble(
  Research = "Overall Effect", Medication = "",
  mean = b_all, lower = ci_lb_all, upper = ci_ub_all,
  is_summary = TRUE
)

# ★ Overall 预测区间行
overall_pi_row <- tibble(
  Research = sprintf("  Prediction Interval: "),
  Medication = "",
  mean = b_all, lower = pi_lb_all, upper = pi_ub_all,
  is_summary = FALSE
)

# --- Overall 异质性信息行 ---
overall_het_row <- tibble(
  Research = het_all_text, Medication = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = FALSE
)

# ============================================================
# 6. 合并所有行
# ============================================================
forest_all <- bind_rows(
  on_label,
  on_data,
  on_subtotal,
  on_pi_row,          # ★ 预测区间
  on_het_row,
  blank_row,
  off_label,
  off_data,
  off_subtotal,
  off_pi_row,         # ★ 预测区间
  off_het_row,
  blank_row,
  overall_row,
  overall_pi_row,     # ★ 预测区间
  overall_het_row
)

# ============================================================
# 7. 构建 tabletext
# ============================================================
ci_text <- ifelse(
  is.na(forest_all$mean), "",
  sprintf("%.2f [%.2f, %.2f]", forest_all$mean, forest_all$lower, forest_all$upper)
)

tabletext <- cbind(
  c("Study",        forest_all$Research),
  c("Medication",   forest_all$Medication),
  c("MD [95% CI]",  ci_text)
)

# ============================================================
# 8. 构建绘图输入
# ============================================================
fn_mean       <- c(NA, forest_all$mean)
fn_lower      <- c(NA, forest_all$lower)
fn_upper      <- c(NA, forest_all$upper)
fn_is_summary <- c(FALSE, forest_all$is_summary)

# ★ 标记预测区间行的位置（用于 shapes_gp）
# 在 fn_mean 中，预测区间行对应的索引
pi_rows <- which(grepl("Prediction Interval", c("header", forest_all$Research)))

# ============================================================
# 9. 绘图
# ============================================================
forestplot(
  labeltext  = tabletext,
  mean       = fn_mean,
  lower      = fn_lower,
  upper      = fn_upper,
  is.summary = fn_is_summary,
  zero       = 0,
  
  col = fpColors(
    box     = "black",
    line    = "black",
    summary = "black",
    zero    = "gray50"
  ),
  
  boxsize  = 0.25,
  lwd.ci   = 1.5,
  lwd.zero = 1,
  
  txt_gp = fpTxtGp(
    label   = gpar(fontfamily = "Arial", cex = 0.8),
    ticks   = gpar(fontfamily = "Arial", cex = 0.7),
    xlab    = gpar(fontfamily = "Arial", cex = 0.85),
    title   = gpar(fontfamily = "Arial", cex = 1.0, fontface = "bold"),
    summary = gpar(fontfamily = "Arial", cex = 0.8, fontface = "bold")
  ),
  
  # ★ 预测区间行用红色虚线，其余保持默认
  shapes_gp = fpShapesGp(
    default = gpar(col = "black", fill = "black"),
    lines = lapply(seq_along(fn_mean), function(i) {
      if (i %in% pi_rows) {
        gpar(col = "red", lty = 2, lwd = 2)
      } else {
        gpar(col = "black", lwd = 1.5)
      }
    }),
    box = lapply(seq_along(fn_mean), function(i) {
      if (i %in% pi_rows) {
        gpar(col = "red", fill = NA)   # 预测区间不画实心box
      } else {
        gpar(col = "black", fill = "black")
      }
    })
  ),
  
  xlab   = "Mean Difference (Favours DBS \u2190    \u2192 Favours Baseline)",
  clip   = c(-35, 15),
  xticks = seq(-40, 10, by = 10),
  graph.pos = 3,
  title  = "FOG_AC: Forest Plot with Subgroup Analysis"
)


