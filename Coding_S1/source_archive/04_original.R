# ============================================================
# 0. 数据准备 & Meta-Analysis
# ============================================================
library(metafor)
library(forestplot)
library(grid)
library(dplyr)
library(tibble)
library(readr)
colnames(dat2)
dat <- read_csv("<ORIGINAL_LOCAL_PATH_REDACTED>")
colnames(dat)
make_effect <- function(row){
  N  <- as.numeric(row["N"])
  MeanCo <- as.numeric(row["Meanco"])
  MeanCtrl <- as.numeric(row["MeanCtrl"])
  SDco <- as.numeric(row["SDco"])
  SDctrl <- as.numeric(row["SDctrl"])
  md_within <- as.numeric(row["MD_co_STN_final"])
  se_within <- as.numeric(row["SE_co_STN_final"])
  
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

dat2 <- dat %>%
  rowwise() %>%
  mutate(tmp = list(make_effect(cur_data()))) %>%
  mutate(yi = tmp[1], sei = tmp[2]) %>%
  ungroup()
updrs <- dat2 %>%
  filter(Assessment == "UPDRS_3") %>% 
  mutate(Group = ifelse(Frequency >= 119, 'High', 'Low')) %>% 
  select(Research, Frequency, Assessment, Medication, N, 
         MD_co_base_final, SE_co_base_final, yi, sei, tmp,Group) %>% slice(-8)
updrs_high <- updrs %>% filter(Group == "High")
updrs_low  <- updrs %>% filter(Group == "Low")

res_high    <- rma.uni(yi = updrs_high$yi, sei = updrs_high$sei, method = "REML", test = "knha")
res_low     <- rma.uni(yi = updrs_low$yi,  sei = updrs_low$sei,  method = "REML", test = "knha")
res_overall <- rma.uni(yi = updrs$yi, sei = updrs$sei, 
                       method = "DL", test = "knha")
updrs
head(updrs)
# ============================================================# ============================================================
# ★★★ Frequency Subgroup Interaction Test
# ============================================================

# --- Z-test ---
diff_b       <- as.numeric(res_high$b) - as.numeric(res_low$b)
se_diff      <- sqrt(res_high$se^2 + res_low$se^2)
z_interact   <- diff_b / se_diff
p_interact_z <- 2 * pnorm(-abs(z_interact))

cat("====== Frequency Subgroup Interaction (Z-test) ======\n")
cat(sprintf("  Difference (High - Low) = %.2f\n", diff_b))
cat(sprintf("  SE of difference        = %.2f\n", se_diff))
cat(sprintf("  z = %.3f, p = %.4f\n\n", z_interact, p_interact_z))

# --- Meta-regression with Knapp-Hartung ---
res_mod <- rma.uni(yi  = updrs$yi,
                   sei = updrs$sei,
                   mods = ~ factor(Group),
                   data = updrs,
                   method = "REML",
                   test = "knha")

cat("====== Frequency Subgroup Interaction (Meta-regression, Knapp-Hartung) ======\n")
cat(sprintf("  Moderator coefficient (Low vs High): %.2f (95%% CI %.2f to %.2f)\n",
            res_mod$b[2], res_mod$ci.lb[2], res_mod$ci.ub[2]))
cat(sprintf("  t = %.3f, p = %.4f\n", res_mod$zval[2], res_mod$pval[2]))
cat(sprintf("  Omnibus test of moderator: QM = %.3f, p = %.4f\n\n",
            res_mod$QM, res_mod$QMp))

# --- Overall P ---
p_overall <- res_overall$pval
cat("====== Overall Pooled Effect ======\n")
cat(sprintf("  Pooled MD = %.2f (95%% CI %.2f to %.2f)\n",
            as.numeric(res_overall$b), res_overall$ci.lb, res_overall$ci.ub))
cat(sprintf("  t = %.3f, %s\n\n", res_overall$zval, fmt_p(p_overall)))

# ============================================================
# 开始绘图
# ============================================================
cairo_pdf("Figure_Forest_UPDRS_co_base_HighLow.pdf",
          width  = 12,
          height = 9,
          family = "Arial")

# ============================================================
# 1. 提取 rma 结果
# ============================================================
b_high     <- as.numeric(res_high$b)
ci_lb_high <- res_high$ci.lb
ci_ub_high <- res_high$ci.ub

b_low     <- as.numeric(res_low$b)
ci_lb_low <- res_low$ci.lb
ci_ub_low <- res_low$ci.ub

b_all     <- as.numeric(res_overall$b)
ci_lb_all <- res_overall$ci.lb
ci_ub_all <- res_overall$ci.ub

# ============================================================
# 2. 提取异质性统计量
# ============================================================
# --- High ---
I2_high   <- round(res_high$I2, 1)
tau2_high <- round(res_high$tau2, 2)
Q_high    <- round(res_high$QE, 2)
Qp_high   <- res_high$QEp

# --- Low ---
I2_low   <- round(res_low$I2, 1)
tau2_low <- round(res_low$tau2, 2)
Q_low    <- round(res_low$QE, 2)
Qp_low   <- res_low$QEp

# --- Overall ---
I2_all   <- round(res_overall$I2, 1)
tau2_all <- round(res_overall$tau2, 2)
Q_all    <- round(res_overall$QE, 2)
Qp_all   <- res_overall$QEp

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
het_high_text <- sprintf(
  "    (I\u00B2 = %.1f%%, \u03C4\u00B2 = %.2f, Q = %.2f, %s)",
  I2_high, tau2_high, Q_high, fmt_p(Qp_high)
)

het_low_text <- sprintf(
  "    (I\u00B2 = %.1f%%, \u03C4\u00B2 = %.2f, Q = %.2f, %s)",
  I2_low, tau2_low, Q_low, fmt_p(Qp_low)
)

het_all_text <- sprintf(
  "    (I\u00B2 = %.1f%%, \u03C4\u00B2 = %.2f, Q = %.2f, %s)",
  I2_all, tau2_all, Q_all, fmt_p(Qp_all)
)

# ============================================================
# 5. 构建数据行（含异质性信息行）
# ============================================================

# --- High 亚组标题 ---
high_label <- tibble(
  Research = "High Frequency (\u2265100 Hz)", Frequency = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = TRUE
)

# --- High 个体研究 ---
high_data <- updrs_high %>%
  transmute(
    Research,
    Frequency = as.character(Frequency),
    mean  = yi,
    lower = yi - 1.96 * sei,
    upper = yi + 1.96 * sei,
    is_summary = FALSE
  )

# --- High subtotal（菱形） ---
high_subtotal <- tibble(
  Research = "  Subtotal (High Frequency)", Frequency = "",
  mean = b_high, lower = ci_lb_high, upper = ci_ub_high,
  is_summary = TRUE
)

# --- High 异质性信息行 ---
high_het_row <- tibble(
  Research = het_high_text, Frequency = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = FALSE
)

# --- 空行 ---
blank_row <- tibble(
  Research = "", Frequency = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = FALSE
)

# --- Low 亚组标题 ---
low_label <- tibble(
  Research = "Low Frequency (<100 Hz)", Frequency = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = TRUE
)

# --- Low 个体研究 ---
low_data <- updrs_low %>%
  transmute(
    Research,
    Frequency = as.character(Frequency),
    mean  = yi,
    lower = yi - 1.96 * sei,
    upper = yi + 1.96 * sei,
    is_summary = FALSE
  )

# --- Low subtotal（菱形） ---
low_subtotal <- tibble(
  Research = "  Subtotal (Low Frequency)", Frequency = "",
  mean = b_low, lower = ci_lb_low, upper = ci_ub_low,
  is_summary = TRUE
)

# --- Low 异质性信息行 ---
low_het_row <- tibble(
  Research = het_low_text, Frequency = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = FALSE
)

# --- Overall（菱形） ---
overall_row <- tibble(
  Research = "Overall Effect", Frequency = "",
  mean = b_all, lower = ci_lb_all, upper = ci_ub_all,
  is_summary = TRUE
)

# --- Overall 异质性信息行 ---
overall_het_row <- tibble(
  Research = het_all_text, Frequency = "",
  mean = NA_real_, lower = NA_real_, upper = NA_real_,
  is_summary = FALSE
)

# ============================================================
# 6. 合并所有行
# ============================================================
forest_all <- bind_rows(
  high_label,
  high_data,
  high_subtotal,
  high_het_row,
  blank_row,
  low_label,
  low_data,
  low_subtotal,
  low_het_row,
  blank_row,
  overall_row,
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
  c("Study",          forest_all$Research),
  c("Frequency (Hz)", forest_all$Frequency),
  c("MD [95% CI]",    ci_text)
)

# ============================================================
# 8. 构建绘图输入
# ============================================================
fn_mean       <- c(NA, forest_all$mean)
fn_lower      <- c(NA, forest_all$lower)
fn_upper      <- c(NA, forest_all$upper)
fn_is_summary <- c(FALSE, forest_all$is_summary)

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
  
  xlab   = "Mean Difference (Favours CO-stimulation \u2190    \u2192 Favours STN-only)",
  clip   = c(-50, 30),
  xticks = seq(-50, 30, by = 10),
  graph.pos = 3,
  title  = "UPDRS-III (CO - STN): Subgroup Analysis by Stimulation Frequency"
)

dev.off()


