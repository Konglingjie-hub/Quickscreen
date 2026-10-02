cairo_pdf("Figure_Forest_FOG_Q_co_base.pdf",
          width  = 12,
          height = 9,
          family = "Arial")
# ============================================================
# 1. 提取 rma 结果
# ============================================================
b_on     <- as.numeric(res_on$b)
ci_lb_on <- res_on$ci.lb
ci_ub_on <- res_on$ci.ub

b_off     <- as.numeric(res_off$b)
ci_lb_off <- res_off$ci.lb
ci_ub_off <- res_off$ci.ub

b_all     <- as.numeric(res_updrs$b)
ci_lb_all <- res_updrs$ci.lb
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

# --- ★ ON 异质性信息行（纯文本，不画图形） ---
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

# --- ★ OFF 异质性信息行 ---
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

# --- ★ Overall 异质性信息行 ---
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
  on_het_row,         # ★ 异质性统计
  blank_row,
  off_label,
  off_data,
  off_subtotal,
  off_het_row,        # ★ 异质性统计
  blank_row,
  overall_row,
  overall_het_row     # ★ 异质性统计
)
#########

forest_all <- bind_rows(
  off_label,
  off_data,
  off_subtotal,
  off_het_row    # ★ 异质性统计
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
  
  xlab   = "Mean Difference (Favours DBS \u2190    \u2192 Favours Baseline)",
  clip   = c(-35, 15),
  xticks = seq(-30, 10, by = 10),
  graph.pos = 3,
  title  = "FOG_AC: Forest Plot with Subgroup Analysis"
)

dev.off()
