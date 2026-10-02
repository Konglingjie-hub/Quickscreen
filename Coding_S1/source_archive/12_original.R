leave_one_out_manual <- function(data, method = "REML", test = "z") {
  n <- nrow(data)
  results <- tibble(
    excluded_study = character(n),
    estimate = numeric(n),
    se = numeric(n),
    ci.lb = numeric(n),
    ci.ub = numeric(n),
    pval = numeric(n),
    I2 = numeric(n),
    tau2 = numeric(n),
    QE = numeric(n),
    QEp = numeric(n)
  )
  
  for (i in 1:n) {
    data_i <- data[-i, ]
    res_i <- rma.uni(yi = data_i$yi, sei = data_i$sei, method = method, test = test)
    
    results$excluded_study[i] <- data$Research[i]
    results$estimate[i] <- as.numeric(res_i$b)
    results$se[i] <- res_i$se
    results$ci.lb[i] <- res_i$ci.lb
    results$ci.ub[i] <- res_i$ci.ub
    results$pval[i] <- res_i$pval
    results$I2[i] <- res_i$I2
    results$tau2[i] <- res_i$tau2
    results$QE[i] <- res_i$QE
    results$QEp[i] <- res_i$QEp
  }
  
  return(results)
}

# --- 执行 Leave-One-Out ---
l1o_overall <- leave_one_out_manual(updrs, method = "REML", test = "knha")
l1o_med_on  <- leave_one_out_manual(updrs_on, method = "REML", test = "z")
l1o_med_off <- leave_one_out_manual(updrs_off, method = "REML", test = "z")

# 查看结果
print(l1o_overall)
print(l1o_med_on)
print(l1o_med_off)


# ============================================================
# 导出 Leave-One-Out 结果表格
# ============================================================
export_l1o_table <- function(l1o_results, original_res, label) {
  
  orig_est <- as.numeric(original_res$b)
  
  l1o_results %>%
    mutate(
      group = label,
      change_from_original = estimate - orig_est,
      ci_text = sprintf("%.2f [%.2f, %.2f]", estimate, ci.lb, ci.ub),
      I2_pct = sprintf("%.1f%%", I2),
      significance = ifelse(ci.lb > 0 | ci.ub < 0, "Significant", "Not significant"),
      direction_change = case_when(
        (orig_est < 0 & ci.ub >= 0) | (orig_est > 0 & ci.lb <= 0) ~ "YES - Changes significance",
        TRUE ~ "No"
      )
    ) %>%
    select(group, excluded_study, estimate, ci.lb, ci.ub, ci_text, 
           change_from_original, I2_pct, significance, direction_change)
}

# 生成表格
table_overall <- export_l1o_table(l1o_overall, res_updrs, "Overall")
table_on      <- export_l1o_table(l1o_med_on, res_on, "Medication ON")
table_off     <- export_l1o_table(l1o_med_off, res_off, "Medication OFF")

# 合并并导出
l1o_full_table <- bind_rows(table_overall, table_on, table_off)
write_csv(l1o_full_table, "Leave_One_Out_Results_UPDRS3.csv")

# 打印概览
cat("\n========== Leave-One-Out Summary ==========\n")
cat("\n--- Overall ---\n")
cat(sprintf("Original estimate: %.2f [%.2f, %.2f]\n", 
            as.numeric(res_updrs$b), res_updrs$ci.lb, res_updrs$ci.ub))
cat(sprintf("Range after exclusion: %.2f to %.2f\n", 
            min(l1o_overall$estimate), max(l1o_overall$estimate)))
cat(sprintf("Any study changes significance? %s\n",
            ifelse(any(table_overall$direction_change == "YES - Changes significance"), "YES", "NO")))

cat("\n--- Medication ON ---\n")
cat(sprintf("Original estimate: %.2f [%.2f, %.2f]\n", 
            as.numeric(res_on$b), res_on$ci.lb, res_on$ci.ub))
cat(sprintf("Range after exclusion: %.2f to %.2f\n", 
            min(l1o_med_on$estimate), max(l1o_med_on$estimate)))

cat("\n--- Medication OFF ---\n")
cat(sprintf("Original estimate: %.2f [%.2f, %.2f]\n", 
            as.numeric(res_off$b), res_off$ci.lb, res_off$ci.ub))
cat(sprintf("Range after exclusion: %.2f to %.2f\n", 
            min(l1o_med_off$estimate), max(l1o_med_off$estimate)))

library(ggplot2)

# ============================================================
# ggplot2 版本的 Leave-One-Out 图
# ============================================================
plot_l1o_ggplot <- function(l1o_results, original_res, title_text) {
  
  orig_est <- as.numeric(original_res$b)
  orig_lb  <- original_res$ci.lb
  orig_ub  <- original_res$ci.ub
  
  # 准备数据
  df <- l1o_results %>%
    mutate(
      excluded_study = factor(excluded_study, levels = rev(excluded_study))
    )
  
  # 绘图
  p <- ggplot(df, aes(x = estimate, y = excluded_study)) +
    # 原始汇总效应的置信区间（浅色背景带）
    annotate("rect", 
             xmin = orig_lb, xmax = orig_ub, 
             ymin = -Inf, ymax = Inf,
             alpha = 0.15, fill = "red") +
    # 原始汇总效应线
    geom_vline(xintercept = orig_est, linetype = "dashed", color = "red", linewidth = 0.8) +
    # 零线
    geom_vline(xintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.5) +
    # 每个研究排除后的置信区间
    geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), height = 0.2, linewidth = 0.6, color = "steelblue") +
    # 每个研究排除后的点估计
    geom_point(size = 3, color = "steelblue", shape = 15) +
    # 标签
    labs(
      title = title_text,
      subtitle = sprintf("Original pooled estimate: %.2f [%.2f, %.2f] (red dashed line & shaded area)",
                         orig_est, orig_lb, orig_ub),
      x = "Mean Difference",
      y = "Study Excluded"
    ) +
    theme_minimal(base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# --- 绘图 ---

library(patchwork)

# --- 生成三个图 ---
p_overall <- plot_l1o_ggplot(l1o_overall, res_updrs, 
                             "Leave-One-Out: UPDRS-III Overall")

p_on <- plot_l1o_ggplot(l1o_med_on, res_on, 
                        "Leave-One-Out: UPDRS-III Medication ON")

p_off <- plot_l1o_ggplot(l1o_med_off, res_off, 
                         "Leave-One-Out: UPDRS-III Medication OFF")

# --- 合并排列并输出为一个 PDF ---
# 方案1: 垂直排列（3行1列）
combined_plot <- p_overall / p_on / p_off +
  plot_annotation(
    title = "Leave-One-Out Sensitivity Analysis: UPDRS-III",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
  )
combined_plot

ggsave("Figure_LeaveOneOut_UPDRS_combined_co_stn.pdf",
       plot = combined_plot,
       width = 10,
       height = 14,
       device = cairo_pdf)



