library(dplyr)
library(readr)
library(metafor)
library(forestplot)

# 读入你的文件
dat <- read_csv("<ORIGINAL_LOCAL_PATH_REDACTED>")

make_effect <- function(row){
  N  <- as.numeric(row["N"])
  MeanCo <- as.numeric(row["Meanco"])
  MeanCtrl <- as.numeric(row["MeanCtrl"])
  SDco <- as.numeric(row["SDco"])
  SDctrl <- as.numeric(row["SDctrl"])
  md_within <- as.numeric(row["Within__mean_co_STN"])
  se_within <- as.numeric(row["Within_SE_co_STN"])
  md_ind <- as.numeric(row["Mean Difference"])
  se_ind <- as.numeric(row["SE"])
  
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

res_updrs <- rma.uni(yi = updrs$yi,
                     sei = updrs$sei,
                     method = "REML",
                     test = "knha")

library(metafor)
res_on  <- rma.uni(yi, sei, data = updrs_on, method = "REML")
res_off <- rma.uni(yi, sei, data = updrs_off, method = "REML")
res_updrs <- rma.uni(yi = updrs$yi,
                     sei = updrs$sei,
                     method = "REML",
                     test = "knha")


###########
forest_df <- bind_rows(
  updrs_on  %>% mutate(Group = "Medication ON"),
  updrs_off %>% mutate(Group = "Medication OFF")
)
summary_on  <- tibble(Research = "Medicaiton ON",  yi = res_on$b,  sei = res_on$se)
summary_off <- tibble(Research = "Medication OFF", yi = res_off$b, sei = res_off$se)
summary_row <- tibble(
  Research = "Overall effect",
  yi = round(res_updrs$b, 2),sei=res_updrs$se
)

forest_df <- bind_rows(forest_df, summary_on, summary_off,summary_row)

tabletext <- cbind(
  c("Study", forest_df$Research),
  c("Medication", as.character(forest_df$Medication)),
  c("Mean Difference [95% CI]",
    sprintf("%.2f [%.2f, %.2f]",
            forest_df$yi,
            forest_df$yi - 1.96*forest_df$sei,
            forest_df$yi + 1.96*forest_df$sei)),
  c('Forest Plot')
)
tabletext
# 在表格下方添加三行 “ON”, "OFF", "BOTH"
tabletext[c(10, 11, 12), 2] <- c("ON", "OFF", "BOTH")


forestplot(
  labeltext = tabletext,
  mean  = c(NA, forest_df$yi),
  lower = c(NA, forest_df$yi - 1.96*forest_df$sei),
  upper = c(NA, forest_df$yi + 1.96*forest_df$sei),
  zero  = 0,
  col = fpColors(box = "#000",
                 line = "#000",
                 summary = "#000"),
  xlab = "Mean Difference",
  title = "UPDRS-III (Medication ON vs OFF)"
)
str(forest_df)
str(tabletext)


library(metafor)
library(forestplot)
library(grid)
library(dplyr)

# ============================================================
# 1. 运行 Meta-Analysis
# ============================================================
# 筛选有效数据（排除 NA）
valid_idx <- which(!is.na(as.numeric(forest_df$yi)) & !is.na(forest_df$sei))

ma_result <- rma(
  yi  = as.numeric(forest_df$yi[valid_idx]),
  sei = forest_df$sei[valid_idx],
  method = "REML"   # 随机效应模型
)

# 查看结果
summary(ma_result)

# 提取关键参数
pooled_est   <- coef(ma_result)
pooled_se    <- ma_result$se
pooled_lower <- ma_result$ci.lb
pooled_upper <- ma_result$ci.ub
pooled_p     <- ma_result$pval
I2           <- ma_result$I2
tau2         <- ma_result$tau2
Q            <- ma_result$QE
Q_p          <- ma_result$QEp

# ============================================================
# 2. 构建 tabletext（含汇总行和异质性信息）
# ============================================================
# 汇总统计文字
summary_text <- sprintf("Pooled estimate: %.2f [%.2f, %.2f], p = %.3f",
                        pooled_est, pooled_lower, pooled_upper, pooled_p)
heterogeneity_text <- sprintf("I² = %.1f%%, τ² = %.3f, Q = %.2f, p = %.3f",
                              I2, tau2, Q, Q_p)

# 新的 tabletext：原有行 + Overall 行 + 空白行（放异质性信息）
tabletext_new <- rbind(
  tabletext[1, ],                                              # 表头
  tabletext[2:nrow(tabletext), ],                              # 原始数据行（去掉已有的Overall行如果有的话）
  c("Overall (RE Model)", 
    sprintf("%.2f [%.2f, %.2f]", pooled_est, pooled_lower, pooled_upper),
    sprintf("p = %.3f", pooled_p),
    ""),
  c(heterogeneity_text, "", "", "")
)

# 如果你的 tabletext 已经包含了 Overall 行，直接使用下面的方式
# 在最后追加异质性统计行
tabletext_final <- rbind(
  tabletext,
  c(heterogeneity_text, "", "", "")
)

# ============================================================
# 3. 构建 mean/lower/upper 向量（含汇总菱形）
# ============================================================
n_studies <- nrow(forest_df)

# 判断 tabletext 最后几行是否已有 Overall
# 假设 tabletext 已有 Overall 行在第 10 行，我们重新构建：
mean_vec  <- c(NA, as.numeric(forest_df$yi))       # 表头NA + 11个研究
lower_vec <- c(NA, as.numeric(forest_df$yi) - 1.96 * forest_df$sei)
upper_vec <- c(NA, as.numeric(forest_df$yi) + 1.96 * forest_df$sei)

# 如果需要在最后加 pooled 汇总行（追加到异质性文字行）
mean_vec  <- c(mean_vec, NA)       # 异质性行无点
lower_vec <- c(lower_vec, NA)
upper_vec <- c(upper_vec, NA)

# ============================================================
# 4. 标记哪些行是汇总行（菱形显示）
# ============================================================
# 找到 "Overall" 行的位置
is_summary_row <- grepl("Overall", tabletext_final[, 1])

# ============================================================
# ============================================================
# 6. 绘制 Forest Plot
# ============================================================
n_rows <- nrow(tabletext_final)

grid.newpage()


# --- 绘制 forestplot ---
forestplot(
  labeltext  = tabletext_final,
  mean       = mean_vec,
  lower      = lower_vec,
  upper      = upper_vec,
  is.summary = is_summary_row,
  zero       = 0,
  col        = fpColors(box     = "black",
                        line    = "black",
                        summary = "black",
                        zero    = "gray50"),
  xlab       = "Mean Difference",
  title      = "UPDRS-III: Forest Plot of STN+SNr vs. STN Stimulation Effect Across Medication Conditions",
  boxsize    = 0.2,
  new_page   = FALSE,
  # 汇总菱形样式
  shapes_gp  = fpShapesGp(
    default = gpar(col = "black", fill = "black"),
    box     = gpar(col = "black", fill = "black"),
    lines   = gpar(col = "black"),
    summary = gpar(col = "black", fill = "black")
  )
)

for (j in 1:nrow(legend_info)) {
  x_pos <- 0.15 + (j - 1) * 0.28
  grid.rect(x = unit(x_pos, "npc"), y = unit(legend_y, "npc"),
            width = unit(0.025, "npc"), height = unit(0.012, "npc"),
            gp = gpar(fill = legend_info$color[j], col = "gray50"))
  grid.text(legend_info$label[j],
            x = unit(x_pos + 0.02, "npc"), y = unit(legend_y, "npc"),
            just = "left", gp = gpar(fontsize = 7))
}

