# ============================================================
# 方向定义:
# CO - STN:      Meanco - MeanCtrl      (MeanCtrl = STN-only)
# CO - baseline: Meanco - Meanbaseline  (Meanbaseline = baseline)
# ============================================================

# 辅助函数
md_se_indep <- function(mean_S, sd_S, n_S, mean_C, sd_C, n_C) {
  md <- mean_C - mean_S
  se <- sqrt(sd_S^2 / n_S + sd_C^2 / n_C)
  list(md = md, se = se)
}

md_se_paired <- function(mean_S, sd_S, mean_C, sd_C, n, rho = 0.5) {
  md <- mean_C - mean_S
  se <- sqrt((sd_S^2 + sd_C^2 - 2 * rho * sd_S * sd_C) / n)
  list(md = md, se = se)
}
colnames(data2)
# ============================================================
# 主函数
# ============================================================
calculate_all_effects <- function(data, rho = 0.5) {
  
  data$MD_co_STN_final   <- NA_real_
  data$SE_co_STN_final   <- NA_real_
  data$MD_co_base_final  <- NA_real_
  data$SE_co_base_final  <- NA_real_
  data$method_co_STN     <- NA_character_
  data$method_co_base    <- NA_character_
  
  for (i in 1:nrow(data)) {
    
    type <- tolower(trimws(as.character(data$Type[i])))
    
    is_paired      <- type %in% c("cross", "crossover", "cross-over", "within", "paired")
    is_independent <- type %in% c("independent", "between", "parallel", "indep","Independent")
    
    # ==============================================================
    # ========== 比较1: CO - STN (Meanco - MeanCtrl) ===============
    # ==============================================================
    
    if (!is.na(data$Within__mean_co_STN[i]) & !is.na(data$Within_SE_co_STN[i])) {
      data$MD_co_STN_final[i] <- data$Within__mean_co_STN[i]
      data$SE_co_STN_final[i] <- data$Within_SE_co_STN[i]
      data$method_co_STN[i]   <- "direct_reported"
      
    } else if (is_paired) {
      if (!is.na(data$Meanco[i]) & !is.na(data$SDco[i]) &
          !is.na(data$MeanCtrl[i]) & !is.na(data$SDctrl[i]) &
          !is.na(data$N[i])) {
        res <- md_se_paired(
          mean_S = data$MeanCtrl[i],       # STN-only 条件
          sd_S   = data$SDctrl[i],
          mean_C = data$Meanco[i],          # CO 条件
          sd_C   = data$SDco[i],
          n      = data$N[i],
          rho    = rho
        )
        data$MD_co_STN_final[i] <- res$md
        data$SE_co_STN_final[i] <- res$se
        data$method_co_STN[i]   <- paste0("paired_rho=", rho)
      }
      
    } else if (is_independent) {
      n_ctrl <- data$Nctrl[i]
      n_co   <- data$Nco[i]
      if (!is.na(data$Meanco[i]) & !is.na(data$SDco[i]) &
          !is.na(data$MeanCtrl[i]) & !is.na(data$SDctrl[i]) &
          !is.na(n_ctrl) & !is.na(n_co)) {
        res <- md_se_indep(
          mean_S = data$MeanCtrl[i],       # STN-only 条件
          sd_S   = data$SDctrl[i],
          n_S    = n_ctrl,
          mean_C = data$Meanco[i],          # CO 条件
          sd_C   = data$SDco[i],
          n_C    = n_co
        )
        data$MD_co_STN_final[i] <- res$md
        data$SE_co_STN_final[i] <- res$se
        data$method_co_STN[i]   <- "independent"
      }
    }
    
    # ==============================================================
    # ========== 比较2: CO - baseline (Meanco - Meanbaseline) ======
    # ==============================================================
    
    if (!is.na(data$Within_mean_co_baseline[i]) & !is.na(data$Within_SE_co_STN[i])) {
      data$MD_co_base_final[i] <- data$Within_mean_co_baseline[i]
      data$SE_co_base_final[i] <- data$Within_SE_co_STN[i]
      data$method_co_base[i]   <- "direct_reported"
      
    } else if (is_paired) {
      if (!is.na(data$Meanco[i]) & !is.na(data$SDco[i]) &
          !is.na(data$Meanbaseline[i]) & !is.na(data$SDbaseline[i]) &
          !is.na(data$N[i])) {
        res <- md_se_paired(
          mean_S = data$Meanbaseline[i],   # baseline 条件
          sd_S   = data$SDbaseline[i],
          mean_C = data$Meanco[i],          # CO 条件
          sd_C   = data$SDco[i],
          n      = data$N[i],
          rho    = rho
        )
        data$MD_co_base_final[i] <- res$md
        data$SE_co_base_final[i] <- res$se
        data$method_co_base[i]   <- paste0("paired_rho=", rho)
      }
      
    } else if (is_independent) {
      n_base <- data$N[i]
      n_co   <- data$Nco[i]library(dplyr)
library(readr)
library(metafor)
library(forestplot)

# 读入你的文件
dat <- read_csv("<ORIGINAL_LOCAL_PATH_REDACTED>")

      if (!is.na(data$Meanco[i]) & !is.na(data$SDco[i]) &
          !is.na(data$Meanbaseline[i]) & !is.na(data$SDbaseline[i]) &
          !is.na(n_base) & !is.na(n_co)) {
        res <- md_se_indep(
          mean_S = data$Meanbaseline[i],   # baseline 条件
          sd_S   = data$SDbaseline[i],
          n_S    = n_base,
          mean_C = data$Meanco[i],          # CO 条件
          sd_C   = data$SDco[i],
          n_C    = n_co
        )
        data$MD_co_base_final[i] <- res$md
        data$SE_co_base_final[i] <- res$se
        data$method_co_base[i]   <- "independent"
      }
    }
  }
  
  return(data)
}

# ============================================================
# 运行
# ============================================================
data <- calculate_all_effects(data2, rho = 0.5)
write.csv(data, "frequency.csv", row.names = FALSE, na = "")
# 查看结果
print(data[, c("Research", "Type", 
               "MD_co_STN_final", "SE_co_STN_final", "method_co_STN",
               "MD_co_base_final", "SE_co_base_final", "method_co_base")])

str(data2)
print(result_check)

