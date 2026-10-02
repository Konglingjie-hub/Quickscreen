library(dplyr)

r <- 0.5  # 配对相关系数（可根据实际调整）
df <- read.csv("<ORIGINAL_LOCAL_PATH_REDACTED>")
df <- df %>%
  mutate(
    # ====== 第一对：Within_mean / Within_SE ======
    Within_mean_co_baseline = ifelse(
      is.na(Within_mean_co_baseline),
      Meanco - Meanbaseline,
      Within_mean_co_baseline
    ),
    
    Within_SE_co_baseline = ifelse(
      is.na(Within_SE_co_baseline),
      case_when(
        Type == "Cross" ~ 
          sqrt(SDco^2 + SDbaseline^2 - 2 * r * SDco * SDbaseline) / sqrt(N),
        Type == "Independent" ~ 
          sqrt(SDco^2 / N + SDbaseline^2 / N)
      ),
      Within_SE_co_baseline
    ),
    
    # ====== 第二对：Mean_Difference / SE ======
    Mean_Difference_CO_base = ifelse(
      is.na(Mean_Difference_CO_base),
      Meanco - Meanbaseline,
      Mean_Difference_CO_base
    ),
    
    SE_CO_base = ifelse(
      is.na(SE_CO_base),
      case_when(
        Type == "Cross" ~ 
          sqrt(SDco^2 + SDbaseline^2 - 2 * r * SDco * SDbaseline) / sqrt(N),
        Type == "Independent" ~ 
          sqrt(SDco^2 / N + SDbaseline^2 / N)
      ),
      SE_CO_base
    )
  )
write.csv(df, file = "output.csv", row.names = FALSE, na = "")

