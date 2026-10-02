# 计算 Low 和 High（95% 置信区间）
library(dplyr)
dat <- read.csv()
colnames(dat)
dat <- dat %>% mutate(Low=Within__mean_co_STN-1.96*Within__mean_co_STN,
                      High=Within__mean_co_STN+1.96*Within__mean_co_STN,
                      Mean_lowHigh_text= paste0(
                        round(Within__mean_co_STN, 2), " (", 
                        round(Low, 2), ", ", 
                        round(High, 2), ")"
                      ))

Low  <- data$Within__mean_co_STN - 1.96 * data$Within_SE_co_STN
High <- data$Within__mean_co_STN + 1.96 * data$Within_SE_co_STN

dat <- read.csv("<ORIGINAL_LOCAL_PATH_REDACTED>")
# 生成文本列：格式为 "Mean (Low, High)"
dat <- dat %>% mutate(Mean_low_high = paste0(
  round(Within_mean_co_baseline, 2), " (", 
  round(Lower, 2), ", ", 
  round(Upper, 2), ")"
) )
Mean_LowHigh_text <- paste0(
  round(Within_mean_co_baseline, 2), " (", 
  round(Lower, 2), ", ", 
  round(Upper, 2), ")"
)
colnames(dat)
# 合并为数据框
result <- data.frame(
  Low              = round(Low, 2),
  High             = round(High, 2),
  Mean_LowHigh     = Mean_LowHigh_text
)

head(result)
write.csv(dat,'forest_co_base_plot.csv')
6.67+1.96*5.28

