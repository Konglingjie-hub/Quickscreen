library(dplyr)
FoG <- read.csv("<ORIGINAL_LOCAL_PATH_REDACTED>")
FoG <- FoG %>% mutate(
  FoG2_1 = FOG_Q2- FOG_base,
  FOG3_2 = FOG_Q3-FOG_base
)
head(FoG)
Fog2 <- FoG %>%
  mutate(
    improve = ifelse(grepl("b$", Patients),
                     FoG2_1 - FOG3_2,   # 以 b 结尾的组
                     FOG3_2 - FoG2_1)   # 以 a 结尾的组
  )
head(Fog2)
mean(Fog2$improve)
SD(Fog2$improve)
mean(Fog2$improve)
sd(Fog2$improve)

Fog2 <- FoG %>%
  mutate(
    all = ifelse(grepl("b$", Patients),
                     FoG2_1 - FOG3_2,   # 以 b 结尾的组
                     FOG3_2 - FoG2_1)   # 以 a 结尾的组
  )
df <- FoG
# 如果 Patients 是因子，先转字符
df$Patients <- as.character(df$Patients)

idx_a <- grepl("a$", df$Patients)
idx_b <- grepl("b$", df$Patients)

# 组合1：a 的 Q2 + b 的 Q3
long1 <- c(df$FOG_Q2[idx_a], df$FOG_Q3[idx_b])
mean1 <- mean(long1, na.rm = TRUE)
sd1   <- sd(long1, na.rm = TRUE)

# 组合2：b 的 Q2 + a 的 Q3
long2 <- c(df$FOG_Q2[idx_b], df$FOG_Q3[idx_a])
mean2 <- mean(long2, na.rm = TRUE)
sd2   <- sd(long2, na.rm = TRUE)

mean1; sd1
mean2; sd2

