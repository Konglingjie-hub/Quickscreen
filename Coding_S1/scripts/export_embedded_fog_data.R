dat <- data.frame(
  study = c(
    "Weiss_2013", "Weiss_2013", "Weiss_2013",
    "Tan", "Artusi", "Cebi", "Horn", "Weiss_2025"
  ),
  assessment = c(
    "axial_UPDRS2_3", "FOG_AC", "FOG_Q",
    "FOG_Q", "NFOG_Q", "FOG_AC", "simple_FOG_AC", "FOG_AC"
  ),
  type = c(
    "Cross", "Cross", "Cross", "Cross",
    "Cross", "Independent", "Cross", "Cross"
  ),
  medication = c("OFF", "OFF", "OFF", "OFF", "ON", "ON", "ON", "OFF"),
  n = c(12, 12, 12, 20, 13, 20, 11, 10),
  mean_ctrl = c(14.25, 14.42, 16.17, 15.75, 16.8, 5.83, 1.727, 20),
  sd_ctrl   = c(5.75, 13.19, 3.83, 4.79, 7.9, 8.37, 2.195, 17.49),
  mean_co   = c(13.42, 8.33, 14.5, 7.4, 16.15, 12.5, 1.545, 21.95),
  sd_co     = c(6.47, 10.91, 4.89, 2.52, 8.76, 6.08, 3.266, 18.52),
  mean_base = c(17.25, 22.17, 14.67, 13, 20.8, 12.17, 4.454, 29.17),
  sd_base   = c(4.31, 11.74, 4.7, 4.96, 4, 6.08, 4.591, 14.45),
  n_co      = c(NA, NA, NA, NA, NA, 10, NA, NA),
  n_ctrl    = c(NA, NA, NA, NA, NA, 10, NA, NA),
  within_mean = c(NA, -13.84, -0.17, NA, NA, NA, -2.91, -10.2),
  within_sd   = c(NA, NA, NA, NA, NA, NA, 3.94, 6.54),
  within_se   = c(NA, 3.27, 1.38, NA, NA, NA, 1.19, 2.068),
  stringsAsFactors = FALSE
)


write.csv(dat, 'data/fog_extracted.csv', row.names=FALSE, na='')
