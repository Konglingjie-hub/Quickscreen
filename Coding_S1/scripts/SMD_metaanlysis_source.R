# Stratified SMD meta-analysis for freezing of gait outcomes
#
# Four prespecified analyses are produced:
#   A. co-stimulation vs STN mono-stimulation, FOG_Q outcomes
#   B. co-stimulation vs STN mono-stimulation, FOG_AC outcomes
#   C. co-stimulation vs baseline, FOG_Q outcomes
#   D. co-stimulation vs baseline, FOG_AC outcomes
#
# A lower score represents less severe freezing of gait. Therefore, a negative
# Hedges' g favours STN-SNr co-stimulation in both comparisons.

# Windows R can start in the C locale, which prevents access to a user library
# when the Windows user name contains Unicode characters. Switching to UTF-8 is
# harmless when already active and makes this script reproducible from Rscript.
if (.Platform$OS.type == "windows") {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE))
}

required_packages <- c("metafor", "forestplot")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them once with install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "), "))."
  )
}

library(metafor)
library(forestplot)
library(grid)

# =============================================================================
# 1. Data
# =============================================================================

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

# Exclude axial_UPDRS2_3 exactly as in the original analysis.
dat2 <- dat[-1L, , drop = FALSE]

# Recover an available within-person change SD from its standard error.
replace_sd <- is.na(dat2$within_sd) & !is.na(dat2$within_se)
dat2$within_sd[replace_sd] <- dat2$within_se[replace_sd] * sqrt(dat2$n[replace_sd])

# "Contains FOG_Q" and "contains FOG_AC" implement the requested split. This
# deliberately classifies NFOG_Q as FOG_Q and simple_FOG_AC as FOG_AC.
dat2$outcome_group <- ifelse(
  grepl("FOG_Q", dat2$assessment, fixed = TRUE), "FOG_Q",
  ifelse(grepl("FOG_AC", dat2$assessment, fixed = TRUE), "FOG_AC", NA_character_)
)
dat2$study_outcome <- paste(dat2$study, dat2$assessment, sep = " - ")

if (anyNA(dat2$outcome_group)) {
  stop(
    "The following assessment(s) were not classified as FOG_Q or FOG_AC: ",
    paste(dat2$assessment[is.na(dat2$outcome_group)], collapse = ", ")
  )
}

# =============================================================================
# 2. Effect-size calculation
# =============================================================================

default_paired_r <- 0.50

calculate_contrast <- function(data, comparison, default_r = default_paired_r) {
  stopifnot(comparison %in% c("co_vs_stn", "co_vs_base"))

  comparator_mean <- if (comparison == "co_vs_stn") data$mean_ctrl else data$mean_base
  comparator_sd   <- if (comparison == "co_vs_stn") data$sd_ctrl else data$sd_base
  comparator_name <- if (comparison == "co_vs_stn") {
    "STN mono-stimulation"
  } else {
    "baseline"
  }

  out <- data
  out$comparison <- comparison
  out$comparator <- comparator_name
  out$yi <- NA_real_
  out$vi <- NA_real_
  out$paired_r <- NA_real_
  out$r_source <- NA_character_

  for (i in seq_len(nrow(out))) {
    if (out$type[i] == "Cross") {
      # For co vs baseline, a reported/derived change SD permits study-specific
      # recovery of the within-person correlation. No corresponding change SD
      # is available for co vs STN, so the prespecified r = 0.50 is used.
      r_i <- default_r
      r_source_i <- "assumed"
      if (comparison == "co_vs_base" && !is.na(out$within_sd[i])) {
        r_candidate <- (
          out$sd_co[i]^2 + comparator_sd[i]^2 - out$within_sd[i]^2
        ) / (2 * out$sd_co[i] * comparator_sd[i])
        if (is.finite(r_candidate) && r_candidate >= -1 && r_candidate <= 1) {
          r_i <- r_candidate
          r_source_i <- "derived from change SD"
        }
      }

      sd_pooled <- sqrt((out$sd_co[i]^2 + comparator_sd[i]^2) / 2)
      d_i <- (out$mean_co[i] - comparator_mean[i]) / sd_pooled
      correction_j <- 1 - 3 / (4 * (out$n[i] - 1) - 1)

      out$yi[i] <- correction_j * d_i
      out$vi[i] <- correction_j^2 * (
        2 * (1 - r_i) / out$n[i] + d_i^2 / (2 * out$n[i])
      )
      out$paired_r[i] <- r_i
      out$r_source[i] <- r_source_i
    } else if (out$type[i] == "Independent") {
      # n_ctrl is also the available comparator-arm n for the baseline contrast
      # in the independent Cebi study (n_co = 10; comparator n = 10).
      esc <- metafor::escalc(
        measure = "SMD",
        m1i = out$mean_co[i], sd1i = out$sd_co[i], n1i = out$n_co[i],
        m2i = comparator_mean[i], sd2i = comparator_sd[i], n2i = out$n_ctrl[i]
      )
      out$yi[i] <- as.numeric(esc$yi)
      out$vi[i] <- as.numeric(esc$vi)
      out$r_source[i] <- "independent groups"
    } else {
      stop("Unknown design type: ", out$type[i])
    }
  }

  out$ci_lb <- out$yi - qnorm(0.975) * sqrt(out$vi)
  out$ci_ub <- out$yi + qnorm(0.975) * sqrt(out$vi)
  out
}

effects <- rbind(
  calculate_contrast(dat2, "co_vs_stn"),
  calculate_contrast(dat2, "co_vs_base")
)
row.names(effects) <- NULL

# =============================================================================
# 3. Separate random-effects models by comparison and outcome group
# =============================================================================

analysis_specs <- data.frame(
  analysis_id = c("A", "B", "C", "D"),
  comparison = c("co_vs_stn", "co_vs_stn", "co_vs_base", "co_vs_base"),
  outcome_group = c("FOG_Q", "FOG_AC", "FOG_Q", "FOG_AC"),
  panel_title = c(
    "A) Co-stimulation vs STN mono-stimulation: FOG_Q",
    "B) Co-stimulation vs STN mono-stimulation: FOG_AC",
    "C) Co-stimulation vs baseline: FOG_Q",
    "D) Co-stimulation vs baseline: FOG_AC"
  ),
  stringsAsFactors = FALSE
)

models <- vector("list", nrow(analysis_specs))
panel_data <- vector("list", nrow(analysis_specs))
names(models) <- names(panel_data) <- analysis_specs$analysis_id

for (j in seq_len(nrow(analysis_specs))) {
  spec <- analysis_specs[j, ]
  keep <- effects$comparison == spec$comparison &
    effects$outcome_group == spec$outcome_group
  d <- effects[keep, , drop = FALSE]

  if (nrow(d) < 2L) {
    stop("Analysis ", spec$analysis_id, " has fewer than two effect sizes.")
  }

  # REML random-effects model with Hartung-Knapp inference. This is the
  # single-level analogue of the original script's t-based multilevel model;
  # after stratification there is only one outcome per study in each model.
  fit <- metafor::rma.uni(
    yi = yi,
    vi = vi,
    data = d,
    method = "REML",
    test = "knha"
  )

  panel_data[[spec$analysis_id]] <- d
  models[[spec$analysis_id]] <- fit
}

# =============================================================================
# 4. Publication-style forest plots
# =============================================================================

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.3f", p)
}

make_panel <- function(data, model, title, comparison) {
  k <- nrow(data)
  pooled_est <- as.numeric(coef(model)[1])
  pooled_lb <- as.numeric(model$ci.lb)
  pooled_ub <- as.numeric(model$ci.ub)
  prediction <- predict(model)
  pi_lb <- as.numeric(prediction$pi.lb)
  pi_ub <- as.numeric(prediction$pi.ub)

  effect_text <- sprintf("%.2f [%.2f, %.2f]", data$yi, data$ci_lb, data$ci_ub)
  pooled_text <- sprintf("%.2f [%.2f, %.2f]", pooled_est, pooled_lb, pooled_ub)
  pi_text <- sprintf("%.2f [%.2f, %.2f]", pooled_est, pi_lb, pi_ub)
  heterogeneity_1 <- sprintf("I\u00B2 = %.1f%%; \u03C4\u00B2 = %.3f", model$I2, model$tau2)
  heterogeneity_2 <- sprintf(
    "Q = %.2f, df = %d, %s",
    model$QE, model$k - model$p, format_p(model$QEp)
  )

  labels <- list(
    c(
      "Study - Outcome", data$study_outcome, NA,
      "Overall (RE Model)", "Prediction Interval",
      heterogeneity_1, heterogeneity_2
    ),
    c("N", as.character(data$n), NA, NA, NA, NA, NA),
    c("Medication", data$medication, NA, NA, NA, NA, NA),
    c("SMD [95% CI]", effect_text, NA, pooled_text, pi_text, NA, NA)
  )

  row_mean  <- c(NA, data$yi, NA, pooled_est, pooled_est, NA, NA)
  row_lower <- c(NA, data$ci_lb, NA, pooled_lb, pi_lb, NA, NA)
  row_upper <- c(NA, data$ci_ub, NA, pooled_ub, pi_ub, NA, NA)
  n_rows <- length(row_mean)
  prediction_row <- k + 4L

  line_styles <- replicate(n_rows, gpar(col = "black", lwd = 1.25), simplify = FALSE)
  box_styles <- replicate(
    n_rows,
    gpar(col = "black", fill = "black"),
    simplify = FALSE
  )
  line_styles[[prediction_row]] <- gpar(col = "#D55E00", lty = 2, lwd = 1.6)
  box_styles[[prediction_row]] <- gpar(col = "#D55E00", fill = "white", lwd = 1.2)

  xlab <- if (comparison == "co_vs_stn") {
    "Hedges' g  (Favours STN-SNr co-stimulation \u2190    \u2192 Favours STN mono-stimulation)"
  } else {
    "Hedges' g  (Favours STN-SNr co-stimulation \u2190    \u2192 Favours baseline)"
  }

  forestplot::forestplot(
    labeltext = labels,
    mean = row_mean,
    lower = row_lower,
    upper = row_upper,
    is.summary = c(
      TRUE, rep(FALSE, k), FALSE, TRUE, FALSE, FALSE, FALSE
    ),
    zero = 0,
    clip = c(-4, 4),
    xticks = seq(-3, 3, by = 1),
    graph.pos = 4,
    graphwidth = unit(44, "mm"),
    colgap = unit(3.5, "mm"),
    boxsize = 0.19,
    lwd.ci = 1.25,
    ci.vertices = TRUE,
    ci.vertices.height = 0.10,
    lwd.zero = 0.9,
    xlab = xlab,
    title = title,
    col = forestplot::fpColors(
      box = "black",
      line = "black",
      summary = "black",
      zero = "gray35"
    ),
    txt_gp = forestplot::fpTxtGp(
      label = gpar(fontfamily = "Arial", cex = 0.69),
      ticks = gpar(fontfamily = "Arial", cex = 0.63),
      xlab = gpar(fontfamily = "Arial", cex = 0.68),
      title = gpar(fontfamily = "Arial", cex = 0.90, fontface = "bold"),
      summary = gpar(fontfamily = "Arial", cex = 0.70, fontface = "bold")
    ),
    shapes_gp = forestplot::fpShapesGp(
      default = gpar(col = "black", fill = "black"),
      lines = line_styles,
      box = box_styles,
      summary = gpar(col = "black", fill = "black")
    ),
    mar = unit(c(5, 3, 5, 3), "mm")
  )
}

panels <- vector("list", nrow(analysis_specs))
names(panels) <- analysis_specs$analysis_id
for (j in seq_len(nrow(analysis_specs))) {
  id <- analysis_specs$analysis_id[j]
  panels[[id]] <- make_panel(
    panel_data[[id]],
    models[[id]],
    analysis_specs$panel_title[j],
    analysis_specs$comparison[j]
  )
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}
output_dir <- file.path(script_dir, "SMD_by_outcome")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

draw_combined_figure <- function(filename, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    cairo_pdf(filename, width = 17, height = 13, family = "Arial")
  } else {
    png(
      filename,
      width = 17,
      height = 13,
      units = "in",
      res = 300,
      type = "cairo",
      bg = "white"
    )
  }
  on.exit(dev.off(), add = TRUE)

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  for (j in seq_len(nrow(analysis_specs))) {
    row_i <- if (j <= 2L) 1L else 2L
    col_i <- if (j %% 2L == 1L) 1L else 2L
    pushViewport(viewport(layout.pos.row = row_i, layout.pos.col = col_i))
    print(panels[[analysis_specs$analysis_id[j]]])
    upViewport()
  }
  upViewport()
}

draw_single_panel <- function(panel, filename) {
  cairo_pdf(filename, width = 10.5, height = 7.2, family = "Arial")
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  print(panel)
}

combined_pdf <- file.path(output_dir, "Figure_SMD_by_outcome_2x2.pdf")
combined_png <- file.path(output_dir, "Figure_SMD_by_outcome_2x2.png")
draw_combined_figure(combined_pdf, "pdf")
draw_combined_figure(combined_png, "png")

for (id in analysis_specs$analysis_id) {
  suffix <- paste(
    panel_data[[id]]$comparison[1],
    panel_data[[id]]$outcome_group[1],
    sep = "_"
  )
  draw_single_panel(
    panels[[id]],
    file.path(output_dir, paste0("Panel_", id, "_", suffix, ".pdf"))
  )
}

# =============================================================================
# 5. Machine-readable results
# =============================================================================

effect_export <- effects[, c(
  "comparison", "outcome_group", "study", "assessment", "study_outcome",
  "type", "medication", "n", "yi", "vi", "ci_lb", "ci_ub",
  "paired_r", "r_source"
)]
write.csv(
  effect_export,
  file.path(output_dir, "SMD_study_effects.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

summary_rows <- lapply(seq_len(nrow(analysis_specs)), function(j) {
  id <- analysis_specs$analysis_id[j]
  fit <- models[[id]]
  pred <- predict(fit)
  data.frame(
    analysis_id = id,
    comparison = analysis_specs$comparison[j],
    outcome_group = analysis_specs$outcome_group[j],
    k = fit$k,
    pooled_g = as.numeric(coef(fit)[1]),
    ci_lb = as.numeric(fit$ci.lb),
    ci_ub = as.numeric(fit$ci.ub),
    p_value = as.numeric(fit$pval),
    tau2 = as.numeric(fit$tau2),
    I2_percent = as.numeric(fit$I2),
    Q = as.numeric(fit$QE),
    Q_df = fit$k - fit$p,
    Q_p_value = as.numeric(fit$QEp),
    prediction_lb = as.numeric(pred$pi.lb),
    prediction_ub = as.numeric(pred$pi.ub),
    stringsAsFactors = FALSE
  )
})
summary_export <- do.call(rbind, summary_rows)
write.csv(
  summary_export,
  file.path(output_dir, "SMD_model_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nCompleted four stratified random-effects SMD meta-analyses.\n")
print(summary_export, row.names = FALSE, digits = 4)
cat("\nOutputs written to:\n", normalizePath(output_dir, winslash = "/"), "\n")

