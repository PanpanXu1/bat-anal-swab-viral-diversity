#!/usr/bin/env Rscript

# Terrain-comparison analysis using a prespecified model formula.

required_packages <- c("mgcv", "emmeans", "ggplot2", "readr", "dplyr", "tidyr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required R packages before running this script: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(mgcv)
  library(emmeans)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

input_filename <- "adjusted_shannon_model_input.csv"
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
workflow_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
input_csv <- file.path(workflow_dir, "input", input_filename)
output_dir <- file.path(workflow_dir, "output")
figures_dir <- file.path(output_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

cleanup_output_paths <- c(
  file.path(output_dir, "primary_model_diagnostics.txt"),
  file.path(output_dir, "sessionInfo.txt"),
  file.path(output_dir, "primary_adjusted_predicted_Shannon_by_terrain_equal_weights.csv"),
  file.path(output_dir, "primary_adjusted_predicted_Shannon_by_terrain_proportional_weights.csv"),
  file.path(output_dir, "mountain_vs_other_terrain_proportional_weights_Tukey.csv"),
  file.path(figures_dir, "raw_Shannon_by_terrain.pdf"),
  file.path(figures_dir, "residuals_vs_fitted.pdf"),
  file.path(figures_dir, "qq_plot.pdf"),
  file.path(figures_dir, "observed_vs_fitted.pdf"),
  file.path(figures_dir, "shannon_vs_log10_pool_size_primary_fit.pdf")
)
unlink(cleanup_output_paths[file.exists(cleanup_output_paths)])

if (!file.exists(input_csv)) {
  stop("Input CSV does not exist: ", input_csv)
}

write_csv_clean <- function(x, path) {
  readr::write_csv(as.data.frame(x), path, na = "")
}

save_plot <- function(plot, stem) {
  times_new_roman_pdf <- function(filename, width, height, ...) {
    grDevices::cairo_pdf(file = filename, width = width, height = height, family = "Times New Roman", ...)
  }
  ggplot2::ggsave(
    file.path(figures_dir, paste0(stem, ".pdf")),
    plot,
    width = 6.5,
    height = 5.2,
    device = times_new_roman_pdf
  )
}

raw <- readr::read_csv(input_csv, show_col_types = FALSE)
required_columns <- c("Terrain", "Shannon index", "Season", "Year", "SampleSize (log10)")
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

dat <- tibble::tibble(
  terrain = factor(raw[["Terrain"]]),
  shannon = suppressWarnings(as.numeric(raw[["Shannon index"]])),
  season = factor(raw[["Season"]]),
  year = factor(raw[["Year"]]),
  # The value is precomputed in the input table. Do not apply log10() again.
  log10_pool_size = suppressWarnings(as.numeric(raw[["SampleSize (log10)"]]))
) |>
  tidyr::drop_na()

if (nrow(dat) == 0) stop("No complete rows remain after filtering.")

# Fit the prespecified adjusted Shannon diversity model.
primary_formula <- shannon ~ terrain + log10_pool_size + season + factor(year)
primary_model <- mgcv::gam(primary_formula, data = dat, method = "REML")
pool_median <- stats::median(dat$log10_pool_size)

emm_equal <- emmeans::emmeans(
  primary_model,
  ~ terrain,
  at = list(log10_pool_size = pool_median),
  weights = "equal",
  data = dat
)
emm_prop <- emmeans::emmeans(
  primary_model,
  ~ terrain,
  at = list(log10_pool_size = pool_median),
  weights = "proportional",
  data = dat
)
emm_equal_df <- as.data.frame(summary(emm_equal, infer = c(TRUE, TRUE)))
emm_prop_df <- as.data.frame(summary(emm_prop, infer = c(TRUE, TRUE)))

emm_combined_df <- dplyr::bind_rows(
  dplyr::mutate(emm_equal_df, weighting = "equal_weights", .before = terrain),
  dplyr::mutate(emm_prop_df, weighting = "proportional_weights", .before = terrain)
)

write_csv_clean(
  emm_combined_df,
  file.path(output_dir, "primary_adjusted_predicted_Shannon_by_terrain.csv")
)

mountain_tukey <- function(emm_object) {
  pairwise_df <- as.data.frame(summary(pairs(emm_object, adjust = "tukey"), infer = c(TRUE, TRUE)))
  mountain_df <- pairwise_df |>
    dplyr::filter(grepl("Mountain", contrast, fixed = TRUE)) |>
    dplyr::mutate(
      reverse_direction = grepl(" - Mountain$", contrast),
      contrast = ifelse(
        reverse_direction,
        paste0("Mountain - ", sub(" - Mountain$", "", contrast)),
        contrast
      ),
      estimate = ifelse(reverse_direction, -estimate, estimate),
      lower.CL_original = lower.CL,
      lower.CL = ifelse(reverse_direction, -upper.CL, lower.CL),
      upper.CL = ifelse(reverse_direction, -lower.CL_original, upper.CL),
      t.ratio = ifelse(reverse_direction, -t.ratio, t.ratio),
      adjustment = "tukey"
    ) |>
    dplyr::select(contrast, estimate, SE, df, lower.CL, upper.CL, t.ratio, p.value, adjustment)
  mountain_df
}

write_csv_clean(
  mountain_tukey(emm_equal),
  file.path(output_dir, "mountain_vs_other_terrain_equal_weights_Tukey.csv")
)

writeLines(
  capture.output(print(summary(primary_model))),
  file.path(output_dir, "primary_model_summary.txt")
)

plot_theme <- ggplot2::theme_bw(base_size = 16, base_family = "Times New Roman") +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 16),
    axis.text = ggplot2::element_text(size = 16)
  )

adjusted_plot <- function(emm_df) {
  ggplot2::ggplot(emm_df, ggplot2::aes(terrain, emmean)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower.CL, ymax = upper.CL),
      width = 0.14, color = "#2C7FB8"
    ) +
    ggplot2::geom_point(size = 3.5, color = "#D95F0E") +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(x = "Terrain", y = "Adjusted predicted Shannon diversity") +
    plot_theme
}

p_equal <- adjusted_plot(emm_equal_df)
p_prop <- adjusted_plot(emm_prop_df)

save_plot(p_equal, "primary_adjusted_predicted_Shannon_equal_weights")
save_plot(p_prop, "primary_adjusted_predicted_Shannon_proportional_weights")

cat("Analysis completed.\n")
cat("Output folder: ", normalizePath(output_dir), "\n", sep = "")
