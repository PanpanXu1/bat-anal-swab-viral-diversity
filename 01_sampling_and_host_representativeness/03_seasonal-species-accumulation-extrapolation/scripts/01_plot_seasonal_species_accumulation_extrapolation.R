# Plot seasonal species accumulation and iNEXT-based extrapolation curves.

suppressPackageStartupMessages(library(iNEXT))
set.seed(42)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(sub("^--file=", "", script_arg[1])) else getwd()
workflow_dir <- if (file.exists(file.path("input", "individual_records_with_season_and_terrain.csv"))) {
  "."
} else {
  file.path(script_dir, "..")
}
input_csv <- file.path(workflow_dir, "input", "individual_records_with_season_and_terrain.csv")
figures_dir <- file.path(workflow_dir, "output", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("Number", "Season", "Species")
if (length(setdiff(required, names(data))) > 0) {
  stop("Input data are missing required columns.")
}
data[] <- lapply(data, function(x) trimws(as.character(x)))
data <- data[data$Number != "" & data$Species != "", ]
data <- data[!duplicated(data$Number), ]

season_order <- c("Spring", "Summer", "Autumn", "Winter")
season_colors <- c(
  Spring = "#6BAF5E",
  Summer = "#F0A27A",
  Autumn = "#7C79B8",
  Winter = "#D9D98A"
)

abundance_vector <- function(subset) as.numeric(table(subset$Species))

get_curve <- function(abund) {
  estimate <- iNEXT(
    abund,
    q = 0,
    datatype = "abundance",
    endpoint = 2 * sum(abund),
    se = TRUE,
    nboot = 200
  )$iNextEst

  if (is.list(estimate) && !is.data.frame(estimate)) estimate <- estimate$size_based
  if ("Method" %in% names(estimate)) names(estimate)[names(estimate) == "Method"] <- "method"
  if ("m" %in% names(estimate)) names(estimate)[names(estimate) == "m"] <- "m_size"
  if ("x" %in% names(estimate)) names(estimate)[names(estimate) == "x"] <- "m_size"
  if ("size" %in% names(estimate)) names(estimate)[names(estimate) == "size"] <- "m_size"
  estimate[estimate$Order.q == 0, ]
}

interpolate_qd <- function(curve, target) {
  approx(curve$m_size, curve$qD, xout = target, rule = 2)$y
}

run_group <- function(group_column, level) {
  subset <- data[data[[group_column]] == level, ]
  abund <- abundance_vector(subset)
  curve <- get_curve(abund)
  n <- sum(abund)
  q_n <- interpolate_qd(curve, n)
  q_2n <- interpolate_qd(curve, 2 * n)
  q_plus_100 <- interpolate_qd(curve, n + 100)
  list(
    curve = curve,
    n = n,
    gain_to_2n = q_2n - q_n,
    slope_per_100 = q_plus_100 - q_n
  )
}

curve_ylim <- function(curves) {
  upper_values <- unlist(lapply(curves, function(result) {
    curve <- result$curve
    if ("qD.UCL" %in% names(curve)) curve$qD.UCL else curve$qD
  }))
  c(0, ceiling(max(upper_values, na.rm = TRUE) * 1.05))
}

plot_season <- function(season, result) {
  curve <- result$curve
  color <- season_colors[[season]]
  y_limits <- curve_ylim(list(result))
  output_file <- file.path(
    figures_dir,
    sprintf("seasonal_species_accumulation_%s.pdf", tolower(season))
  )

  pdf(output_file, width = 4.2, height = 2.6, family = "Times")
  par(mar = c(3.4, 3.2, 0.3, 0.4), mgp = c(1.9, 0.45, 0), tcl = -0.25)
  plot(
    curve$m_size, curve$qD,
    type = "n",
    xlab = "Number of Samples",
    ylab = "Species Richness",
    ylim = y_limits,
    cex.lab = 0.75,
    cex.axis = 0.65,
    las = 1
  )
  grid(col = "#E6E6E6", lwd = 0.8)
  if (all(c("qD.LCL", "qD.UCL") %in% names(curve))) {
    polygon(
      c(curve$m_size, rev(curve$m_size)),
      c(curve$qD.LCL, rev(curve$qD.UCL)),
      col = adjustcolor(color, alpha.f = 0.38),
      border = NA
    )
  }
  lines(curve$m_size, curve$qD, lwd = 1.4, col = "black")
  box(lwd = 0.8)
  label_text <- sprintf(
    "Season: %s\nGain (2N) = %.1f | Slope (per100) = %.1f",
    season,
    result$gain_to_2n,
    result$slope_per_100
  )
  text(
    x = max(curve$m_size, na.rm = TRUE) * 0.33,
    y = y_limits[2] * 0.13,
    labels = label_text,
    adj = c(0, 0),
    cex = 0.62
  )
  dev.off()
}

for (season in season_order) {
  plot_season(season, run_group("Season", season))
}
