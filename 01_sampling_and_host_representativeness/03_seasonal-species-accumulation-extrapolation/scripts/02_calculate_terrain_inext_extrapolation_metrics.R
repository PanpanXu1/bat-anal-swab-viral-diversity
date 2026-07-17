# Calculate terrain-level iNEXT extrapolation metrics used to annotate the terrain accumulation panel.

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
tables_dir <- file.path(workflow_dir, "output", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("Number", "Terrain", "Species")
if (length(setdiff(required, names(data))) > 0) {
  stop("Input data are missing required columns.")
}
data[] <- lapply(data, function(x) trimws(as.character(x)))
data <- data[data$Number != "" & data$Species != "", ]
data <- data[!duplicated(data$Number), ]

terrain_order <- c("Plain", "Mountain", "Hill", "Mesa")

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

rows <- lapply(terrain_order, function(terrain) {
  subset <- data[data$Terrain == terrain, ]
  abund <- abundance_vector(subset)
  curve <- get_curve(abund)
  n <- sum(abund)
  q_n <- interpolate_qd(curve, n)
  q_2n <- interpolate_qd(curve, 2 * n)
  q_plus_100 <- interpolate_qd(curve, n + 100)
  data.frame(
    Terrain = terrain,
    N = n,
    Observed_richness = length(abund),
    Gain_to_2N = q_2n - q_n,
    Slope_per_100 = q_plus_100 - q_n,
    stringsAsFactors = FALSE
  )
})

metrics <- do.call(rbind, rows)
write.csv(
  metrics,
  file.path(tables_dir, "terrain_inext_extrapolation_metrics.csv"),
  row.names = FALSE
)
print(metrics)
