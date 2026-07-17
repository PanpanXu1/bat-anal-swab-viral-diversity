options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260712)
if (!requireNamespace("mgcv", quietly = TRUE)) stop("Required package 'mgcv' is unavailable.")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
workflow_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
input_csv <- file.path(workflow_dir, "input", "environmental_variables_selected_by_corr0.8_VIF10.csv")
table_dir <- file.path(workflow_dir, "output", "tables")
figure_dir <- file.path(workflow_dir, "output", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Deliberately fixed in the order used by preliminary_environmental_contribution.csv.
env_labels <- c(
  "Global Linear Hydrography (GLH)", "Fractional Vegetation Cover (FVC)",
  "Normalized Difference Vegetation Index (NDVI)", "Mean Temperature of Wettest Quarter (Bio8)",
  "Temperature Seasonality (standard deviation x100) (Bio4)", "Global Railway (GR)",
  "Precipitation of Driest Quarter (Bio17)", "Global Mammal Richness (GMR)",
  "China High-Resolution Ecological Environment Quality (CHEQ)", "Human Footprint (HFT)",
  "Mean Diurnal Range (Mean of monthly max temp - min temp) (Bio2)",
  "China Population Spatial Distribution (PSD)", "Max Temperature of Warmest Month (Bio5)",
  "Precipitation of Warmest Quarter (Bio18)"
)
env_vars <- paste0("env_", seq_along(env_labels))
landuse_cols <- c(
  "landuseName Forest land closed canopy", "landuseName High coverage grassland",
  "landuseName Low coverage grassland", "landuseName Medium coverage grassland",
  "landuseName Other construction land", "landuseName Other forest land",
  "landuseName Paddy field", "landuseName Reservoirs and ponds",
  "landuseName Rivers and canals", "landuseName Rural residential area",
  "landuseName Shrubland", "landuseName Sparse woodland", "landuseName Urban built up area"
)
required <- c("Shannon index", "SampleSize (log10)", "Season", "Year", "Terrain", "Host genus",
              "longitude", "latitude", landuse_cols, env_labels)
raw <- utils::read.csv(input_csv, check.names = FALSE)
missing_columns <- setdiff(required, names(raw))
if (length(missing_columns)) stop("Missing required columns: ", paste(missing_columns, collapse = ", "))

strict_numeric <- function(x, column, allow_missing = TRUE) {
  original_missing <- is.na(x) | (is.character(x) & trimws(x) == "")
  converted <- suppressWarnings(as.numeric(x))
  bad_conversion <- !original_missing & is.na(converted)
  if (any(bad_conversion)) {
    examples <- unique(as.character(x[bad_conversion]))
    stop("Numeric conversion failed for column '", column, "'; example value(s): ",
         paste(utils::head(examples, 3), collapse = ", "))
  }
  if (!allow_missing && any(is.na(converted))) stop("Column '", column, "' contains missing values.")
  if (any(!is.na(converted) & !is.finite(converted))) stop("Column '", column, "' contains non-finite numeric values.")
  converted
}
raw$shannon <- strict_numeric(raw[["Shannon index"]], "Shannon index")
raw$log10_pool_size <- strict_numeric(raw[["SampleSize (log10)"]], "SampleSize (log10)")
raw$season <- factor(raw[["Season"]]); raw$year <- factor(raw[["Year"]])
raw$terrain <- factor(raw[["Terrain"]]); raw$host_genus <- factor(raw[["Host genus"]])
raw$longitude <- strict_numeric(raw[["longitude"]], "longitude")
raw$latitude <- strict_numeric(raw[["latitude"]], "latitude")
for (i in seq_along(env_vars)) raw[[env_vars[i]]] <- strict_numeric(raw[[env_labels[i]]], env_labels[i])

landuse_mapping <- data.frame(
  original_landuse = c("landuseName Forest land closed canopy", "landuseName Other forest land",
    "landuseName Sparse woodland", "landuseName Shrubland", "landuseName High coverage grassland",
    "landuseName Medium coverage grassland", "landuseName Low coverage grassland", "landuseName Paddy field",
    "landuseName Other construction land", "landuseName Urban built up area", "landuseName Rural residential area",
    "landuseName Reservoirs and ponds", "landuseName Rivers and canals"),
  landuse_broad = c("Closed_or_other_forest", "Closed_or_other_forest", "Sparse_woodland", rep("Open_vegetation", 4),
                    "Cropland", rep("Built_up", 3), rep("Other_or_unclassified", 2)))
land_mat <- matrix(0, nrow(raw), nrow(landuse_mapping), dimnames = list(NULL, landuse_mapping$original_landuse))
for (v in landuse_mapping$original_landuse) land_mat[, v] <- strict_numeric(raw[[v]], v, allow_missing = FALSE)
if (any(!land_mat %in% c(0, 1))) {
  where <- which(!land_mat %in% c(0, 1), arr.ind = TRUE)[1, ]
  stop("Land-use dummy column '", colnames(land_mat)[where[2]], "' must contain only 0/1; found ", land_mat[where[1], where[2]])
}
dummy_sum <- rowSums(land_mat); landuse_broad <- rep("Other_or_unclassified", nrow(raw))
singles <- which(dummy_sum == 1)
if (length(singles)) {
  hit <- apply(land_mat[singles, , drop = FALSE], 1, function(z) colnames(land_mat)[which(z == 1)[1]])
  landuse_broad[singles] <- landuse_mapping$landuse_broad[match(hit, landuse_mapping$original_landuse)]
}
multiples <- which(dummy_sum > 1)
if (length(multiples) >= 5) landuse_broad[multiples] <- "Multiple_landuse"
# With fewer than five multi-hot rows, retain the final-model rule: Other_or_unclassified.
if (any(dummy_sum < 0 | dummy_sum > ncol(land_mat))) stop("Invalid land-use row encoding detected.")
raw$landuse_broad <- factor(landuse_broad)
model_fields <- c("shannon", "log10_pool_size", "season", "year", "terrain", "host_genus", "longitude", "latitude",
                  "landuse_broad", env_vars)
dat <- droplevels(raw[stats::complete.cases(raw[, model_fields]), , drop = FALSE])
if (nrow(dat) < 20) stop("Too few common complete cases for model assessment: ", nrow(dat))
cat(sprintf("Rows: original_n=%d, complete_n=%d, deleted_n=%d\n", nrow(raw), nrow(dat), nrow(raw) - nrow(dat)))

base_terms <- c("log10_pool_size", "season", "year", "landuse_broad", "s(terrain,bs='re')",
                "s(host_genus,bs='re')", "s(longitude,latitude,bs='tp',k=10)")
fit_model <- function(d, focus, smooth_k = NULL) {
  other <- setdiff(env_vars, focus)
  focus_term <- if (is.null(smooth_k)) focus else sprintf("s(%s,bs='tp',k=%d)", focus, smooth_k)
  f <- stats::as.formula(paste("shannon ~", paste(c(base_terms, other, focus_term), collapse = "+")))
  mgcv::gam(f, data = d, method = "ML")
}
capture <- function(expr) tryCatch(list(value = force(expr), error = NULL),
                                   error = function(e) list(value = NULL, error = conditionMessage(e)))
clean_error <- function(x) gsub("[\r\n]+", " ", x)
require_core_scalar <- function(value, variable, stage) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value))
    stop("Variable '", variable, "' stage ", stage, " produced a missing, non-finite, or non-scalar statistic.")
  as.numeric(value)
}
normalize_curve_result <- function(result, grid) {
  if (!is.null(result$error) || length(result$value) != length(grid)) {
    if (is.null(result$error)) result$error <- "prediction returned a vector with unexpected length"
    result$value <- rep(NA_real_, length(grid))
  }
  result
}
smooth_stats <- function(model, focus) {
  if (is.null(model)) return(c(edf = NA_real_, p = NA_real_))
  tab <- summary(model)$s.table; rn <- gsub(" ", "", rownames(tab)); key <- paste0("s(", focus, ")")
  j <- match(key, rn)
  if (is.na(j)) return(c(edf = NA_real_, p = NA_real_))
  c(edf = unname(tab[j, "edf"]), p = unname(tab[j, ncol(tab)]))
}
partial_curve <- function(model, focus, grid, reference_data = dat) {
  nd <- reference_data[rep(1, length(grid)), , drop = FALSE]
  for (v in env_vars) nd[[v]] <- stats::median(reference_data[[v]])
  nd[[focus]] <- grid
  pr <- stats::predict(model, newdata = nd, type = "terms")
  target <- if (any(colnames(pr) == focus)) focus else grep(paste0("^s\\(", focus, "\\)"), colnames(pr), value = TRUE)[1]
  if (is.na(target) || !length(target)) stop("Focal term absent from prediction matrix")
  z <- as.numeric(pr[, target]); z - mean(z, na.rm = TRUE)
}
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || stats::sd(x[ok]) < 1e-10 || stats::sd(y[ok]) < 1e-10) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok]))
}
pass <- function(x) isTRUE(is.finite(x) && x)

AIC_THRESHOLD <- -2; EDF_THRESHOLD <- 1.5; P_THRESHOLD <- 0.05
K_P_THRESHOLD <- 0.05; EDF_MARGIN_THRESHOLD <- 0.5
CURVE_COR_THRESHOLD <- 0.98; CURVE_DIFF_THRESHOLD <- 0.20
results <- vector("list", length(env_vars)); plot_data <- vector("list", length(env_vars))
for (i in seq_along(env_vars)) {
  focus <- env_vars[i]; x <- dat[[focus]]; unique_n <- length(unique(x))
  failures <- character()
  kvals <- if (unique_n < 8) c(4L, 5L, 6L) else c(4L, 6L, 8L); primary_k <- 6L
  linear_result <- capture(fit_model(dat, focus))
  if (!is.null(linear_result$error)) stop("Variable '", env_labels[i], "' stage primary_linear_fit failed: ", linear_result$error)
  linear_fit <- linear_result$value
  smooth_results <- setNames(lapply(kvals, function(k) capture(fit_model(dat, focus, k))), kvals)
  primary_result <- smooth_results[[as.character(primary_k)]]
  if (!is.null(primary_result$error)) stop("Variable '", env_labels[i], "' stage primary_k6_smooth_fit failed: ", primary_result$error)
  for (k in kvals[kvals != primary_k]) if (!is.null(smooth_results[[as.character(k)]]$error))
    failures <- c(failures, paste0("smooth_k", k, "_fit:", clean_error(smooth_results[[as.character(k)]]$error)))
  smooth_fits <- lapply(smooth_results, `[[`, "value")
  primary_fit <- smooth_fits[[as.character(primary_k)]]
  ss <- smooth_stats(primary_fit, focus)
  primary_edf <- require_core_scalar(unname(ss["edf"]), env_labels[i], "primary_smooth_edf")
  primary_p <- require_core_scalar(unname(ss["p"]), env_labels[i], "primary_smooth_p_value")
  linear_aic_result <- capture(AIC(linear_fit)); smooth_aic_result <- capture(AIC(primary_fit))
  if (!is.null(linear_aic_result$error)) stop("Variable '", env_labels[i], "' stage primary_linear_AIC failed: ", linear_aic_result$error)
  if (!is.null(smooth_aic_result$error)) stop("Variable '", env_labels[i], "' stage primary_k6_smooth_AIC failed: ", smooth_aic_result$error)
  linear_aic <- require_core_scalar(linear_aic_result$value, env_labels[i], "primary_linear_AIC")
  smooth_aic <- require_core_scalar(smooth_aic_result$value, env_labels[i], "primary_k6_smooth_AIC")
  delta_aic <- smooth_aic - linear_aic
  set.seed(20260712 + i)
  kc_result <- capture(mgcv::k.check(primary_fit))
  kc <- kc_result$value
  if (!is.null(kc_result$error)) failures <- c(failures, paste0("k_check:", clean_error(kc_result$error)))
  kr <- if (is.null(kc)) integer() else grep(paste0("s\\(", focus, "\\)"), rownames(kc))
  if (is.null(kc_result$error) && !length(kr)) failures <- c(failures, "k_check:focal smooth absent from k.check result")
  k_index <- if (length(kr)) unname(kc[kr[1], "k-index"]) else NA_real_
  k_p <- if (length(kr)) unname(kc[kr[1], "p-value"]) else NA_real_
  max_edf <- primary_k - 1; edf_margin <- max_edf - primary_edf
  qs <- stats::quantile(x, c(.025, .975), names = FALSE, type = 7)
  grid <- seq(qs[1], qs[2], length.out = 100)
  curve_results <- lapply(names(smooth_fits), function(k) {
    if (is.null(smooth_fits[[k]])) return(list(value = rep(NA_real_, length(grid)), error = "smooth fit unavailable"))
    normalize_curve_result(capture(partial_curve(smooth_fits[[k]], focus, grid)), grid)
  })
  for (j in seq_along(curve_results)) if (!is.null(curve_results[[j]]$error))
    failures <- c(failures, paste0("smooth_k", names(smooth_fits)[j], "_prediction:", clean_error(curve_results[[j]]$error)))
  curves <- do.call(cbind, lapply(curve_results, `[[`, "value"))
  cors <- if (ncol(curves) >= 2) stats::cor(curves, use = "pairwise.complete.obs") else matrix(NA_real_, 1, 1)
  min_cor <- if (ncol(curves) >= 2 && all(is.finite(cors[upper.tri(cors)]))) min(cors[upper.tri(cors)]) else NA_real_
  primary_curve <- curves[, match(primary_k, kvals)]
  curve_range <- diff(range(primary_curve, finite = TRUE))
  max_rel_diff <- if (is.finite(curve_range) && curve_range > 1e-10 && all(is.finite(curves))) {
    max(sapply(seq_len(ncol(curves)), function(j) max(abs(curves[, j] - primary_curve)))) / curve_range
  } else NA_real_
  trim <- droplevels(dat[x >= qs[1] & x <= qs[2], , drop = FALSE])
  trim_linear_result <- capture(fit_model(trim, focus)); trim_smooth_result <- capture(fit_model(trim, focus, primary_k))
  if (!is.null(trim_linear_result$error)) failures <- c(failures, paste0("trim_linear_fit:", clean_error(trim_linear_result$error)))
  if (!is.null(trim_smooth_result$error)) failures <- c(failures, paste0("trim_smooth_fit:", clean_error(trim_smooth_result$error)))
  trim_linear <- trim_linear_result$value; trim_smooth <- trim_smooth_result$value
  trim_ss <- smooth_stats(trim_smooth, focus)
  if (!is.null(trim_smooth) && (length(trim_ss["edf"]) != 1L || !is.finite(trim_ss["edf"])))
    failures <- c(failures, "trim_smooth_edf:missing or non-finite statistic")
  if (!is.null(trim_smooth) && (length(trim_ss["p"]) != 1L || !is.finite(trim_ss["p"])))
    failures <- c(failures, "trim_smooth_p_value:missing or non-finite statistic")
  trim_aic_result <- if (is.null(trim_linear) || is.null(trim_smooth)) list(value = NA_real_, error = NULL) else
    capture(AIC(trim_smooth) - AIC(trim_linear))
  if (!is.null(trim_aic_result$error)) failures <- c(failures, paste0("trim_AIC:", clean_error(trim_aic_result$error)))
  trim_delta <- trim_aic_result$value
  if (length(trim_delta) != 1L || !is.finite(trim_delta)) failures <- c(failures, "trim_delta_AIC:missing or non-finite statistic")
  trim_curve_result <- if (is.null(trim_smooth)) list(value = rep(NA_real_, length(grid)), error = "trim smooth fit unavailable") else
    normalize_curve_result(capture(partial_curve(trim_smooth, focus, grid, reference_data = trim)), grid)
  if (!is.null(trim_curve_result$error)) failures <- c(failures, paste0("trim_smooth_prediction:", clean_error(trim_curve_result$error)))
  trim_curve <- trim_curve_result$value
  trim_cor <- safe_cor(primary_curve, trim_curve)
  criteria <- c(delta_aic <= AIC_THRESHOLD, primary_edf > EDF_THRESHOLD, primary_p < P_THRESHOLD,
                k_p >= K_P_THRESHOLD, edf_margin >= EDF_MARGIN_THRESHOLD,
                min_cor > CURVE_COR_THRESHOLD && max_rel_diff < CURVE_DIFF_THRESHOLD,
                trim_delta <= AIC_THRESHOLD && trim_ss["edf"] > EDF_THRESHOLD && trim_ss["p"] < P_THRESHOLD && trim_cor > CURVE_COR_THRESHOLD)
  criteria <- vapply(criteria, pass, logical(1)); selected <- if (all(criteria)) "nonlinear" else "linear"
  failed <- which(!criteria); reason <- if (!length(failed)) "All seven nonlinear-form criteria passed." else paste0("Linear retained; failed C", paste(failed, collapse = ", C"), ".")
  if (length(failures)) reason <- paste(reason, paste0("computation_failure:", failures, collapse = "; "))
  results[[i]] <- data.frame(variable = env_labels[i], n = nrow(dat), unique_n = unique_n,
    k_values_tested = paste(kvals, collapse = "/"), primary_k = primary_k,
    linear_AIC = linear_aic, smooth_AIC = smooth_aic, delta_AIC = delta_aic,
    smooth_edf = primary_edf, smooth_p_value = primary_p, k_index = k_index, k_check_p_value = k_p,
    maximum_available_edf = max_edf, edf_ceiling_margin = edf_margin,
    minimum_curve_correlation = min_cor, maximum_relative_curve_difference = max_rel_diff,
    trimmed_n = nrow(trim), trimmed_delta_AIC = trim_delta, trimmed_edf = trim_ss["edf"],
    trimmed_p_value = trim_ss["p"], trimmed_curve_correlation = trim_cor,
    criterion_1_aic = criteria[1], criterion_2_edf = criteria[2], criterion_3_p = criteria[3],
    criterion_4_k_check = criteria[4], criterion_5_edf_ceiling = criteria[5],
    criterion_6_curve_stability = criteria[6], criterion_7_boundary_robustness = criteria[7],
    selected_form = selected, selected_k = if (selected == "nonlinear") primary_k else NA_integer_,
    decision_reason = reason, check.names = FALSE)
  linear_curve_result <- normalize_curve_result(capture(partial_curve(linear_fit, focus, grid)), grid)
  if (!is.null(linear_curve_result$error)) {
    out_reason <- paste0("computation_failure:linear_prediction:", clean_error(linear_curve_result$error))
    results[[i]]$decision_reason <- paste(results[[i]]$decision_reason, out_reason)
  }
  plot_data[[i]] <- list(x = x, grid = grid, linear = linear_curve_result$value, curves = curves,
                         kvals = kvals, criteria = criteria, selected = selected)
  message(sprintf("Assessed %d/14: %s", i, env_labels[i]))
}
out <- do.call(rbind, results)
main_csv <- file.path(table_dir, "environmental_predictor_model_form_assessment.csv")
utils::write.csv(out, main_csv, row.names = FALSE, na = "")
criterion_cols <- grep("^criterion_", names(out), value = TRUE)
summary_table <- rbind(
  data.frame(summary_type = "criterion", category = paste0("C", 1:7),
             count = vapply(out[criterion_cols], sum, integer(1)), total = nrow(out)),
  data.frame(summary_type = "selected_form", category = c("linear", "nonlinear"),
             count = c(sum(out$selected_form == "linear"), sum(out$selected_form == "nonlinear")), total = nrow(out)))
utils::write.csv(summary_table, file.path(table_dir, "environmental_predictor_model_form_criteria_summary.csv"), row.names = FALSE)

write_diagnostics_pdf <- function(path, plot_data, env_labels) {
  grDevices::pdf(path, width = 8.5, height = 7, family = "serif")
  device <- grDevices::dev.cur()
  on.exit(if (!is.null(grDevices::dev.list()) && device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
  cols <- c("#1B9E77", "#D95F02", "#7570B3")
  for (i in seq_along(plot_data)) {
  z <- plot_data[[i]]; yr <- range(c(z$linear, z$curves), finite = TRUE)
  if (!all(is.finite(yr)) || diff(yr) < 1e-10) yr <- c(-1, 1)
  graphics::plot(z$grid, z$linear, type = "n", ylim = yr, xlim = range(z$x, finite = TRUE),
    xlab = env_labels[i], ylab = "Centered partial effect",
    main = paste0(env_labels[i], "\n", paste0("C", 1:7, "=", ifelse(z$criteria, "PASS", "FAIL"), collapse = " | "),
                  " | decision: ", z$selected))
  if (any(is.finite(z$linear))) graphics::lines(z$grid, z$linear, lwd = 2, lty = 2)
  for (j in seq_along(z$kvals)) if (any(is.finite(z$curves[, j])))
    graphics::lines(z$grid, z$curves[, j], col = cols[j], lwd = 2)
  graphics::rug(z$x, col = grDevices::adjustcolor("black", .35))
  graphics::abline(v = stats::quantile(z$x, c(.025, .975)), lty = 3, col = "grey45")
  graphics::legend("topleft", legend = c("linear", paste0("smooth k=", z$kvals)),
                   col = c("black", cols), lty = c(2, 1, 1, 1), lwd = 2, bty = "n")
  }
}
write_diagnostics_pdf(file.path(figure_dir, "environmental_predictor_model_form_diagnostics.pdf"), plot_data, env_labels)
message("Common complete cases: ", nrow(dat), "; outputs written to ", normalizePath(file.path(workflow_dir, "output")))
