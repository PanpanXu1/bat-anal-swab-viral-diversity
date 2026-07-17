options(stringsAsFactors = FALSE)
pdf_family <- "Times"

if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Required package 'mgcv' is unavailable.")
}

write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
abbr <- function(x) {
  out <- sub("^.*\\(([^()]*)\\)$", "\\1", x)
  ifelse(out == x, x, out)
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
workflow_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
module_dir <- normalizePath(file.path(workflow_dir, ".."), mustWork = FALSE)
input_csv <- file.path(
  workflow_dir,
  "input",
  "environmental_variables_selected_by_corr0.8_VIF10.csv"
)
predictor_csv <- file.path(
  workflow_dir,
  "input",
  "parsimonious_environmental_predictors.csv"
)
form_assessment_csv <- file.path(
  workflow_dir,
  "input",
  "environmental_predictor_model_form_assessment.csv"
)
output_dir <- file.path(workflow_dir, "output")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

raw <- utils::read.csv(input_csv, check.names = FALSE)
predictors <- utils::read.csv(predictor_csv, check.names = FALSE, colClasses = "character", na.strings = "NA")
assessment <- utils::read.csv(form_assessment_csv, check.names = FALSE, colClasses = "character", na.strings = "NA")
criterion_cols <- paste0("criterion_", 1:7, c("_aic", "_edf", "_p", "_k_check", "_edf_ceiling", "_curve_stability", "_boundary_robustness"))
screened_environmental_predictors <- c(
  "Mean Diurnal Range (Mean of monthly max temp - min temp) (Bio2)",
  "Temperature Seasonality (standard deviation x100) (Bio4)",
  "Max Temperature of Warmest Month (Bio5)",
  "Mean Temperature of Wettest Quarter (Bio8)",
  "Precipitation of Driest Quarter (Bio17)",
  "Precipitation of Warmest Quarter (Bio18)",
  "Human Footprint (HFT)",
  "Global Mammal Richness (GMR)",
  "China High-Resolution Ecological Environment Quality (CHEQ)",
  "Normalized Difference Vegetation Index (NDVI)",
  "Fractional Vegetation Cover (FVC)",
  "China Population Spatial Distribution (PSD)",
  "Global Railway (GR)",
  "Global Linear Hydrography (GLH)"
)
required_predictor_cols <- c("variable", "selected_form", "selected_k", "model_form_source")
required_assessment_cols <- c("variable", "selected_form", "selected_k", criterion_cols)
validate_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
}
validate_columns(predictors, required_predictor_cols, "Parsimonious predictor table")
validate_columns(assessment, required_assessment_cols, "Model-form assessment table")
if (!nrow(predictors)) stop("Parsimonious predictor table contains no predictors.")
if (anyNA(predictors$variable) || any(!nzchar(trimws(predictors$variable)))) stop("Parsimonious predictor table has missing variable names.")
if (anyDuplicated(predictors$variable)) stop("Parsimonious predictor table has duplicate variables: ", paste(unique(predictors$variable[duplicated(predictors$variable)]), collapse = ", "))
if (anyNA(assessment$variable) || any(!nzchar(trimws(assessment$variable)))) stop("Model-form assessment table has missing variable names.")
if (anyDuplicated(assessment$variable)) stop("Model-form assessment table has duplicate variables: ", paste(unique(assessment$variable[duplicated(assessment$variable)]), collapse = ", "))
missing_screened_assessments <- setdiff(screened_environmental_predictors, assessment$variable)
unknown_assessments <- setdiff(assessment$variable, screened_environmental_predictors)
if (length(missing_screened_assessments) || length(unknown_assessments) || nrow(assessment) != length(screened_environmental_predictors)) {
  details <- c(
    if (length(missing_screened_assessments)) paste0("missing: ", paste(missing_screened_assessments, collapse = ", ")),
    if (length(unknown_assessments)) paste0("unknown: ", paste(unknown_assessments, collapse = ", "))
  )
  stop(
    "Model-form assessment table must contain exactly one row for each of the 14 screened environmental predictors",
    if (length(details)) paste0(" (", paste(details, collapse = "; "), ")") else "",
    "."
  )
}
missing_assessments <- setdiff(predictors$variable, assessment$variable)
if (length(missing_assessments)) stop("Parsimonious predictors are missing from model-form assessment: ", paste(missing_assessments, collapse = ", "))

parse_bool <- function(x, column) {
  if (anyNA(x) || any(!x %in% c("TRUE", "FALSE"))) stop("Model-form assessment has invalid boolean values in ", column, "; only exact TRUE/FALSE are allowed.")
  x == "TRUE"
}
for (nm in criterion_cols) assessment[[nm]] <- parse_bool(assessment[[nm]], nm)
if (anyNA(predictors$selected_form) || any(!predictors$selected_form %in% c("linear", "nonlinear"))) stop("Parsimonious predictor table has invalid selected_form; only exact linear/nonlinear are allowed.")
if (anyNA(assessment$selected_form) || any(!assessment$selected_form %in% c("linear", "nonlinear"))) stop("Model-form assessment table has invalid selected_form; only exact linear/nonlinear are allowed.")
all_seven_true <- apply(assessment[, criterion_cols, drop = FALSE], 1, all)
all_true_but_linear <- assessment$variable[all_seven_true & assessment$selected_form != "nonlinear"]
if (length(all_true_but_linear)) stop("Model-form assessment all-seven inconsistency: C1-C7 are all TRUE but selected_form is not nonlinear for: ", paste(all_true_but_linear, collapse = ", "))
nonlinear_but_failed <- assessment$variable[!all_seven_true & assessment$selected_form == "nonlinear"]
if (length(nonlinear_but_failed)) stop("Model-form assessment all-seven inconsistency: selected_form is nonlinear but C1-C7 are not all TRUE for: ", paste(nonlinear_but_failed, collapse = ", "))
expected_source <- "workflow_04_seven_criterion_assessment"
if (any(is.na(predictors$model_form_source) | predictors$model_form_source != expected_source)) stop("Parsimonious predictor model_form_source must be '", expected_source, "'.")
selected_assessment <- assessment[match(predictors$variable, assessment$variable), , drop = FALSE]
parse_selected_k <- function(raw_k, selected_form, label) {
  is_empty <- is.na(raw_k) | raw_k == ""
  if (any(selected_form == "linear" & !is_empty)) stop(label, ": linear predictors must have truly empty/NA selected_k values.")
  nonlinear_idx <- which(selected_form == "nonlinear")
  if (any(is_empty[nonlinear_idx])) stop(label, ": nonlinear predictors require non-empty selected_k values.")
  numeric_pattern <- "^[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$"
  if (length(nonlinear_idx) && any(!grepl(numeric_pattern, raw_k[nonlinear_idx], perl = TRUE))) stop(label, ": nonlinear selected_k must be strictly numeric without surrounding whitespace.")
  out <- rep(NA_real_, length(raw_k))
  out[nonlinear_idx] <- suppressWarnings(as.numeric(raw_k[nonlinear_idx]))
  if (length(nonlinear_idx) && any(!is.finite(out[nonlinear_idx]) | out[nonlinear_idx] < 4 | out[nonlinear_idx] != floor(out[nonlinear_idx]))) stop(label, ": nonlinear selected_k must be a finite integer >= 4.")
  out
}
predictor_k <- parse_selected_k(predictors$selected_k, predictors$selected_form, "Parsimonious predictor table")
assessment_all_k <- parse_selected_k(assessment$selected_k, assessment$selected_form, "Model-form assessment table")
assessment_k <- assessment_all_k[match(predictors$variable, assessment$variable)]
nonlinear <- predictors$selected_form == "nonlinear"
same_k <- (is.na(predictor_k) & is.na(assessment_k)) | (!is.na(predictor_k) & !is.na(assessment_k) & predictor_k == assessment_k)
if (any(predictors$selected_form != selected_assessment$selected_form | !same_k)) stop("Parsimonious predictor and model-form assessment selected_form/selected_k values are inconsistent.")
if (any(nonlinear & !apply(selected_assessment[, criterion_cols, drop = FALSE], 1, all))) stop("Nonlinear predictors require all seven assessment criteria (C1-C7) to be TRUE.")

env_vars <- paste0("env_", seq_len(nrow(predictors)))
environmental_mapping <- data.frame(
  internal = env_vars,
  label = predictors$variable,
  selected_form = predictors$selected_form,
  selected_k = predictor_k,
  stringsAsFactors = FALSE
)
environmental_mapping$formula_term <- ifelse(
  environmental_mapping$selected_form == "linear",
  environmental_mapping$internal,
  sprintf("s(%s, bs='tp', k=%d)", environmental_mapping$internal, environmental_mapping$selected_k)
)
environmental_mapping$summary_term <- ifelse(environmental_mapping$selected_form == "linear", environmental_mapping$internal, paste0("s(", environmental_mapping$internal, ")"))
env_labels <- environmental_mapping$label
env_label_by_summary <- stats::setNames(environmental_mapping$label, environmental_mapping$summary_term)
form_check <- data.frame(
  variable = environmental_mapping$label,
  selected_assessment[, criterion_cols, drop = FALSE],
  upstream_selected_form = environmental_mapping$selected_form,
  upstream_selected_k = environmental_mapping$selected_k,
  final_formula_term = environmental_mapping$formula_term,
  final_model_form = ifelse(nonlinear, "environmental nonlinear smooth", "environmental linear fixed effect"),
  action = ifelse(nonlinear, "validated and implemented as nonlinear smooth", "validated and implemented as linear fixed effect"),
  check.names = FALSE
)
names(form_check)[2:8] <- paste0("C", 1:7)
write_csv(form_check, file.path(table_dir, "environmental_predictor_model_form_check.csv"))
needed <- c(
  "Shannon index", "SampleSize (log10)", "Season", "Year", "Terrain", "Host genus",
  "longitude", "latitude", env_labels
)
missing_needed <- setdiff(needed, names(raw))
if (length(missing_needed)) stop("Missing required columns: ", paste(missing_needed, collapse = ", "))

strict_numeric <- function(x, column) {
  original_nonmissing <- !is.na(x)
  converted <- suppressWarnings(as.numeric(x))
  invalid <- original_nonmissing & (is.na(converted) | !is.finite(converted))
  if (any(invalid)) stop("Column '", column, "' contains non-numeric or non-finite values.")
  converted
}
validate_nonempty_factor <- function(x, column) {
  invalid <- !is.na(x) & !nzchar(trimws(as.character(x)))
  if (any(invalid)) stop("Column '", column, "' contains empty non-missing factor values.")
  factor(x)
}

raw$shannon <- strict_numeric(raw[["Shannon index"]], "Shannon index")
raw$log10_pool_size <- strict_numeric(raw[["SampleSize (log10)"]], "SampleSize (log10)")
raw$season <- validate_nonempty_factor(raw[["Season"]], "Season")
raw$year <- validate_nonempty_factor(raw[["Year"]], "Year")
raw$terrain <- validate_nonempty_factor(raw[["Terrain"]], "Terrain")
raw$host_genus <- validate_nonempty_factor(raw[["Host genus"]], "Host genus")
raw$longitude <- strict_numeric(raw[["longitude"]], "longitude")
raw$latitude <- strict_numeric(raw[["latitude"]], "latitude")
for (i in seq_along(env_labels)) raw[[env_vars[i]]] <- strict_numeric(raw[[env_labels[i]]], env_labels[i])

landuse_mapping <- data.frame(
  original_landuse = c(
    "landuseName Forest land closed canopy", "landuseName Other forest land",
    "landuseName Sparse woodland", "landuseName Shrubland",
    "landuseName High coverage grassland", "landuseName Medium coverage grassland",
    "landuseName Low coverage grassland", "landuseName Paddy field",
    "landuseName Other construction land", "landuseName Urban built up area",
    "landuseName Rural residential area", "landuseName Reservoirs and ponds",
    "landuseName Rivers and canals"
  ),
  landuse_broad = c(
    "Closed_or_other_forest", "Closed_or_other_forest", "Sparse_woodland",
    "Open_vegetation", "Open_vegetation", "Open_vegetation", "Open_vegetation",
    "Cropland", "Built_up", "Built_up", "Built_up", "Other_or_unclassified", "Other_or_unclassified"
  )
)
missing_landuse <- setdiff(landuse_mapping$original_landuse, names(raw))
if (length(missing_landuse)) stop("Missing required land-use dummy columns: ", paste(missing_landuse, collapse = ", "))
land_mat <- vapply(
  landuse_mapping$original_landuse,
  function(v) strict_numeric(raw[[v]], v),
  numeric(nrow(raw))
)
colnames(land_mat) <- landuse_mapping$original_landuse
if (anyNA(land_mat) || any(!is.finite(land_mat))) stop("Land-use dummy columns must not contain missing or non-finite values.")
if (any(!land_mat %in% c(0, 1))) stop("Land-use dummy columns must contain only 0 or 1.")
dummy_sum <- rowSums(land_mat)
landuse_broad <- rep("Other_or_unclassified", nrow(raw))
singles <- which(dummy_sum == 1)
if (length(singles)) {
  hit <- apply(land_mat[singles, , drop = FALSE], 1, function(x) colnames(land_mat)[which(x == 1)[1]])
  landuse_broad[singles] <- landuse_mapping$landuse_broad[match(hit, landuse_mapping$original_landuse)]
}
multiples <- which(dummy_sum > 1)
if (length(multiples) >= 5) landuse_broad[multiples] <- "Multiple_landuse"
if (length(multiples) > 0 && length(multiples) < 5) landuse_broad[multiples] <- "Other_or_unclassified"
raw$landuse_broad <- factor(landuse_broad)

model_vars <- c(
  "shannon", "log10_pool_size", "season", "year", "terrain", "host_genus",
  "longitude", "latitude", "landuse_broad", env_vars
)
dat <- droplevels(raw[stats::complete.cases(raw[, model_vars]), , drop = FALSE])
cat("Input rows:", nrow(raw), "; complete model rows:", nrow(dat), "; removed:", nrow(raw) - nrow(dat), "\n")
if (!nrow(dat)) stop("No complete rows remain for the final model.")
for (v in c("season", "year", "terrain", "host_genus", "landuse_broad")) {
  if (nlevels(dat[[v]]) < 2) stop("Complete-case model data require at least two levels for '", v, "'.")
}
spatial_k <- 10L
if (nrow(unique(dat[, c("longitude", "latitude")])) < spatial_k) {
  stop("Complete-case model data require at least ", spatial_k, " unique longitude/latitude combinations for the spatial smooth.")
}
parametric_terms <- c("log10_pool_size", "season", "factor(year)", "landuse_broad",
                      environmental_mapping$internal[environmental_mapping$selected_form == "linear"])
parametric_formula <- stats::as.formula(paste("~", paste(parametric_terms, collapse = " + ")))
parametric_rank <- qr(stats::model.matrix(parametric_formula, data = dat))$rank
minimum_complete_n <- max(50L, parametric_rank + 11L)
if (nrow(dat) < minimum_complete_n) {
  stop("Insufficient complete-case rows: need at least ", minimum_complete_n,
       " (maximum of 50 and parametric design rank + 11), found ", nrow(dat), ".")
}

spatial_term <- sprintf("s(longitude, latitude, bs='tp', k=%d)", spatial_k)
final_terms <- c(
  "log10_pool_size", "season", "factor(year)", "landuse_broad", environmental_mapping$formula_term,
  "s(terrain, bs='re')", "s(host_genus, bs='re')", spatial_term
)
final_formula <- paste("shannon ~", paste(final_terms, collapse = " + "))
fit_gam <- function(formula_text, method) mgcv::gam(stats::as.formula(formula_text), data = dat, method = method)
final_model <- fit_gam(final_formula, method = "REML")
final_summary <- summary(final_model)

extract_results <- function(model) {
  su <- summary(model)
  p <- as.data.frame(su$p.table, check.names = FALSE)
  p$raw_term <- rownames(p)
  p$term <- ifelse(p$raw_term %in% names(env_label_by_summary), env_label_by_summary[p$raw_term], p$raw_term)
  rownames(p) <- NULL
  names(p)[1:4] <- c("estimate", "std_error", "statistic", "p_value")
  p$term_type <- "parametric"
  p$edf <- NA_real_
  p <- p[, c("raw_term", "term", "term_type", "estimate", "std_error", "edf", "statistic", "p_value")]
  s <- as.data.frame(su$s.table, check.names = FALSE)
  names(s)[1] <- "edf"
  names(s)[ncol(s) - 1] <- "statistic"
  names(s)[ncol(s)] <- "p_value"
  s$raw_term <- rownames(s)
  s$term <- ifelse(s$raw_term %in% names(env_label_by_summary), env_label_by_summary[s$raw_term], s$raw_term)
  rownames(s) <- NULL
  s$term_type <- "smooth"
  s$estimate <- NA_real_
  s$std_error <- NA_real_
  s <- s[, c("raw_term", "term", "term_type", "estimate", "std_error", "edf", "statistic", "p_value")]
  rbind(p, s)
}

term_results <- extract_results(final_model)
env_results <- term_results[term_results$raw_term %in% environmental_mapping$summary_term, , drop = FALSE]
significant_env <- env_results[env_results$p_value < 0.05, , drop = FALSE]
term_role <- function(raw_term) {
  if (raw_term == "(Intercept)") return("intercept")
  if (raw_term == "log10_pool_size") return("pool sample size")
  if (grepl("^season", raw_term)) return("sampling season")
  if (grepl("^factor\\(year\\)", raw_term)) return("sampling year")
  if (grepl("^landuse_broad", raw_term)) return("land-use category")
  if (raw_term %in% environmental_mapping$summary_term) return("environmental predictor")
  if (raw_term == "s(terrain)") return("terrain background")
  if (raw_term == "s(host_genus)") return("host-genus background")
  if (raw_term == "s(longitude,latitude)") return("spatial structure")
  "other"
}
model_form <- function(raw_term, term_type) {
  if (raw_term == "(Intercept)") return("intercept")
  if (raw_term == "log10_pool_size") return("linear fixed effect")
  env_idx <- match(raw_term, environmental_mapping$summary_term)
  if (!is.na(env_idx)) return(if (environmental_mapping$selected_form[env_idx] == "nonlinear") "environmental nonlinear smooth" else "environmental linear fixed effect")
  if (grepl("^season|^factor\\(year\\)|^landuse_broad", raw_term)) return("categorical fixed effect")
  if (raw_term %in% c("s(terrain)", "s(host_genus)")) return("random-effect smooth")
  if (raw_term == "s(longitude,latitude)") return("two-dimensional smooth")
  term_type
}
model_summary <- data.frame(
  term = term_results$term,
  role = vapply(term_results$raw_term, term_role, character(1)),
  model_form = mapply(model_form, term_results$raw_term, term_results$term_type, USE.NAMES = FALSE),
  term_type = term_results$term_type,
  estimate = term_results$estimate,
  std_error = term_results$std_error,
  edf = term_results$edf,
  statistic = term_results$statistic,
  p_value = term_results$p_value
)
write_csv(model_summary, file.path(table_dir, "model_fit_summary.csv"))

metrics <- function(model) c(AIC = AIC(model), adjusted_R2 = summary(model)$r.sq, deviance_explained = summary(model)$dev.expl)
comparison_full_model <- fit_gam(final_formula, method = "ML")
full_metrics <- metrics(comparison_full_model)
term_specs <- data.frame(
  term_group = c("pool_sample_size", "temporal", "temporal", "landuse", rep("environmental", length(env_vars)),
                 "random_effect", "random_effect", "spatial"),
  term_id = c("log10_pool_size", "season", "factor(year)", "landuse_broad", environmental_mapping$summary_term,
              "terrain_random_effect", "host_genus_random_effect", "spatial_smooth"),
  term_name = c("log10_pool_size", "season", "factor(year)", "landuse_broad", env_labels,
                "terrain_random_effect", "host_genus_random_effect", "spatial_smooth"),
  removed_term_text = c("log10_pool_size", "season", "factor(year)", "landuse_broad", environmental_mapping$formula_term,
                        "s(terrain, bs='re')", "s(host_genus, bs='re')", spatial_term)
)

term_p <- function(name) {
  if (name %in% environmental_mapping$summary_term || name == "log10_pool_size") {
    return(term_results$p_value[match(name, term_results$raw_term)])
  }
  if (name == "host_genus_random_effect") return(term_results$p_value[match("s(host_genus)", term_results$term)])
  if (name == "terrain_random_effect") return(term_results$p_value[match("s(terrain)", term_results$term)])
  if (name == "spatial_smooth") return(term_results$p_value[match("s(longitude,latitude)", term_results$term)])
  prefix <- switch(name, season = "^season", `factor(year)` = "^factor\\(year\\)", landuse_broad = "^landuse_broad")
  idx <- grep(prefix, rownames(final_summary$p.table))
  if (length(idx)) min(final_summary$p.table[idx, ncol(final_summary$p.table)]) else NA_real_
}

drop_rows <- lapply(seq_len(nrow(term_specs)), function(i) {
  spec <- term_specs[i, ]
  reduced_terms <- final_terms[seq_along(final_terms) != match(spec$removed_term_text, final_terms)]
  reduced <- fit_gam(paste("shannon ~", paste(reduced_terms, collapse = " + ")), method = "ML")
  reduced_metrics <- metrics(reduced)
  data.frame(
    term_group = spec$term_group,
    term_name = spec$term_name,
    model_comparison_method = "ML",
    delta_deviance_explained = unname(full_metrics["deviance_explained"] - reduced_metrics["deviance_explained"]),
    delta_adjusted_R2 = unname(full_metrics["adjusted_R2"] - reduced_metrics["adjusted_R2"]),
    delta_AIC = unname(reduced_metrics["AIC"] - full_metrics["AIC"]),
    full_model_term_p_value_if_available = term_p(spec$term_id)
  )
})
contribution <- do.call(rbind, drop_rows)
contribution <- contribution[order(-contribution$delta_deviance_explained, -contribution$delta_AIC), ]
contribution$rank <- seq_len(nrow(contribution))
contribution <- contribution[, c("rank", "term_group", "term_name", "model_comparison_method", "delta_deviance_explained",
                                 "delta_adjusted_R2", "delta_AIC", "full_model_term_p_value_if_available")]
write_csv(contribution, file.path(table_dir, "gamm_drop_one_contribution_summary.csv"))

plot_pdf <- function(path, draw, width = 9, height = 6.5) {
  grDevices::pdf(path, width = width, height = height, family = pdf_family)
  old <- graphics::par(family = pdf_family)
  on.exit({ graphics::par(old); grDevices::dev.off() }, add = TRUE)
  draw()
}

terms_matrix <- stats::predict(final_model, type = "terms")
key_plot_idx <- unique(c(
  match(significant_env$raw_term, environmental_mapping$summary_term),
  which(environmental_mapping$selected_form == "nonlinear")
))
key_plot_idx <- key_plot_idx[!is.na(key_plot_idx)]
key_plot_vars <- environmental_mapping$internal[key_plot_idx]
key_plot_labels <- environmental_mapping$label[key_plot_idx]
key_plot_forms <- environmental_mapping$selected_form[key_plot_idx]
key_plot_k <- environmental_mapping$selected_k[key_plot_idx]
if (!length(key_plot_vars)) {
  key_plot_idx <- seq_len(min(2, nrow(environmental_mapping)))
  key_plot_vars <- environmental_mapping$internal[key_plot_idx]
  key_plot_labels <- environmental_mapping$label[key_plot_idx]
  key_plot_forms <- environmental_mapping$selected_form[key_plot_idx]
  key_plot_k <- environmental_mapping$selected_k[key_plot_idx]
}
keep <- !is.na(key_plot_vars)
key_plot_vars <- key_plot_vars[keep]
key_plot_labels <- key_plot_labels[keep]
key_plot_forms <- key_plot_forms[keep]
key_plot_k <- key_plot_k[keep]

plot_observed_association <- function(x, y, xlab, selected_form, selected_k) {
  grid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 200)
  if (selected_form == "nonlinear") {
    fit <- mgcv::gam(y ~ s(x, bs = "tp", k = selected_k), method = "REML")
    pred_raw <- stats::predict(fit, newdata = data.frame(x = grid), se.fit = TRUE)
    pred <- cbind(fit = pred_raw$fit, lwr = pred_raw$fit - 1.96 * pred_raw$se.fit, upr = pred_raw$fit + 1.96 * pred_raw$se.fit)
    cat("Descriptive observed association plot:", xlab, "-> nonlinear smooth (k=", selected_k, ")\n")
  } else {
    fit <- stats::lm(y ~ x)
    pred <- stats::predict(fit, newdata = data.frame(x = grid), interval = "confidence")
    cat("Descriptive observed association plot:", xlab, "-> linear lm\n")
  }
  graphics::plot(x, y, pch = 16, col = "#386CB080", xlab = xlab, ylab = "Shannon index",
                 cex = 0.75, las = 1)
  graphics::grid(col = "#E5E5E5", lwd = 1)
  graphics::polygon(c(grid, rev(grid)), c(pred[, "lwr"], rev(pred[, "upr"])),
                    col = "#F4A34055", border = NA)
  graphics::points(x, y, pch = 16, col = "#386CB080", cex = 0.75)
  graphics::lines(grid, pred[, "fit"], col = "#E31A1C", lwd = 2.5)
  graphics::box()
}

plot_pdf(file.path(figure_dir, "environmental_observed_associations_key_variables.pdf"), function() {
  old <- graphics::par(mfrow = c(1, length(key_plot_vars)), mar = c(5.5, 5.5, 1, 1), cex.lab = 1.25,
                       cex.axis = 1.05)
  on.exit(graphics::par(old), add = TRUE)
  for (i in seq_along(key_plot_vars)) {
    v <- key_plot_vars[i]
    label <- key_plot_labels[i]
    plot_observed_association(dat[[v]], dat$shannon, paste0(label, " (standardized)"), key_plot_forms[i], key_plot_k[i])
  }
}, width = 5.2 * max(1, length(key_plot_vars)), height = 5.8)

plot_pdf(file.path(figure_dir, "ML_drop_one_contribution.pdf"), function() {
  z <- contribution[contribution$term_group == "environmental", , drop = FALSE]
  z$term_abbreviation <- abbr(z$term_name)
  z <- z[order(z$delta_deviance_explained), ]
  values <- 100 * z$delta_deviance_explained
  positive_max <- max(0, values, na.rm = TRUE)
  tick_candidates <- base::pretty(c(0, positive_max * 1.10), n = 5)
  x_ticks <- tick_candidates[tick_candidates >= 0]
  x_upper <- max(x_ticks)
  if (x_upper <= positive_max) {
    tick_step <- if (length(x_ticks) >= 2L) min(diff(x_ticks)) else max(0.5, positive_max * 0.25)
    x_upper <- x_upper + tick_step
    x_ticks <- c(x_ticks, x_upper)
  }
  negative_min <- min(0, values, na.rm = TRUE)
  x_lower <- if (negative_min < 0) min(negative_min * 1.10, -0.03 * x_upper) else 0
  graphics::par(mar = c(5.5, 5.5, 1, 2), cex.lab = 1.25, cex.axis = 1.05)
  graphics::barplot(values, names.arg = z$term_abbreviation, horiz = TRUE, las = 1,
                    col = "#2C7FB8", border = NA, cex.names = 0.95,
                    xlim = c(x_lower, x_upper), xaxt = "n", xaxs = "i",
                    xlab = "Delta deviance explained (ML drop-one, %)")
  graphics::axis(1, at = x_ticks)
  graphics::abline(v = 0, lty = 2, col = "grey40")
}, width = 11, height = 8)

cat("Parsimonious GAMM completed.\n")
cat("Adjusted R2:", final_summary$r.sq, "\n")
cat("Deviance explained (%):", 100 * final_summary$dev.expl, "\n")
cat("Drop-one model comparison method: ML\n")
