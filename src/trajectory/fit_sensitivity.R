#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(gbmt)
})

root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
tables <- file.path(root, "results", "trajectory", "tables")
models <- file.path(root, "results", "trajectory", "models", "sensitivity", "gbmtm")
reports <- file.path(root, "results", "trajectory", "reports")
logs <- file.path(root, "results", "trajectory", "logs")
for (directory in c(tables, models, reports, logs)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

classes <- 1:6
z_vars <- c("log2_easix_interp_z", "log2_siri_interp_z")
model_data <- as.data.frame(fread(file.path(tables, "all_icu_gbmtm_prepared_training.csv")))

specifications <- rbindlist(lapply(1:4, function(degree) {
  data.table(
    spec_id = c(
      paste0("level_d", degree, "_pruning"),
      paste0("level_d", degree, "_no_pruning")
    ),
    degree = degree,
    scaling = 0L,
    pruning = c(TRUE, FALSE),
    description = c(
      paste0("d=", degree, " with automatic pruning of nonsignificant highest-order terms"),
      paste0("d=", degree, " retaining all degree 0-", degree, " time terms")
    ),
    description_en = c(
      paste0("Degree ", degree, " with automatic pruning"),
      paste0("Degree ", degree, " without automatic pruning")
    )
  )
}))

requested <- trimws(Sys.getenv("GBMTM_SPEC_IDS", unset = ""))
if (nzchar(requested)) {
  keep <- trimws(strsplit(requested, ",", fixed = TRUE)[[1]])
  specifications <- specifications[spec_id %chin% keep]
}
if (!nrow(specifications)) stop("No GBMTM sensitivity specifications selected.")

legacy_map <- c(
  level_d1_pruning = "level_linear",
  level_d2_pruning = "level_quadratic",
  level_d3_pruning = "level_cubic",
  level_d2_no_pruning = "level_quadratic_no_pruning"
)

calc_entropy <- function(post) {
  post <- as.matrix(post)
  if (ncol(post) <= 1) return(NA_real_)
  p <- pmax(post, .Machine$double.eps)
  1 + sum(p * log(p)) / (nrow(p) * log(ncol(p)))
}

profile_rms <- function(assignments) {
  k <- length(unique(assignments))
  if (k <= 1) return(NA_real_)
  id_map <- data.table(id = unique(model_data$id), trajectory_class = as.integer(assignments))
  d <- merge(as.data.table(model_data), id_map, by = "id")
  daily <- d[, .(
    easix = mean(log2_easix_interp_z),
    siri = mean(log2_siri_interp_z)
  ), by = .(trajectory_class, time)]
  wide <- dcast(daily, trajectory_class ~ time, value.var = c("easix", "siri"))
  distance_matrix <- as.matrix(dist(as.matrix(wide[, -1]))) / sqrt(ncol(wide) - 1)
  min(distance_matrix[upper.tri(distance_matrix)])
}

fit_one <- function(spec, k) {
  path <- file.path(models, paste0("gbmtm_", spec$spec_id, "_k", k, ".rds"))
  if (file.exists(path)) return(readRDS(path))

  legacy_id <- unname(legacy_map[spec$spec_id])
  if (length(legacy_id) && !is.na(legacy_id)) {
    legacy_path <- file.path(models, paste0("gbmtm_", legacy_id, "_k", k, ".rds"))
    if (file.exists(legacy_path)) {
      fit <- readRDS(legacy_path)
      saveRDS(fit, path)
      return(fit)
    }
  }

  set.seed(20260729 + 100L * spec$degree + 10L * as.integer(!spec$pruning) + k)
  fit <- gbmt(
    x.names = z_vars, unit = "id", time = "time", ng = k,
    d = spec$degree, data = model_data, scaling = 0L,
    pruning = spec$pruning, tol = 1e-4, maxit = 1000, quiet = TRUE
  )
  saveRDS(fit, path)
  fit
}

all_rows <- list()
for (i in seq_len(nrow(specifications))) {
  spec <- specifications[i]
  message("GBMTM specification: ", spec$spec_id)
  rows <- list()
  for (k in classes) {
    message("  fitting K=", k)
    fit <- tryCatch(fit_one(spec, k), error = function(e) e)
    if (inherits(fit, "error")) {
      rows[[as.character(k)]] <- data.table(
        spec_id = spec$spec_id, k = k, status = "fit_failed", error = conditionMessage(fit)
      )
      next
    }
    props <- as.numeric(prop.table(table(factor(fit$assign, levels = seq_len(k)))))
    avepp <- if (k == 1) 1 else fit$appa
    converged <- if (!is.null(fit$em) && ncol(fit$em)) {
      any(fit$em[, ncol(fit$em)] == 1)
    } else {
      NA
    }
    rows[[as.character(k)]] <- data.table(
      spec_id = spec$spec_id,
      k = k,
      AIC = unname(fit$ic["aic"]),
      BIC = unname(fit$ic["bic"]),
      entropy = calc_entropy(posterior(fit)),
      min_avepp = min(avepp, na.rm = TRUE),
      mean_avepp = mean(avepp, na.rm = TRUE),
      min_class_pct = min(props) * 100,
      min_pairwise_profile_rms = profile_rms(fit$assign),
      converged = converged,
      class_proportions = paste0(
        "C", seq_along(props), "=", sprintf("%.1f%%", props * 100), collapse = "; "
      ),
      status = "fit_ok",
      error = ""
    )
  }
  result <- rbindlist(rows, fill = TRUE)
  result <- merge(result, spec, by = "spec_id", all.x = TRUE)
  fwrite(
    result,
    file.path(tables, paste0("all_icu_gbmtm_parameter_sensitivity_", spec$spec_id, ".csv"))
  )
  all_rows[[spec$spec_id]] <- result
}

# Build the active grid only from this run so archived scaling analyses cannot leak back in.
candidate_dt <- rbindlist(all_rows, fill = TRUE, use.names = TRUE)
setorder(candidate_dt, degree, -pruning, k)
candidate_dt[, diagnostic_eligible_10pct :=
  status == "fit_ok" & k >= 2 & entropy >= 0.80 & min_class_pct >= 10 & min_avepp >= 0.70]
candidate_dt[, profile_separated_10pct :=
  diagnostic_eligible_10pct & min_pairwise_profile_rms >= 0.50]
fwrite(candidate_dt, file.path(tables, "all_icu_gbmtm_parameter_sensitivity_all_candidates.csv"))

summaries <- candidate_dt[status == "fit_ok", {
  non_null <- .SD[k >= 2]
  diagnostic <- non_null[diagnostic_eligible_10pct %in% TRUE]
  separated <- non_null[profile_separated_10pct %in% TRUE]
  list(
    diagnostic_eligible_k_10pct = paste(diagnostic$k, collapse = ","),
    profile_separated_k_10pct = paste(separated$k, collapse = ","),
    candidate_fit_strategy = paste0(
      "Single fit per K/specification; K=1 is the no-group benchmark; ",
      "all active models use scaling=0"
    )
  )
}, by = .(spec_id, degree, scaling, pruning, description, description_en)]
setorder(summaries, degree, -pruning)
fwrite(summaries, file.path(tables, "all_icu_gbmtm_parameter_sensitivity_summary.csv"))

report <- c(
  "# GBMTM degree and pruning sensitivity",
  "",
  "Candidate K=1-6 includes K=1 as the ungrouped reference. Compare degree 1-4 with and without pruning (eight specifications), all with scaling=0.",
  "",
  "Diagnostic thresholds: minimum class proportion >=10%, entropy >=0.80, minimum class AvePP >=0.70; minimum pairwise standardized trajectory RMS >=0.50 assesses separation without outcomes.",
  "",
  paste(capture.output(print(summaries)), collapse = "\n")
)
writeLines(report, file.path(reports, "all_icu_gbmtm_parameter_sensitivity_report.md"))
