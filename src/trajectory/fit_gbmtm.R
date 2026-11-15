#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(gbmt)
})

root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
tables <- file.path(root, "results", "trajectory", "tables")
models <- file.path(root, "results", "trajectory", "models")
reports <- file.path(root, "results", "trajectory", "reports")
figures <- file.path(root, "results", "trajectory", "figures")
logs <- file.path(root, "results", "trajectory", "logs")
for (d in c(tables, models, reports, figures, logs)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

train_file <- file.path(tables, "all_icu_d1d7_daily_mimiciv.csv")
valid_file <- file.path(tables, "all_icu_d1d7_daily_mimiciii.csv")
stopifnot(file.exists(train_file), file.exists(valid_file))

tag <- "all_icu_mimiciv_d1d7_easix_siri"
classes <- 1:6
degree <- 2L
minimum_class_pct <- 10
set.seed(20260727)

log_file <- file.path(logs, paste0("gbmtm_", tag, ".log"))
log_conn <- file(log_file, open = "wt")
on.exit(close(log_conn), add = TRUE)
log_msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  writeLines(x, log_conn)
  flush(log_conn)
}

clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

eligible_ids <- function(d) {
  d[, .(
    n_easix_days = sum(is.finite(easix_day) & easix_day > 0),
    n_siri_days = sum(is.finite(siri_day) & siri_day > 0),
    n_paired_days = sum(
      is.finite(easix_day) & easix_day > 0 & is.finite(siri_day) & siri_day > 0
    )
  ), by = icu_stay_id][n_easix_days >= 2 & n_siri_days >= 2]
}

interpolate_one <- function(day, value) {
  ok <- is.finite(value)
  if (sum(ok) < 2) return(rep(NA_real_, length(day)))
  stats::approx(day[ok], value[ok], xout = day, rule = 2, ties = mean)$y
}

prepare_training <- function(d) {
  ids <- eligible_ids(d)
  d <- d[icu_stay_id %chin% ids$icu_stay_id]
  d[, log2_easix := fifelse(is.finite(easix_day) & easix_day > 0, log2(easix_day), NA_real_)]
  d[, log2_siri := fifelse(is.finite(siri_day) & siri_day > 0, log2(siri_day), NA_real_)]
  bounds <- rbindlist(lapply(c("log2_easix", "log2_siri"), function(v) {
    q <- quantile(d[[v]], c(0.01, 0.99), na.rm = TRUE, names = FALSE)
    data.table(variable = v, winsor_lower = q[1], winsor_upper = q[2])
  }))
  for (v in bounds$variable) {
    p <- bounds[variable == v]
    d[[v]] <- clip(d[[v]], p$winsor_lower, p$winsor_upper)
    d[, paste0(v, "_interp") := interpolate_one(day, get(v)), by = icu_stay_id]
  }
  scaling <- rbindlist(lapply(c("log2_easix_interp", "log2_siri_interp"), function(v) {
    data.table(variable = v, mean = mean(d[[v]]), sd = sd(d[[v]]))
  }))
  for (v in scaling$variable) {
    p <- scaling[variable == v]
    d[[paste0(v, "_z")]] <- (d[[v]] - p$mean) / p$sd
  }
  list(data = d, eligibility = ids, bounds = bounds, scaling = scaling)
}

prepare_validation <- function(d, bounds, scaling) {
  ids <- eligible_ids(d)
  d <- d[icu_stay_id %chin% ids$icu_stay_id]
  d[, log2_easix := fifelse(is.finite(easix_day) & easix_day > 0, log2(easix_day), NA_real_)]
  d[, log2_siri := fifelse(is.finite(siri_day) & siri_day > 0, log2(siri_day), NA_real_)]
  for (v in bounds$variable) {
    p <- bounds[variable == v]
    d[[v]] <- clip(d[[v]], p$winsor_lower, p$winsor_upper)
    d[, paste0(v, "_interp") := interpolate_one(day, get(v)), by = icu_stay_id]
  }
  for (v in scaling$variable) {
    p <- scaling[variable == v]
    d[[paste0(v, "_z")]] <- (d[[v]] - p$mean) / p$sd
  }
  list(data = d, eligibility = ids)
}

train_raw <- fread(train_file)
valid_raw <- fread(valid_file)
train <- prepare_training(train_raw)
valid <- prepare_validation(valid_raw, train$bounds, train$scaling)

preprocessing <- merge(
  train$bounds,
  data.table(
    variable = sub("_interp$", "", train$scaling$variable),
    interpolated_variable = train$scaling$variable,
    standardization_mean = train$scaling$mean,
    standardization_sd = train$scaling$sd
  ),
  by = "variable"
)
preprocessing[, `:=`(
  development_dataset = "MIMIC-IV",
  time_window = "D1-D7",
  time_unit = "day",
  interpolation = "linear_between_observed_days_with_nearest_edge_carry",
  minimum_observed_days_per_index = 2L
)]
fwrite(preprocessing, file.path(tables, "all_icu_gbmtm_preprocessing.csv"))

eligibility <- rbind(
  train$eligibility[, dataset := "mimiciv"],
  valid$eligibility[, dataset := "mimiciii"]
)
setcolorder(eligibility, c("dataset", "icu_stay_id", "n_easix_days", "n_siri_days", "n_paired_days"))
fwrite(eligibility, file.path(tables, "all_icu_d1d7_dual_index_eligibility.csv"))

z_vars <- c("log2_easix_interp_z", "log2_siri_interp_z")
model_data <- as.data.frame(train$data[, .(
  id = paste0("mimiciv:", icu_stay_id),
  time = as.integer(day),
  log2_easix_interp_z,
  log2_siri_interp_z
)])
fwrite(as.data.table(model_data), file.path(tables, "all_icu_gbmtm_prepared_training.csv"))

calc_entropy <- function(post) {
  post <- as.matrix(post)
  if (ncol(post) <= 1) return(NA_real_)
  p <- pmax(post, .Machine$double.eps)
  1 + sum(p * log(p)) / (nrow(p) * log(ncol(p)))
}

metrics <- list()
fits <- list()
clinical_rows <- list()
level_label <- function(x) {
  if (x < -0.25) return("low")
  if (x > 0.25) return("high")
  "intermediate"
}
trend_label <- function(x) {
  if (x < -0.2) return("decreasing")
  if (x > 0.2) return("increasing")
  "stable"
}
clinical_profile_summary <- function(assignments, k) {
  id_map <- data.table(id = unique(model_data$id), trajectory_class = as.integer(assignments))
  d <- merge(as.data.table(model_data), id_map, by = "id")
  profiles <- d[, .(
    easix_mean = mean(log2_easix_interp_z),
    siri_mean = mean(log2_siri_interp_z),
    easix_day1 = mean(log2_easix_interp_z[time == 1]),
    easix_day7 = mean(log2_easix_interp_z[time == 7]),
    siri_day1 = mean(log2_siri_interp_z[time == 1]),
    siri_day7 = mean(log2_siri_interp_z[time == 7]),
    n = uniqueN(id)
  ), by = trajectory_class]
  daily <- d[, .(
    easix = mean(log2_easix_interp_z),
    siri = mean(log2_siri_interp_z)
  ), by = .(trajectory_class, time)]
  wide <- dcast(daily, trajectory_class ~ time, value.var = c("easix", "siri"))
  min_pairwise_rms <- if (k == 1) {
    NA_real_
  } else {
    distance_matrix <- as.matrix(dist(as.matrix(wide[, -1]))) / sqrt(ncol(wide) - 1)
    min(distance_matrix[upper.tri(distance_matrix)])
  }
  profiles[, `:=`(
    k = k,
    min_pairwise_profile_rms = min_pairwise_rms,
    profile_label = paste0(
      "EASIX ", vapply(easix_mean, level_label, character(1)), "-",
      vapply(easix_day7 - easix_day1, trend_label, character(1)), "; SIRI ",
      vapply(siri_mean, level_label, character(1)), "-",
      vapply(siri_day7 - siri_day1, trend_label, character(1))
    )
  )]
  profiles[]
}
for (k in classes) {
  model_path <- file.path(models, paste0("gbmtm_", tag, "_k", k, ".rds"))
  if (file.exists(model_path)) {
    log_msg("Loading existing GBMTM K=", k)
    fit <- readRDS(model_path)
  } else {
    log_msg("Fitting GBMTM K=", k)
    fit <- tryCatch(
      gbmt(
        x.names = z_vars,
        unit = "id",
        time = "time",
        d = degree,
        ng = k,
        data = model_data,
        scaling = 0,
        pruning = TRUE,
        quiet = TRUE
      ),
      error = function(e) e
    )
  }
  if (inherits(fit, "error")) {
    metrics[[as.character(k)]] <- data.table(
      k = k, AIC = NA_real_, BIC = NA_real_, entropy = NA_real_,
      min_avepp = NA_real_, mean_avepp = NA_real_, min_class_pct = NA_real_,
      class_proportions = NA_character_, status = "fit_failed",
      error = conditionMessage(fit)
    )
    log_msg("K=", k, " failed: ", conditionMessage(fit))
    next
  }
  fits[[as.character(k)]] <- fit
  saveRDS(fit, model_path)
  props <- as.numeric(prop.table(table(factor(fit$assign, levels = seq_len(k)))))
  post <- posterior(fit)
  avepp <- if (k == 1) 1 else fit$appa
  clinical <- clinical_profile_summary(fit$assign, k)
  clinical_rows[[as.character(k)]] <- clinical
  min_pairwise_profile_rms <- unique(clinical$min_pairwise_profile_rms)
  clinically_interpretable <- k >= 2 && is.finite(min_pairwise_profile_rms) && min_pairwise_profile_rms >= 0.50
  metrics[[as.character(k)]] <- data.table(
    k = k,
    AIC = unname(fit$ic["aic"]),
    BIC = unname(fit$ic["bic"]),
    entropy = calc_entropy(post),
    min_avepp = min(avepp, na.rm = TRUE),
    mean_avepp = mean(avepp, na.rm = TRUE),
    min_class_pct = min(props) * 100,
    min_pairwise_profile_rms = min_pairwise_profile_rms,
    clinically_interpretable = clinically_interpretable,
    class_proportions = paste0("C", seq_along(props), "=", sprintf("%.1f%%", props * 100), collapse = "; "),
    status = "fit_ok",
    error = ""
  )
}

existing_metric_file <- file.path(tables, "all_icu_gbmtm_candidate_metrics_k1_k6.csv")
metric_dt <- if (!length(metrics) && file.exists(existing_metric_file)) {
  fread(existing_metric_file)
} else {
  rbindlist(metrics, fill = TRUE)
}
setorder(metric_dt, k)
valid_candidates <- metric_dt[
  status == "fit_ok" & k >= 2 & min_class_pct >= minimum_class_pct &
    min_avepp >= 0.70 & entropy >= 0.80 & clinically_interpretable
]
geometric_elbow_k <- function(d) {
  d <- d[status == "fit_ok" & k %between% c(1, 6)][order(k)]
  if (nrow(d) != 6) return(NA_integer_)
  x <- (d$k - min(d$k)) / diff(range(d$k))
  y <- (d$BIC - min(d$BIC)) / diff(range(d$BIC))
  start <- c(x[1], y[1])
  end <- c(x[length(x)], y[length(y)])
  line <- end - start
  distance <- abs(line[1] * (start[2] - y) - (start[1] - x) * line[2]) /
    sqrt(sum(line^2))
  distance[c(1, length(distance))] <- NA_real_
  d$k[which.max(distance)]
}
elbow_k <- geometric_elbow_k(metric_dt)
selection_status <- "selected_geometric_bic_elbow_with_10pct_diagnostic_threshold"
if (!nrow(valid_candidates)) {
  valid_candidates <- metric_dt[status == "fit_ok" & k >= 2]
  selection_status <- "exploratory_no_candidate_met_full_diagnostic_criteria"
}
if (!nrow(valid_candidates)) stop("No GBMTM candidate with K=2-6 fitted successfully.")
if (is.finite(elbow_k) && elbow_k %in% valid_candidates$k) {
  selected_k <- elbow_k
} else {
  selected_k <- valid_candidates[order(abs(k - elbow_k), BIC)][1, k]
  selection_status <- "selected_nearest_diagnostically_eligible_k_to_geometric_elbow"
}
metric_dt[, selected := k == selected_k]
metric_dt[, selection_status := fifelse(selected, selection_status, "not_selected")]
fwrite(metric_dt, file.path(tables, "all_icu_gbmtm_candidate_metrics_k1_k6.csv"))
fwrite(rbindlist(clinical_rows), file.path(tables, "all_icu_gbmtm_clinical_interpretability_k1_k6.csv"))

selected_model <- fits[[as.character(selected_k)]]
if (is.null(selected_model)) {
  selected_model <- readRDS(file.path(models, paste0("gbmtm_", tag, "_k", selected_k, ".rds")))
}
saveRDS(selected_model, file.path(models, "gbmtm_all_icu_mimiciv_d1d7_easix_siri_selected.rds"))

fitted_profiles <- rbindlist(lapply(seq_along(selected_model$fitted), function(class_id) {
  fitted <- as.data.table(selected_model$fitted[[class_id]])
  fitted[, `:=`(gbmtm_class = class_id, day = seq_len(.N))]
  melt(
    fitted,
    id.vars = c("gbmtm_class", "day"),
    variable.name = "index",
    value.name = "fitted_z"
  )
}))
fitted_profiles[, `:=`(
  selected_k = selected_k,
  degree = as.integer(selected_model$call$d),
  pruning = as.logical(selected_model$call$pruning)
)]
fwrite(fitted_profiles, file.path(tables, "all_icu_gbmtm_selected_fitted_profiles.csv"))

prediction_grid <- seq(1, 7, length.out = 100)
prediction_bands <- rbindlist(lapply(seq_along(selected_model$reg), function(class_id) {
  rbindlist(lapply(z_vars, function(index_name) {
    prediction <- suppressWarnings(predict.lm(
      selected_model$reg[[class_id]][[index_name]],
      newdata = data.frame(t = prediction_grid),
      interval = "prediction",
      level = 0.95
    ))
    data.table(
      gbmtm_class = class_id,
      day = prediction_grid,
      index = index_name,
      fitted_z = prediction[, "fit"],
      prediction_low = prediction[, "lwr"],
      prediction_high = prediction[, "upr"],
      selected_k = selected_k,
      degree = as.integer(selected_model$call$d),
      pruning = as.logical(selected_model$call$pruning)
    )
  }))
}))
fwrite(prediction_bands, file.path(tables, "all_icu_gbmtm_selected_prediction_bands.csv"))

fixed_posterior <- function(model, newdata) {
  x_names <- model$call$x.names
  ids <- as.character(newdata[[model$call$unit]])
  tval <- as.numeric(newdata[[model$call$time]])
  x <- as.matrix(newdata[, x_names, drop = FALSE])
  kmax <- length(model$prior)
  row_ll <- matrix(NA_real_, nrow(x), kmax)
  for (k in seq_len(kmax)) {
    mu <- matrix(NA_real_, nrow(x), length(x_names))
    for (j in seq_along(x_names)) {
      mu[, j] <- predict(model$reg[[k]][[x_names[j]]], newdata = data.frame(t = tval, time = tval))
    }
    sigma <- as.matrix(model$Sigma[[k]]) + diag(1e-8, length(x_names))
    inv <- solve(sigma)
    log_det <- as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
    residual <- x - mu
    mahal <- rowSums((residual %*% inv) * residual)
    row_ll[, k] <- -0.5 * (ncol(x) * log(2 * pi) + log_det + mahal)
  }
  ll <- data.table(id = ids)
  for (k in seq_len(kmax)) ll[[paste0("ll_", k)]] <- row_ll[, k]
  ll_cols <- paste0("ll_", seq_len(kmax))
  summed <- ll[, lapply(.SD, sum), by = id, .SDcols = ll_cols]
  log_post <- as.matrix(summed[, ..ll_cols])
  log_post <- sweep(log_post, 2, log(model$prior), "+")
  max_log <- apply(log_post, 1, max)
  post <- exp(sweep(log_post, 1, max_log, "-"))
  post <- post / rowSums(post)
  out <- data.table(id = summed$id, gbmtm_class = max.col(post, ties.method = "first"))
  for (k in seq_len(kmax)) out[[paste0("post_prob_class_", k)]] <- post[, k]
  out[, max_post_prob := apply(post, 1, max)]
  out
}

train_projection <- fixed_posterior(selected_model, model_data)
train_projection[, `:=`(
  dataset = "mimiciv",
  icu_stay_id = sub("^mimiciv:", "", id),
  projection_role = "development_reconstruction"
)]
internal_reference <- data.table(
  id = unique(model_data$id),
  fitted_class = as.integer(selected_model$assign)
)
internal_reconstruction <- merge(
  train_projection[, .(id, projected_class = gbmtm_class)],
  internal_reference,
  by = "id"
)[, .(
  n = .N,
  exact_agreement = mean(projected_class == fitted_class)
)]
fwrite(internal_reconstruction, file.path(tables, "all_icu_gbmtm_internal_reconstruction.csv"))

valid_model_data <- as.data.frame(valid$data[, .(
  id = paste0("mimiciii:", icu_stay_id),
  time = as.integer(day),
  log2_easix_interp_z,
  log2_siri_interp_z
)])
valid_projection <- fixed_posterior(selected_model, valid_model_data)
valid_projection[, `:=`(
  dataset = "mimiciii",
  icu_stay_id = sub("^mimiciii:", "", id),
  projection_role = "frozen_mimiciv_cross_version_transport"
)]
assignments <- rbind(train_projection, valid_projection, fill = TRUE)
assignments[, selected_k := selected_k]
assignments[, selection_status := selection_status]
fwrite(assignments, file.path(tables, "all_icu_gbmtm_assignments_mimiciv_mimiciii.csv"))

distribution <- assignments[, .(
  n = .N,
  pct = 100 * .N / nrow(.SD),
  mean_max_post_prob = mean(max_post_prob)
), by = .(dataset, gbmtm_class)]
distribution[, pct := 100 * n / sum(n), by = dataset]
fwrite(distribution, file.path(tables, "all_icu_gbmtm_class_distribution.csv"))

profile_data <- rbind(
  train$data[, .(dataset = "mimiciv", icu_stay_id, day, log2_easix_interp_z, log2_siri_interp_z)],
  valid$data[, .(dataset = "mimiciii", icu_stay_id, day, log2_easix_interp_z, log2_siri_interp_z)]
)
profile_data[, icu_stay_id := as.character(icu_stay_id)]
profile_data <- merge(
  profile_data,
  assignments[, .(dataset, icu_stay_id, gbmtm_class)],
  by = c("dataset", "icu_stay_id")
)
profiles <- melt(
  profile_data,
  id.vars = c("dataset", "icu_stay_id", "day", "gbmtm_class"),
  variable.name = "index",
  value.name = "z"
)[, .(
  mean_z = mean(z),
  se_z = sd(z) / sqrt(.N),
  n_stays = uniqueN(icu_stay_id)
), by = .(dataset, gbmtm_class, day, index)]
profiles[, `:=`(ci_low = mean_z - 1.96 * se_z, ci_high = mean_z + 1.96 * se_z)]
fwrite(profiles, file.path(tables, "all_icu_gbmtm_trajectory_profiles.csv"))

selected_row <- metric_dt[k == selected_k]
report <- c(
  "# All-ICU D1-D7 EASIX/SIRI GBMTM",
  "",
  "- Development: MIMIC-IV; historical cross-version transport: MIMIC-III.",
  "- Inclusion: adult, first ICU stay per patient, ICU LOS >=24h, hospital LOS >=24h.",
  "- Inputs: daily log2(EASIX) and log2(SIRI), D1-D7; at least two observed days for each index.",
  "- Missing daily values: linear interpolation between observations with nearest observed edge carry.",
  paste0("- Candidate classes: K=1-6; selected K=", selected_k, "."),
  paste0("- Selection status: `", selection_status, "`."),
  paste0("- Selected BIC: ", selected_row$BIC, "."),
  paste0("- Selected entropy: ", sprintf("%.3f", selected_row$entropy), "."),
  paste0("- Minimum class proportion: ", sprintf("%.2f%%", selected_row$min_class_pct), "."),
  paste0("- Minimum AvePP: ", sprintf("%.3f", selected_row$min_avepp), "."),
  paste0("- Minimum pairwise profile RMS: ", sprintf("%.3f", selected_row$min_pairwise_profile_rms), "."),
  paste0("- Geometric BIC elbow over K=1-6: K=", elbow_k, "."),
  "- Selection hierarchy: post hoc geometric BIC elbow, then diagnostic confirmation (entropy >=0.80, minimum class >=10%, minimum AvePP >=0.70, and quantitative profile RMS >=0.50).",
  "- MIMIC-III assignments use the frozen MIMIC-IV preprocessing and complete GBMTM parameters; no refitting.",
  "",
  "Mortality outcomes were excluded from model fitting and class assignment."
)
writeLines(report, file.path(reports, "all_icu_gbmtm_d1d7_report.md"))
log_msg("Selected K=", selected_k, "; status=", selection_status)
