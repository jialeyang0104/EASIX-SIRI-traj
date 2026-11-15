#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
root <- normalizePath(Sys.getenv("V5_ENDPOINT_ROOT", unset = file.path(project_root, "src", "endpoints")), mustWork = TRUE)
source(file.path(root, "scripts", "05_70_v01r_covariate_contract.R"))
tables <- file.path(root, "aggregate")
reports <- file.path(root, "aggregate")
logs <- file.path(root, "logs")
for (directory in c(tables, reports, logs)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

log_con <- file(file.path(logs, "all_icu_survival_baseline_extended.log"), open = "wt")
on.exit(close(log_con), add = TRUE)
log_msg <- function(...) {
  message(...)
  writeLines(paste0(...), log_con)
  flush(log_con)
}

collapse_race <- function(x) {
  x <- toupper(trimws(as.character(x)))
  fcase(
    grepl("WHITE", x), "White",
    grepl("BLACK|AFRICAN", x), "Black",
    grepl("HISPANIC|LATINO", x), "Hispanic",
    grepl("ASIAN", x), "Asian",
    is.na(x) | x == "" | grepl("UNKNOWN|UNABLE|DECLINED", x), "Unknown",
    default = "Other"
  )
}

prepare_day1 <- function(dataset_name) {
  d <- fread(file.path(root, "data_private", paste0("all_icu_d1d7_daily_", dataset_name, ".csv")))[day == 1]
  d[, `:=`(
    dataset = dataset_name,
    icu_stay_id = as.character(icu_stay_id),
    sex = factor(fcase(toupper(sex) == "M", "Male", toupper(sex) == "F", "Female", default = "Other"),
                 levels = c("Female", "Male", "Other")),
    race_group = factor(collapse_race(race_ethnicity),
                        levels = c("White", "Black", "Hispanic", "Asian", "Other", "Unknown")),
    log2_easix = fifelse(is.finite(easix_day) & easix_day > 0, log2(easix_day), NA_real_),
    log2_siri = fifelse(is.finite(siri_day) & siri_day > 0, log2(siri_day), NA_real_),
    aki_day1 = as.integer(aki_stage_day1 >= 1),
    age_ge65 = as.integer(age >= 65)
  )]
  d[, male := as.integer(sex == "Male")]
  d
}

dev <- prepare_day1("mimiciv")
val <- prepare_day1("mimiciii")

prepare_observation_flags <- function(dataset_name) {
  daily <- fread(file.path(root, "data_private", paste0("all_icu_d1d7_daily_", dataset_name, ".csv")))
  daily[, icu_stay_id := as.character(icu_stay_id)]
  flags <- daily[, .(
    late_easix_day5_7 = any(day %between% c(5, 7) & is.finite(easix_day) & easix_day > 0),
    late_siri_day5_7 = any(day %between% c(5, 7) & is.finite(siri_day) & siri_day > 0),
    late_paired_day5_7 = any(
      day %between% c(5, 7) & is.finite(easix_day) & easix_day > 0 &
        is.finite(siri_day) & siri_day > 0
    )
  ), by = icu_stay_id]
  calendar <- rbindlist(lapply(c("EASIX", "SIRI"), function(index_name) {
    observed <- if (index_name == "EASIX") {
      is.finite(daily$easix_day) & daily$easix_day > 0
    } else {
      is.finite(daily$siri_day) & daily$siri_day > 0
    }
    data.table(day = daily$day, observed = observed)[, .(
      n_patients_observed = sum(observed),
      n_parent_cohort = uniqueN(daily$icu_stay_id)
    ), by = day][, `:=`(
      dataset = dataset_name,
      index = index_name,
      observed_pct = 100 * n_patients_observed / n_parent_cohort
    )]
  }))
  list(flags = flags, calendar = calendar)
}

observation_dev <- prepare_observation_flags("mimiciv")
observation_val <- prepare_observation_flags("mimiciii")
calendar_observation <- rbind(observation_dev$calendar, observation_val$calendar, fill = TRUE)
setcolorder(calendar_observation, c(
  "dataset", "index", "day", "n_patients_observed", "n_parent_cohort", "observed_pct"
))
fwrite(calendar_observation, file.path(tables, "all_icu_calendar_day_observation_counts.csv"))

fit_scale <- function(x) {
  q <- quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  clipped <- pmin(pmax(x, q[1]), q[2])
  list(lower = q[1], upper = q[2], mean = mean(clipped, na.rm = TRUE), sd = sd(clipped, na.rm = TRUE))
}
apply_scale <- function(x, p) (pmin(pmax(x, p$lower), p$upper) - p$mean) / p$sd
scales <- list(easix = fit_scale(dev$log2_easix), siri = fit_scale(dev$log2_siri))
for (d in list(dev, val)) {
  d[, `:=`(
    log2_easix_z = apply_scale(log2_easix, scales$easix),
    log2_siri_z = apply_scale(log2_siri, scales$siri)
  )]
}

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
format_continuous <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  q <- quantile(x, c(0.25, 0.50, 0.75), names = FALSE)
  sprintf("%.1f (%.1f-%.1f)", q[2], q[1], q[3])
}
format_binary <- function(x) {
  ok <- !is.na(x)
  if (!any(ok)) return("NA")
  sprintf("%s (%.1f%%)", format(sum(x[ok] == 1), big.mark = ",", trim = TRUE), 100 * mean(x[ok] == 1))
}
safe_p_continuous <- function(x, group) {
  ok <- is.finite(x) & !is.na(group)
  if (sum(ok) < 3 || uniqueN(group[ok]) < 2) return(NA_real_)
  tryCatch(kruskal.test(x[ok], group[ok])$p.value, error = function(e) NA_real_)
}
safe_p_binary <- function(x, group) {
  ok <- !is.na(x) & !is.na(group)
  tab <- table(x[ok], group[ok])
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  tryCatch(suppressWarnings(chisq.test(tab, correct = FALSE)$p.value), error = function(e) NA_real_)
}
smd_continuous <- function(x, group) {
  lev <- levels(droplevels(factor(group)))
  if (length(lev) != 2) return(NA_real_)
  x1 <- x[group == lev[1] & is.finite(x)]
  x2 <- x[group == lev[2] & is.finite(x)]
  pooled <- sqrt((var(x1) + var(x2)) / 2)
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  abs(mean(x1) - mean(x2)) / pooled
}
smd_binary <- function(x, group) {
  lev <- levels(droplevels(factor(group)))
  if (length(lev) != 2) return(NA_real_)
  p1 <- mean(x[group == lev[1]] == 1, na.rm = TRUE)
  p2 <- mean(x[group == lev[2]] == 1, na.rm = TRUE)
  denom <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  abs(p1 - p2) / denom
}

baseline_spec <- list(
  list("Age, years", "age", "continuous"),
  list("Male sex", "male", "binary"),
  list("White race/ethnicity", "white", "binary"),
  list("Sepsis", "sepsis", "binary"),
  list("OASIS", "oasis_day1", "continuous"),
  list("SOFA", "sofa_day1", "continuous"),
  list("SIRS score", "sirs_day1", "continuous"),
  list("Mechanical ventilation, Day 1", "mechanical_ventilation_day1", "binary"),
  list("Vasopressor, Day 1", "vasopressor_day1", "binary"),
  list("CRRT, Day 1", "crrt_day1", "binary"),
  list("AKI, Day 1", "aki_day1", "binary"),
  list("ICU length of stay, days", "icu_los_days", "continuous"),
  list("Hospital length of stay, days", "hospital_los_days", "continuous"),
  list("Day 1 EASIX", "easix_day", "continuous"),
  list("Day 1 SIRI", "siri_day", "continuous"),
  list("28-day mortality", "mortality_28d", "binary"),
  list("90-day mortality", "mortality_90d", "binary")
)

make_baseline_table <- function(d, group_col, group_labels, include_smd = FALSE) {
  d <- copy(d)
  d[, white := as.integer(race_group == "White")]
  d[, group_value := factor(get(group_col), levels = names(group_labels), labels = unname(group_labels))]
  rows <- list(data.table(Characteristic = "Patients", Type = "count"))
  for (label in unname(group_labels)) {
    rows[[1]][[label]] <- format(sum(d$group_value == label, na.rm = TRUE), big.mark = ",", trim = TRUE)
  }
  rows[[1]][["P value"]] <- ""
  if (include_smd) rows[[1]][["SMD"]] <- ""
  for (spec in baseline_spec) {
    label <- spec[[1]]; variable <- spec[[2]]; type <- spec[[3]]
    row <- data.table(Characteristic = label, Type = type)
    for (group_label in unname(group_labels)) {
      values <- d[group_value == group_label][[variable]]
      row[[group_label]] <- if (type == "continuous") format_continuous(values) else format_binary(values)
    }
    p <- if (type == "continuous") safe_p_continuous(d[[variable]], d$group_value) else safe_p_binary(d[[variable]], d$group_value)
    row[["P value"]] <- format_p(p)
    if (include_smd) {
      smd <- if (type == "continuous") smd_continuous(d[[variable]], d$group_value) else smd_binary(d[[variable]], d$group_value)
      row[["SMD"]] <- ifelse(is.finite(smd), sprintf("%.3f", smd), "NA")
    }
    rows[[length(rows) + 1L]] <- row
  }
  rbindlist(rows, fill = TRUE)
}

overall <- rbind(copy(dev), copy(val), fill = TRUE)
overall[, cohort := factor(dataset, levels = c("mimiciv", "mimiciii"))]
table1 <- make_baseline_table(
  overall, "cohort", c(mimiciv = "MIMIC-IV", mimiciii = "MIMIC-III"), include_smd = FALSE
)
fwrite(table1, file.path(tables, "all_icu_main_table1_cohort_baseline_comparison.csv"))

gb_assignments <- fread(file.path(root, "data_private", "all_icu_gbmtm_assignments_mimiciv_mimiciii.csv"))
gb_assignments[, icu_stay_id := as.character(icu_stay_id)]
eligibility <- fread(file.path(root, "data_private", "all_icu_d1d7_dual_index_eligibility.csv"))
eligibility[, icu_stay_id := as.character(icu_stay_id)]
attach_gbmtm <- function(d, dataset_name) {
  out <- merge(
    copy(d),
    gb_assignments[dataset == dataset_name, .(icu_stay_id, trajectory_class = gbmtm_class, max_post_prob)],
    by = "icu_stay_id"
  )
  out <- merge(
    out,
    eligibility[dataset == dataset_name, .(icu_stay_id, n_easix_days, n_siri_days, n_paired_days)],
    by = "icu_stay_id",
    all.x = TRUE
  )
  observation_flags <- if (dataset_name == "mimiciv") observation_dev$flags else observation_val$flags
  out <- merge(out, observation_flags, by = "icu_stay_id", all.x = TRUE)
  out[, trajectory_class := factor(trajectory_class, levels = 1:3)]
  out
}
traj_dev <- attach_gbmtm(dev, "mimiciv")
traj_val <- attach_gbmtm(val, "mimiciii")
table2_dev <- make_baseline_table(
  traj_dev, "trajectory_class", c(`1` = "Class 1", `2` = "Class 2", `3` = "Class 3")
)
table2_dev[, Dataset := "MIMIC-IV"]
table2_val <- make_baseline_table(
  traj_val, "trajectory_class", c(`1` = "Class 1", `2` = "Class 2", `3` = "Class 3")
)
table2_val[, Dataset := "MIMIC-III"]
fwrite(rbind(table2_dev, table2_val, fill = TRUE), file.path(tables, "all_icu_main_table2_gbmtm_class_baseline_comparison.csv"))

missingness_rows <- rbindlist(lapply(list(mimiciv = dev, mimiciii = val), function(d) {
  rbindlist(lapply(c(
    "easix_day", "siri_day", "sirs_day1", "oasis_day1", "sofa_day1",
    "race_group", "mechanical_ventilation_day1", "vasopressor_day1", "crrt_day1", "aki_stage_day1"
  ), function(variable) {
    data.table(
      variable = variable,
      n_total = nrow(d),
      n_available = sum(!is.na(d[[variable]])),
      available_pct = 100 * mean(!is.na(d[[variable]]))
    )
  }))
}), idcol = "dataset")
fwrite(missingness_rows, file.path(tables, "all_icu_supp_table_baseline_missingness.csv"))

base_covariates <- v01r_outcome_covariates()
static_specs <- list(
  static_easix = c("log2_easix_z"),
  static_siri = c("log2_siri_z"),
  static_joint = c("log2_easix_z", "log2_siri_z"),
  static_joint_interaction = c("log2_easix_z * log2_siri_z"),
  static_sirs_secondary = c("sirs_day1")
)

cox_terms <- list()
ph_rows <- list()
cox_counter <- 0L
ph_counter <- 0L

extract_cox <- function(fit, metadata) {
  s <- summary(fit)
  term_names <- rownames(s$coefficients)
  result <- data.table(
    term = term_names,
    estimate = s$coefficients[term_names, "coef"],
    se = s$coefficients[term_names, if ("robust se" %in% colnames(s$coefficients)) "robust se" else "se(coef)"],
    hazard_ratio = s$conf.int[term_names, "exp(coef)"],
    ci_low = s$conf.int[term_names, "lower .95"],
    ci_high = s$conf.int[term_names, "upper .95"],
    p_value = s$coefficients[term_names, "Pr(>|z|)"]
  )
  for (name in names(metadata)) result[[name]] <- metadata[[name]]
  setcolorder(result, c(names(metadata), "term", "estimate", "se", "hazard_ratio", "ci_low", "ci_high", "p_value"))
  result
}

extract_ph <- function(fit, metadata) {
  z <- tryCatch(cox.zph(fit, transform = "km"), error = function(e) NULL)
  if (is.null(z)) return(NULL)
  tab <- as.data.table(z$table, keep.rownames = "term")
  if ("p" %in% names(tab)) setnames(tab, "p", "p_value")
  for (name in names(metadata)) tab[[name]] <- metadata[[name]]
  setcolorder(tab, c(names(metadata), setdiff(names(tab), names(metadata))))
  tab
}

fit_cox_model <- function(d, time_col, event_col, exposure_rhs, metadata) {
  rhs <- paste(c(exposure_rhs, base_covariates), collapse = " + ")
  formula <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", rhs))
  required <- all.vars(formula)
  dd <- d[complete.cases(d[, ..required]) & get(time_col) > 0]
  if (nrow(dd) < 100 || sum(dd[[event_col]]) < 20) return(NULL)
  fit <- tryCatch(
    suppressWarnings(coxph(formula, data = dd, ties = "efron", x = TRUE, model = TRUE)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  metadata$formula <- paste(deparse(formula), collapse = " ")
  metadata$n_eligible <- nrow(d)
  metadata$n_excluded_missing_or_nonpositive <- nrow(d)-nrow(dd)
  metadata$n <- nrow(dd)
  metadata$events <- sum(dd[[event_col]])
  list(terms = extract_cox(fit, metadata), ph = extract_ph(fit, metadata), fit = fit, data = dd)
}

all_full_terms <- list()
extract_cox_original <- extract_cox
extract_cox <- function(fit,metadata) {
 z <- extract_cox_original(fit,metadata)
 all_full_terms[[length(all_full_terms)+1L]] <<- z
 z
}
horizons <- list(
  `28-day` = c(time = "followup_28d_days", event = "event_28d"),
  `90-day` = c(time = "followup_90d_days", event = "event_90d")
)
for (dataset_name in c("mimiciv", "mimiciii")) {
  d <- if (dataset_name == "mimiciv") dev else val
  for (horizon_name in names(horizons)) {
    horizon <- horizons[[horizon_name]]
    for (model_name in names(static_specs)) {
      admission_result <- fit_cox_model(
        d, horizon[["time"]], horizon[["event"]], static_specs[[model_name]],
        list(
          dataset = dataset_name, horizon = horizon_name,
          analysis = "static_day1_admission_inception", model = model_name
        )
      )
      if (!is.null(admission_result)) {
        cox_counter <- cox_counter + 1L
        cox_terms[[cox_counter]] <- admission_result$terms
        if (!is.null(admission_result$ph)) {
          ph_counter <- ph_counter + 1L
          ph_rows[[ph_counter]] <- admission_result$ph
        }
      }
      landmark_data <- copy(d)[get(horizon[["time"]]) > 1]
      landmark_data[, static_landmark_time := get(horizon[["time"]]) - 1]
      landmark_result <- fit_cox_model(
        landmark_data, "static_landmark_time", horizon[["event"]], static_specs[[model_name]],
        list(
          dataset = dataset_name, horizon = horizon_name,
          analysis = "static_day1_24h_landmark", model = model_name
        )
      )
      if (is.null(landmark_result)) next
      cox_counter <- cox_counter + 1L
      cox_terms[[cox_counter]] <- landmark_result$terms
      if (!is.null(landmark_result$ph)) {
        ph_counter <- ph_counter + 1L
        ph_rows[[ph_counter]] <- landmark_result$ph
      }
    }
  }
}

piecewise_rows <- list()
piecewise_counter <- 0L
piecewise_covariate_tv_rows <- list()
piecewise_covariate_tv_counter <- 0L
piecewise_specs <- list(
  static_day1_admission_inception_piecewise = list(
    time_col = "followup_90d_days", minimum_time = 0, cuts = c(7, 28),
    labels = c("ICU admission to Day 7", "Admission Day 7 to Day 28", "Admission Day 28 to Day 90")
  ),
  static_day1_24h_landmark_piecewise = list(
    time_col = "static_landmark_time", minimum_time = 1, cuts = c(6, 27),
    labels = c("24-h landmark to admission Day 7", "Admission Day 7 to Day 28", "Admission Day 28 to Day 90")
  )
)
for (dataset_name in c("mimiciv", "mimiciii")) {
  d <- if (dataset_name == "mimiciv") copy(dev) else copy(val)
  required <- c(
    "icu_stay_id", "followup_90d_days", "event_90d", "log2_easix_z", "log2_siri_z",
    base_covariates
  )
  for (analysis_name in names(piecewise_specs)) {
    specification <- piecewise_specs[[analysis_name]]
    dd <- d[complete.cases(d[, ..required]) & followup_90d_days > specification$minimum_time]
    if (analysis_name == "static_day1_24h_landmark_piecewise") {
      dd[, static_landmark_time := followup_90d_days - 1]
    }
    split <- survSplit(
      as.formula(paste0("Surv(", specification$time_col, ", event_90d) ~ .")),
      data = as.data.frame(dd), cut = specification$cuts,
      start = "tstart", end = "tstop", event = "event_split", episode = "interval_id"
    )
    split$interval <- factor(
      split$interval_id, levels = 1:3,
      labels = specification$labels
    )
  for (interval_id in 1:3) {
    split[[paste0("easix_interval_", interval_id)]] <- split$log2_easix_z * as.integer(split$interval_id == interval_id)
    split[[paste0("siri_interval_", interval_id)]] <- split$log2_siri_z * as.integer(split$interval_id == interval_id)
    split[[paste0("age_interval_", interval_id)]] <- split$age * as.integer(split$interval_id == interval_id)
    split[[paste0("oasis_interval_", interval_id)]] <- split$oasis_day1 * as.integer(split$interval_id == interval_id)
  }
  exposure_terms <- c(paste0("easix_interval_", 1:3), paste0("siri_interval_", 1:3))
  formula <- as.formula(paste0(
    "Surv(tstart, tstop, event_split) ~ ",
    paste(c(exposure_terms, base_covariates, "strata(interval)", "cluster(icu_stay_id)"), collapse = " + ")
  ))
    fit <- coxph(formula, data = split, ties = "efron", x = TRUE)
    extracted <- extract_cox(
      fit,
      list(
        dataset = dataset_name, horizon = "90-day",
        analysis = analysis_name, model = "static_joint"
      )
    )[term %chin% exposure_terms]
    extracted[, `:=`(
      exposure = fifelse(grepl("^easix", term), "EASIX", "SIRI"),
      interval = specification$labels[as.integer(sub(".*_", "", term))],
      n = nrow(dd),
      events = sum(dd$event_90d)
    )]
    piecewise_counter <- piecewise_counter + 1L
    piecewise_rows[[piecewise_counter]] <- extracted

  covariate_tv_terms <- c(
    paste0("age_interval_", 1:3), paste0("oasis_interval_", 1:3)
  )
  sensitivity_formula <- as.formula(paste0(
    "Surv(tstart, tstop, event_split) ~ ",
    paste(c(
      exposure_terms, covariate_tv_terms, "sex", "race_group",
      "strata(interval)", "cluster(icu_stay_id)"
    ), collapse = " + ")
  ))
    sensitivity_fit <- coxph(sensitivity_formula, data = split, ties = "efron", x = TRUE)
    sensitivity_extracted <- extract_cox(
      sensitivity_fit,
      list(
        dataset = dataset_name, horizon = "90-day",
        analysis = paste0(analysis_name, "_covariate_tv_sensitivity"),
        model = "static_joint"
      )
    )[term %chin% exposure_terms]
    sensitivity_extracted[, `:=`(
      exposure = fifelse(grepl("^easix", term), "EASIX", "SIRI"),
      interval = specification$labels[as.integer(sub(".*_", "", term))],
      n = nrow(dd),
      events = sum(dd$event_90d)
    )]
    piecewise_covariate_tv_counter <- piecewise_covariate_tv_counter + 1L
    piecewise_covariate_tv_rows[[piecewise_covariate_tv_counter]] <- sensitivity_extracted
  }
}
fwrite(rbindlist(piecewise_rows), file.path(tables, "all_icu_cox_static_piecewise_effects.csv"))
fwrite(
  rbindlist(piecewise_covariate_tv_rows),
  file.path(tables, "all_icu_cox_static_piecewise_covariate_tv_sensitivity.csv")
)

trajectory_models <- list()
trajectory_ph <- list()
trajectory_counter <- 0L
trajectory_ph_counter <- 0L
trajectory_sensitivity <- list()
trajectory_sensitivity_counter <- 0L

fit_trajectory_cox <- function(d, dataset_name, horizon_name, analysis_name, restriction = NULL) {
  horizon <- horizons[[horizon_name]]
  dd <- copy(d)
  time_col <- horizon[["time"]]
  event_col <- horizon[["event"]]
  if (analysis_name == "day7_landmark") {
    dd <- dd[icu_los_days >= 7 & mortality_7d == 0 & followup_90d_days > 7]
    dd[, analysis_time := get(time_col) - 7]
    time_col <- "analysis_time"
  }
  if (!is.null(restriction) && restriction == "paired_days_ge2") dd <- dd[n_paired_days >= 2]
  if (!is.null(restriction) && restriction == "posterior_ge_0.90") dd <- dd[max_post_prob >= 0.90]
  if (!is.null(restriction) && restriction == "late_observation_day5_7") {
    dd <- dd[late_easix_day5_7 & late_siri_day5_7]
  }
  result <- fit_cox_model(
    dd, time_col, event_col, "trajectory_class",
    list(
      dataset = dataset_name, horizon = horizon_name,
      analysis = analysis_name, model = "gbmtm_k3_d2_pruning",
      restriction = ifelse(is.null(restriction), "none", restriction)
    )
  )
  result
}

for (dataset_name in c("mimiciv", "mimiciii")) {
  d <- if (dataset_name == "mimiciv") traj_dev else traj_val
  for (horizon_name in names(horizons)) {
    for (analysis_name in c("icu_admission_inception", "day7_landmark")) {
      result <- fit_trajectory_cox(d, dataset_name, horizon_name, analysis_name)
      if (!is.null(result)) {
        trajectory_counter <- trajectory_counter + 1L
        trajectory_models[[trajectory_counter]] <- result$terms
        if (!is.null(result$ph)) {
          trajectory_ph_counter <- trajectory_ph_counter + 1L
          trajectory_ph[[trajectory_ph_counter]] <- result$ph
        }
      }
    }
    for (restriction in c("paired_days_ge2", "posterior_ge_0.90", "late_observation_day5_7")) {
      for (analysis_name in c("icu_admission_inception", "day7_landmark")) {
        result <- fit_trajectory_cox(d, dataset_name, horizon_name, analysis_name, restriction)
        if (!is.null(result)) {
          trajectory_sensitivity_counter <- trajectory_sensitivity_counter + 1L
          trajectory_sensitivity[[trajectory_sensitivity_counter]] <- result$terms
        }
      }
    }
  }
}

static_cox <- rbindlist(cox_terms, fill = TRUE)
static_ph <- rbindlist(ph_rows, fill = TRUE)
trajectory_cox <- rbindlist(trajectory_models, fill = TRUE)
trajectory_ph_dt <- rbindlist(trajectory_ph, fill = TRUE)
trajectory_sensitivity_dt <- rbindlist(trajectory_sensitivity, fill = TRUE)
fwrite(static_cox, file.path(tables, "all_icu_cox_static_effects.csv"))
fwrite(static_ph, file.path(tables, "all_icu_cox_static_ph_assumptions.csv"))
fwrite(trajectory_cox, file.path(tables, "all_icu_cox_gbmtm_trajectory_effects.csv"))
fwrite(trajectory_ph_dt, file.path(tables, "all_icu_cox_gbmtm_ph_assumptions.csv"))
fwrite(trajectory_sensitivity_dt, file.path(tables, "all_icu_cox_gbmtm_sensitivity_effects.csv"))

km_curves <- list()
km_risk <- list()
km_tests <- list()
km_counter <- 0L
for (dataset_name in c("mimiciv", "mimiciii")) {
  d0 <- if (dataset_name == "mimiciv") traj_dev else traj_val
  for (horizon_name in names(horizons)) {
    horizon <- horizons[[horizon_name]]
    for (analysis_name in c("icu_admission_inception", "day7_landmark")) {
      d <- copy(d0)
      time_col <- horizon[["time"]]
      event_col <- horizon[["event"]]
      max_time <- if (horizon_name == "28-day") 28 else 90
      if (analysis_name == "day7_landmark") {
        d <- d[icu_los_days >= 7 & mortality_7d == 0 & followup_90d_days > 7]
        d[, analysis_time := get(time_col) - 7]
        time_col <- "analysis_time"
        max_time <- max_time - 7
      }
      d <- d[complete.cases(d[, c(time_col, event_col, "trajectory_class"), with = FALSE])]
      fit <- survfit(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ trajectory_class")), data = d, conf.type = "log-log")
      times <- seq(0, max_time, by = 1)
      s <- summary(fit, times = times, extend = TRUE)
      km_counter <- km_counter + 1L
      km_curves[[km_counter]] <- data.table(
        dataset = dataset_name, horizon = horizon_name, analysis = analysis_name,
        time = s$time,
        trajectory_class = as.integer(sub("trajectory_class=", "", s$strata)),
        survival = s$surv, ci_low = s$lower, ci_high = s$upper,
        n_risk = s$n.risk, n_event = s$n.event
      )
      risk_times <- if (horizon_name == "28-day") c(0, 7, 14, max_time) else c(0, 30, 60, max_time)
      risk_times <- sort(unique(risk_times[risk_times <= max_time]))
      rs <- summary(fit, times = risk_times, extend = TRUE)
      km_risk[[km_counter]] <- data.table(
        dataset = dataset_name, horizon = horizon_name, analysis = analysis_name,
        time = rs$time,
        trajectory_class = as.integer(sub("trajectory_class=", "", rs$strata)),
        n_risk = rs$n.risk
      )
      sd <- survdiff(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ trajectory_class")), data = d)
      km_tests[[km_counter]] <- data.table(
        dataset = dataset_name, horizon = horizon_name, analysis = analysis_name,
        n = nrow(d), events = sum(d[[event_col]]),
        logrank_chisq = sd$chisq, df = length(sd$n) - 1,
        logrank_p = pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
      )
    }
  }
}
fwrite(rbindlist(km_curves), file.path(tables, "all_icu_km_gbmtm_curve_data.csv"))
fwrite(rbindlist(km_risk), file.path(tables, "all_icu_km_gbmtm_numbers_at_risk.csv"))
fwrite(rbindlist(km_tests), file.path(tables, "all_icu_km_gbmtm_logrank_tests.csv"))

source(file.path(root, "scripts", "subgroup_models.R"))
subgroup_full_terms <- fit_subgroups(traj_dev, traj_val, tables)
all_full_terms[[length(all_full_terms)+1L]] <- subgroup_full_terms

flow_rows <- list()
for (dataset_name in c("mimiciv", "mimiciii")) {
  d <- if (dataset_name == "mimiciv") dev else val
  trajectory <- if (dataset_name == "mimiciv") traj_dev else traj_val
  raw_flow <- fread(file.path(root, "aggregate", paste0("all_icu_flow_counts_", dataset_name, ".csv")))
  raw_flow[, `:=`(stage_type = "sequential", parent_stage = shift(step, type = "lag"))]
  raw_flow[step_order == 1, parent_stage := NA_character_]
  raw_flow[, stage_id := paste0("cohort_", step_order)]
  additional <- data.table(
    dataset = dataset_name,
    step_order = max(raw_flow$step_order) + 1:5,
    step = c(
      "Day 1 EASIX available", "Day 1 SIRI available", "Day 1 EASIX and SIRI jointly available",
      "D1-D7 trajectory eligible", "GBMTM classes analysed"
    ),
    n = c(
      sum(is.finite(d$log2_easix)), sum(is.finite(d$log2_siri)),
      sum(is.finite(d$log2_easix) & is.finite(d$log2_siri)),
      nrow(trajectory), nrow(trajectory)
    ),
    excluded_from_previous = NA_integer_,
    stage_type = c("branch", "branch", "branch", "branch", "analysis"),
    parent_stage = c(
      tail(raw_flow$step, 1), tail(raw_flow$step, 1), tail(raw_flow$step, 1),
      tail(raw_flow$step, 1), "D1-D7 trajectory eligible"
    ),
    stage_id = c(
      "static_easix", "static_siri", "static_joint", "trajectory", "gbmtm"
    )
  )
  flow_rows[[dataset_name]] <- rbind(raw_flow, additional, fill = TRUE)
}
fwrite(rbindlist(flow_rows), file.path(tables, "all_icu_flow_diagram_counts.csv"))

main_static <- static_cox[
  term %chin% c("log2_easix_z", "log2_siri_z", "sirs_day1") &
    model %chin% c("static_easix", "static_siri", "static_joint", "static_sirs_secondary") &
    analysis == "static_day1_24h_landmark"
]
main_static[, `HR (95% CI)` := sprintf("%.2f (%.2f-%.2f)", hazard_ratio, ci_low, ci_high)]
main_static[, `P value` := vapply(p_value, format_p, character(1))]
fwrite(
  main_static[, .(Dataset = dataset, Horizon = horizon, Model = model, Term = term, N = n, Events = events, `HR (95% CI)`, `P value`)],
  file.path(tables, "all_icu_main_table3_static_cox.csv")
)

main_dynamic <- trajectory_cox[
  grepl("^trajectory_class", term) & restriction == "none" &
    analysis == "day7_landmark"
]
main_dynamic[, `HR (95% CI)` := sprintf("%.2f (%.2f-%.2f)", hazard_ratio, ci_low, ci_high)]
main_dynamic[, `P value` := vapply(p_value, format_p, character(1))]
fwrite(
  main_dynamic[, .(Dataset = dataset, Horizon = horizon, Analysis = analysis, Class = term, N = n, Events = events, `HR (95% CI)`, `P value`)],
  file.path(tables, "all_icu_main_table4_dynamic_cox.csv")
)


# Explicit primary/legacy time-origin manifest from fitted models.
manifest <- unique(rbindlist(list(static_cox, trajectory_cox), fill=TRUE)[, .(dataset,horizon,analysis,model,formula,n_eligible,n_excluded_missing_or_nonpositive,n,events)])
manifest[, `:=`(primary = analysis %chin% c("static_day1_24h_landmark", "day7_landmark"), origin_days = fifelse(analysis=="day7_landmark",7,fifelse(analysis=="static_day1_24h_landmark",1,0)))]
manifest[, remaining_horizon_days := as.numeric(sub("-day", "", horizon))-origin_days]
fwrite(manifest,file.path(tables,"model_manifest.csv"))
fwrite(rbindlist(lapply(names(scales),function(nm) data.table(exposure=nm,as.data.table(scales[[nm]])))),file.path(tables,"legacy_scale_parameters.csv"))
# Day 7 trajectory piecewise effects: fixed admission Day 28 boundary.
pw <- list()
for (dataset_name in c("mimiciv","mimiciii")) {
 d <- copy(if(dataset_name=="mimiciv") traj_dev else traj_val)
 d <- d[icu_los_days>=7 & mortality_7d==0 & followup_90d_days>7]
 d[, analysis_time:=followup_90d_days-7]
 required <- c("analysis_time","event_90d","trajectory_class",base_covariates)
 dd <- d[complete.cases(d[,..required])]
 sp <- survSplit(Surv(analysis_time,event_90d)~.,data=as.data.frame(dd),cut=21,start="tstart",end="tstop",event="event_split",episode="interval_id")
 for (j in 1:2) for(k in 2:3) sp[[paste0("class",k,"_interval",j)]] <- as.integer(sp$trajectory_class==k & sp$interval_id==j)
 et <- c("class2_interval1","class3_interval1","class2_interval2","class3_interval2")
 formula <- as.formula(paste("Surv(tstart,tstop,event_split) ~",paste(c(et,base_covariates,"strata(interval_id)","cluster(icu_stay_id)"),collapse=" + ")))
 fit <- coxph(formula,data=sp,ties="efron",x=TRUE)
 z <- extract_cox(fit,list(dataset=dataset_name,horizon="90-day",analysis="day7_landmark_piecewise",model="gbmtm_k3_d2_pruning",n=nrow(dd),events=sum(dd$event_90d),formula=paste(deparse(formula),collapse=" ")))
 z[, interval:=fifelse(grepl("interval1$",term),"Admission Day 7 to Day 28",fifelse(grepl("interval2$",term),"Admission Day 28 to Day 90","common adjustment"))]
 pw[[dataset_name]]<-z
}
fwrite(rbindlist(pw),file.path(tables,"all_icu_cox_gbmtm_piecewise_effects.csv"))
fwrite(rbindlist(all_full_terms,fill=TRUE),file.path(tables,"all_models_full_terms.csv"))
log_msg("Landmark core refits completed.")
