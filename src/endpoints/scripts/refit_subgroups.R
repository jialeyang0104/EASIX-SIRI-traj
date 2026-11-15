#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(survival)})
project_root <- normalizePath(".", mustWork=TRUE)
root <- file.path(project_root, "src", "endpoints")
tables <- file.path(root, "aggregate")
source(file.path(root, "scripts", "subgroup_models.R"))

# Same Day 1 transformations and fixed assignments as 05_58_v5_landmark.R.
assignments <- fread(file.path(root, "data_private", "all_icu_gbmtm_assignments_mimiciv_mimiciii.csv"))
assignments[, icu_stay_id := as.character(icu_stay_id)]
prepare <- function(db) {
  d <- fread(file.path(root, "data_private", paste0("all_icu_d1d7_daily_", db, ".csv")))[day == 1]
  d[, icu_stay_id := as.character(icu_stay_id)]
  x <- toupper(trimws(as.character(d$race_ethnicity)))
  d[, race_group := factor(fcase(
    grepl("WHITE", x), "White", grepl("BLACK|AFRICAN", x), "Black",
    grepl("HISPANIC|LATINO", x), "Hispanic", grepl("ASIAN", x), "Asian",
    is.na(x) | x == "" | grepl("UNKNOWN|UNABLE|DECLINED", x), "Unknown", default="Other"),
    levels=c("White", "Black", "Hispanic", "Asian", "Other", "Unknown"))]
  d[, sex := factor(fcase(toupper(sex)=="M", "Male", toupper(sex)=="F", "Female", default="Other"), levels=c("Female", "Male", "Other"))]
  d[, `:=`(dataset=db, male=as.integer(sex=="Male"), age_ge65=as.integer(age >= 65), aki_day1=as.integer(aki_stage_day1 >= 1))]
  d <- merge(d, assignments[dataset==db, .(icu_stay_id, trajectory_class=gbmtm_class)], by="icu_stay_id")
  d[, trajectory_class := factor(trajectory_class, levels=1:3)]
  stopifnot(!anyDuplicated(d$icu_stay_id))
  d
}
old <- fread(file.path(tables, "all_icu_cox_gbmtm_subgroup_interactions.csv"))
new_terms <- fit_subgroups(prepare("mimiciv"), prepare("mimiciii"), tables)
new <- fread(file.path(tables, "all_icu_cox_gbmtm_subgroup_interactions.csv"))
comparison <- merge(old, new, by=c("dataset", "analysis", "subgroup"), suffixes=c("_old_age_omitted", "_continuous_age"))
stopifnot(nrow(comparison)==28L,
  all(comparison$n_old_age_omitted==comparison$n_continuous_age),
  all(comparison$events_old_age_omitted==comparison$events_continuous_age),
  all(abs(comparison[subgroup!="age", interaction_p_old_age_omitted - interaction_p_continuous_age]) < 1e-10))
fwrite(comparison, file.path(tables, "subgroup_revision_comparison.csv"))
all_terms <- fread(file.path(tables, "all_models_full_terms.csv"))
all_terms <- rbindlist(list(all_terms[is.na(subgroup) | subgroup==""], new_terms), fill=TRUE)
fwrite(all_terms, file.path(tables, "all_models_full_terms.csv"))
cat("Focused subgroup refit completed:", format(Sys.time(), tz="UTC"), "UTC\n")
cat("All 28 interaction samples and event counts matched; six non-age raw tests unchanged per database/origin.\n")
print(new[subgroup=="age"])
print(fread(file.path(tables, "all_icu_cox_gbmtm_continuous_age_interactions.csv"))[, .(dataset, analysis, n, events, interaction_p)])
print(sessionInfo())
