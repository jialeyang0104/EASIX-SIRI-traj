#!/usr/bin/env Rscript

# Version 01R single source of truth for outcome-association adjustment.
# Prediction models and selection/weight models remain purpose-specific and do
# not use this contract unless they subsequently fit an outcome model.

v01r_outcome_covariates <- function(sex_variable = "sex", race_variable = "race_group") {
  c("age", sex_variable, race_variable, "oasis_day1")
}

v01r_adjustment_set_id <- "v01r_core4_age_sex_race_oasis"
v01r_adjustment_label <- "age + sex + race/ethnicity + Day 1 OASIS"
v01r_treatment_multiplicity_policy <- paste(
  "Exact raw global interaction P values are reported for four main and eight",
  "supplementary exploratory tests; no treatment-interaction multiplicity correction is applied."
)

v01r_assert_core4 <- function(formula, sex_variable = "sex", race_variable = "race_group",
                              allowed_extra = character()) {
  variables <- all.vars(formula)
  required <- v01r_outcome_covariates(sex_variable, race_variable)
  missing <- setdiff(required, variables)
  if (length(missing)) {
    stop("Version 01R outcome formula is missing core covariates: ", paste(missing, collapse = ", "))
  }
  forbidden <- setdiff(
    c("sepsis", "sofa_day1", "apsiii_day1", "mechanical_ventilation_day1",
      "crrt_day1", "aki_day1", "log2_easix_day1", "log2_siri_day1"),
    allowed_extra
  )
  present_forbidden <- intersect(forbidden, variables)
  if (length(present_forbidden)) {
    stop("Version 01R outcome formula contains undeclared additional adjustment terms: ",
         paste(present_forbidden, collapse = ", "))
  }
  invisible(TRUE)
}
