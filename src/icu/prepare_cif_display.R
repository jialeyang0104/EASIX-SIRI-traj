#!/usr/bin/env Rscript
# Prepare the requested CIF display without refitting association models.
suppressPackageStartupMessages({library(data.table); library(cmprsk); library(jsonlite)})
root <- file.path(Sys.getenv("PROJECT_ROOT", unset=getwd()), "src/icu")
out <- file.path(root, "aggregate")
settings <- fromJSON(file.path(dirname(dirname(root)), "config/icu_analysis_contract.json"))$cif_display
origin <- settings$origin_admission_day
last_day <- settings$maximum_admission_day
patients <- fread(file.path(root, "data_private/icu_mortality_analysis.csv"),
    select=c("dataset", "gbmtm_class", "trajectory_eligible", "icu_status", "icu_endpoint_time_days"))[
        trajectory_eligible==TRUE & !is.na(icu_status)]
original <- fread(file.path(out, "cif_curve_data.csv"))[event==settings$event]
curves <- tests <- risks <- summary_rows <- list()
for(db in c("mimiciv", "mimiciii")) {
  x <- copy(patients[dataset==db])
  x[, `:=`(display_time=pmin(icu_endpoint_time_days, last_day)-origin,
           display_status=fifelse(icu_endpoint_time_days<=last_day, icu_status, 0L))]
  fit <- cuminc(x$display_time, x$display_status, group=x$gbmtm_class, cencode=0)
  tests[[length(tests)+1L]] <- data.table(dataset=db, event=settings$event,
      origin_admission_day=origin, maximum_admission_day=last_day,
      n=nrow(x), icu_deaths=sum(x$display_status==1L),
      live_discharges=sum(x$display_status==2L), right_censored=sum(x$display_status==0L),
      gray_chisq=fit$Tests[as.character(settings$event), "stat"],
      df=fit$Tests[as.character(settings$event), "df"],
      p_value=fit$Tests[as.character(settings$event), "pv"])
  for(cls in 1:3) {
    g <- x[gbmtm_class==cls]
    source <- original[dataset==db & trajectory_class==cls]
    # Every class has observed ICU follow-up through Day 60. Carry only the
    # latest step value to the exact plotting boundary, within this support.
    stopifnot(max(source$admission_day)>=last_day, all(g$icu_status %in% 1:2))
    visible <- copy(source[admission_day<=last_day])
    if(tail(visible$admission_day, 1)<last_day) {
      boundary <- copy(tail(visible, 1))
      boundary[, `:=`(admission_day=last_day, days_since_landmark=last_day-origin)]
      visible <- rbind(visible, boundary)
    }
    direct <- fit[[paste(cls, settings$event)]]
    stopifnot(abs(tail(visible$cif, 1)-tail(direct$est, 1))<1e-10,
              abs(tail(visible$variance, 1)-tail(direct$var, 1))<1e-10)
    curves[[length(curves)+1L]] <- visible
    summary_rows[[length(summary_rows)+1L]] <- data.table(dataset=db, trajectory_class=cls,
        n=nrow(g), original_curve_end_day=max(source$admission_day), plotted_end_day=last_day,
        cif_at_day60=tail(visible$cif, 1), ci_low=tail(visible$ci_low, 1), ci_high=tail(visible$ci_high, 1),
        icu_deaths_by_day60=sum(g$display_status==1L),
        live_discharges_by_day60=sum(g$display_status==2L), still_in_icu_at_day60=sum(g$display_status==0L),
        numerical_check="Day-60 CIF and variance match independently administratively censored cuminc")
    for(day in settings$risk_days) risks[[length(risks)+1L]] <- data.table(dataset=db,
        trajectory_class=cls, admission_day=day, event_free_icu_risk_n=sum(g$icu_endpoint_time_days>=day), cohort_n=nrow(g))
  }
}
outputs <- list(cif_display_curve_data=rbindlist(curves), cif_display_gray_tests=rbindlist(tests),
    cif_display_numbers_at_risk=rbindlist(risks), cif_display_endpoint_summary=rbindlist(summary_rows))
for(nm in names(outputs)) fwrite(outputs[[nm]], file.path(out, paste0(nm, ".csv.tmp")))
for(nm in names(outputs)) stopifnot(file.rename(file.path(out, paste0(nm, ".csv.tmp")), file.path(out, paste0(nm, ".csv"))))
print(outputs$cif_display_gray_tests)
print(outputs$cif_display_endpoint_summary)
