#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(survival); library(cmprsk); library(jsonlite)})

root <- file.path(Sys.getenv("PROJECT_ROOT", unset=getwd()), "src/icu")
out <- file.path(root, "aggregate")
source(file.path(dirname(root), "endpoints/scripts/05_70_v01r_covariate_contract.R"))
d <- fread(file.path(root, "data_private/icu_mortality_analysis.csv"))
stopifnot(uniqueN(d, by=c("dataset", "icu_stay_id")) == nrow(d))
d[, trajectory_class := factor(gbmtm_class, levels=1:3)]
d[, sex_group := factor(sex_group, levels=c("Female", "Male", "Other"))]
d[, race_group := fifelse(race_group %chin% c("Other", "Unknown"), "Other/Unknown", race_group)]
d[, race_group := factor(race_group, levels=c("White", "Black", "Hispanic", "Asian", "Other/Unknown"))]
base <- v01r_outcome_covariates(sex_variable="sex_group")
all_terms <- manifest <- diagnostics <- list()

sample_data <- function(data, rhs, origin) {
  x <- copy(data)
  x[, time := icu_endpoint_time_days-origin]
  required <- unique(c(all.vars(as.formula(paste("~", rhs))), "time", "icu_status", "icu_stay_id"))
  x <- x[complete.cases(x[, ..required]) & time > 0]
  for (v in names(x)[vapply(x, is.factor, logical(1))]) x[, (v) := droplevels(get(v))]
  x
}

fit_model <- function(data, rhs, metadata, kind="Fine-Gray", ph=FALSE) {
  x <- sample_data(data, rhs, metadata$origin_days)
  warning_text <- character()
  if (kind == "Fine-Gray") {
    x[, event_state := factor(icu_status, levels=0:2,
                              labels=c("Censored", "ICU death", "Live ICU discharge"))]
    expansion_formula <- as.formula(paste("Surv(time, event_state) ~", rhs, "+ icu_stay_id"))
    expanded <- finegray(expansion_formula, data=x, etype="ICU death", na.action=na.fail,
                         id=icu_stay_id)
    formula <- as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", rhs))
    fit <- withCallingHandlers(
      coxph(formula, data=expanded, weights=fgwt, cluster=icu_stay_id,
            ties="breslow", singular.ok=FALSE, x=TRUE, y=TRUE, model=TRUE,
            na.action=na.fail, control=coxph.control(iter.max=50)),
      warning=function(w) {warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning")})
    effect_name <- "sHR"
  } else {
    formula <- as.formula(paste("Surv(time, icu_status == 1L) ~", rhs))
    fit <- withCallingHandlers(
      coxph(formula, data=x, cluster=icu_stay_id, ties="efron", singular.ok=FALSE,
            x=TRUE, y=TRUE, model=TRUE, na.action=na.fail, control=coxph.control(iter.max=50)),
      warning=function(w) {warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning")})
    effect_name <- "csHR"
  }
  if(any(grepl("infinite|converge", warning_text, ignore.case=TRUE)))
    stop(paste(metadata$dataset, metadata$analysis, metadata$model, paste(warning_text, collapse="; ")))
  stopifnot(all(is.finite(coef(fit))), all(is.finite(vcov(fit))))
  s <- summary(fit)
  meta <- c(metadata, list(estimator=kind, effect_measure=effect_name,
      n_eligible=nrow(data), unknown_icu_outcome=sum(is.na(data$icu_status)),
      n=nrow(x), excluded_unknown_or_missing=nrow(data)-nrow(x),
      icu_deaths=sum(x$icu_status==1L), live_discharges=sum(x$icu_status==2L),
      formula=paste(deparse(formula, width.cutoff=500), collapse=" "),
      warning=paste(warning_text, collapse="; ")))
  rows <- data.table(term=names(coef(fit)), log_effect=as.numeric(coef(fit)),
      se=sqrt(diag(vcov(fit))), effect=exp(coef(fit)),
      ci_low=exp(coef(fit)-qnorm(.975)*sqrt(diag(vcov(fit)))),
      ci_high=exp(coef(fit)+qnorm(.975)*sqrt(diag(vcov(fit)))),
      p_value=2*pnorm(-abs(coef(fit)/sqrt(diag(vcov(fit))))))
  for (nm in names(meta)) rows[, (nm) := meta[[nm]]]
  all_terms[[length(all_terms)+1L]] <<- rows
  manifest[[length(manifest)+1L]] <<- as.data.table(meta)
  if (ph) {
    z <- cox.zph(fit, transform=function(t) log1p(t), terms=TRUE, singledf=FALSE)$table
    z <- data.table(term=rownames(z), chisq=z[, "chisq"], df=z[, "df"], p_value=z[, "p"])
    for (nm in c("dataset", "analysis", "model")) z[, (nm) := metadata[[nm]]]
    z[, method := if(kind=="Fine-Gray")
        "Exploratory log(1+t) Schoenfeld diagnostic on weighted Fine-Gray pseudo-data" else
        "Cause-specific Cox log(1+t) Schoenfeld diagnostic"]
    diagnostics[[length(diagnostics)+1L]] <<- z
  }
  cat(kind, metadata$dataset, metadata$analysis, metadata$model,
      "N", nrow(x), "deaths", sum(x$icu_status==1L), "\n")
  list(fit=fit, data=x, terms=rows)
}

primary <- cause_specific <- cross_checks <- time_varying <- list()
specifications <- list(static_easix="log2_easix_z", static_siri="log2_siri_z",
                       static_joint=c("log2_easix_z", "log2_siri_z"), gbmtm_k3="trajectory_class")
for (db in c("mimiciv", "mimiciii")) for (model_name in names(specifications)) {
  trajectory <- model_name == "gbmtm_k3"
  origin <- if(trajectory) 7L else 1L
  flag <- if(trajectory) "trajectory_eligible" else "static_eligible"
  rhs <- paste(c(specifications[[model_name]], base), collapse=" + ")
  x <- copy(d[dataset==db & get(flag)==TRUE])
  meta <- list(dataset=db, analysis=if(trajectory) "day7_landmark" else "static_day1_24h_landmark",
               model=model_name, origin_days=origin, role="primary_association")
  result <- fit_model(x, rhs, meta, ph=TRUE)
  primary[[length(primary)+1L]] <- result$terms[grepl("^trajectory_class|^log2_", term)]
  cs <- fit_model(x, rhs, meta, kind="Cause-specific Cox", ph=TRUE)
  cause_specific[[length(cause_specific)+1L]] <- cs$terms[grepl("^trajectory_class|^log2_", term)]
  if (trajectory) {
    # A second established Fine-Gray implementation checks this candidate on
    # the identical trajectory sample; Breslow ties match cmprsk::crr.
    xx <- result$data
    cov1 <- model.matrix(as.formula(paste("~", rhs)), data=xx)[, -1, drop=FALSE]
    direct <- crr(xx$time, xx$icu_status, cov1=cov1, failcode=1, cencode=0,
                  gtol=1e-8, maxiter=50)
    stopifnot(direct$converged)
    delta <- max(abs(coef(result$fit)-direct$coef))
    stopifnot(delta < .001)
    cross_checks[[length(cross_checks)+1L]] <- data.table(dataset=db, n=nrow(xx),
        check="Fine-Gray coefficients: survival expansion versus cmprsk::crr",
        max_absolute_log_coefficient_difference=delta, tolerance=.001, status="PASS")
    class_columns <- grep("^trajectory_class", colnames(cov1))
    tv <- crr(xx$time, xx$icu_status, cov1=cov1, cov2=cov1[, class_columns, drop=FALSE],
              tf=function(t) cbind(log1p(t), log1p(t)), failcode=1, cencode=0,
              gtol=1e-7, maxiter=50)
    stopifnot(tv$converged, all(is.finite(tv$coef)), all(is.finite(tv$var)))
    index <- (length(tv$coef)-1L):length(tv$coef)
    statistic <- as.numeric(crossprod(tv$coef[index], solve(tv$var[index, index], tv$coef[index])))
    time_varying[[length(time_varying)+1L]] <- data.table(dataset=db, analysis="day7_landmark",
        model="Fine-Gray with both class-by-log(1+days since Day7) terms", n=nrow(xx),
        icu_deaths=sum(xx$icu_status==1L), chisq=statistic, df=2L,
        p_value=pchisq(statistic, 2L, lower.tail=FALSE),
        class2_baseline_log_shr=tv$coef["trajectory_class2"],
        class3_baseline_log_shr=tv$coef["trajectory_class3"],
        class2_logtime_coefficient=tv$coef[index[1]], class3_logtime_coefficient=tv$coef[index[2]])
  }
}

# CIFs estimate the actual competing-risk probabilities, not 1-KM.
cif_curves <- cif_tests <- cif_risk <- cif_estimates <- list()
cif_interval <- function(estimate, variance) {
  lower <- upper <- estimate
  middle <- estimate > 0 & estimate < 1
  z <- log(-log1p(-estimate[middle]))
  se <- sqrt(pmax(variance[middle], 0))/((1-estimate[middle]) * -log1p(-estimate[middle]))
  lower[middle] <- 1-exp(-exp(z-qnorm(.975)*se))
  upper[middle] <- 1-exp(-exp(z+qnorm(.975)*se))
  list(low=lower, high=upper)
}
for (db in c("mimiciv", "mimiciii")) {
  x <- copy(d[dataset==db & trajectory_eligible==TRUE & !is.na(icu_status)])
  x[, time := icu_endpoint_time_days-7]
  z <- cuminc(x$time, x$icu_status, group=x$trajectory_class, cencode=0)
  cif_tests[[length(cif_tests)+1L]] <- data.table(dataset=db, event=as.integer(rownames(z$Tests)),
      gray_chisq=z$Tests[, "stat"], df=z$Tests[, "df"], p_value=z$Tests[, "pv"],
      n=nrow(x), icu_deaths=sum(x$icu_status==1L), live_discharges=sum(x$icu_status==2L))
  for (cls in 1:3) {
    g <- x[trajectory_class==as.character(cls)]
    for (event in 1:2) {
      component <- z[[paste(cls, event)]]
      ci <- cif_interval(component$est, component$var)
      curve <- data.table(dataset=db, trajectory_class=cls, event=event,
          days_since_landmark=component$time, admission_day=component$time+7,
          cif=component$est, variance=component$var, ci_low=ci$low, ci_high=ci$high, n=nrow(g))
      cif_curves[[length(cif_curves)+1L]] <- curve
      for (day in c(14, 28, 90, Inf)) {
        index <- tail(which(component$time <= day-7), 1)
        cif_estimates[[length(cif_estimates)+1L]] <- data.table(dataset=db, trajectory_class=cls,
            event=event, horizon=if(is.infinite(day)) "Complete index ICU stay" else paste0("Admission Day ", day),
            cif=component$est[index], ci_low=ci$low[index], ci_high=ci$high[index],
            n=nrow(g), observed_events=sum(g$icu_status==event))
      }
      stopifnot(abs(tail(component$est, 1)-mean(g$icu_status==event)) < 1e-10)
    }
    for (day in c(7, 14, 28, 60, 90)) cif_risk[[length(cif_risk)+1L]] <- data.table(
      dataset=db, trajectory_class=cls, admission_day=day,
      event_free_icu_risk_n=sum(g$time >= day-7), cohort_n=nrow(g))
  }
}

# The seven existing clinical subgroups use the same Day 7 landmark and core
# adjustment. Joint robust Wald tests target the two class interaction terms.
definitions <- list(
  sepsis=list(variable="sepsis", labels=c("No sepsis", "Sepsis")),
  ventilation=list(variable="mechanical_ventilation_day1", labels=c("No ventilation", "Ventilation")),
  vasopressor=list(variable="vasopressor_day1", labels=c("No vasopressor", "Vasopressor")),
  crrt=list(variable="crrt_day1", labels=c("No CRRT", "CRRT")),
  aki=list(variable="aki_day1", labels=c("No AKI", "AKI")),
  age=list(variable="age_ge65", labels=c("Age <65", "Age >=65")),
  sex=list(variable="male", labels=c("Female", "Male")))
subgroup_effects <- subgroup_tests <- subgroup_status <- list()
for (db in c("mimiciv", "mimiciii")) for (sg in names(definitions)) {
  definition <- definitions[[sg]]
  data <- copy(d[dataset==db & trajectory_eligible==TRUE])
  data[, subgroup_value := get(definition$variable)]
  covariates <- if(sg=="sex") setdiff(base, "sex_group") else base
  design_rhs <- paste(c("trajectory_class", "subgroup_value", covariates), collapse=" + ")
  design_data <- sample_data(data, design_rhs, 7L)
  support <- merge(CJ(trajectory_class=as.character(1:3), subgroup_value=0:1),
      design_data[, .N, by=.(trajectory_class, subgroup_value)],
      by=c("trajectory_class", "subgroup_value"), all.x=TRUE)
  support[is.na(N), N:=0L]
  active <- as.integer(support[, .(both=all(N>0)), by=trajectory_class][both==TRUE, trajectory_class])
  stopifnot(length(active)>=2L)
  for(cls in active) design_data[, (paste0("subgroup_class", cls)) := subgroup_value*as.integer(trajectory_class==as.character(cls))]
  rhs <- paste(c("trajectory_class", paste0("subgroup_class", active), covariates), collapse=" + ")
  meta <- list(dataset=db, analysis="day7_landmark", model=paste0("subgroup_interaction_", sg),
               origin_days=7L, role="subgroup_interaction")
  r <- fit_model(design_data, rhs, meta)
  beta <- coef(r$fit); V <- vcov(r$fit)
  contrast <- matrix(0, nrow=length(active)-1L, ncol=length(beta), dimnames=list(NULL, names(beta)))
  for(j in seq_len(nrow(contrast))) {
    contrast[j, paste0("subgroup_class", active[j+1L])] <- 1
    contrast[j, paste0("subgroup_class", active[1L])] <- -1
  }
  delta <- contrast %*% beta; covariance <- contrast %*% V %*% t(contrast)
  stat <- as.numeric(crossprod(delta, solve(covariance, delta)))
  subgroup_tests[[length(subgroup_tests)+1L]] <- data.table(dataset=db, subgroup=sg,
      n=nrow(r$data), icu_deaths=sum(r$data$icu_status==1L), wald_chisq=stat, df=length(active)-1L,
      estimable_classes=paste(active, collapse=";"),
      interaction_p=pchisq(stat, length(active)-1L, lower.tail=FALSE))
  for (level in 0:1) {
    x <- data[subgroup_value==level]
    rhs <- paste(c("trajectory_class", covariates), collapse=" + ")
    xx <- sample_data(x, rhs, 7L)
    reason <- if(nrow(xx)<100L || sum(xx$icu_status==1L)<20L)
        "Below existing subgroup support threshold (N>=100, ICU deaths>=20)" else
        if(!all(as.character(1:3) %in% as.character(xx$trajectory_class)))
        "A trajectory class or Class 1 reference is absent" else ""
    factor_covariates <- covariates[vapply(xx[, ..covariates], is.factor, logical(1))]
    zero_event_levels <- unlist(lapply(factor_covariates, function(v)
      xx[, .(events=sum(icu_status==1L)), by=v][events==0L][[v]]), use.names=FALSE)
    if(reason=="" && length(zero_event_levels))
      reason <- paste("Zero ICU deaths in adjustment factor levels:", paste(zero_event_levels, collapse="; "))
    subgroup_status[[length(subgroup_status)+1L]] <- data.table(dataset=db, subgroup=sg,
        stratum=definition$labels[level+1L], n=nrow(xx), icu_deaths=sum(xx$icu_status==1L),
        status=if(reason=="") "estimated" else "not_estimable", reason=reason)
    if(reason!="") next
    meta$model <- paste0("subgroup_", sg, "_", level)
    meta$role <- "subgroup_stratum"
    result <- fit_model(x, rhs, meta)
    rows <- result$terms[grepl("^trajectory_class", term)]
    rows[, `:=`(subgroup=sg, stratum=definition$labels[level+1L])]
    subgroup_effects[[length(subgroup_effects)+1L]] <- rows
  }
}
subgroup_tests <- rbindlist(subgroup_tests)
subgroup_tests[, interaction_p_bh := p.adjust(interaction_p, method="BH"), by=dataset]

# Match the existing fixed-cutoff supplementary treatment comparisons. These
# retain the complete assigned 24 h cohort and its retrospective class scope.
specs <- data.table(exposure=c("nee_any_24h", "nee_high", "fluid_high", "fluid_kg_high",
                              "documented_mv_24h", "documented_crrt_24h"),
    contrast=c("Recorded vs no recorded NEE", ">0.3 vs 0<NEE<=0.3 microgram/kg/min",
               ">2000 vs <=2000 mL", ">30 vs <=30 mL/kg", "Recorded vs no recorded MV", "Recorded vs no recorded CRRT"))
treatment_effects <- treatment_tests <- treatment_cells <- list()
for (db in c("mimiciv", "mimiciii")) for (i in seq_len(nrow(specs))) {
  spec <- specs[i]
  x <- copy(d[dataset==db & treatment_eligible==TRUE])
  x[, nee_high := fifelse(nee_any_24h==1L & is.finite(nee_max_24h) & nee_max_24h>0,
                          as.integer(nee_max_24h>.3), NA_integer_)]
  io <- x$input_documented==1L & x$output_documented==1L
  x[, fluid_high := fifelse(io & is.finite(fluid_balance_ml_24h),
                             as.integer(fluid_balance_ml_24h>2000), NA_integer_)]
  x[, fluid_kg_high := fifelse(io & between(weight_kg, 20, 300) & is.finite(fluid_balance_ml_kg_24h),
                                as.integer(fluid_balance_ml_kg_24h>30), NA_integer_)]
  x[, exposure_value := as.integer(get(spec$exposure))]
  cov <- c("age", "sex_group", "race_group", "oasis_without_mv")
  xx <- sample_data(x, paste(c("trajectory_class", "exposure_value", cov), collapse=" + "), 1L)
  support <- merge(CJ(trajectory_class=as.character(1:3), exposure_value=0:1),
      xx[, .(n=.N, icu_deaths=sum(icu_status==1L), live_discharges=sum(icu_status==2L)),
         by=.(trajectory_class, exposure_value)], by=c("trajectory_class", "exposure_value"), all.x=TRUE)
  support[is.na(n), `:=`(n=0L, icu_deaths=0L, live_discharges=0L)]
  support[, `:=`(dataset=db, exposure=spec$exposure)]
  treatment_cells[[length(treatment_cells)+1L]] <- support
  active <- as.integer(support[, .(both_arms=all(n>0)), by=trajectory_class][both_arms==TRUE, trajectory_class])
  stopifnot(length(active)>=2L)
  for (cls in active) xx[, (paste0("exposure_class", cls)) := exposure_value*as.integer(trajectory_class==as.character(cls))]
  rhs <- paste(c("trajectory_class", paste0("exposure_class", active), cov), collapse=" + ")
  meta <- list(dataset=db, analysis="assigned_day1_treatment", model=spec$exposure, origin_days=1L,
               role="retrospective_treatment_association", source_n=nrow(x),
               source_unknown_icu_outcome=sum(is.na(x$icu_status)), excluded_before_design=nrow(x)-nrow(xx))
  r <- fit_model(xx, rhs, meta)
  beta <- coef(r$fit); V <- vcov(r$fit)
  contrast <- matrix(0, nrow=length(active)-1L, ncol=length(beta), dimnames=list(NULL, names(beta)))
  for (j in seq_len(nrow(contrast))) {
    contrast[j, paste0("exposure_class", active[j+1L])] <- 1
    contrast[j, paste0("exposure_class", active[1L])] <- -1
  }
  delta <- contrast %*% beta; variance <- contrast %*% V %*% t(contrast)
  statistic <- as.numeric(crossprod(delta, solve(variance, delta)))
  interaction <- pchisq(statistic, length(active)-1L, lower.tail=FALSE)
  treatment_tests[[length(treatment_tests)+1L]] <- data.table(dataset=db, exposure=spec$exposure,
      contrast=spec$contrast, n=nrow(xx), icu_deaths=sum(xx$icu_status==1L),
      wald_chisq=statistic, df=length(active)-1L, interaction_p=interaction,
      estimable_classes=paste(active, collapse=";"), adjustment="OASIS without MV points + age + sex + race/ethnicity")
  for (cls in 1:3) {
    rows <- if(cls %in% active) r$terms[term==paste0("exposure_class", cls)] else
        data.table(effect=NA_real_, ci_low=NA_real_, ci_high=NA_real_, p_value=NA_real_)
    cell <- support[trajectory_class==as.character(cls)]
    rows[, `:=`(dataset=db, exposure=spec$exposure, contrast=spec$contrast, trajectory_class=cls,
        n=nrow(xx), icu_deaths=sum(xx$icu_status==1L), interaction_p=interaction,
        low_n=cell[exposure_value==0, n], high_n=cell[exposure_value==1, n],
        low_deaths=cell[exposure_value==0, icu_deaths], high_deaths=cell[exposure_value==1, icu_deaths],
        status=if(cls %in% active) "estimated" else "not_estimable", sparse_cell=any(cell$n<10L | cell$icu_deaths<5L),
        reason=if(cls %in% active) "" else "No recorded exposure in this class")]
    treatment_effects[[length(treatment_effects)+1L]] <- rows
  }
}

outputs <- list(
  fine_gray_estimates=rbindlist(primary, fill=TRUE),
  cause_specific_cox_estimates=rbindlist(cause_specific, fill=TRUE),
  model_manifest=rbindlist(manifest, fill=TRUE), full_model_terms=rbindlist(all_terms, fill=TRUE),
  proportionality_diagnostics=rbindlist(diagnostics, fill=TRUE),
  fine_gray_implementation_checks=rbindlist(cross_checks), class_time_varying_diagnostics=rbindlist(time_varying),
  cif_curve_data=rbindlist(cif_curves), cif_gray_tests=rbindlist(cif_tests),
  cif_numbers_at_risk=rbindlist(cif_risk), cif_landmark_estimates=rbindlist(cif_estimates),
  subgroup_estimates=rbindlist(subgroup_effects, fill=TRUE), subgroup_interaction_tests=subgroup_tests,
  subgroup_model_status=rbindlist(subgroup_status),
  treatment_estimates=rbindlist(treatment_effects, fill=TRUE), treatment_interaction_tests=rbindlist(treatment_tests),
  treatment_support_cells=rbindlist(treatment_cells))
stopifnot(nrow(outputs$fine_gray_estimates)==12L, nrow(outputs$cause_specific_cox_estimates)==12L,
          nrow(outputs$subgroup_interaction_tests)==14L, nrow(outputs$treatment_interaction_tests)==12L,
          nrow(outputs$treatment_estimates)==36L)
for (nm in names(outputs)) fwrite(outputs[[nm]], file.path(out, paste0(nm, ".csv.tmp")))
for (nm in names(outputs)) stopifnot(file.rename(file.path(out, paste0(nm, ".csv.tmp")), file.path(out, paste0(nm, ".csv"))))
write_json(list(status="PASS", completed=as.character(Sys.time()),
    R=as.character(getRversion()), survival=as.character(packageVersion("survival")),
    cmprsk=as.character(packageVersion("cmprsk")),
    primary_models=8L, cause_specific_models=8L, subgroup_interaction_models=14L,
    treatment_models=12L, article_outputs_updated=FALSE),
    file.path(out, "model_execution.json"), pretty=TRUE, auto_unbox=TRUE)
print(outputs$fine_gray_estimates[, .(dataset, analysis, model, term, n, icu_deaths, effect, ci_low, ci_high, p_value)])
