# Shared subgroup implementation for the full analysis and the focused refit.
# Continuous age is retained in age strata and binary age interactions.
fit_subgroups <- function(traj_dev, traj_val, tables) {
  base_covariates <- c("age", "sex", "race_group", "oasis_day1")
  subgroup_defs <- list(
    sepsis = list(variable="sepsis", labels=c("No sepsis", "Sepsis"), remove=character()),
    ventilation = list(variable="mechanical_ventilation_day1", labels=c("No ventilation", "Ventilation"), remove=character()),
    vasopressor = list(variable="vasopressor_day1", labels=c("No vasopressor", "Vasopressor"), remove=character()),
    crrt = list(variable="crrt_day1", labels=c("No CRRT", "CRRT"), remove=character()),
    aki = list(variable="aki_day1", labels=c("No AKI", "AKI"), remove=character()),
    age = list(variable="age_ge65", labels=c("Age <65", "Age >=65"), remove=character()),
    sex = list(variable="male", labels=c("Female", "Male"), remove="sex")
  )
  terms <- function(fit, metadata) {
    z <- summary(fit)
    out <- data.table(term=rownames(z$coefficients), estimate=z$coefficients[, "coef"],
      se=z$coefficients[, "se(coef)"], hazard_ratio=z$conf.int[, "exp(coef)"],
      ci_low=z$conf.int[, "lower .95"], ci_high=z$conf.int[, "upper .95"],
      p_value=z$coefficients[, "Pr(>|z|)"])
    for (nm in names(metadata)) out[[nm]] <- metadata[[nm]]
    setcolorder(out, c(names(metadata), setdiff(names(out), names(metadata))))
    out
  }
  formula_text <- function(f) paste(deparse(f, width.cutoff=500L), collapse=" ")
  effects <- interactions <- full_terms <- fit_manifest <- continuous <- list()
  formulas <- rbindlist(lapply(names(subgroup_defs), function(nm) {
    cov <- setdiff(base_covariates, subgroup_defs[[nm]]$remove)
    data.table(subgroup=nm, omitted_covariate=paste(subgroup_defs[[nm]]$remove, collapse=""),
      stratum_formula=paste("Surv(analysis_time, event_horizon) ~", paste(c("trajectory_class", cov), collapse=" + ")),
      interaction_formula=paste("Surv(analysis_time, event_horizon) ~ trajectory_class * subgroup_value +", paste(cov, collapse=" + ")),
      bh_family="seven tests per database and origin")
  }))
  for (dataset_name in c("mimiciv", "mimiciii")) {
    for (horizon in c(28L,90L)) {
    for (analysis_name in c("icu_admission_inception", "day7_landmark")) {
      d0 <- copy(if (dataset_name == "mimiciv") traj_dev else traj_val)
      d0[, event_horizon := get(paste0("event_",horizon,"d"))]
      d0[, followup_horizon := get(paste0("followup_",horizon,"d_days"))]
      if (analysis_name == "day7_landmark") {
        d0 <- d0[icu_los_days >= 7 & mortality_7d == 0 & followup_horizon > 7]
        d0[, analysis_time := followup_horizon - 7]
      } else d0[, analysis_time := followup_horizon]
      for (nm in names(subgroup_defs)) {
        definition <- subgroup_defs[[nm]]
        d <- copy(d0)
        d[, subgroup_value := as.integer(get(definition$variable))]
        cov <- setdiff(base_covariates, definition$remove)
        reduced <- as.formula(paste("Surv(analysis_time, event_horizon) ~", paste(c("trajectory_class", "subgroup_value", cov), collapse=" + ")))
        full <- as.formula(paste("Surv(analysis_time, event_horizon) ~ trajectory_class * subgroup_value +", paste(cov, collapse=" + ")))
        required <- unique(c(all.vars(reduced), all.vars(full)))
        dd <- d[complete.cases(d[, ..required]) & analysis_time > 0]
        stopifnot(nrow(dd) >= 200, uniqueN(dd$subgroup_value) == 2, sum(dd$event_horizon) >= 30)
        reduced_fit <- coxph(reduced, data=dd, ties="efron", na.action=na.fail)
        full_fit <- coxph(full, data=dd, ties="efron", na.action=na.fail)
        lrt <- anova(reduced_fit, full_fit, test="LRT")
        interactions[[length(interactions)+1L]] <- data.table(dataset=dataset_name,
          analysis=analysis_name, horizon=paste0(horizon,"-day"), subgroup=nm, n=nrow(dd), events=sum(dd$event_horizon),
          interaction_p=lrt$`Pr(>|Chi|)`[2])
        fit_manifest[[length(fit_manifest)+1L]] <- data.table(dataset=dataset_name,
          analysis=analysis_name, horizon=paste0(horizon,"-day"), subgroup=nm, n_eligible=nrow(d),
          n_excluded_missing_or_nonpositive=nrow(d)-nrow(dd), n=nrow(dd), events=sum(dd$event_horizon),
          reduced_n=reduced_fit$n, full_n=full_fit$n, reduced_formula=formula_text(reduced),
          full_formula=formula_text(full), reduced_loglik=as.numeric(logLik(reduced_fit)),
          full_loglik=as.numeric(logLik(full_fit)), lrt_df=lrt$Df[2],
          interaction_p=lrt$`Pr(>|Chi|)`[2])
        for (level in 0:1) {
          stratum <- dd[subgroup_value == level]
          if (nrow(stratum) < 100 || sum(stratum$event_horizon) < 20) next
          formula <- as.formula(formulas[subgroup == nm, stratum_formula])
          fit <- coxph(formula, data=stratum, ties="efron", na.action=na.fail)
          z <- terms(fit, list(dataset=dataset_name, analysis=analysis_name, horizon=paste0(horizon,"-day"), subgroup=nm,
            stratum=definition$labels[level+1], n=nrow(stratum), events=sum(stratum$event_horizon),
            formula=formula_text(formula)))
          full_terms[[length(full_terms)+1L]] <- z
          effects[[length(effects)+1L]] <- z[grepl("^trajectory_class", term)]
        }
      }
      # Separate exploratory continuous-age interaction; outside the seven binary tests.
      reduced <- Surv(analysis_time, event_horizon) ~ trajectory_class + age + sex + race_group + oasis_day1
      full <- Surv(analysis_time, event_horizon) ~ trajectory_class * age + sex + race_group + oasis_day1
      required <- all.vars(full)
      dd <- d0[complete.cases(d0[, ..required]) & analysis_time > 0]
      fit0 <- coxph(reduced, data=dd, ties="efron", na.action=na.fail)
      fit1 <- coxph(full, data=dd, ties="efron", na.action=na.fail)
      lrt <- anova(fit0, fit1, test="LRT")
      continuous[[length(continuous)+1L]] <- data.table(dataset=dataset_name, analysis=analysis_name, horizon=paste0(horizon,"-day"),
        n=nrow(dd), events=sum(dd$event_horizon), reduced_formula=formula_text(reduced),
        full_formula=formula_text(full), lrt_df=lrt$Df[2], interaction_p=lrt$`Pr(>|Chi|)`[2])
    }
  }
  }
  interactions <- rbindlist(interactions)
  interactions[, interaction_p_bh := p.adjust(interaction_p, method="BH"), by=.(dataset, analysis, horizon)]
  stopifnot(all(interactions[, .N, by=.(dataset, analysis, horizon)]$N == 7L),
    all(is.finite(interactions$interaction_p)), all(is.finite(interactions$interaction_p_bh)))
  fit_manifest <- merge(rbindlist(fit_manifest),
    interactions[, .(dataset, analysis, horizon, subgroup, interaction_p_bh)],
    by=c("dataset", "analysis", "horizon", "subgroup"), sort=FALSE)
  outputs <- list(
    all_icu_cox_gbmtm_subgroup_interactions.csv=interactions,
    all_icu_cox_gbmtm_subgroup_stratified_effects.csv=rbindlist(effects),
    subgroup_formula_manifest.csv=formulas,
    subgroup_fit_manifest.csv=fit_manifest,
    all_icu_cox_gbmtm_continuous_age_interactions.csv=rbindlist(continuous))
  stage <- tempfile("subgroup-results-")
  dir.create(stage)
  on.exit(unlink(stage, recursive=TRUE))
  for (nm in names(outputs)) fwrite(outputs[[nm]], file.path(stage, nm))
  for (nm in names(outputs)) stopifnot(file.copy(file.path(stage, nm), file.path(tables, nm), overwrite=TRUE))
  rbindlist(full_terms, fill=TRUE)
}
