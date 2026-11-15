suppressPackageStartupMessages(library(data.table))
project_root <- normalizePath(getwd())
root <- file.path(project_root,"src","endpoints")
fixed <- file.path(project_root,"results","trajectory")
model <- readRDS(file.path(fixed,"models","gbmtm_all_icu_mimiciv_d1d7_easix_siri_selected.rds"))
params <- fread(file.path(fixed,"tables","all_icu_gbmtm_preprocessing.csv"))
bounds <- params[,.(variable,winsor_lower,winsor_upper)]
scaling <- params[,.(variable=interpolated_variable,mean=standardization_mean,sd=standardization_sd)]
old <- fread(file.path(fixed,"tables","all_icu_gbmtm_assignments_mimiciv_mimiciii.csv"))
old[,icu_stay_id:=as.character(icu_stay_id)]
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


assignments <- eligibility <- profiles <- comparisons <- list()
for (db in c("mimiciv","mimiciii")) {
  raw <- fread(file.path(root,"data_private",paste0("all_icu_d1d7_daily_",db,".csv")))
  raw[,icu_stay_id:=as.character(icu_stay_id)]
  d <- prepare_validation(raw,bounds,scaling)
  md <- as.data.frame(d$data[,.(id=paste0(db,":",icu_stay_id),time=as.integer(day),log2_easix_interp_z,log2_siri_interp_z)])
  z <- fixed_posterior(model,md)
  z[,`:=`(dataset=db,icu_stay_id=sub("^[^:]+:","",id),selected_k=3L,selection_status="fixed K3 quadratic pruned")]
  previous <- old[dataset==db]
  common <- merge(z,previous,by=c("dataset","icu_stay_id"),suffixes=c("_new","_old"))
  stopifnot(nrow(common)==nrow(previous),all(common$gbmtm_class_new==common$gbmtm_class_old))
  comparisons[[db]] <- data.table(dataset=db,previous=nrow(previous),current=nrow(z),added=nrow(z)-nrow(previous),same_class_fraction=mean(common$gbmtm_class_new==common$gbmtm_class_old),degree=as.integer(model$call$d),pruning=as.logical(model$call$pruning))
  assignments[[db]] <- z
  eligibility[[db]] <- d$eligibility[,dataset:=db]
  wide <- merge(d$data,z[,.(icu_stay_id,gbmtm_class)],by="icu_stay_id")
  profile <- melt(wide,id.vars=c("icu_stay_id","day","gbmtm_class"),measure.vars=c("log2_easix_interp_z","log2_siri_interp_z"),variable.name="index",value.name="z")
  profiles[[db]] <- profile[,.(mean_z=mean(z),se_z=sd(z)/sqrt(.N),n_stays=uniqueN(icu_stay_id)),by=.(gbmtm_class,day,index)][,dataset:=db]
}
z <- rbindlist(assignments,fill=TRUE)
fwrite(z,file.path(root,"data_private","all_icu_gbmtm_assignments_mimiciv_mimiciii.csv"))
fwrite(rbindlist(eligibility),file.path(root,"data_private","all_icu_d1d7_dual_index_eligibility.csv"))
fwrite(z[,.(n=.N,mean_max_post_prob=mean(max_post_prob)),by=.(dataset,gbmtm_class)][,pct:=100*n/sum(n),by=dataset],file.path(root,"aggregate","all_icu_gbmtm_class_distribution.csv"))
fwrite(rbindlist(profiles)[,`:=`(ci_low=mean_z-1.96*se_z,ci_high=mean_z+1.96*se_z)],file.path(root,"aggregate","all_icu_gbmtm_trajectory_profiles.csv"))
fwrite(rbindlist(comparisons),file.path(root,"aggregate","fixed_projection_reconciliation.csv"))
print(rbindlist(comparisons))
