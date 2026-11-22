#!/usr/bin/env Rscript
# Refit the Figure S6 fixed-class and Figure S7 joint probability models.
suppressPackageStartupMessages({library(data.table); library(survival)})
root <- file.path(Sys.getenv("PROJECT_ROOT", unset=getwd()), "src")
out <- file.path(root, "treatment/aggregate")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
key <- c("dataset","icu_stay_id")
cov <- c("age","sex_group","race_group","oasis_without_mv")
d <- fread(file.path(root,"icu/data_private/icu_mortality_analysis.csv"))[treatment_eligible==TRUE]
p <- fread(file.path(root,"treatment/data_private/day1_probabilities.csv"), select=c(key,"prob_class_2","prob_class_3"))
stopifnot(nrow(d)==4025L,nrow(p)==4025L,uniqueN(d,by=key)==4025L,uniqueN(p,by=key)==4025L)
d <- merge(d,p,by=key)
stopifnot(nrow(d)==4025L,all(is.finite(d$prob_class_2)),all(is.finite(d$prob_class_3)),
          all(d$prob_class_2>=0 & d$prob_class_3>=0 & d$prob_class_2+d$prob_class_3<=1+1e-6))
d[, fluid_ok := input_documented %in% 1 & output_documented %in% 1 &
    is.finite(fluid_balance_ml_kg_24h) & between(weight_kg,20,300)]
d[, nee_ok := nee_any_24h %in% 1 & is.finite(nee_max_24h) & nee_max_24h>0]
d[, exposure_any_nee := as.integer(nee_any_24h)]
d[, exposure_nee_dose := fifelse(nee_ok,as.integer(nee_max_24h>.30),NA_integer_)]
d[, exposure_net_fluid := fifelse(fluid_ok,as.integer(fluid_balance_ml_kg_24h>35),NA_integer_)]
d[, exposure_crrt := as.integer(documented_crrt_24h)]
d[, c("p2","p3") := list(10*prob_class_2,10*prob_class_3)]
# Match the published-source reconstruction: original 28-day race/sex coding;
# collapsed Other/Unknown and competing-risk time/status for the ICU model.
race <- toupper(trimws(d$race_ethnicity))
d[, race_28 := fcase(grepl("WHITE",race),"White",grepl("BLACK|AFRICAN",race),"Black",
   grepl("HISPANIC|LATINO",race),"Hispanic",grepl("ASIAN",race),"Asian",
   is.na(race)|race==""|grepl("UNKNOWN|UNABLE|DECLINED",race),"Unknown",default="Other")]
sample_at <- function(db, ep, name) {
  x <- copy(d[dataset==db & !is.na(get(name))])
  x[, exposure := get(name)]
  stopifnot(all(x$exposure %in% 0:1))
  if (ep=="ICU") {
    x[race_group %in% c("Other","Unknown"), race_group := "Other/Unknown"]
    x[, c("time","status") := list(icu_endpoint_time_days-1,icu_status)]
  } else x[, c("race_group","sex_group","time","status") := list(race_28,sex,followup_28d_days-1,event_28d)]
  needed <- c(cov,"time","status","gbmtm_class")
  x <- x[complete.cases(x[,..needed]) & is.finite(time) & time>0 &
         status %in% (if(ep=="ICU") 0:2 else 0:1)]
  x[, c("sex_group","race_group","state") :=
      list(factor(sex_group),factor(race_group),
           factor(status,levels=0:2,labels=c("Censored","ICU death","Live discharge")))]
  x
}

fit_endpoint <- function(x, ep, terms) {
  rhs <- paste(c(terms,cov),collapse="+")
  opts <- coxph.control(iter.max=50)
  withCallingHandlers({
    if (ep=="ICU") {
      fg <- finegray(as.formula(paste("Surv(time,state)~",rhs,"+icu_stay_id")),x,
                     id=icu_stay_id,etype="ICU death",na.action=na.fail)
      fit <- coxph(as.formula(paste("Surv(fgstart,fgstop,fgstatus)~",rhs)),fg,
                   weights=fgwt,cluster=icu_stay_id,ties="breslow",na.action=na.fail,singular.ok=FALSE,control=opts)
    } else fit <- coxph(as.formula(paste("Surv(time,status)~",rhs)),x,
                         cluster=icu_stay_id,ties="efron",na.action=na.fail,singular.ok=FALSE,control=opts)
    if (!all(is.finite(coef(fit))) || !all(is.finite(vcov(fit)))) stop("Non-finite model estimates")
    fit
  }, warning=function(w) {
    if (grepl("converg|infinite|singular",conditionMessage(w),ignore.case=TRUE)) stop(conditionMessage(w))
  })
}

fit_fixed <- function(x, ep, exposure_name, cutoff) {
  # Encode one exposure coefficient per supported class. A class with
  # no exposed patients still contributes to the shared adjustment fit.
  x[, klass := factor(gbmtm_class, levels=1:3)]
  arms <- lapply(1:3, function(k) {
    g <- x[gbmtm_class==k]
    c(tabulate(g$exposure+1L, nbins=2L),
      tabulate(g$exposure[g$status==1]+1L, nbins=2L))
  })
  active <- which(vapply(arms, function(a) all(a[1:2]>0L), logical(1)))
  terms <- paste0("arm_class", active)
  for (j in seq_along(active))
    x[, (terms[j]) := exposure*as.integer(gbmtm_class==active[j])]
  fit <- fit_endpoint(x, ep, c("klass",terms))
  b <- coef(fit)[terms]
  V <- vcov(fit)[terms,terms,drop=FALSE]
  L <- cbind(-1,diag(length(active)-1L))
  delta <- as.vector(L %*% b)
  pv <- pchisq(as.numeric(crossprod(delta,solve(L %*% V %*% t(L),delta))),
               df=nrow(L), lower.tail=FALSE)
  se <- sqrt(diag(V))
  rbindlist(lapply(1:3,function(k) {
    j <- match(k,active); a <- arms[[k]]
    data.table(exposure=exposure_name,cutoff=cutoff,trajectory_class=k,
               ratio=exp(b[j]),ci_low=exp(b[j]-qnorm(.975)*se[j]),
               ci_high=exp(b[j]+qnorm(.975)*se[j]),p_value=2*pnorm(-abs(b[j]/se[j])),
               interaction_p=pv,low_n=a[1],high_n=a[2],low_events=a[3],high_events=a[4],
               status=if(is.na(j)) "not_estimable" else "estimated")
  }))
}

exps <- data.table(
  exposure=c("exposure_any_nee","exposure_nee_dose","exposure_net_fluid","exposure_crrt"),
  exposure_label=c("Any recorded NEE","NEE >0.30 vs 0<NEE≤0.30 µg/kg/min",
                   "Net fluid balance >35 vs ≤35 mL/kg","Recorded CRRT"),
  cutoff=c(NA_real_,.30,35,NA_real_))
curves <- list(); tests <- list(); histograms <- list(); fixed <- list()
z <- qnorm(.975)
edges <- seq(0,1,length.out=21)
for (i in seq_len(nrow(exps))) for (db in c("mimiciv","mimiciii")) for (ep in c("28-day","ICU")) {
  e <- exps[i]
  x <- sample_at(db,ep,e$exposure)
  fixed_name <- c("nee_any_24h","nee_dose","net_fluid","documented_crrt_24h")[i]
  fixed[[length(fixed)+1L]] <- fit_fixed(copy(x), ep, fixed_name, e$cutoff)[, c("dataset","outcome") := list(db,ep)]
  n <- tabulate(x$exposure+1L,nbins=2L); ev <- tabulate(x$exposure[x$status==1]+1L,nbins=2L)
  id <- c(as.list(e),list(dataset=db,outcome=ep))
  cell <- as.data.table(c(id,list(n=nrow(x),events=sum(x$status==1),low_n=n[1],high_n=n[2],
    low_events=ev[1],high_events=ev[2],status="estimated",reason="",interaction_statistic=NA_real_,
    interaction_df=2L,interaction_p=NA_real_,interaction_p2=NA_real_,interaction_p3=NA_real_,
    ratio=NA_real_,ci_low=NA_real_,ci_high=NA_real_,slope_difference_p=NA_real_)))

  # Keep observed marginal distributions even when model support is insufficient.
  for (cl in 2:3) for (arm in 0:1) {
    pp <- x[[paste0("prob_class_",cl)]][x$exposure==arm]
    count <- tabulate(findInterval(pp,edges,rightmost.closed=TRUE),nbins=length(edges)-1L)
    stopifnot(sum(count)==n[arm+1L])
    histograms[[length(histograms)+1L]] <- as.data.table(c(id,list(predicted_class=cl,exposure_arm=arm,
      bin_lower=head(edges,-1),bin_upper=tail(edges,-1),count=count,arm_n=n[arm+1L])))
  }

  stats <- tryCatch({
    if (any(n<30L | ev<10L | n-ev<10L)) stop("Arm support below 30 participants, 10 deaths or 10 non-deaths")
    fit <- fit_endpoint(x,ep,"exposure*(p2+p3)")
    keys <- c("exposure","exposure:p2","exposure:p3")
    b <- coef(fit)[keys]; V <- vcov(fit)[keys,keys,drop=FALSE]
    if (!all(is.finite(b)) || !all(is.finite(V)) ||
        min(eigen(V,symmetric=TRUE,only.values=TRUE)$values)<=0)
      stop("Singular or non-finite robust contrast covariance")
    bb <- b[2:3]; VV <- V[2:3,2:3,drop=FALSE]
    joint <- as.numeric(crossprod(bb,solve(VV,bb)))
    L <- matrix(c(0,-1,1),nrow=1)
    shift <- as.numeric(L%*%b); se <- sqrt(as.numeric(L%*%V%*%t(L)))
    list(b=b,V=V,joint=joint,component_p=2*pnorm(-abs(bb/sqrt(diag(VV)))),
         shift=shift,se=se,shift_p=2*pnorm(-abs(shift/se)))
  },error=identity)
  if (inherits(stats,"error")) {
    cell[, c("status","reason") := list("not_estimable",conditionMessage(stats))]
  } else {
    s <- stats
    cell[, c("interaction_statistic","interaction_p","interaction_p2","interaction_p3",
             "ratio","ci_low","ci_high","slope_difference_p") :=
      list(s$joint,pchisq(s$joint,df=2,lower.tail=FALSE),s$component_p[1],s$component_p[2],
           exp(s$shift),exp(s$shift-z*s$se),exp(s$shift+z*s$se),s$shift_p)]
    for (cl in 2:3) {
      pp <- x[[paste0("prob_class_",cl)]]
      support <- sapply(0:1,function(arm) quantile(pp[x$exposure==arm],c(.025,.975)))
      lo <- max(support[1,]); hi <- min(support[2,])
      if (lo>=hi) next
      grid <- seq(lo,hi,length.out=101)
      # Vary the selected class on its Class-1 edge; fix the other class to zero.
      X <- cbind(1,if(cl==2L) 10*grid else 0,if(cl==3L) 10*grid else 0)
      est <- as.vector(X%*%s$b); se <- sqrt(rowSums((X%*%s$V)*X))
      curves[[length(curves)+1L]] <- as.data.table(c(id,list(predicted_class=cl,probability=grid,
        ratio=exp(est),ci_low=exp(est-z*se),ci_high=exp(est+z*se),interaction_p=cell$interaction_p,
        slope_difference_p=cell$slope_difference_p,n=cell$n,events=cell$events)))
    }
  }
  tests[[length(tests)+1L]] <- cell
}
tests <- rbindlist(tests); curves <- rbindlist(curves); histograms <- rbindlist(histograms)
stopifnot(nrow(tests)==16L,all(tests$low_n+tests$high_n==tests$n),
          all(tests$low_events+tests$high_events==tests$events),nrow(histograms)==16L*2L*2L*20L,
          all(is.finite(curves$ratio) & is.finite(curves$ci_low) & is.finite(curves$ci_high)),
          all(curves$ci_low<=curves$ratio & curves$ratio<=curves$ci_high))
fixed <- rbindlist(fixed)
stopifnot(nrow(fixed)==48L,uniqueN(fixed,by=c("dataset","outcome","exposure","trajectory_class"))==48L)
tables <- list(figure_s6_fixed=fixed,figure_s7_four_treatment_curves=curves,
               figure_s7_four_treatment_tests=tests,figure_s7_four_treatment_histograms=histograms)
for (name in names(tables)) {
  tmp <- tempfile(tmpdir=out)
  fwrite(tables[[name]],tmp)
  stopifnot(file.rename(tmp,file.path(out,paste0(name,".csv"))))
}
cat("Wrote 48 fixed-class rows and 16 probability-model tests to the local treatment aggregates.\n")
