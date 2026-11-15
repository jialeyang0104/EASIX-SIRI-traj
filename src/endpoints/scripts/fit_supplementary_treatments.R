#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(survival)})
root <- file.path(Sys.getenv("PROJECT_ROOT", unset=getwd()), "src/endpoints")
out <- file.path(root,"aggregate")
labels <- fread(file.path(root,"data_private/all_icu_gbmtm_assignments_mimiciv_mimiciii.csv"))[selected_k==3,
  .(dataset,icu_stay_id,gbmtm_class)]
treatment <- fread(file.path(root,"data_private/current_day1_treatments.csv"))
daily <- rbindlist(lapply(c("mimiciv","mimiciii"),function(db)
  fread(file.path(root,"data_private",paste0("all_icu_d1d7_daily_",db,".csv")))[day==1]))
d <- merge(merge(labels,daily,by=c("dataset","icu_stay_id")),treatment,by=c("dataset","icu_stay_id"))
stopifnot(nrow(d)==4025L,uniqueN(d,by=c("dataset","icu_stay_id"))==4025L)
d[,trajectory_class:=factor(gbmtm_class,levels=1:3)]
d[,sex_group:=factor(sex)]
race<-toupper(trimws(d$race_ethnicity))
d[,race_group:=factor(fcase(grepl("WHITE",race),"White",grepl("BLACK|AFRICAN",race),"Black",
  grepl("HISPANIC|LATINO",race),"Hispanic",grepl("ASIAN",race),"Asian",
  is.na(race)|race==""|grepl("UNKNOWN|UNABLE|DECLINED",race),"Unknown",default="Other"))]
d[,nee_high:=fifelse(nee_any_24h==1 & is.finite(nee_max_24h) & nee_max_24h>0,as.integer(nee_max_24h>0.3),NA_integer_)]
io<-d$input_documented==1 & d$output_documented==1
d[,fluid_high:=fifelse(io & is.finite(fluid_balance_ml_24h),as.integer(fluid_balance_ml_24h>2000),NA_integer_)]
d[,fluid_kg_high:=fifelse(io & between(weight_kg,20,300) & is.finite(fluid_balance_ml_kg_24h),
  as.integer(fluid_balance_ml_kg_24h>30),NA_integer_)]
specs<-data.table(exposure=c("nee_any_24h","nee_high","fluid_high","fluid_kg_high","documented_mv_24h","documented_crrt_24h"),
  exposure_order=0:5,exposure_label=c("Any NEE","Maximum NEE among exposed","Net fluid balance","Net fluid balance per kg","Mechanical ventilation","CRRT"),
  contrast_label=c("Recorded vs no recorded exposure",">0.3 vs 0<NEE≤0.3 µg/kg/min",">2,000 vs ≤2,000 mL",">30 vs ≤30 mL/kg",
    "Recorded vs no recorded receipt","Recorded vs no recorded receipt"),cutoff=c(NA,.3,2000,30,NA,NA))
estimates<-tests<-cells<-ph_rows<-list()
for(db in c("mimiciv","mimiciii")) for(h in c(28L,90L)) for(i in seq_len(nrow(specs))) {
  spec<-specs[i];x<-copy(d[dataset==db]);x[,exposure_value:=as.integer(get(spec$exposure))]
  x[,time:=get(paste0("followup_",h,"d_days"))-1];x[,event:=get(paste0("event_",h,"d"))]
  required<-c("time","event","exposure_value","age","sex_group","race_group","oasis_without_mv","trajectory_class")
  x<-x[complete.cases(x[,..required]) & time>0]
  for(v in c("sex_group","race_group","trajectory_class")) x[,(v):=droplevels(get(v))]
  support<-merge(CJ(trajectory_class=as.character(1:3),exposure_value=0:1),
    x[,.(n=.N,events=sum(event)),by=.(trajectory_class,exposure_value)],by=c("trajectory_class","exposure_value"),all.x=TRUE)
  support[is.na(n),`:=`(n=0L,events=0L)]
  support[,`:=`(dataset=db,horizon=h,exposure=spec$exposure)]
  cells[[length(cells)+1L]]<-support
  meta<-c(list(dataset=db,horizon=h,model_n=nrow(x),model_events=sum(x$event),time_origin_days=1L,
    encoding="cutoff_binary",adjustment="age + sex + race/ethnicity + OASIS without mechanical-ventilation points",
    source_analysis="current_24h_landmark_cutoffs",interval="95% ICU-stay-clustered robust confidence interval"),as.list(spec))
  active<-as.integer(support[,.(both_arms=all(n>0)),by=trajectory_class][both_arms==TRUE,trajectory_class])
  for(cls in active)x[,(paste0("exposure_class",cls)):=exposure_value*as.integer(trajectory_class==as.character(cls))]
  formula<-as.formula(paste("Surv(time,event) ~ trajectory_class +",paste(c(paste0("exposure_class",active),
    "age","sex_group","race_group","oasis_without_mv","cluster(icu_stay_id)"),collapse=" + ")))
  warning_text<-character()
  cat("Fitting",db,h,spec$exposure,"N",nrow(x),"\n")
  fit<-withCallingHandlers(coxph(formula,data=x,ties="efron",x=TRUE,y=TRUE,model=TRUE,na.action=na.fail,singular.ok=FALSE),
    warning=function(w){warning_text<<-c(warning_text,conditionMessage(w));invokeRestart("muffleWarning")})
  beta<-coef(fit);V<-vcov(fit)
  stopifnot(all(is.finite(beta)),all(is.finite(V)),!any(grepl("infinite|converge",warning_text,ignore.case=TRUE)))
  stopifnot(length(active)>=2L)
  contrast<-matrix(0,nrow=length(active)-1L,ncol=length(beta),dimnames=list(NULL,names(beta)))
  for(j in seq_len(nrow(contrast))){contrast[j,paste0("exposure_class",active[j+1L])]<-1;contrast[j,paste0("exposure_class",active[1])]<--1}
  difference<-contrast%*%beta;covariance<-contrast%*%V%*%t(contrast)
  chi<-as.numeric(crossprod(difference,solve(covariance,difference)))
  df<-length(active)-1L;global_p<-pchisq(chi,df,lower.tail=FALSE)
  tests[[length(tests)+1L]]<-as.data.table(c(meta,list(global_interaction_p=global_p,interaction_df=df,interaction_chisq=chi,
    estimable_classes=paste(active,collapse=";"),
    source_path="src/endpoints/aggregate/supp_treatment_interaction_tests.csv",
    formula=paste(deparse(formula,width.cutoff=500),collapse=" "),warning=paste(warning_text,collapse="; "))))
  for(cls in 1:3) {
    loghr<-se<-NA_real_
    if(cls %in% active){term<-paste0("exposure_class",cls);loghr<-beta[term];se<-sqrt(V[term,term])}
    arm<-support[trajectory_class==as.character(cls)]
    estimates[[length(estimates)+1L]]<-as.data.table(c(meta,list(global_interaction_p=global_p,interaction_df=df,
      trajectory_class=cls,hazard_ratio=exp(loghr),ci_low=exp(loghr-qnorm(.975)*se),ci_high=exp(loghr+qnorm(.975)*se),
      p_value=2*pnorm(-abs(loghr/se)),low_n=arm[exposure_value==0,n],high_n=arm[exposure_value==1,n],
      low_events=arm[exposure_value==0,events],high_events=arm[exposure_value==1,events],
      model_status=if(cls %in% active) "estimated" else "not_estimable",reason=if(cls %in% active) "" else "No recorded exposure in this class",
      sparse_cell=any(arm$n<10L | arm$events<5L),
      source_path="src/endpoints/aggregate/supp_treatment_forest_estimates.csv")))
  }
  ph<-cox.zph(fit,transform="km",terms=TRUE,singledf=FALSE)$table
  ph_rows[[length(ph_rows)+1L]]<-data.table(dataset=db,horizon=h,exposure=spec$exposure,term=rownames(ph),
    chisq=ph[,"chisq"],df=ph[,"df"],p_value=ph[,"p"])
}
outputs<-list(supp_treatment_forest_estimates=rbindlist(estimates,fill=TRUE),
  supp_treatment_interaction_tests=rbindlist(tests,fill=TRUE),supp_treatment_support_cells=rbindlist(cells),
  supp_treatment_ph_diagnostics=rbindlist(ph_rows))
stopifnot(nrow(outputs$supp_treatment_forest_estimates)==72L,nrow(outputs$supp_treatment_interaction_tests)==24L)
for(name in names(outputs))fwrite(outputs[[name]],file.path(out,paste0(name,".csv.tmp")))
for(name in names(outputs))stopifnot(file.rename(file.path(out,paste0(name,".csv.tmp")),file.path(out,paste0(name,".csv"))))
print(outputs$supp_treatment_interaction_tests[,.(dataset,horizon,exposure,model_n,model_events,global_interaction_p)])
