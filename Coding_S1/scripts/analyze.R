# Adaptable analysis examples: replace datasets and model settings for your study.
# No correlation-value sensitivity analyses are added.
library(metafor)
read_data <- function(f) read.csv(file.path("data",f), na.strings=c("","NA"), check.names=FALSE)
summaries <- list()
save_fit <- function(fit,id,scale) {
  pr <- predict(fit)
  tau <- if(inherits(fit,"rma.mv")) sum(fit$sigma2) else fit$tau2
  summaries[[id]] <<- data.frame(analysis=id,scale=scale,k=fit$k,estimate=as.numeric(coef(fit)[1]),
    ci_lb=fit$ci.lb[1],ci_ub=fit$ci.ub[1],p=fit$pval[1],tau2_total=tau,
    QE=fit$QE,QEp=fit$QEp,prediction_lb=pr$pi.lb[1],prediction_ub=pr$pi.ub[1])
  pdf(file.path("outputs",paste0(id,"_forest.pdf")),width=10,height=7)
  forest(fit, xlab=scale)
  dev.off()
}
motor_effect <- function(d,comparison) {
  md <- if(comparison=="stn") "Within_mean_co_STN" else "Within_mean_co_baseline"
  se <- if(comparison=="stn") "Within_SE_co_STN" else "Within_SE_co_baseline"
  md_alt <- if(comparison=="stn") "Mean_Difference_CO_STN" else "Mean_Difference_CO_base"
  se_alt <- if(comparison=="stn") "SE_CO_STN" else "SE_CO_base"
  paired <- !is.na(d[[md]]) & !is.na(d[[se]])
  d$yi <- ifelse(paired,d[[md]],d[[md_alt]])
  d$sei <- ifelse(paired,d[[se]],d[[se_alt]])
  d$effect_source <- ifelse(paired,"retained within-subject estimate","retained alternative estimate")
  d[d$Assessment=="UPDRS_3" & is.finite(d$yi) & is.finite(d$sei) & d$sei>0,,drop=FALSE]
}
loo <- function(d,id,test="knha") {
  studies <- unique(sub("_[12]$","",d$Research))
  do.call(rbind,lapply(studies,function(s) {
    dd <- d[sub("_[12]$","",d$Research)!=s,,drop=FALSE]
    f <- rma.uni(yi=dd$yi,sei=dd$sei,method="REML",test=test)
    data.frame(analysis=id,excluded_study=s,k=f$k,estimate=as.numeric(f$b),ci_lb=f$ci.lb,ci_ub=f$ci.ub,p=f$pval)
  }))
}
loo_rows <- list()
for(cmp in c("stn","base")) {
  d <- motor_effect(read_data("motor_overall.csv"),cmp)
  id <- paste0("motor_overall_",cmp)
  write.csv(d,file.path("outputs",paste0(id,"_effects.csv")),row.names=FALSE)
  save_fit(rma.uni(yi=d$yi,sei=d$sei,method="REML",test="knha"),id,"Mean difference")
  loo_rows[[id]] <- loo(d,id)
  ds <- motor_effect(read_data("motor_medication.csv"),cmp)
  ds <- ds[ds$Research!="Valldeoriola" & ds$Medication %in% c("ON","OFF"),]
  for(med in c("ON","OFF")) {
    dm <- ds[ds$Medication==med,]
    sid <- paste0("motor_",med,"_",cmp)
    save_fit(rma.uni(yi=dm$yi,sei=dm$sei,method="REML",test="knha"),sid,"Mean difference")
    loo_rows[[sid]] <- loo(dm,sid,"z")
  }
  mod <- rma.uni(yi=ds$yi,sei=ds$sei,mods=~factor(Medication),data=ds,method="REML",test="knha")
  writeLines(capture.output(print(mod)),file.path("outputs",paste0("motor_medication_moderator_",cmp,".txt")))
  dh <- motor_effect(read_data("motor_high_risk_excluded.csv"),cmp)
  save_fit(rma.uni(yi=dh$yi,sei=dh$sei,method="REML",test="z"),paste0("motor_high_risk_excluded_",cmp),"Mean difference")
}
write.csv(do.call(rbind,loo_rows),"outputs/motor_leave_one_study_out.csv",row.names=FALSE)
fog <- read_data("fog_extracted.csv")
replace <- is.na(fog$within_sd) & !is.na(fog$within_se)
fog$within_sd[replace] <- fog$within_se[replace]*sqrt(fog$n[replace])
fog$outcome_group <- ifelse(grepl("FOG_Q",fog$assessment,fixed=TRUE),"FOG_Q",
  ifelse(grepl("FOG_AC",fog$assessment,fixed=TRUE),"FOG_AC","axial"))
fog$effect_id <- seq_len(nrow(fog))
fog_effects <- list()
for(cmp in c("stn","base")) {
  d <- fog
  cm <- if(cmp=="stn") d$mean_ctrl else d$mean_base
  cs <- if(cmp=="stn") d$sd_ctrl else d$sd_base
  d$yi <- d$vi <- d$paired_r <- NA_real_
  d$r_source <- NA_character_
  for(i in seq_len(nrow(d))) {
    if(d$type[i]=="Cross") {
      r <- 0.5
      src <- "assumed r=0.5"
      if(cmp=="base" & !is.na(d$within_sd[i])) {
        candidate <- (d$sd_co[i]^2+cs[i]^2-d$within_sd[i]^2)/(2*d$sd_co[i]*cs[i])
        if(is.finite(candidate) & abs(candidate)<=1) {r <- candidate; src <- "derived from retained change SD"}
      }
      raw_d <- (d$mean_co[i]-cm[i])/sqrt((d$sd_co[i]^2+cs[i]^2)/2)
      J <- 1-3/(4*(d$n[i]-1)-1)
      d$yi[i] <- J*raw_d
      d$vi[i] <- J^2*(2*(1-r)/d$n[i]+raw_d^2/(2*d$n[i]))
      d$paired_r[i] <- r
      d$r_source[i] <- src
    } else {
      e <- escalc(measure="SMD",m1i=d$mean_co[i],sd1i=d$sd_co[i],n1i=d$n_co[i],
        m2i=cm[i],sd2i=cs[i],n2i=d$n_ctrl[i])
      d$yi[i] <- e$yi; d$vi[i] <- e$vi
      d$r_source[i] <- "independent groups"
    }
  }
  d$comparison <- cmp
  fog_effects[[cmp]] <- d
  for(gr in c("FOG_Q","FOG_AC")) {
    dg <- d[d$outcome_group==gr,]
    save_fit(rma.uni(yi=dg$yi,vi=dg$vi,method="REML",test="knha"),
      paste0("fog_",gr,"_",cmp),"Hedges' g")
  }
}
write.csv(do.call(rbind,fog_effects),"outputs/fog_effects_with_correlations.csv",row.names=FALSE)
f <- read_data("frequency_extracted.csv")
for(cmp in c("STN","base")) {
  yi_col <- paste0("MD_co_",cmp,"_final"); se_col <- paste0("SE_co_",cmp,"_final")
  d <- f[f$Assessment=="UPDRS_3" & f$Frequency>=60 & !is.na(f[[yi_col]]) & !is.na(f[[se_col]]),]
  d$yi <- d[[yi_col]]; d$sei <- d[[se_col]]
  d$frequency_group <- ifelse(d$Frequency>=119,"High","Low")
  for(gr in c("High","Low")) {
    dg <- d[d$frequency_group==gr,]
    if(nrow(dg)>=2) save_fit(rma.uni(yi=dg$yi,sei=dg$sei,method="REML",test="knha"),
      paste0("frequency_",gr,"_",cmp),"Mean difference")
  }
  mod <- rma.uni(yi=yi,sei=sei,mods=~factor(frequency_group),data=d,method="REML",test="knha")
  writeLines(capture.output(print(mod)),file.path("outputs",paste0("frequency_moderator_",cmp,".txt")))
}
summary_table <- do.call(rbind,summaries)
write.csv(summary_table,"outputs/model_summary.csv",row.names=FALSE)

