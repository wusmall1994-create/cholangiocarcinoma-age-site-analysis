Sys.setenv(CCA_CLINICAL_MODE='audit')
source(file.path(Sys.getenv('CCA_CLINICAL_SCRIPTS_DIR'),'competing_outcomes.R'))
x<-lapply(list.files(file.path(outnew,'auxiliary_bootstrap'),pattern='[.]rds$',full.names=TRUE),readRDS)
logs<-rbindlist(lapply(x,function(z)if(length(z$diagnostics))rbindlist(z$diagnostics)[,replicate:=z$replicate]))
pick<-logs[,.(replicate=replicate[1]),by=message]
ids<-split(seq_len(nrow(binary_data)),binary_data$patient_id)
report<-list()
for(b in pick$replicate){
 set.seed(202610150L+b);z<-binary_data[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]
 for(ev in c('cs_event','other_event')){
  ww<-character()
  f<-withCallingHandlers(coxph(as.formula(paste('Surv(survival_time_months,',ev,') ~',ns_term,'* site_group +',paste(cov_binary,collapse='+'))),data=z,ties='efron'),warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')})
  if(length(ww))for(message in ww){
   idx<-as.integer(unlist(regmatches(message,gregexpr('[0-9]+',message))))
   for(j in idx){
    term<-names(coef(f))[j]
    covar<-cov_binary[vapply(cov_binary,function(v)startsWith(term,v),logical(1))]
    if(length(covar)){lev<-substring(term,nchar(covar)+1);group<-z[as.character(get(covar))==lev]}else{group<-z[0]}
    report[[length(report)+1L]]<-data.table(replicate=b,event=ev,term=term,coefficient=coef(f)[j],category_records=nrow(group),category_events=sum(group[[ev]]),message=message)
   }
  }
 }
}
ans<-rbindlist(report);fwrite(ans,file.path(outnew,'original_bootstrap_warning_mechanism.csv'));print(ans)
numerical<-unique(ans[!grepl('Unknown',term),replicate])
checks<-list()
for(b in numerical){
 set.seed(202610150L+b);z<-binary_data[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]
 rhs<-paste(ns_term,'* site_group +',paste(cov_binary,collapse='+'))
 fitone<-function(ev)suppressWarnings(coxph(as.formula(paste('Surv(survival_time_months,',ev,') ~',rhs)),data=z,ties='efron',control=coxph.control(iter.max=100,eps=1e-11)))
 a<-fitone('cs_event');o<-fitone('other_event')
 f<-list(cancer=a,other=o,bh_cancer=basehaz(a,centered=TRUE),bh_other=basehaz(o,centered=TRUE))
 p<-predict_states(z,f,'site_group',cov_binary,ridge=FALSE)
 before<-readRDS(file.path(outnew,'auxiliary_bootstrap',sprintf('rep_%04d.rds',b)))$value$original
 q<-merge(p,before,by=c('horizon_months','age','exposure_level'),suffixes=c('_tight','_original'))
 checks[[length(checks)+1L]]<-data.table(replicate=b,max_probability_difference=max(abs(as.matrix(q[,paste0(c('cancer','other','alive','all_death'),'_tight'),with=FALSE])-as.matrix(q[,paste0(c('cancer','other','alive','all_death'),'_original'),with=FALSE]))))
}
if(length(checks)){cc<-rbindlist(checks);fwrite(cc,file.path(outnew,'numerical_warning_tolerance_check.csv'));print(cc)}
