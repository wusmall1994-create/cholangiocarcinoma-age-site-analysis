s<-readLines('11_horizon_anatomy_absolute_risks.R')
eval(parse(text=s[seq_len(grep('^binary_point_long <-',s)[1]-1L)]))
setDTthreads(1)
out<-file.path(results_dir,'composition_decomposition','warning_audit')
logs<-fread(file.path(out,'warnings.csv'))
ids<-split(seq_len(nrow(schema_data)),schema_data$patient_id)
set.seed(20261011);counts<-list();checks<-list()
for(b in 1:200){
 z<-schema_data[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]
 if(b%in%logs$replicate){
  counts[[length(counts)+1L]]<-data.table(replicate=b,race_unknown_N=sum(z$race_factor=='Non-Hispanic Unknown Race'),race_unknown_other=sum(z$other_event[z$race_factor=='Non-Hispanic Unknown Race']),rural_unknown_N=sum(z$rural3=='Unknown'),rural_unknown_other=sum(z$other_event[z$rural3=='Unknown']))
 }
 if(b==128){
  f<-as.formula(paste('Surv(survival_time_months,cs_event) ~',paste(c(paste0(ns_term,' * schema3'),cov_schema),collapse=' + ')))
  for(ep in c(1e-9,1e-11)){
   ww<-character();fit<-withCallingHandlers(coxph(f,data=z,ties='efron',control=coxph.control(eps=ep,iter.max=100)),warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')})
   checks[[length(checks)+1L]]<-data.table(eps=ep,term=names(coef(fit)),coefficient=coef(fit),SE=sqrt(diag(vcov(fit))),iterations=fit$iter,loglik=tail(fit$loglik,1),warnings=paste(ww,collapse=';'))
  }
 }
}
fwrite(rbindlist(counts),file.path(out,'warning_replicate_event_counts.csv'))
fwrite(rbindlist(checks),file.path(out,'replicate128_tolerance_check.csv'))
print(rbindlist(checks)[,.(eps,iterations,loglik,warnings)][!duplicated(eps)])
