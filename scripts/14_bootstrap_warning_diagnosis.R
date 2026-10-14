options(stringsAsFactors=FALSE)
s<-readLines('11_horizon_anatomy_absolute_risks.R')
eval(parse(text=s[seq_len(grep('^binary_point_long <-',s)[1]-1L)]))
out<-file.path(results_dir,'composition_decomposition','warning_audit')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
target<-copy(schema_data)
logs<-list()
set.seed(20261011)
ids<-split(seq_len(nrow(target)),target$patient_id)
for(b in 1:200){
 z<-target[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]
 for(endpoint in c('cs_event','other_event')){
  rhs<-paste(c(paste0(ns_term,' * schema3'),cov_schema),collapse=' + ')
  fit<-withCallingHandlers(coxph(as.formula(paste0('Surv(survival_time_months,',endpoint,') ~ ',rhs)),data=z,ties='efron'),warning=function(w){
   logs[[length(logs)+1L]]<<-data.table(replicate=b,endpoint=endpoint,message=conditionMessage(w))
   invokeRestart('muffleWarning')
  })
  if(any(abs(coef(fit))>10,na.rm=TRUE)){
   fwrite(data.table(replicate=b,endpoint=endpoint,term=names(coef(fit)),coefficient=coef(fit)),file.path(out,paste0('coefficients_',b,'_',endpoint,'.csv')))
  }
 }
 if(length(logs)) fwrite(rbindlist(logs),file.path(out,'warnings.csv'))
 if(b%%20==0)cat('Diagnosed',b,'/200\n')
}
fwrite(target[,.(N=.N,cancer=sum(cs_event),other=sum(other_event)),by=race_factor],file.path(out,'race_event_counts.csv'))
cat('Warning count:',length(logs),'\n')
