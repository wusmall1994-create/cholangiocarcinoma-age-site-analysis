Sys.setenv(CCA_CLINICAL_MODE='audit')
source(file.path(Sys.getenv('CCA_CLINICAL_SCRIPTS_DIR'),'competing_outcomes.R'))
worker<-as.integer(Sys.getenv('CCA_WORKER','1'));workers<-as.integer(Sys.getenv('CCA_WORKERS','4'))
B<-as.integer(Sys.getenv('CCA_CLINICAL_B','1000'))
bd<-file.path(outnew,'auxiliary_bootstrap');dir.create(bd,showWarnings=FALSE)
original_ids<-split(seq_len(nrow(binary_data)),binary_data$patient_id)
harmonized_ids<-split(seq_len(nrow(base)),base$patient_id)
for(b in seq.int(worker,B,by=workers)){
 path<-file.path(bd,sprintf('rep_%04d.rds',b));if(file.exists(path))next
 warnings_log<-list();start<-proc.time()[3]
 val<-tryCatch(capture({
  set.seed(202610150L+b)
  z<-binary_data[unlist(original_ids[sample.int(length(original_ids),length(original_ids),replace=TRUE)],use.names=FALSE)]
  # Preserve levels; absent original categories have zero standardization weight.
  f<-fit_competing_models(z,'site_group',cov_binary)
  original<-predict_states(z,f,'site_group',cov_binary,ridge=FALSE)
  set.seed(202610140L+b)
  h<-base[unlist(harmonized_ids[sample.int(length(harmonized_ids),length(harmonized_ids),replace=TRUE)],use.names=FALSE)][!is.na(schema3)]
  stages<-rbindlist(lapply(c('Localized','Regional','Distant'),function(st){
   d<-h[as.character(stage4)==st]
   f<-fit_ridge(d,'schema3',.1)
   p<-predict_states(d,f,'schema3',cov_schema,age_values=c(50,75));p[,stage:=st];p
  }))
  list(original=original,stages=stages)
 },paste0('aux_',b)),error=function(e){warnings_log[[length(warnings_log)+1L]]<<-data.table(label='ERROR',message=conditionMessage(e));NULL})
 saveRDS(list(replicate=b,value=val,diagnostics=warnings_log,seconds=proc.time()[3]-start),path)
 cat('Replicate',b,'seconds',round(proc.time()[3]-start,1),'diagnostics',length(warnings_log),'\n');flush.console()
}
