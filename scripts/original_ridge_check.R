Sys.setenv(CCA_CLINICAL_MODE='audit')
source(file.path(Sys.getenv('CCA_CLINICAL_SCRIPTS_DIR'),'competing_outcomes.R'))
cov_schema<-cov_binary
worker<-as.integer(Sys.getenv('CCA_WORKER','1'));workers<-as.integer(Sys.getenv('CCA_WORKERS','4'))
bd<-file.path(outnew,'original_ridge_bootstrap');dir.create(bd,showWarnings=FALSE)
if(worker==1){p<-predict_states(binary_data,fit_ridge(binary_data,'site_group',.1),'site_group',cov_binary);fwrite(p,file.path(outnew,'original_ridge_states_POINT.csv'))}
ids<-split(seq_len(nrow(binary_data)),binary_data$patient_id)
for(b in seq.int(worker,1000,by=workers)){
 path<-file.path(bd,sprintf('rep_%04d.rds',b));if(file.exists(path))next
 set.seed(202610150L+b);z<-binary_data[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]
 warnings_log<-list()
 p<-tryCatch(capture(predict_states(z,fit_ridge(z,'site_group',.1),'site_group',cov_binary),'ridge_binary'),error=function(e){warnings_log[[length(warnings_log)+1L]]<<-data.table(label='ERROR',message=conditionMessage(e));NULL})
 saveRDS(list(replicate=b,value=p,diagnostics=warnings_log),path)
 if(b%%20%in%1:4)cat('Ridge check',b,'\n')
}
