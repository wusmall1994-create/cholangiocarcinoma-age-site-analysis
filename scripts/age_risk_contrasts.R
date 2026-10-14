# Summarize existing replicates only; no new model fitting.
library(data.table)
root <- Sys.getenv('CCA_CLINICAL_OUT')
paths <- list.files(file.path(root,'auxiliary_bootstrap'),pattern='^rep_[0-9]+[.]rds$',full.names=TRUE)
stopifnot(length(paths)==1000)
b <- rbindlist(lapply(paths,function(p){z<-readRDS(p);copy(z$value$original)[,replicate:=z$replicate]}))
p <- fread(file.path(root,'original_binary_states.csv'))
contrast <- function(d,boot=FALSE){
 id <- c('horizon_months','age',if(boot)'replicate')
 z <- melt(d,id.vars=c(id,'exposure_level'),measure.vars=c('cancer','other','alive','all_death'),variable.name='state',value.name='risk')
 w <- dcast(z,as.formula(paste(paste(c(id,'state'),collapse='+'),'~exposure_level')),value.var='risk')
 w[,rd:=eCCA-iCCA]
 rbindlist(lapply(c(65,75),function(a){
  v<-dcast(w[age%in%c(50,a)],as.formula(paste(paste(c('horizon_months','state',if(boot)'replicate'),collapse='+'),'~age')),value.var='rd')
  v[,`:=`(older_age=a,estimate=get(as.character(a))-`50`)];setnames(v,c('50',as.character(a)),c('RD50','RD_older'));v
 }))
}
bb<-contrast(b,TRUE); pp<-contrast(p)
ci<-bb[,.(lower95=quantile(estimate,.025,type=7),upper95=quantile(estimate,.975,type=7),replicates=.N),by=.(horizon_months,state,older_age)]
fwrite(merge(pp,ci,by=c('horizon_months','state','older_age')),file.path(root,'original_age65_age75_contrasts_with_CI.csv'))
