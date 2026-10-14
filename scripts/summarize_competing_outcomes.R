options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table))
root<-Sys.getenv('CCA_CLINICAL_OUT')
read_boot<-function(folder){
 paths<-list.files(file.path(root,folder),pattern='^rep_[0-9]+[.]rds$',full.names=TRUE)
 x<-lapply(paths,readRDS)
 logs<-rbindlist(lapply(x,function(r){if(length(r$diagnostics))rbindlist(r$diagnostics)[,replicate:=r$replicate]}),fill=TRUE)
 if(nrow(logs))fwrite(logs,file.path(root,paste0(folder,'_diagnostics.csv')))
 status<-data.table(requested=1000,completed=length(x),successful=sum(vapply(x,function(r)!is.null(r$value),logical(1))),with_diagnostics=sum(vapply(x,function(r)length(r$diagnostics)>0,logical(1))))
 fwrite(status,file.path(root,paste0(folder,'_status.csv')));print(status)
 if(length(x)!=1000 || status$successful!=1000)stop('Incomplete or failed replicates; no final intervals released')
 stopifnot(identical(sort(vapply(x,function(z)as.integer(z$replicate),integer(1))),1:1000))
 x
}
summarize<-function(d,keys,value='estimate'){
 d[,.(estimate_mean=mean(get(value)),lower95=quantile(get(value),.025,type=7),upper95=quantile(get(value),.975,type=7),replicates=uniqueN(replicate)),by=keys]
}
longstates<-function(d,extra=character())melt(d,id.vars=c('horizon_months','age','exposure_level',extra),measure.vars=c('cancer','other','alive','all_death'),variable.name='state',value.name='estimate')
contrasts<-function(d,extra=character()){
 keys<-c('horizon_months','age','state',extra)
 w<-dcast(d,as.formula(paste(paste(keys,collapse='+'),'~ exposure_level')),value.var='estimate')
 rbindlist(lapply(setdiff(unique(d$exposure_level),'iCCA'),function(ex){q<-copy(w);q[,`:=`(comparison=paste0(ex,' minus iCCA'),estimate=get(ex)-iCCA)];q[,c(keys,'comparison','estimate'),with=FALSE]}))
}
core<-read_boot('bootstrap')
states<-rbindlist(lapply(core,function(r)copy(r$value$states)[,replicate:=r$replicate]))
dec<-rbindlist(lapply(core,function(r)copy(r$value$decomposition)[,replicate:=r$replicate]))
point<-fread(file.path(root,'harmonized_decomposition_states.csv'))
keys<-c('horizon_months','component','state','model')
fwrite(merge(point,summarize(dec,keys),by=keys),file.path(root,'harmonized_decomposition_with_CI.csv'))
paired<-dcast(dec,horizon_months+component+state+replicate~model,value.var='estimate');paired[,estimate:=interval-PH]
pp<-dcast(point,horizon_months+component+state~model,value.var='estimate');pp[,estimate:=interval-PH]
fwrite(merge(pp,summarize(paired,c('horizon_months','component','state')),by=c('horizon_months','component','state')),file.path(root,'interval_minus_PH_paired_CI.csv'))
state_long<-longstates(states,c('model','replicate'))
state_point<-longstates(fread(file.path(root,'harmonized_EOD_states.csv')),c('model'))
ks<-c('horizon_months','age','exposure_level','state','model')
fwrite(merge(state_point,summarize(state_long,ks),by=ks),file.path(root,'harmonized_subtype_states_with_CI.csv'))
cr<-contrasts(state_long,c('model','replicate'));cp<-contrasts(state_point,'model')
kr<-c('horizon_months','age','state','model','comparison')
fwrite(merge(cp,summarize(cr,kr),by=kr),file.path(root,'harmonized_subtype_RDs_with_CI.csv'))
aux<-read_boot('auxiliary_bootstrap')
for(kind in c('original','stages')){
 bs<-rbindlist(lapply(aux,function(r)copy(r$value[[kind]])[,replicate:=r$replicate]))
 file<-if(kind=='original')'original_binary_states.csv' else 'stage_conditional_states_POINT_ONLY.csv'
 extra<-if(kind=='original')character() else 'stage'
 bp<-longstates(bs,c(extra,'replicate'));pt<-longstates(fread(file.path(root,file)),extra)
 keys<-c('horizon_months','age','exposure_level','state',extra)
 fwrite(merge(pt,summarize(bp,keys),by=keys),file.path(root,paste0(kind,'_states_with_CI.csv')))
 bc<-contrasts(bp,c(extra,'replicate'));pc<-contrasts(pt,extra)
 keys<-c('horizon_months','age','state',extra,'comparison')
 fwrite(merge(pc,summarize(bc,keys),by=keys),file.path(root,paste0(kind,'_RDs_with_CI.csv')))
 # Pair age contrasts inside each bootstrap replicate, never subtract CI endpoints.
 ages<-dcast(bc[age%in%c(50,75)],as.formula(paste(paste(c('horizon_months','state',extra,'comparison','replicate'),collapse='+'),'~ age')),value.var='estimate');ages[,estimate:=`75`-`50`]
 ap<-dcast(pc[age%in%c(50,75)],as.formula(paste(paste(c('horizon_months','state',extra,'comparison'),collapse='+'),'~ age')),value.var='estimate');ap[,estimate:=`75`-`50`]
 keys<-c('horizon_months','state',extra,'comparison')
 fwrite(merge(ap,summarize(ages,keys),by=keys),file.path(root,paste0(kind,'_age_contrasts_with_CI.csv')))
}
cat('Core and auxiliary summaries completed; percentile type 7, 1000 paired replicates.\n')
