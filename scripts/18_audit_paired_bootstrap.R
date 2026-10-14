options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table))
root<-file.path(Sys.getenv('CCA_RESULTS_DIR'),'composition_decomposition','paired_1000')
files<-list.files(root,pattern='^rep_[0-9]+\\.rds$',full.names=TRUE)
objects<-lapply(files,function(p)tryCatch(readRDS(p),error=function(e)NULL))
objects<-Filter(Negate(is.null),objects)
bad<-Filter(function(x)is.null(x$value)||length(x$diagnostics)>0,objects)
cat('Saved:',length(objects),'Warned/failed:',length(bad),'\n')
if(length(bad))print(rbindlist(lapply(bad,function(x)data.table(replicate=x$replicate,message=paste(x$diagnostics,collapse=' | ')))))
if(Sys.getenv('CCA_AUDIT_MODE','status')=='final'){
 stopifnot(length(objects)==1000L,length(bad)==0L,identical(sort(vapply(objects,`[[`,integer(1),'replicate')),1:1000))
 bb<-rbindlist(lapply(objects,function(x){q<-copy(x$value);q[,replicate:=x$replicate];q}))
 stopifnot(nrow(bb)==78000L,all(is.finite(bb$estimate)),!anyDuplicated(bb[,.(replicate,analysis,horizon_months,component)]))
 sums<-dcast(bb[analysis%in%c('EOD','EOD_8160')],replicate+analysis+horizon_months~component,value.var='estimate')
 error<-max(abs(sums$total-sums$composition-sums$perihilar_risk-sums$distal_risk-sums$iCCA_reference))
 stopifnot(error<1e-10,max(abs(sums$total-(sums$RD75-sums$RD50)))<1e-10)
 point<-fread(file.path(root,'penalty_sensitivity_points.csv'))[theta==.1]
 summarise<-function(b,p){
  ci<-b[,.(lower95=quantile(estimate,.025,type=7),upper95=quantile(estimate,.975,type=7),successful=uniqueN(replicate)),by=.(analysis,horizon_months,component)]
  merge(p,ci,by=c('analysis','horizon_months','component'))
 }
 paired<-function(d,hasrep=FALSE){
  keys<-c(if(hasrep)'replicate','horizon_months','component')
  ans<-list()
  for(a in c('EOD_8160','D2_exclude_discordant','D3_8160_only','D4_reallocate_8162')){
   ref<-if(a=='EOD_8160')'EOD' else 'D1_topography'
   m<-merge(d[analysis==a,c(keys,'estimate'),with=FALSE],d[analysis==ref,c(keys,'estimate'),with=FALSE],by=keys)
   m[,`:=`(analysis=paste0(a,'_minus_',ref),estimate=estimate.x-estimate.y)]
   ans[[a]]<-m[,c(keys,'analysis','estimate'),with=FALSE]
  }
  rbindlist(ans)
 }
 ci<-summarise(bb,point)
 pp<-paired(point);pb<-paired(bb,TRUE);pc<-summarise(pb,pp)
 fwrite(bb,file.path(root,'bootstrap_estimates.csv'))
 fwrite(ci,file.path(root,'estimates_with_ci.csv'))
 fwrite(pc,file.path(root,'paired_differences_with_ci.csv'))
 fwrite(pb,file.path(root,'paired_bootstrap_differences.csv'))
 # Combined non-composition component: calculate within replicate, not by summing CI limits.
 nc<-sums[,.(replicate,analysis,horizon_months,estimate=total-composition)]
 nci<-nc[,.(lower95=quantile(estimate,.025),upper95=quantile(estimate,.975)),by=.(analysis,horizon_months)]
 fwrite(nci,file.path(root,'noncomposition_ci.csv'))
 fwrite(data.table(status='PASS',replicates=1000,rows=nrow(bb),warnings=0,errors=0,max_additivity_error=error),file.path(root,'final_audit.csv'))
 print(ci[horizon_months==60,.(analysis,component,estimate=100*estimate,lower95=100*lower95,upper95=100*upper95)])
 print(pc[horizon_months==60,.(analysis,component,estimate=100*estimate,lower95=100*lower95,upper95=100*upper95)])
 capture.output(sessionInfo(),file=file.path(root,'audit_session_info.txt'))
}
