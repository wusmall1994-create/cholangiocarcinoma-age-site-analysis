library(data.table)
p <- Sys.getenv('CCA_RESULTS_DIR'); stopifnot(nzchar(p))
files<-list.files(file.path(p,'bootstrap'),pattern='^rep_[0-9]+[.]rds$',full.names=TRUE)
stopifnot(length(files)==2000L)
reps<-lapply(files,readRDS)
stopifnot(identical(sort(vapply(reps,function(x)x$replicate,1L)),1:2000))
diag<-rbindlist(lapply(reps,function(x)data.table(replicate=x$replicate,seed=x$seed,seconds=x$seconds,n_estimates=if(is.null(x$value))0L else nrow(x$value),warnings=length(x$warnings))))
fwrite(diag,file.path(p,'bootstrap_diagnostics.csv'))
stopifnot(all(diag$n_estimates==diag$n_estimates[1]),all(diag$n_estimates>0))
w<-rbindlist(lapply(reps,function(x)if(length(x$warnings))rbindlist(x$warnings)[,replicate:=x$replicate]else NULL))
if(nrow(w))fwrite(w,file.path(p,'bootstrap_warning_log.csv'))
b<-rbindlist(lapply(reps,function(x)copy(x$value)[,replicate:=x$replicate]))
stopifnot(all(is.finite(b$estimate)))
key<-c('sequence_group','definition','state','contrast')
ci<-b[,.(lower95=quantile(estimate,.025,type=7),upper95=quantile(estimate,.975,type=7),B=.N),by=key]
ans<-merge(fread(file.path(p,'ridge_point_contrasts.csv')),ci,by=key)
stopifnot(all(ans$B==2000),nrow(ans)==nrow(reps[[1]]$value))
fwrite(ans,file.path(p,'ridge_contrasts_with_CI.csv'))
fwrite(b,file.path(p,'aggregate_bootstrap_contrasts.csv'))
saveRDS(b,file.path(p,'aggregate_bootstrap_contrasts.rds'))
main<-ans[state=='cancer' & contrast%chin%c('RD_age50','RD_age75','RD_change_75_minus_50') | definition=='Direct_allcause_including_unknown']
main[,c('estimate','lower95','upper95'):=lapply(.SD,function(x)100*x),.SDcols=c('estimate','lower95','upper95')]
fwrite(main,file.path(p,'key_results_percentage_points.csv'))
# Replicate paired identities, distinct from interval endpoint subtraction.
checks<-list()
for(g in unique(b$sequence_group)){
 for(de in c('Common_biliary','Common_hepatobiliary','Scenario621')){
  x<-b[sequence_group==g & definition%chin%c('SEER_original',de,paste0(de,'_minus_SEER'))]
  if(uniqueN(x$definition)!=3L)next
  x<-dcast(x,replicate+state+contrast~definition,value.var='estimate')
  err<-max(abs(x[[paste0(de,'_minus_SEER')]]-(x[[de]]-x$SEER_original)),na.rm=TRUE)
  checks[[length(checks)+1L]]<-data.table(check=paste(g,de,'paired identity'),max_error=err,pass=err<1e-12)
 }
}
fwrite(rbindlist(checks),file.path(p,'paired_identity_QA.csv'))
# Independent random-half CI diagnostic (not another inferential CI).
half<-b[state=='cancer' & sequence_group=='All' & contrast%chin%c('RD_age75','RD_change_75_minus_50')]
half[,half:=ifelse(replicate<=1000,'first1000','second1000')]
hc<-half[,.(lower95=100*quantile(estimate,.025,type=7),upper95=100*quantile(estimate,.975,type=7)),by=.(definition,state,contrast,half)]
fwrite(hc,file.path(p,'half_sample_interval_stability.csv'))
# Quantile Monte Carlo SE via resampling the already-estimated bootstrap distribution.
set.seed(202610199)
mc<-rbindlist(lapply(split(half,by=c('definition','state','contrast')),function(x){
 q<-replicate(400,quantile(sample(x$estimate,2000,replace=TRUE),c(.025,.975),type=7))
 data.table(definition=x$definition[1],state=x$state[1],contrast=x$contrast[1],lower_MCSE_pp=100*sd(q[1,]),upper_MCSE_pp=100*sd(q[2,]))
}))
fwrite(mc,file.path(p,'quantile_MCSE.csv'))
cat('B=',nrow(diag),'warnings=',sum(diag$warnings),'contrasts=',nrow(ans),'\n')
print(main[sequence_group=='All' & contrast=='RD_change_75_minus_50'])
