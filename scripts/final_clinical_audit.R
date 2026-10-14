suppressPackageStartupMessages(library(data.table))
p<-Sys.getenv('CCA_CLINICAL_OUT')
checks<-list()
for(dir in c('bootstrap','auxiliary_bootstrap','original_ridge_bootstrap')){
 files<-list.files(file.path(p,dir),pattern='^rep_[0-9]+[.]rds$',full.names=TRUE)
 stopifnot(length(files)==1000)
 maxsum<-0;maxadd<-0;nwarn<-0
 for(f in files){
  z<-readRDS(f);stopifnot(!is.null(z$value));nwarn<-nwarn+as.integer(length(z$diagnostics)>0)
  tables<-if(dir=='bootstrap')list(z$value$states) else if(dir=='auxiliary_bootstrap')list(z$value$original,z$value$stages) else list(z$value)
  for(t in tables){
   m<-as.matrix(t[,.(cancer,other,alive,all_death)])
   stopifnot(all(is.finite(m)),all(m>=0 & m<=1))
   maxsum<-max(maxsum,max(abs(t$cancer+t$other+t$alive-1)),max(abs(t$all_death-t$cancer-t$other)))
  }
  if(dir=='bootstrap'){
   d<-dcast(z$value$decomposition,horizon_months+state+model~component,value.var='estimate')
   maxadd<-max(maxadd,max(abs(d$total-d$composition-d$perihilar_risk-d$distal_risk-d$iCCA_reference)),max(abs(d$total-(d$RD75-d$RD50))))
  }
 }
 stopifnot(maxsum<1e-10,maxadd<1e-10)
 checks[[dir]]<-data.table(analysis=dir,replicates=length(files),replicates_with_diagnostics=nwarn,max_probability_identity_error=maxsum,max_decomposition_identity_error=maxadd,status='PASS')
}
for(f in list.files(p,pattern='with_CI[.]csv$|paired.*CI[.]csv$',full.names=TRUE)){
 d<-fread(f);stopifnot(all(is.finite(d$estimate)),all(is.finite(d$lower95)),all(is.finite(d$upper95)),all(d$lower95<=d$upper95),all(d$replicates==1000))
}
result<-rbindlist(checks);fwrite(result,file.path(p,'FINAL_NUMERICAL_AUDIT.csv'));print(result)
capture.output(sessionInfo(),file=file.path(p,'R_SESSION_INFO.txt'))
