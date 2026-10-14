Sys.setenv(CCA_CLINICAL_MODE='audit')
source(file.path(Sys.getenv('CCA_CLINICAL_SCRIPTS_DIR'),'competing_outcomes.R'))
original<-fread(file.path(outnew,'original_binary_states.csv'))
old<-fread(file.path(Sys.getenv('CCA_AGGREGATE_RESULTS'),'reporting_extensions/horizon_standardized_cif.csv'))
print(names(old))
compare<-merge(original,old,by=c('horizon_months','age','exposure_level'))
stopifnot(nrow(compare)==nrow(original),max(abs(compare$cancer-compare$risk))<1e-10)
cat('Original cancer-risk replication PASS\n')
stopifnot(all(binary_data$cs_event+binary_data$other_event==binary_data$os_event),all(target$cs_event+target$other_event==target$os_event))
minus<-fread(file.path(outnew,'original_single_row_exclusion_sensitivity.csv'))
z<-merge(original,minus,by=c('age','horizon_months','exposure_level'),suffixes=c('_full','_excluded'))
for(v in c('cancer','other','alive','all_death'))cat('single-row maximum absolute change',v,max(abs(z[[paste0(v,'_full')]]-z[[paste0(v,'_excluded')]])),'\n')
rd<-function(d,by){
 z<-melt(d,id.vars=by,measure.vars=c('cancer','other','alive','all_death'),variable.name='state',value.name='risk')
 w<-dcast(z,as.formula(paste(paste(setdiff(c(by,'state'),'exposure_level'),collapse='+'),'~ exposure_level')),value.var='risk')
 for(v in setdiff(unique(d$exposure_level),'iCCA'))w[,(paste0('RD_',v)):=get(v)-iCCA]
 w
}
fwrite(rd(original,c('horizon_months','age','exposure_level')),file.path(outnew,'original_binary_state_RDs.csv'))
stage<-fread(file.path(outnew,'stage_conditional_states_POINT_ONLY.csv'))
fwrite(rd(stage,c('horizon_months','age','exposure_level','stage')),file.path(outnew,'stage_conditional_RDs_POINT_ONLY.csv'))
cat('Stage-specific age 50/75 cancer contrasts\n')
print(rd(stage,c('horizon_months','age','exposure_level','stage'))[horizon_months==60 & age%in%c(50,75) & state=='cancer'])
# Individual-level full-target integration independently checks exact aggregation.
fit<-readRDS(file.path(outnew,'local_models_DO_NOT_UPLOAD.rds'))$interval
nd<-copy(target);nd[,`:=`(age_years=75,schema3=factor('Distal eCCA',levels=levels(target$schema3)))]
ss<-rep(1,nrow(nd));cc<-oo<-rep(0,nrow(nd))
for(k in 1:3){
 nd[,period:=k];X<-interval_design(nd,'schema3')
 rc<-exp(drop(X%*%coef(fit$cancer))-sum(fit$cancer$means*coef(fit$cancer)))
 ro<-exp(drop(X%*%coef(fit$other))-sum(fit$other$means*coef(fit$other)))
 hc<-subset(fit$bh_cancer,as.integer(factor(strata))==k);ho<-subset(fit$bh_other,as.integer(factor(strata))==k)
 tt<-sort(unique(c(hc$time,ho$time)))
 dc<-diff(c(0,cumhaz_at(hc,tt)));do<-diff(c(0,cumhaz_at(ho,tt)))
 for(j in seq_along(tt)){
  ac<-dc[j]*rc;ao<-do[j]*ro;a<-ac+ao
  cc<-cc+ss*(-expm1(-a))*ifelse(a>0,ac/a,0)
  oo<-oo+ss*(-expm1(-a))*ifelse(a>0,ao/a,0)
  ss<-ss*exp(-a)
 }
}
reference<-fread(file.path(outnew,'harmonized_EOD_states.csv'))[model=='interval' & age==75 & horizon_months==60 & exposure_level=='Distal eCCA']
err<-max(abs(c(mean(cc),mean(oo),mean(ss))-as.numeric(reference[,.(cancer,other,alive)])))
stopifnot(err<1e-10)
cat('Full-row interval prediction equivalence PASS; maximum error',err,'\n')
