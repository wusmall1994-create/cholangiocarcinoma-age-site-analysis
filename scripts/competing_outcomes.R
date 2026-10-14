# Joint competing-outcome, interval-specific, and stage-conditional models.
options(stringsAsFactors=FALSE)
stopifnot(nzchar(Sys.getenv('CCA_DERIVED_DIR')), nzchar(Sys.getenv('CCA_CLINICAL_OUT')), nzchar(Sys.getenv('CCA_CLINICAL_SCRIPTS_DIR')))
Sys.setenv(CCA_RESULTS_DIR=file.path(Sys.getenv('CCA_CLINICAL_OUT'),'source_runtime'))
setwd(Sys.getenv('CCA_CLINICAL_SCRIPTS_DIR'))
s<-readLines('16_paired_bootstrap.R')
eval(parse(text=s[seq_len(grep('^mode<-',s)[1]-1L)]))
outnew<-Sys.getenv('CCA_CLINICAL_OUT')
dir.create(outnew,recursive=TRUE,showWarnings=FALSE)
setDTthreads(1)
warnings_log<-list()
capture<-function(expr,label){withCallingHandlers(expr,warning=function(w){warnings_log[[length(warnings_log)+1L]]<<-data.table(label=label,message=conditionMessage(w));invokeRestart('muffleWarning')})}

# Joint interval-specific model: age-by-anatomy block varies by period;
# demographic/context coefficients shared; separate baseline hazards.
fit_interval<-function(d,exposure){
 z<-copy(d);z[,time60:=pmin(survival_time_months,60)]
 z[,evt:=as.integer(cs_event==1 & survival_time_months<=60)]
 z<-as.data.table(survSplit(Surv(time60,evt)~.,as.data.frame(z),cut=c(12,36),start='start',end='stop',event='evt',episode='period'))
 z[,other_evt:=as.integer(other_event==1 & stop==survival_time_months & survival_time_months<=60)]
 X<-interval_design(z,exposure)
 one<-function(ev){
  dat<-list(y=Surv(z$start,z$stop,z[[ev]]),X=X,period=z$period)
  f<-coxph(y~ridge(X,theta=.1,scale=FALSE)+strata(period),data=dat,ties='efron',x=TRUE,model=TRUE,control=coxph.control(iter.max=100,eps=1e-11))
  stopifnot(all(is.finite(coef(f))))
  f
 }
 a<-one('evt');b<-one('other_evt')
 list(cancer=a,other=b,bh_cancer=basehaz(a,centered=TRUE),bh_other=basehaz(b,centered=TRUE))
}
interval_design<-function(d,exposure){
 X<-design(d,exposure)
 varying<-grepl('splines::|^schema3|^site_group',colnames(X))
 ans<-cbind(X[,!varying,drop=FALSE],do.call(cbind,lapply(1:3,function(k){q<-X[,varying,drop=FALSE]*(d$period==k);colnames(q)<-paste0(colnames(q),'_period',k);q})))
 ans
}

predict_states<-function(targ,fit,exposure,covars,age_values=c(40,50,65,75),interval=FALSE,ridge=TRUE){
 u<-targ[,.(weight=.N),by=covars]; w<-u$weight/sum(u$weight)
 lev<-levels(targ[[exposure]]);res<-list()
 for(age in age_values)for(level in lev){
  nd<-copy(u);nd[,age_years:=age];set(nd,j=exposure,value=factor(level,levels=lev))
  ss<-rep(1,nrow(nd));cc<-oo<-rep(0,nrow(nd))
  for(period_value in if(interval)1:3 else 1){
   nd[,period:=period_value]
   bhc<-fit$bh_cancer;bho<-fit$bh_other
   if(interval){bhc<-bhc[as.integer(factor(bhc$strata))==period_value,];bho<-bho[as.integer(factor(bho$strata))==period_value,]}
   limit<-if(interval)c(12,36,60)[period_value] else 60
   times<-sort(unique(c(bhc$time[bhc$time<=limit],bho$time[bho$time<=limit],if(interval)limit else c(12,36,60))))
   dc<-diff(c(0,cumhaz_at(bhc,times)));do<-diff(c(0,cumhaz_at(bho,times)))
   if(ridge){
    X<-if(interval)interval_design(nd,exposure) else design(nd,exposure)
    rc<-exp(drop(X%*%coef(fit$cancer))-sum(fit$cancer$means*coef(fit$cancer)))
    ro<-exp(drop(X%*%coef(fit$other))-sum(fit$other$means*coef(fit$other)))
   }else{
    rc<-exp(predict(fit$cancer,newdata=nd,type='lp',reference='sample'))
    ro<-exp(predict(fit$other,newdata=nd,type='lp',reference='sample'))
   }
   for(k in seq_along(times)){
    hc<-dc[k]*rc;ho<-do[k]*ro;hh<-hc+ho;event<-ss*(-expm1(-hh));fraction<-ifelse(hh>0,hc/hh,0)
    cc<-cc+event*fraction;oo<-oo+event*(1-fraction);ss<-ss*exp(-hh)
    if(times[k]%in%c(12,36,60))res[[length(res)+1L]]<-data.table(horizon_months=times[k],age=age,exposure_level=level,cancer=sum(cc*w),other=sum(oo*w),alive=sum(ss*w),all_death=sum((cc+oo)*w))
   }
  }
 }
 ans<-rbindlist(res)
 stopifnot(max(abs(ans$cancer+ans$other+ans$alive-1))<1e-10,all(is.finite(as.matrix(ans[,.(cancer,other,alive)]))))
 ans
}
decompose_states<-function(p,z,label){
 e<-copy(z[schema3!='iCCA']);e[,distal:=as.integer(schema3=='Distal eCCA')]
 cm<-glm(as.formula(paste('distal ~',ns_term)),data=e,family=binomial())
 weights<-predict(cm,newdata=data.frame(age_years=c(50,75)),type='response')
 rbindlist(lapply(c('cancer','other','all_death','alive'),function(state){
  r<-dcast(p[age%in%c(50,75)],horizon_months+age~exposure_level,value.var=state)
  r[,weight:=weights[match(age,c(50,75))]]
  r[,RD:=(1-weight)*`Perihilar eCCA`+weight*`Distal eCCA`-iCCA]
  a<-r[,{
   x<-.SD[age==50];y<-.SD[age==75];wm<-(x$weight+y$weight)/2
   v<-c((y$weight-x$weight)*(x$`Distal eCCA`+y$`Distal eCCA`-x$`Perihilar eCCA`-y$`Perihilar eCCA`)/2,(1-wm)*(y$`Perihilar eCCA`-x$`Perihilar eCCA`),wm*(y$`Distal eCCA`-x$`Distal eCCA`),-(y$iCCA-x$iCCA))
   stopifnot(abs(sum(v)-(y$RD-x$RD))<1e-10)
   list(component=c('composition','perihilar_risk','distal_risk','iCCA_reference','total','RD50','RD75'),estimate=c(v,y$RD-x$RD,x$RD,y$RD))
  },by=horizon_months];a[,`:=`(state=state,model=label)];a
 }))
}

cat('Cohorts',nrow(binary_data),nrow(base),nrow(target),'\n')
fwrite(data.table(analysis=c('original_binary','harmonized_D1','harmonized_EOD'),N=c(nrow(binary_data),nrow(base),nrow(target)),patients=c(uniqueN(binary_data$patient_id),uniqueN(base$patient_id),uniqueN(target$patient_id))),file.path(outnew,'cohorts.csv'))
mode_new<-Sys.getenv('CCA_CLINICAL_MODE','point')
if(mode_new=='point'){
 originalfit<-capture(fit_competing_models(binary_data,'site_group',cov_binary),'original_PH')
 original<-predict_states(binary_data,originalfit,'site_group',cov_binary,ridge=FALSE)
 fwrite(original,file.path(outnew,'original_binary_states.csv'))
 ph<-capture(fit_ridge(target,'schema3',.1),'EOD_PH')
 ip<-capture(fit_interval(target,'schema3'),'EOD_interval')
 states<-rbindlist(lapply(c('PH','interval'),function(m){p<-predict_states(target,if(m=='PH')ph else ip,'schema3',cov_schema,interval=m=='interval');p[,model:=m];p}))
 fwrite(states,file.path(outnew,'harmonized_EOD_states.csv'))
 dd<-rbindlist(lapply(c('PH','interval'),function(m)decompose_states(states[model==m],target,m)))
 fwrite(dd,file.path(outnew,'harmonized_decomposition_states.csv'))
 # Exact replication of previously validated PH cancer estimates.
 old<-decomp(target,target,.1,'EOD')
 check<-merge(dd[model=='PH' & state=='cancer'],old,by=c('horizon_months','component'),suffixes=c('_new','_old'))
 stopifnot(max(abs(check$estimate_new-check$estimate_old))<1e-10)
 fwrite(check,file.path(outnew,'PH_replication_audit.csv'))
 saveRDS(list(original=originalfit,ph=ph,interval=ip),file.path(outnew,'local_models_DO_NOT_UPLOAD.rds'))
 print(dd[horizon_months==60 & component%in%c('composition','total','RD50','RD75')])
 cat('Columns for descriptive extensions\n');print(names(target))
 if(length(warnings_log))fwrite(rbindlist(warnings_log),file.path(outnew,'warnings.csv'))
 cat('POINT COMPLETED; warnings',length(warnings_log),'\n')
}
if(mode_new=='bootstrap'){
 B<-as.integer(Sys.getenv('CCA_CLINICAL_B','1000'))
 worker<-as.integer(Sys.getenv('CCA_WORKER','1'));workers<-as.integer(Sys.getenv('CCA_WORKERS','4'))
 bd<-file.path(outnew,'bootstrap');dir.create(bd,showWarnings=FALSE)
 ids<-split(seq_len(nrow(base)),base$patient_id)
 for(b in seq.int(worker,B,by=workers)){
  path<-file.path(bd,sprintf('rep_%04d.rds',b));if(file.exists(path))next
  set.seed(202610140L+b)
  z<-base[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]
  targ<-z[!is.na(schema3)];warnings_log<-list();started<-proc.time()[3]
  value<-tryCatch(capture({
   pp<-predict_states(targ,fit_ridge(targ,'schema3',.1),'schema3',cov_schema,age_values=c(50,75))
   pi<-predict_states(targ,fit_interval(targ,'schema3'),'schema3',cov_schema,age_values=c(50,75),interval=TRUE)
   ds<-rbindlist(list(decompose_states(pp,targ,'PH'),decompose_states(pi,targ,'interval')))
   pp[,model:='PH'];pi[,model:='interval']
   list(states=rbindlist(list(pp,pi)),decomposition=ds)
  },paste0('replicate_',b)),error=function(e){warnings_log[[length(warnings_log)+1L]]<<-data.table(label='ERROR',message=conditionMessage(e));NULL})
  saveRDS(list(replicate=b,value=value,diagnostics=warnings_log,seconds=proc.time()[3]-started),path)
  cat('Replicate',b,'/',B,'seconds',round(proc.time()[3]-started,1),'diagnostics',length(warnings_log),'\n');flush.console()
 }
}
if(mode_new=='descriptive'){
 for(label in c('original_binary','harmonized_EOD')){
  z<-copy(if(label=='original_binary')binary_data else target);ex<-if(label=='original_binary')'site_group' else 'schema3'
  z[,age_band:=cut(age_years,c(15,40,50,65,75,85),right=FALSE)]
  for(v in c('stage4','surgery_primary','Chemotherapy recode (yes, no/unk)')){
   z[,category:=as.character(get(v))]
   tab<-z[,.(N=.N),by=c('age_band',ex,'category')];tab[,denominator:=sum(N),by=c('age_band',ex)];tab[,proportion:=N/denominator]
   fwrite(tab,file.path(outnew,paste0(label,'_',if(v=='stage4')'stage' else if(v=='surgery_primary')'surgery' else 'recorded_chemotherapy','_distribution.csv')))
  }
 }
 # Stage-restricted fits and stage-specific empirical targets; not mediation.
 stage_results<-list();stage_counts<-list()
 for(st in c('Localized','Regional','Distant')){
  z<-copy(target[as.character(stage4)==st]);cat('Stage',st,'N',nrow(z),'\n');stopifnot(nrow(z)>100)
  stage_counts[[st]]<-z[,.(N=.N,cancer=sum(cs_event),other=sum(other_event)),by=schema3][,stage:=st]
  f<-capture(fit_ridge(z,'schema3',.1),paste0('stage_',st))
  p<-predict_states(z,f,'schema3',cov_schema);p[,stage:=st];stage_results[[st]]<-p
 }
 fwrite(rbindlist(stage_results),file.path(outnew,'stage_conditional_states_POINT_ONLY.csv'))
 fwrite(rbindlist(stage_counts),file.path(outnew,'stage_cohort_counts.csv'))
 if(length(warnings_log))fwrite(rbindlist(warnings_log),file.path(outnew,'stage_warnings.csv'))
 # Audit original model warning and compare dropping the single unknown-income row.
 m<-readRDS(file.path(outnew,'local_models_DO_NOT_UPLOAD.rds'))$original
 cat('Original other model coefficient 17:',names(coef(m$other))[17],'\n')
 print(binary_data[,.(N=.N,other=sum(other_event)),by=income3])
 sub<-droplevels(binary_data[income3!='Unknown'])
 f<-capture(fit_competing_models(sub,'site_group',cov_binary),'original_without_unknown_income')
 p<-predict_states(sub,f,'site_group',cov_binary,ridge=FALSE)
 fwrite(p,file.path(outnew,'original_single_row_exclusion_sensitivity.csv'))
}
