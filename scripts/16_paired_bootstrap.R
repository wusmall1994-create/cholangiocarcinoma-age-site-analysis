# Full paired patient bootstrap. No patient-level files are exported.
options(stringsAsFactors=FALSE)
s<-readLines('11_horizon_anatomy_absolute_risks.R')
eval(parse(text=s[seq_len(grep('^binary_point_long <-',s)[1]-1L)]))
setDTthreads(1L)
out<-file.path(results_dir,'composition_decomposition','paired_1000')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
base<-copy(binary_data[diagnosis_year>=2010L]);base[,hist:=as.integer(`Histologic Type ICD-O-3`)]
target<-copy(schema_data)
defs<-c('D1_topography','D2_exclude_discordant','D3_8160_only','D4_reallocate_8162')
make_definition<-function(x,def){
 z<-copy(x)
 if(def==defs[2])z<-z[!(site_group=='iCCA' & hist==8162L)]
 if(def==defs[3])z<-z[hist==8160L]
 if(def==defs[4])z[hist==8162L,site_group:=factor('eCCA',levels=c('iCCA','eCCA'))]
 z
}
design<-function(d,exposure){
 model.matrix(as.formula(paste('~',paste(c(paste0(ns_term,' * ',exposure),cov_schema),collapse=' + '))),data=d)[,-1,drop=FALSE]
}
fit_ridge<-function(d,exposure,theta){
 X<-design(d,exposure)
 fit_one<-function(event){
  dat<-list(y=Surv(d$survival_time_months,d[[event]]),X=X)
  f<-coxph(y~ridge(X,theta=theta,scale=FALSE),data=dat,ties='efron',model=TRUE,x=TRUE,control=coxph.control(iter.max=60,eps=1e-11))
  stopifnot(all(is.finite(coef(f))),all(is.finite(f$loglik)))
  f
 }
 a<-fit_one('cs_event');b<-fit_one('other_event')
 list(cancer=a,other=b,bh_cancer=basehaz(a,centered=TRUE),bh_other=basehaz(b,centered=TRUE))
}
fast_predict<-function(d,fit,exposure,levels){
 # Exact aggregation of identical adjustment covariates, not subsampling.
 u<-d[,.(weight=.N),by=cov_schema]
 w<-u$weight/sum(u$weight)
 ts<-sort(unique(c(fit$bh_cancer$time[fit$bh_cancer$time<=60],fit$bh_other$time[fit$bh_other$time<=60],12,36,60)))
 dc<-diff(c(0,cumhaz_at(fit$bh_cancer,ts)));do<-diff(c(0,cumhaz_at(fit$bh_other,ts)))
 res<-list()
 for(age in c(50,75))for(level in levels){
  nd<-copy(u);nd[,age_years:=age];set(nd,j=exposure,value=factor(level,levels=levels))
  X<-design(nd,exposure)
  # The same centering is used by centered baseline hazards.
  rc<-exp(drop(X%*%coef(fit$cancer))-sum(fit$cancer$means*coef(fit$cancer)))
  ro<-exp(drop(X%*%coef(fit$other))-sum(fit$other$means*coef(fit$other)))
  ss<-rep(1,nrow(nd));cc<-rep(0,nrow(nd))
  for(k in seq_along(ts)){
   tt<-dc[k]*rc+do[k]*ro
   frac<-ifelse(tt>0,dc[k]*rc/tt,0)
   cc<-cc+ss*(-expm1(-tt))*frac;ss<-ss*exp(-tt)
   if(ts[k]%in%c(12,36,60))res[[length(res)+1L]]<-data.table(horizon_months=ts[k],age=age,exposure_level=level,risk=sum(cc*w))
  }
 }
 ans<-rbindlist(res);stopifnot(all(is.finite(ans$risk)),all(ans$risk>=0 & ans$risk<=1));ans
}
decomp<-function(z,targ,theta,label){
 fit<-fit_ridge(z,'schema3',theta)
 p<-fast_predict(targ,fit,'schema3',levels(target$schema3))
 e<-copy(z[schema3!='iCCA']);e[,distal:=as.integer(schema3=='Distal eCCA')]
 cm<-glm(as.formula(paste('distal ~',ns_term)),data=e,family=binomial())
 stopifnot(cm$converged)
 w<-as.numeric(predict(cm,newdata=data.frame(age_years=c(50,75)),type='response'))
 r<-dcast(p,horizon_months+age~exposure_level,value.var='risk')
 r[,weight:=w[match(age,c(50,75))]]
 r[,RD:=(1-weight)*`Perihilar eCCA`+weight*`Distal eCCA`-iCCA]
 ans<-r[,{
 a<-.SD[age==50];b<-.SD[age==75];wm<-(a$weight+b$weight)/2
 v<-c((b$weight-a$weight)*(a$`Distal eCCA`+b$`Distal eCCA`-a$`Perihilar eCCA`-b$`Perihilar eCCA`)/2,
 (1-wm)*(b$`Perihilar eCCA`-a$`Perihilar eCCA`),wm*(b$`Distal eCCA`-a$`Distal eCCA`),-(b$iCCA-a$iCCA))
 total<-b$RD-a$RD;stopifnot(abs(sum(v)-total)<1e-10)
 list(component=c('composition','perihilar_risk','distal_risk','iCCA_reference','total','RD50','RD75'),estimate=c(v,total,a$RD,b$RD))
 },by=horizon_months]
 ans[,analysis:=label];ans
}
analyse<-function(z,theta=.1){
 targ<-z[!is.na(schema3)]
 rr<-list(decomp(targ,targ,theta,'EOD'),decomp(targ[hist==8160L],targ,theta,'EOD_8160'))
 for(def in defs){
  d<-make_definition(z,def)
  p<-binary_differences(fast_predict(targ,fit_ridge(d,'site_group',theta),'site_group',c('iCCA','eCCA')))
  q<-p[,.(component=c('RD50','RD75','total'),estimate=c(risk_difference_eCCA_minus_iCCA[age==50],risk_difference_eCCA_minus_iCCA[age==75],risk_difference_eCCA_minus_iCCA[age==75]-risk_difference_eCCA_minus_iCCA[age==50])),by=horizon_months]
  q[,analysis:=def];rr[[length(rr)+1L]]<-q
 }
 rbindlist(rr)
}
mode<-Sys.getenv('CCA_BOOT_MODE','preflight')
if(mode=='preflight'){
 points<-rbindlist(lapply(c(.01,.1,1),function(th){cat('Theta',th,'\n');q<-analyse(base,th);q[,theta:=th];q}))
 fwrite(points,file.path(out,'penalty_sensitivity_points.csv'))
 # Numerical equivalence test against original unpenalized full-row predictions.
 fit0<-fit_competing_models(target,'schema3',cov_schema)
 old<-predict_grid(target,fit0,'schema3',levels(target$schema3),c(12,36,60),c(50,75))
 new<-fast_predict(target,fit0,'schema3',levels(target$schema3))
 compare<-merge(old,new,by=c('horizon_months','age','exposure_level'))
 stopifnot(max(abs(compare$risk.x-compare$risk.y))<1e-10)
 fwrite(compare,file.path(out,'prediction_equivalence.csv'))
 cat('Exact prediction equivalence PASS\n')
}else{
 worker<-as.integer(Sys.getenv('CCA_WORKER','1'));workers<-as.integer(Sys.getenv('CCA_WORKERS','4'))
 B<-as.integer(Sys.getenv('CCA_PAIRED_B','1000'))
 ids<-split(seq_len(nrow(base)),base$patient_id)
 for(b in seq.int(worker,B,by=workers)){
  path<-file.path(out,sprintf('rep_%04d.rds',b))
  if(file.exists(path))next
  set.seed(202610110L+b)
  z<-base[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]
  logs<-character()
  value<-tryCatch(withCallingHandlers(analyse(z,.1),warning=function(w){logs<<-c(logs,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){logs<<-c(logs,paste('ERROR',conditionMessage(e)));NULL})
  saveRDS(list(replicate=b,value=value,diagnostics=logs),path)
  if(b%%20%in%seq_len(workers))cat('Worker',worker,'replicate',b,'/',B,'warnings',length(logs),'\n')
 }
 cat('Worker',worker,'completed\n')
}
