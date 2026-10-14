s<-readLines('16_paired_bootstrap.R')
eval(parse(text=s[seq_len(grep('^mode<-',s)[1]-1L)]))
fit<-fit_ridge(target,'schema3',.1)
X<-design(target,'schema3')
lp_errors<-sapply(c('cancer','other'),function(k){
 f<-fit[[k]]
 manual<-drop(X%*%coef(f))-sum(f$means*coef(f))
 direct<-as.numeric(predict(f,newdata=list(X=X),type='lp',reference='sample'))
 max(abs(manual-direct))
})
stopifnot(max(lp_errors)<1e-10)
# Full patient rows, independent implementation of the original integration loop.
full<-rbindlist(lapply(c(50,75),function(age)rbindlist(lapply(levels(target$schema3),function(level){
 nd<-copy(target);nd[,age_years:=age];nd[,schema3:=factor(level,levels=levels(target$schema3))]
 X<-design(nd,'schema3')
 rc<-exp(predict(fit$cancer,newdata=list(X=X),type='lp',reference='sample'))
 ro<-exp(predict(fit$other,newdata=list(X=X),type='lp',reference='sample'))
 rbindlist(lapply(c(12,36,60),function(h){
 ts<-sort(unique(c(fit$bh_cancer$time[fit$bh_cancer$time<=h],fit$bh_other$time[fit$bh_other$time<=h],h)))
 dc<-diff(c(0,cumhaz_at(fit$bh_cancer,ts)));do<-diff(c(0,cumhaz_at(fit$bh_other,ts)))
 ss<-rep(1,nrow(nd));cc<-rep(0,nrow(nd))
 for(k in seq_along(ts)){tt<-dc[k]*rc+do[k]*ro;cc<-cc+ss*(1-exp(-tt))*ifelse(tt>0,dc[k]*rc/tt,0);ss<-ss*exp(-tt)}
 data.table(horizon_months=h,age=age,exposure_level=level,risk=mean(cc))
 }))
}))))
fast<-fast_predict(target,fit,'schema3',levels(target$schema3))
cmp<-merge(full,fast,by=c('horizon_months','age','exposure_level'))
stopifnot(max(abs(cmp$risk.x-cmp$risk.y))<1e-10)
fwrite(cmp,file.path(out,'ridge_prediction_equivalence.csv'))
print(lp_errors)
cat('Ridge prediction full-row validation PASS; maximum risk error',max(abs(cmp$risk.x-cmp$risk.y)),'\n')
