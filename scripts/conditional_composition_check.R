# Point-estimate sensitivity, not a replacement for the original estimand.
Sys.setenv(CCA_CLINICAL_MODE='audit')
source(file.path(Sys.getenv('CCA_CLINICAL_SCRIPTS_DIR'),'competing_outcomes.R'))
fit<-readRDS(file.path(outnew,'local_models_DO_NOT_UPLOAD.rds'))$ph
u<-target[,.(weight=.N),by=cov_schema];u[,weight:=weight/sum(weight)]
e<-copy(target[schema3!='iCCA']);e[,distal:=as.integer(schema3=='Distal eCCA')]
cm<-glm(as.formula(paste('distal ~',ns_term,'+ sex_factor + era_schema_risk')),data=e,family=binomial())
stopifnot(cm$converged,all(is.finite(coef(cm))))
individual<-function(age,level){
 nd<-copy(u);nd[,age_years:=age];nd[,schema3:=factor(level,levels=levels(target$schema3))]
 X<-design(nd,'schema3')
 rc<-exp(drop(X%*%coef(fit$cancer))-sum(fit$cancer$means*coef(fit$cancer)))
 ro<-exp(drop(X%*%coef(fit$other))-sum(fit$other$means*coef(fit$other)))
 times<-sort(unique(c(fit$bh_cancer$time[fit$bh_cancer$time<=60],fit$bh_other$time[fit$bh_other$time<=60],60)))
 dc<-diff(c(0,cumhaz_at(fit$bh_cancer,times)));do<-diff(c(0,cumhaz_at(fit$bh_other,times)))
 ss<-rep(1,nrow(nd));cc<-oo<-rep(0,nrow(nd))
 for(k in seq_along(times)){a<-dc[k]*rc;b<-do[k]*ro;h<-a+b;v<-ss*(-expm1(-h));fr<-ifelse(h>0,a/h,0);cc<-cc+v*fr;oo<-oo+v*(1-fr);ss<-ss*exp(-h)}
 data.table(cancer=cc,other=oo,all_death=cc+oo,alive=ss)
}
r<-lapply(c(50,75),function(age){nd<-copy(u);nd[,age_years:=age];list(w=as.numeric(predict(cm,nd,type='response')),I=individual(age,'iCCA'),P=individual(age,'Perihilar eCCA'),D=individual(age,'Distal eCCA'))})
a<-r[[1]];b<-r[[2]];wm<-(a$w+b$w)/2
result<-rbindlist(lapply(c('cancer','other','all_death','alive'),function(st){
 comp<-(b$w-a$w)*(a$D[[st]]+b$D[[st]]-a$P[[st]]-b$P[[st]])/2
 peri<-(1-wm)*(b$P[[st]]-a$P[[st]])
 distal<-wm*(b$D[[st]]-a$D[[st]])
 intra<- -(b$I[[st]]-a$I[[st]])
 rd50<-(1-a$w)*a$P[[st]]+a$w*a$D[[st]]-a$I[[st]]
 rd75<-(1-b$w)*b$P[[st]]+b$w*b$D[[st]]-b$I[[st]]
 stopifnot(max(abs(comp+peri+distal+intra-(rd75-rd50)))<1e-10)
 data.table(state=st,component=c('composition','perihilar_risk','distal_risk','iCCA_reference','total','RD50','RD75'),estimate=vapply(list(comp,peri,distal,intra,rd75-rd50,rd50,rd75),function(v)sum(v*u$weight),numeric(1)))
}))
result[,`:=`(horizon_months=60,model='conditional_age_sex_era_POINT_ONLY')]
fwrite(result,file.path(outnew,'conditional_composition_POINT_ONLY.csv'))
print(result[state=='cancer'])
