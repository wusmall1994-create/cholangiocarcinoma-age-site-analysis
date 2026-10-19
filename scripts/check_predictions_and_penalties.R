s<-readLines(file.path(Sys.getenv('CCA_ANALYSIS_CODE_DIR','.'),'run_unified_ridge.R'));eval(parse(text=s[seq_len(grep('^mode <-',s)[1]-1L)]))
warning_log<-list()
checks<-list()
# Independent survival-form identity and direct predict(type='risk') agreement.
for(g in c('All','Single_primary','Multiple_primary_sequence'))for(de in defs){
 x<-if(g=='All')copy(dat) else copy(dat[sequence_group==g])
 f<-ridge_pair(endpoint(x,de),paste(g,de));st<-pred(f,x,de)
 for(age in c(50,75))for(site in c('iCCA','eCCA')){
  nd<-copy(x);nd[,age_years:=age];nd[,site_group:=factor(site,levels=c('iCCA','eCCA'))]
  X<-model.matrix(rhs,data=nd)[,-1,drop=FALSE]
  rc<-exp(drop(X%*%coef(f$cancer))-sum(f$cancer$means*coef(f$cancer)))
  ro<-exp(drop(X%*%coef(f$other))-sum(f$other$means*coef(f$other)))
  # Full rows, independent of grouped target and iterative CIF integration.
  os<-mean(-expm1(-(cumhaz_at(f$bh_cancer,60)*rc+cumhaz_at(f$bh_other,60)*ro)))
  ix<-st$age==age & st$site==site & st$state=='all_death';val<-st$risk[ix]
  checks[[length(checks)+1L]]<-data.table(group=g,definition=de,age=age,site=site,check='full-row direct survival versus grouped joint integration',error=abs(os-val))
 }
}
fwrite(rbindlist(checks),file.path(endpoint_out,'independent_prediction_QA.csv'))
stopifnot(max(rbindlist(checks)$error)<1e-10)

source_text<-paste(deparse(ridge_pair),collapse='\n')
ridge_pair<-eval(parse(text=gsub('theta = 0.1','theta = theta_revision',source_text,fixed=TRUE)))
pen<-rbindlist(lapply(c(.01,.1,1),function(th){theta_revision<<-th;warning_log<<-list();a<-all_contrasts(full);a[,theta:=th];a}))
fwrite(pen,file.path(endpoint_out,'alternative_penalty_point_checks.csv'))
cat('Independent prediction checks passed\n')
