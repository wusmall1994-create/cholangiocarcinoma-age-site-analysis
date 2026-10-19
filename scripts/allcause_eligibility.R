s<-readLines(file.path(Sys.getenv('CCA_ANALYSIS_CODE_DIR','.'),'run_unified_ridge.R'));eval(parse(text=s[seq_len(grep('^mode <-',s)[1]-1L)]))
expanded<-copy(full);known<-copy(dat);out<-Sys.getenv('CCA_ALLCAUSE_RESULTS_DIR');stopifnot(nzchar(out));dir.create(out,recursive=TRUE,showWarnings=FALSE)
age_values <- c(40,50,65,75)
covars <- cov_binary
capture_warnings <- function(expr){
    w<-character()
    value<-withCallingHandlers(expr,warning=function(x){w<<-c(w,conditionMessage(x));invokeRestart('muffleWarning')})
    list(value=value,warnings=w)
}
grouped_target <- function(d){
    u<-d[,.(weight=.N),by=covars];u[,weight:=weight/sum(weight)];u
}
predict_os <- function(fit,targ,label){
    u<-grouped_target(targ)
    bh<-basehaz(fit,centered=TRUE)
    H<-cumhaz_at(bh,60)
    z<-rbindlist(lapply(age_values,function(age)rbindlist(lapply(c('iCCA','eCCA'),function(site){
        nd<-copy(u);nd[,age_years:=age];nd[,site_group:=factor(site,levels=c('iCCA','eCCA'))]
        lp<-predict(fit,newdata=nd,type='lp',reference='sample')
        risk<-sum((-expm1(-H*exp(lp)))*nd$weight)
        stopifnot(is.finite(risk),risk>=0,risk<=1)
        data.table(specification=label,age=age,site=site,risk=risk)
    }))))
    z
}
fit_os <- function(d){
    coxph(as.formula(paste('Surv(survival_time_months,os_event) ~',paste(c(paste0(ns_term,' * site_group'),covars),collapse=' + '))),data=d,ties='efron',control=coxph.control(iter.max=100))
}
joint_allcause <- function(fits,targ,label){
    u<-grouped_target(targ)
    hc<-cumhaz_at(fits$bh_cancer,60);ho<-cumhaz_at(fits$bh_other,60)
    rbindlist(lapply(age_values,function(age)rbindlist(lapply(c('iCCA','eCCA'),function(site){
        nd<-copy(u);nd[,age_years:=age];nd[,site_group:=factor(site,levels=c('iCCA','eCCA'))]
        rc<-exp(predict(fits$cancer,newdata=nd,type='lp',reference='sample'))
        ro<-exp(predict(fits$other,newdata=nd,type='lp',reference='sample'))
        risk<-sum((-expm1(-(hc*rc+ho*ro)))*nd$weight)
        stopifnot(is.finite(risk),risk>=0,risk<=1)
        data.table(specification=label,age=age,site=site,risk=risk)
    }))))
}
contrasts <- function(states){
    z<-dcast(states,specification+age~site,value.var='risk');z[,RD:=eCCA-iCCA]
    a<-melt(z,id.vars=c('specification','age'),measure.vars=c('iCCA','eCCA','RD'),variable.name='metric',value.name='estimate')
    a[,contrast:=paste0(as.character(metric),'_age',age)]
    changes<-z[,.(contrast='RD_change_75_minus_50',estimate=RD[age==75]-RD[age==50]),by=specification]
    ans<-rbindlist(list(a[,.(specification,contrast,estimate)],changes))
    # Paired specifications all use the same expanded resample; target effects separated.
    wide<-dcast(ans,contrast~specification,value.var='estimate')
    pairs<-list(direct_model_minus_joint=c('A_direct_known_own','J_joint_known_own'),inclusion_fixed_target=c('B_direct_expanded_known_target','A_direct_known_own'),target_distribution=c('C_direct_expanded_own','B_direct_expanded_known_target'),full_expansion=c('C_direct_expanded_own','A_direct_known_own'))
    pairrows<-rbindlist(lapply(names(pairs),function(n){v<-pairs[[n]];data.table(specification=n,contrast=wide$contrast,estimate=wide[[v[1]]]-wide[[v[2]]])}))
    rbindlist(list(ans,pairrows))
}
analyse <- function(d){
    k<-d[known==TRUE]
    fk<-fit_os(k);fe<-fit_os(d);fj<-fit_competing_models(k,'site_group',covars)
    states<-rbindlist(list(predict_os(fk,k,'A_direct_known_own'),predict_os(fe,k,'B_direct_expanded_known_target'),predict_os(fe,d,'C_direct_expanded_own'),joint_allcause(fj,k,'J_joint_known_own')))
    list(states=states,contrasts=contrasts(states),aliased_coefficients=c(known=sum(is.na(coef(fk))),expanded=sum(is.na(coef(fe))),cancer=sum(is.na(coef(fj$cancer))),other=sum(is.na(coef(fj$other)))))
}
mode <- Sys.getenv('CCA_ALLCAUSE_MODE','point')
if(mode=='point'){res<-capture_warnings(analyse(expanded));fwrite(res$value$contrasts,file.path(out,'allcause_point_contrasts.csv'));fwrite(res$value$states,file.path(out,'allcause_point_states.csv'))}
if(mode=='bootstrap'){
    B<-as.integer(Sys.getenv('CCA_ALLCAUSE_B','1000'));worker<-as.integer(Sys.getenv('CCA_WORKER','1'));workers<-as.integer(Sys.getenv('CCA_WORKERS','1'))
    bd<-file.path(out,'bootstrap');dir.create(bd,showWarnings=FALSE)
    ids<-split(seq_len(nrow(expanded)),expanded$patient_id)
    for(b in seq.int(worker,B,by=workers)){
        p<-file.path(bd,sprintf('rep_%04d.rds',b));if(file.exists(p))next
        set.seed(202610180L+b)
        z<-expanded[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]
        tic<-proc.time()[3]
        ans<-tryCatch(capture_warnings(analyse(z)),error=function(e)list(value=NULL,warnings=paste('ERROR:',conditionMessage(e))))
        saveRDS(list(replicate=b,seed=202610180L+b,value=ans$value,warnings=ans$warnings,seconds=proc.time()[3]-tic),p)
        if(b%%20%in%seq_len(workers))cat('Replicate',b,'/',B,'sec',round(proc.time()[3]-tic,2),'warnings',length(ans$warnings),'\n')
        flush.console()
    }
}
