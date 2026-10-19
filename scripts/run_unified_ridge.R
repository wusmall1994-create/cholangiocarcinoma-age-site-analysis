options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(splines)
})

derived_dir <- Sys.getenv("CCA_DERIVED_DIR"); stopifnot(nzchar(derived_dir), dir.exists(derived_dir))

prepare_covariates <- function(x) {
  x <- copy(as.data.table(x))
  x[, patient_id := as.character(`Patient ID`)]
  x[, site_group := factor(site_group, levels = c("iCCA", "eCCA"))]
  x[, sex_factor := factor(Sex, levels = c("Female", "Male"))]
  x[, race_factor := factor(
    `Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)`,
    levels = c(
      "Non-Hispanic White", "Non-Hispanic Black",
      "Non-Hispanic Asian or Pacific Islander",
      "Non-Hispanic American Indian/Alaska Native",
      "Hispanic (All Races)", "Non-Hispanic Unknown Race"
    )
  )]
  x[, marital4 := factor(fcase(
    `Marital status at diagnosis` == "Married (including common law)", "Married",
    `Marital status at diagnosis` == "Widowed", "Widowed",
    `Marital status at diagnosis` == "Unknown", "Unknown",
    default = "Unmarried"
  ), levels = c("Married", "Unmarried", "Widowed", "Unknown"))]
  income_label <- x[["Median household income inflation adj to 2024"]]
  income_first <- suppressWarnings(as.numeric(gsub(",", "", sub(".*?([0-9][0-9,]*).*", "\\1", income_label))))
  income_first[grepl("Unknown", income_label, ignore.case = TRUE)] <- NA_real_
  x[, income3 := factor(fcase(
    is.na(income_first), "Unknown",
    income_first < 70000, "<70,000",
    income_first < 100000, "70,000-99,999",
    default = ">=100,000"
  ), levels = c("<70,000", "70,000-99,999", ">=100,000", "Unknown"))]
  rural_label <- x[["Rural-Urban Continuum Code"]]
  x[, rural3 := factor(fcase(
    grepl("^Counties in metropolitan", rural_label), "Metropolitan",
    grepl("^Nonmetropolitan", rural_label), "Nonmetropolitan",
    default = "Unknown"
  ), levels = c("Metropolitan", "Nonmetropolitan", "Unknown"))]
  x[, era4 := factor(fcase(
    diagnosis_year <= 2009L, "2004-2009",
    diagnosis_year <= 2015L, "2010-2015",
    diagnosis_year <= 2019L, "2016-2019",
    default = "2020-2023"
  ), levels = c("2004-2009", "2010-2015", "2016-2019", "2020-2023"))]
  x[, era_schema_risk := factor(
    fifelse(diagnosis_year <= 2015L, "2010-2015", "2016-2018"),
    levels = c("2010-2015", "2016-2018")
  )]
  x[, schema3 := factor(fcase(
    `EOD Schema ID Recode (2010+)` == "Bile Ducts Intrahepat", "iCCA",
    `EOD Schema ID Recode (2010+)` == "Bile Ducts Perihilar", "Perihilar eCCA",
    `EOD Schema ID Recode (2010+)` == "Bile Duct Distal", "Distal eCCA",
    default = NA_character_
  ), levels = c("iCCA", "Perihilar eCCA", "Distal eCCA"))]
  x[, cause_known := !is.na(competing_event)]
  x[, cs_event := as.integer(competing_event == 1L)]
  x[, other_event := as.integer(competing_event == 2L)]
  x
}

d_hazard <- prepare_covariates(readRDS(file.path(derived_dir, "cca_strict_hazard_2004_2023.rds")))
d_cif <- prepare_covariates(readRDS(file.path(derived_dir, "cca_strict_cif5_2004_2018.rds")))

age_reference <- d_hazard[age_topcoded_85plus == FALSE]
age_knots <- as.numeric(quantile(age_reference$age_years, c(0.10, 0.50, 0.90), na.rm = TRUE, type = 2))
age_bounds <- range(age_reference$age_years, na.rm = TRUE)
numvec <- function(x) paste(format(x, scientific = FALSE, trim = TRUE), collapse = ", ")
ns_term <- paste0(
  "splines::ns(age_years, knots = c(", numvec(age_knots),
  "), Boundary.knots = c(", numvec(age_bounds), "))"
)

cov_binary <- c("sex_factor", "race_factor", "marital4", "income3", "rural3", "era4")
cov_schema <- c("sex_factor", "race_factor", "marital4", "income3", "rural3", "era_schema_risk")
ages <- c(40, 50, 65, 75)
horizons <- c(12, 36, 60)

fit_competing_models <- function(data, exposure, covars) {
  rhs <- paste(c(paste0(ns_term, " * ", exposure), covars), collapse = " + ")
  cancer <- coxph(as.formula(paste0("Surv(survival_time_months, cs_event) ~ ", rhs)),
                   data = data, ties = "efron")
  other <- coxph(as.formula(paste0("Surv(survival_time_months, other_event) ~ ", rhs)),
                  data = data, ties = "efron")
  list(
    cancer = cancer,
    other = other,
    bh_cancer = suppressWarnings(basehaz(cancer, centered = TRUE)),
    bh_other = suppressWarnings(basehaz(other, centered = TRUE))
  )
}

cumhaz_at <- function(base_hazard, times) {
  idx <- findInterval(times, base_hazard$time)
  c(0, base_hazard$hazard)[idx + 1L]
}

standardized_cif <- function(fits, newdata, horizon) {
  event_times <- sort(unique(c(
    fits$bh_cancer$time[fits$bh_cancer$time <= horizon],
    fits$bh_other$time[fits$bh_other$time <= horizon], horizon
  )))
  hc <- cumhaz_at(fits$bh_cancer, event_times)
  ho <- cumhaz_at(fits$bh_other, event_times)
  dhc <- diff(c(0, hc))
  dho <- diff(c(0, ho))
  rc <- exp(predict(fits$cancer, newdata = newdata, type = "lp", reference = "sample"))
  ro <- exp(predict(fits$other, newdata = newdata, type = "lp", reference = "sample"))
  surv <- rep(1, nrow(newdata))
  cif <- rep(0, nrow(newdata))
  for (k in seq_along(event_times)) {
    total <- dhc[k] * rc + dho[k] * ro
    fraction <- ifelse(total > 0, dhc[k] * rc / total, 0)
    cif <- cif + surv * (1 - exp(-total)) * fraction
    surv <- surv * exp(-total)
  }
  mean(cif)
}

predict_grid <- function(data, fits, exposure, exposure_levels, horizons, ages) {
  rbindlist(lapply(horizons, function(horizon) {
    rbindlist(lapply(ages, function(age) {
      rbindlist(lapply(exposure_levels, function(level) {
        nd <- copy(data)
        nd[, age_years := age]
        set(nd, j = exposure, value = factor(level, levels = exposure_levels))
        data.table(
          horizon_months = horizon,
          age = age,
          exposure_level = level,
          risk = standardized_cif(fits, as.data.frame(nd), horizon)
        )
      }))
    }))
  }))
}

binary_differences <- function(long) {
  wide <- dcast(long, horizon_months + age ~ exposure_level, value.var = "risk")
  wide[, risk_difference_eCCA_minus_iCCA := eCCA - iCCA]
  wide
}

schema_differences <- function(long) {
  wide <- dcast(long, horizon_months + age ~ exposure_level, value.var = "risk")
  wide[, `:=`(
    rd_perihilar_minus_iCCA = `Perihilar eCCA` - iCCA,
    rd_distal_minus_iCCA = `Distal eCCA` - iCCA,
    contrast_distal_minus_perihilar = `Distal eCCA` - `Perihilar eCCA`
  )]
  wide
}

binary_data <- droplevels(d_cif[age_topcoded_85plus == FALSE & cause_known == TRUE])
schema_data <- droplevels(binary_data[diagnosis_year >= 2010L & !is.na(schema3)])


expanded <- droplevels(copy(d_cif[age_topcoded_85plus==FALSE]))
expanded[,known:=!is.na(competing_event)]
stopifnot(all(is.finite(expanded$survival_time_months)),all(expanded$survival_time_months>0),all(expanded$os_event %in% 0:1))
known <- copy(expanded[known==TRUE])
stopifnot(nrow(known)==15578L,identical(age_knots,c(52,68,80)))
stopifnot(all(known$os_event==as.integer(known$competing_event!=0L)))
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

setDTthreads(1)
newraw<-fread(Sys.getenv('CCA_COD_CSV'),colClasses='character',check.names=FALSE)
makekey<-function(x)paste(x[['Patient ID']],x[['Record number recode']],sep=':')
stopifnot(!anyDuplicated(makekey(newraw)))
attach_cod<-function(x){
    x<-copy(x);j<-match(makekey(x),makekey(newraw));stopifnot(!anyNA(j))
    x[,cod:=newraw[['COD to site recode']][j]]
    x[,cod_expanded:=newraw[['COD to site recode ICD-O-3 2023 Revision Expanded (1999+)']][j]]
    x[,updated_race:=newraw[['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)']][j]]
    x[,raw_cod_known:=!cod%chin%c('State DC not available or state DC available but no COD','Unknown/missing/invalid COD')]
    x
}
full<-attach_cod(expanded)
dat<-copy(full[known==TRUE]);stopifnot(nrow(dat)==15578L,all(dat$raw_cod_known))
stopifnot(all((full$cod=='Alive')==(full$os_event==0L)))
stopifnot(all((full$cod=='State DC not available or state DC available but no COD')==(!full$known)))
endpoint<-function(x,label){
    z<-copy(x)
    if(label=='Common_biliary')z[,cs_event:=as.integer(os_event==1L & cod%chin%c('Intrahepatic Bile Duct','Other Biliary'))]
    if(label=='Common_hepatobiliary')z[,cs_event:=as.integer(os_event==1L & cod%chin%c('Liver','Intrahepatic Bile Duct','Gallbladder','Other Biliary'))]
    z[,other_event:=os_event-cs_event]
    stopifnot(all(z$cs_event%in%0:1),all(z$other_event%in%0:1),all(z$cs_event+z$other_event==z$os_event))
    z
}
defs<-c('SEER_original','Common_biliary','Common_hepatobiliary')
capture_stage<-function(expr,label){
    withCallingHandlers(expr,warning=function(w){warning_log[[length(warning_log)+1L]]<<-data.table(stage=label,message=conditionMessage(w));invokeRestart('muffleWarning')})
}
fit_pair<-function(z,label){
    rhs<-paste(c(paste0(ns_term,' * site_group'),covars),collapse=' + ')
    a<-capture_stage(coxph(as.formula(paste('Surv(survival_time_months,cs_event) ~',rhs)),data=z,ties='efron',control=coxph.control(iter.max=100)),paste(label,'endpoint_fit'))
    b<-capture_stage(coxph(as.formula(paste('Surv(survival_time_months,other_event) ~',rhs)),data=z,ties='efron',control=coxph.control(iter.max=100)),paste(label,'other_fit'))
    list(cancer=a,other=b,bh_cancer=suppressWarnings(basehaz(a,centered=TRUE)),bh_other=suppressWarnings(basehaz(b,centered=TRUE)))
}
predict_states_cod<-function(fits,targ,age_values,label){
    u<-grouped_target(targ);w<-u$weight
    times<-sort(unique(c(fits$bh_cancer$time[fits$bh_cancer$time<=60],fits$bh_other$time[fits$bh_other$time<=60],60)))
    dc<-diff(c(0,cumhaz_at(fits$bh_cancer,times)));do<-diff(c(0,cumhaz_at(fits$bh_other,times)))
    rbindlist(lapply(age_values,function(age)rbindlist(lapply(c('iCCA','eCCA'),function(site){
        nd<-copy(u);nd[,age_years:=age];nd[,site_group:=factor(site,levels=c('iCCA','eCCA'))]
        rc<-exp(predict(fits$cancer,newdata=nd,type='lp',reference='sample'));ro<-exp(predict(fits$other,newdata=nd,type='lp',reference='sample'))
        surv<-rep(1,nrow(nd));cancer<-other<-rep(0,nrow(nd))
        for(k in seq_along(times)){
            hc<-dc[k]*rc;ho<-do[k]*ro;hh<-hc+ho;mass<-surv*(-expm1(-hh));share<-ifelse(hh>0,hc/hh,0)
            cancer<-cancer+mass*share;other<-other+mass*(1-share);surv<-surv*exp(-hh)
        }
        vals<-c(cancer=sum(cancer*w),other=sum(other*w),alive=sum(surv*w),all_death=sum((cancer+other)*w))
        stopifnot(all(is.finite(vals)),all(vals>=0 & vals<=1),abs(sum(vals[1:3])-1)<1e-10)
        direct_total<-sum((-expm1(-(cumhaz_at(fits$bh_cancer,60)*rc+cumhaz_at(fits$bh_other,60)*ro)))*w)
        stopifnot(abs(vals['all_death']-direct_total)<1e-10)
        data.table(definition=label,age=age,site=site,state=names(vals),risk=as.numeric(vals))
    }))))
}
contrasts_cod<-function(states){
    z<-dcast(states,definition+state+age~site,value.var='risk');z[,RD:=eCCA-iCCA]
    p<-melt(z,id.vars=c('definition','state','age'),measure.vars=c('iCCA','eCCA','RD'),variable.name='metric',value.name='estimate')
    p[,contrast:=paste0(metric,'_age',age)]
    ans<-rbindlist(list(p[,.(definition,state,contrast,estimate)],z[,.(contrast='RD_change_75_minus_50',estimate=RD[age==75]-RD[age==50]),by=.(definition,state)]))
    wide<-dcast(ans,state+contrast~definition,value.var='estimate')
    paired<-rbindlist(lapply(defs[-1],function(def)data.table(definition=paste0(def,'_minus_SEER'),state=wide$state,contrast=wide$contrast,estimate=wide[[def]]-wide$SEER_original)))
    rbindlist(list(ans,paired))
}
analyse_cod<-function(z,ages){
    states<-rbindlist(lapply(defs,function(label){x<-endpoint(z,label);predict_states_cod(fit_pair(x,label),z,ages,label)}))
    list(states=states,contrasts=contrasts_cod(states))
}

rhs<-as.formula(paste('~',paste(c(paste0(ns_term,' * site_group'),covars),collapse=' + ')))
ridge_pair<-function(z,label){
    X<-model.matrix(rhs,data=z)[,-1,drop=FALSE]
    fit<-function(ev){
        zz<-list(y=Surv(z$survival_time_months,z[[ev]]),X=X)
        f<-capture_stage(coxph(y~ridge(X,theta=.1,scale=FALSE),data=zz,ties='efron',x=TRUE,model=TRUE,control=coxph.control(iter.max=100,eps=1e-11)),paste(label,ev,'ridge'))
        stopifnot(all(is.finite(coef(f))))
        f
    }
    a<-fit('cs_event');b<-fit('other_event')
    list(cancer=a,other=b,bh_cancer=basehaz(a,centered=TRUE),bh_other=basehaz(b,centered=TRUE))
}
pred<-function(f,targ,label){
    u<-grouped_target(targ);w<-u$weight
    times<-sort(unique(c(f$bh_cancer$time[f$bh_cancer$time<=60],f$bh_other$time[f$bh_other$time<=60],60)))
    dc<-diff(c(0,cumhaz_at(f$bh_cancer,times)));do<-diff(c(0,cumhaz_at(f$bh_other,times)))
    rbindlist(lapply(c(50,75),function(age)rbindlist(lapply(c('iCCA','eCCA'),function(site){
        nd<-copy(u);nd[,age_years:=age];nd[,site_group:=factor(site,levels=c('iCCA','eCCA'))]
        X<-model.matrix(rhs,data=nd)[,-1,drop=FALSE]
        rc<-exp(drop(X%*%coef(f$cancer))-sum(f$cancer$means*coef(f$cancer)))
        ro<-exp(drop(X%*%coef(f$other))-sum(f$other$means*coef(f$other)))
        sv<-rep(1,nrow(nd));cc<-oo<-rep(0,nrow(nd))
        for(k in seq_along(times)){
            hc<-dc[k]*rc;ho<-do[k]*ro;hh<-hc+ho;mass<-sv*(-expm1(-hh));fr<-ifelse(hh>0,hc/hh,0)
            cc<-cc+mass*fr;oo<-oo+mass*(1-fr);sv<-sv*exp(-hh)
        }
        v<-c(cancer=sum(cc*w),other=sum(oo*w),alive=sum(sv*w),all_death=sum((cc+oo)*w))
        stopifnot(all(is.finite(v)),abs(sum(v[1:3])-1)<1e-10)
        data.table(definition=label,age=age,site=site,state=names(v),risk=as.numeric(v))
    }))))
}

endpoint_out <- Sys.getenv("CCA_RESULTS_DIR"); stopifnot(nzchar(endpoint_out)); dir.create(endpoint_out,recursive=TRUE,showWarnings=FALSE)
full[,sequence_group:=ifelse(`Sequence number`=='One primary only','Single_primary','Multiple_primary_sequence')]
dat <- copy(full[known==TRUE])
full[,selected621:=known & site_group=='eCCA' & other_event==1L & cod=='Intrahepatic Bile Duct']
full[is.na(selected621),selected621:=FALSE]
stopifnot(sum(full$selected621)==621L)
dat <- copy(full[known==TRUE])
all_contrasts <- function(z){
 k <- z[known==TRUE]
 ans <- rbindlist(lapply(c('All','Single_primary','Multiple_primary_sequence'),function(g){
  x <- if(g=='All')copy(k) else copy(k[sequence_group==g])
  st <- rbindlist(lapply(defs,function(label)pred(ridge_pair(endpoint(x,label),paste(g,label)),x,label)))
  # Scenario changes classifications only; same records and empirical target.
  if(g=='All'){
   sc <- endpoint(x,'SEER_original');sc[selected621==TRUE,`:=`(cs_event=1L,other_event=0L)]
   ss <- pred(ridge_pair(sc,'Scenario621'),x,'Scenario621')
  }
  ct <- contrasts_cod(st)
  if(g=='All'){
   w <- dcast(ss,state+age~site,value.var='risk');w[,RD:=eCCA-iCCA]
   scct <- rbindlist(list(w[,.(state,contrast=paste0('RD_age',age),estimate=RD)],w[,.(contrast='RD_change_75_minus_50',estimate=RD[age==75]-RD[age==50]),by=state]))
   scct[,definition:='Scenario621'];ref<-ct[definition=='SEER_original']
   delta<-merge(scct,ref,by=c('state','contrast'),suffixes=c('_s','_r'))
   scct<-rbindlist(list(scct,delta[,.(state,contrast,estimate=estimate_s-estimate_r,definition='Scenario621_minus_SEER')]),use.names=TRUE)
   ct<-rbindlist(list(ct,scct),use.names=TRUE)
  }
  ct[,sequence_group:=g];ct
 }))
 # A direct all-cause model uses all eligible records, including unknown COD.
 X<-model.matrix(rhs,data=z)[,-1,drop=FALSE]
 f<-capture_stage(coxph(y~ridge(X,theta=.1,scale=FALSE),data=list(y=Surv(z$survival_time_months,z$os_event),X=X),ties='efron',x=TRUE,model=TRUE,control=coxph.control(iter.max=100,eps=1e-11)),'Direct_allcause')
 u<-grouped_target(z);H<-cumhaz_at(basehaz(f,centered=TRUE),60)
 st<-rbindlist(lapply(c(50,75),function(age)rbindlist(lapply(c('iCCA','eCCA'),function(site){
  nd<-copy(u);nd[,age_years:=age];nd[,site_group:=factor(site,levels=c('iCCA','eCCA'))]
  X<-model.matrix(rhs,data=nd)[,-1,drop=FALSE]
  lp<-drop(X%*%coef(f))-sum(f$means*coef(f))
  data.table(age=age,site=site,risk=sum((-expm1(-H*exp(lp)))*u$weight))
 }))))
 w<-dcast(st,age~site,value.var='risk');w[,RD:=eCCA-iCCA]
 ct<-rbindlist(list(w[,.(contrast=paste0('RD_age',age),estimate=RD)],data.table(contrast='RD_change_75_minus_50',estimate=w[age==75,RD]-w[age==50,RD])))
 ct[,`:=`(definition='Direct_allcause_including_unknown',state='all_death',sequence_group='Expanded')]
 ans<-rbindlist(list(ans,ct),use.names=TRUE)
 stopifnot(all(is.finite(ans$estimate)),!anyDuplicated(ans[,.(sequence_group,definition,state,contrast)]))
 ans
}
mode <- Sys.getenv('CCA_ENDPOINT_MODE','point')
if(mode=='point'){
 warning_log<-list();tic<-proc.time()[3];ans<-all_contrasts(full)
 fwrite(ans,file.path(endpoint_out,'ridge_point_contrasts.csv'))
 fwrite(if(length(warning_log))rbindlist(warning_log) else data.table(stage=character(),message=character()),file.path(endpoint_out,'ridge_point_warnings.csv'))
 print(ans[state%chin%c('cancer','all_death') & contrast%chin%c('RD_age75','RD_change_75_minus_50')]);cat('SECONDS',proc.time()[3]-tic,'WARNINGS',length(warning_log),'\n')
 # Bidirectional transitions at recorded follow-up, and restricted to deaths <=60 months.
 tables<-rbindlist(lapply(c('Recorded_followup','Within_60_months'),function(scope){
  x<-if(scope=='Recorded_followup')copy(dat) else copy(dat[os_event==0L | survival_time_months<=60])
  x[,age_band:=cut(age_years,c(15,40,50,65,75,85),right=FALSE,labels=c('15-39','40-49','50-64','65-74','75-84'))]
  rbindlist(lapply(defs[-1],function(de){v<-endpoint(x,de);v[,seer_event:=x$cs_event]
   v[,.(N=.N,mean_age=mean(age_years),median_age=as.numeric(median(age_years)),age75plus=sum(age_years>=75)),by=.(site_group,sequence_group,age_band,os_event,seer_event,common_event=cs_event)][,`:=`(definition=de,scope=scope)]
  }))
 }))
 fwrite(tables,file.path(endpoint_out,'bidirectional_transitions_by_age.csv'))
 fwrite(dat[selected621==TRUE,.(N=.N,mean_age=mean(age_years),median_age=as.numeric(median(age_years)),min_month=min(survival_time_months),max_month=max(survival_time_months),deaths_within60=sum(survival_time_months<=60)),by=.(sequence_group)],file.path(endpoint_out,'scenario621_audit.csv'))
 # Map warning coefficient indexes to design columns, with event support by model.
 diagnostic<-list();ii<-0L
 for(g in c('All','Single_primary','Multiple_primary_sequence'))for(de in defs){
  x<-if(g=='All')copy(dat) else copy(dat[sequence_group==g]);x<-endpoint(x,de);X<-model.matrix(rhs,data=x)[,-1,drop=FALSE]
  for(ev in c('cs_event','other_event')){
   warnings_now<-character();f<-withCallingHandlers(coxph(as.formula(paste('Surv(survival_time_months,',ev,') ~',paste(c(paste0(ns_term,' * site_group'),covars),collapse=' + '))),data=x,ties='efron',control=coxph.control(iter.max=100)),warning=function(w){warnings_now<<-c(warnings_now,conditionMessage(w));invokeRestart('muffleWarning')})
   ii<-ii+1L;diagnostic[[ii]]<-data.table(group=g,definition=de,event=ev,index=seq_len(ncol(X)),term=colnames(X),n_nonzero=colSums(X!=0),events_nonzero=colSums((X!=0)*x[[ev]]),coefficient=coef(f),warning=paste(warnings_now,collapse=' | '))
  }
 }
 fwrite(rbindlist(diagnostic),file.path(endpoint_out,'unpenalized_coefficient_support_diagnostics.csv'))
 capture.output(sessionInfo(),file=file.path(endpoint_out,'sessionInfo.txt'))
}
if(mode=='bootstrap'){
 worker<-as.integer(Sys.getenv('CCA_ENDPOINT_WORKER','1'));workers<-as.integer(Sys.getenv('CCA_ENDPOINT_WORKERS','8'));B<-as.integer(Sys.getenv('CCA_ENDPOINT_B','2000'))
 bd<-file.path(endpoint_out,'bootstrap');dir.create(bd,showWarnings=FALSE)
 ids<-split(seq_len(nrow(full)),full$patient_id)
 for(b in seq.int(worker,B,by=workers)){
  p<-file.path(bd,sprintf('rep_%04d.rds',b));if(file.exists(p))next
  set.seed(202610192L+b);z<-full[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]
  warning_log<-list();tic<-proc.time()[3]
  value<-tryCatch(all_contrasts(z),error=function(e){warning_log[[length(warning_log)+1L]]<<-data.table(stage='ERROR',message=conditionMessage(e));NULL})
  saveRDS(list(replicate=b,seed=202610192L+b,value=value,warnings=warning_log,seconds=proc.time()[3]-tic),p)
  if(b%%40<workers)cat('rep',b,'seconds',round(proc.time()[3]-tic,2),'warnings',length(warning_log),'\n');flush.console()
 }
}
