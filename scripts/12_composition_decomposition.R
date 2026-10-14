options(stringsAsFactors=FALSE)
# Reuse validated preparation and competing-risk prediction functions only.
src <- readLines('11_horizon_anatomy_absolute_risks.R')
stop_line <- grep('^binary_point_long <-', src)[1]-1L
eval(parse(text=src[seq_len(stop_line)]))
out_dir <- file.path(results_dir,'composition_decomposition')
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
# All definitions share the 2010-2018 source cohort and standardization target.
target <- copy(schema_data)
base <- copy(binary_data[diagnosis_year>=2010L])
base[,hist := as.integer(`Histologic Type ICD-O-3`)]
make_definition <- function(x,def) {
  z<-copy(x)
  if(def=='D2_exclude_discordant') z<-z[!(site_group=='iCCA' & hist==8162L)]
  if(def=='D3_8160_only') z<-z[hist==8160L]
  if(def=='D4_reallocate_8162') z[hist==8162L,site_group:=factor('eCCA',levels=c('iCCA','eCCA'))]
  z<-droplevels(z)
  z
}
defs<-c('D1_topography','D2_exclude_discordant','D3_8160_only','D4_reallocate_8162')
definition_results<-rbindlist(lapply(defs,function(def){
 z<-make_definition(base,def)
 fit<-fit_competing_models(z,'site_group',cov_schema)
 p<-predict_grid(target,fit,'site_group',c('iCCA','eCCA'),c(12,36,60),c(50,75))
 r<-binary_differences(p)
 r[,`:=`(definition=def,N=nrow(z))]
 r
}))
fwrite(definition_results,file.path(out_dir,'definition_risk_point_estimates.csv'))

# Marginal composition-weighted standardized risks define a descriptive mixture.
# This mixture is not assumed equal to a separately fitted pooled binary model.
# Symmetric Kitagawa allocation exactly partitions the age-75 minus age-50 RD.
decompose <- function(z,standardize_to=z,def='EOD') {
 fit<-fit_competing_models(z,'schema3',cov_schema)
 p<-predict_grid(standardize_to,fit,'schema3',levels(z$schema3),c(12,36,60),c(50,75))
 e<-copy(z[schema3!='iCCA'])
 e[,distal:=as.integer(schema3=='Distal eCCA')]
 composition<-glm(as.formula(paste('distal ~',ns_term)),data=e,family=binomial())
 w<-as.numeric(predict(composition,newdata=data.frame(age_years=c(50,75)),type='response'))
 r<-dcast(p,horizon_months+age~exposure_level,value.var='risk')
 r[,distal_weight:=w[match(age,c(50,75))]]
 r[,mixture_RD:=(1-distal_weight)*`Perihilar eCCA`+distal_weight*`Distal eCCA`-iCCA]
 ans<-r[,{
  a<-.SD[age==50];b<-.SD[age==75]
  wm<-(a$distal_weight+b$distal_weight)/2
  comp<-(b$distal_weight-a$distal_weight)*((a$`Distal eCCA`+b$`Distal eCCA`)-(a$`Perihilar eCCA`+b$`Perihilar eCCA`))/2
  peri<-(1-wm)*(b$`Perihilar eCCA`-a$`Perihilar eCCA`)
  dist<-wm*(b$`Distal eCCA`-a$`Distal eCCA`)
  intra<- -(b$iCCA-a$iCCA)
  total<-b$mixture_RD-a$mixture_RD
  stopifnot(abs(total-comp-peri-dist-intra)<1e-10)
  list(component=c('composition','perihilar_risk','distal_risk','iCCA_reference','total'),estimate=c(comp,peri,dist,intra,total),weight50=a$distal_weight,weight75=b$distal_weight,RD50=a$mixture_RD,RD75=b$mixture_RD)
 },by=horizon_months]
 ans[,definition:=def]
 ans
}
point<-rbindlist(list(decompose(target,def='EOD'),decompose(target[as.integer(`Histologic Type ICD-O-3`)==8160],target,def='EOD_8160')))
fwrite(point,file.path(out_dir,'decomposition_point_estimates.csv'))
print(point)
B<-as.integer(Sys.getenv('CCA_DECOMPOSITION_B','200'))
set.seed(20261011)
ids<-split(seq_len(nrow(target)),target$patient_id)
boot<-vector('list',B)
diagnostics<-list()
for(b in seq_len(B)) {
 z<-target[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]
 boot[[b]]<-tryCatch(withCallingHandlers({q<-decompose(z);q[,replicate:=b];q},warning=function(w){
   diagnostics[[length(diagnostics)+1L]]<<-data.table(replicate=b,type='warning',message=conditionMessage(w))
   invokeRestart('muffleWarning')
 }),error=function(e){
   diagnostics[[length(diagnostics)+1L]]<<-data.table(replicate=b,type='error',message=conditionMessage(e))
   NULL
 })
 if(length(diagnostics)) fwrite(rbindlist(diagnostics),file.path(out_dir,'bootstrap_diagnostics.csv'))
 if(b%%20==0) {cat('Bootstrap',b,'/',B,'\n');flush.console()}
}
bb<-rbindlist(boot)
fwrite(bb,file.path(out_dir,'decomposition_bootstrap.csv'))
ci<-bb[,.(lower95=quantile(estimate,.025),upper95=quantile(estimate,.975),successful=uniqueN(replicate)),by=.(horizon_months,component)]
fwrite(merge(point[definition=='EOD'],ci,by=c('horizon_months','component')),file.path(out_dir,'decomposition_with_ci.csv'))
cat('Completed',uniqueN(bb$replicate),'of',B,'replicates\n')
