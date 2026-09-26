suppressPackageStartupMessages({library(data.table);library(survival);library(splines);library(ggplot2)})
source("00_config.R")
ext_dir <- file.path(results_dir,"reporting_extensions")
dir.create(ext_dir,recursive=TRUE,showWarnings=FALSE)
models <- readRDS(file.path(models_dir,"prepublication_models.rds"))
d <- as.data.table(readRDS(file.path(derived_dir,"cca_strict_hazard_2004_2023.rds")))
d5 <- as.data.table(readRDS(file.path(derived_dir,"cca_strict_cif5_2004_2018.rds")))

# Histology-8160-only three-category EOD analysis.
x <- d[`SEER cause-specific death classification`!="Dead (missing/unknown COD)" & age_years<=84]
x[,schema3:=factor(fcase(`EOD Schema ID Recode (2010+)`=="Bile Ducts Intrahepat","iCCA",
 `EOD Schema ID Recode (2010+)`=="Bile Ducts Perihilar","Perihilar eCCA",
 `EOD Schema ID Recode (2010+)`=="Bile Duct Distal","Distal eCCA",default=NA_character_),levels=c("iCCA","Perihilar eCCA","Distal eCCA"))]
x <- x[diagnosis_year>=2010 & !is.na(schema3)]
mf1 <- model.frame(models$schema_fits[["Model 1"]]$fit)
model2_name <- if ("Model 2 without grade" %in% names(models$schema_fits)) "Model 2 without grade" else "Model 2"
mf2 <- model.frame(models$schema_fits[[model2_name]]$fit)
stopifnot(nrow(x)==nrow(mf1),nrow(x)==nrow(mf2))
z1 <- data.frame(time=x$survival_time_months,event=as.integer(x$competing_event==1L),age_years=x$age_years,
 histology=as.integer(x[["Histologic Type ICD-O-3"]]),patient_cluster=mf1[["(cluster)"]],schema3=mf1$schema3,
 sex_factor=mf1$sex_factor,race_factor=mf1$race_factor,marital4=mf1$marital4,income3=mf1$income3,rural3=mf1$rural3,era_schema=mf1$era_schema)
z2 <- cbind(z1,stage4=mf2$stage4,surgery_primary=mf2$surgery_primary,radiation_any=mf2$radiation_any,chemo_binary=mf2$chemo_binary)
age_term <- "ns(age_years,knots=c(52,68,80),Boundary.knots=c(15,84))"
fit8160 <- function(z,covars){
 z<-droplevels(z[z$histology==8160,]); f0<-as.formula(paste0("Surv(time,event)~",age_term,"+schema3+",covars)); f1<-as.formula(paste0("Surv(time,event)~",age_term,"*schema3+",covars))
 a<-coxph(f0,z,ties="efron"); b<-coxph(f1,z,ties="efron"); stat<-2*(as.numeric(logLik(b))-as.numeric(logLik(a))); df<-attr(logLik(b),"df")-attr(logLik(a),"df")
 data.table(n=nrow(z),cancer_deaths=sum(z$event),chisq=stat,df=df,p=pchisq(stat,df,lower.tail=FALSE))
}
tests<-rbind(cbind(model="Model 1",fit8160(z1,"sex_factor+race_factor+marital4+income3+rural3+era_schema")),
 cbind(model="Model 2",fit8160(z2,"sex_factor+race_factor+marital4+income3+rural3+era_schema+stage4+surgery_primary+radiation_any+chemo_binary")))
fwrite(tests,file.path(ext_dir,"eod_histology8160_tests.csv"))

# Nonparametric CIF and cancer-death restricted mean time lost through 60 months.
d5<-d5[age_years<=84 & !is.na(competing_event)]
d5[,age_group3:=factor(fcase(age_years<=39,"15-39",age_years<=64,"40-64",default="65-84"),levels=c("15-39","40-64","65-84"))]
d5[,site_group:=factor(site_group,levels=c("iCCA","eCCA"))]
aj_curve<-function(x,cause,horizon=60){
 z<-as.data.table(x)[,.(n=.N,d_all=sum(competing_event%in%c(1,2)),d_cause=sum(competing_event==cause)),by=survival_time_months]
 setorder(z,survival_time_months); z[,risk:=rev(cumsum(rev(n)))]; z<-z[survival_time_months<=horizon]
 surv<-1;cif<-0;ans<-list(data.table(month=0,cif=0));for(i in seq_len(nrow(z))){if(z$d_all[i]>0){cif<-cif+surv*z$d_cause[i]/z$risk[i];surv<-surv*(1-z$d_all[i]/z$risk[i])};ans[[length(ans)+1L]]<-data.table(month=z$survival_time_months[i],cif=cif)}
 unique(rbindlist(c(ans,list(data.table(month=horizon,cif=cif)))),by="month")
}
rmtl<-function(x,cause=1){q<-aj_curve(x,cause);sum(head(q$cif,-1)*diff(q$month))}
curves<-d5[,rbind(aj_curve(.SD,1)[,cause:="Cancer death"],aj_curve(.SD,2)[,cause:="Other-cause death"]),by=.(age_group3,site_group)]
point<-d5[,.(cancer_rmtl_months=rmtl(.SD,1),other_rmtl_months=rmtl(.SD,2),N=.N,cancer_deaths=sum(competing_event==1),other_deaths=sum(competing_event==2)),by=.(age_group3,site_group)]
pw<-dcast(point,age_group3~site_group,value.var="cancer_rmtl_months");pw[,difference_eCCA_minus_iCCA:=eCCA-iCCA]
B<-as.integer(Sys.getenv("CCA_RMTL_BOOTSTRAP_B",unset="1000"));set.seed(20260926);clusters<-split(seq_len(nrow(d5)),d5[["Patient ID"]])
boot<-rbindlist(lapply(seq_len(B),function(b){ids<-sample.int(length(clusters),length(clusters),replace=TRUE);q<-d5[unlist(clusters[ids],use.names=FALSE),.(rmtl=rmtl(.SD,1)),by=.(age_group3,site_group)];q<-dcast(q,age_group3~site_group,value.var="rmtl");q[,`:=`(replicate=b,difference_eCCA_minus_iCCA=eCCA-iCCA)]}))
ci<-boot[,.(lower95=quantile(difference_eCCA_minus_iCCA,.025,type=6),upper95=quantile(difference_eCCA_minus_iCCA,.975,type=6)),by=age_group3]
fwrite(curves,file.path(ext_dir,"competing_risk_curves.csv"));fwrite(point,file.path(ext_dir,"age_group_competing_rmtl.csv"));fwrite(merge(pw,ci,by="age_group3"),file.path(ext_dir,"cancer_rmtl_differences.csv"))
p<-ggplot(curves[cause=="Other-cause death"],aes(month,100*cif,color=site_group))+geom_step(linewidth=.75)+facet_wrap(~age_group3,nrow=1)+scale_color_manual(values=c(iCCA="#6B6B6B",eCCA="#1F78B4"))+labs(x="Months since diagnosis",y="Other-cause cumulative incidence (%)",color=NULL)+theme_classic(base_size=9)+theme(legend.position="bottom")
ggsave(file.path(ext_dir,"Figure_S2_other_cause_CIF.png"),p,width=183,height=70,units="mm",dpi=600,bg="white")

# Additive-scale contrasts between standardized fixed-age risk differences.
rd_boot_file <- file.path(source_dir,"adjusted_cif_rd_bootstrap_replicates.csv.gz")
rd_point_file <- file.path(tables_dir,"prepub_adjusted_cif5_risk_differences.csv")
if (file.exists(rd_boot_file) && file.exists(rd_point_file)) {
 rb<-fread(rd_boot_file)[target_type=="Fixed age" & target%in%c("50","65","75")]
 rw<-dcast(rb,replicate~target,value.var="risk_difference_eCCA_minus_iCCA")
 rp<-fread(rd_point_file)[target_type=="Fixed age" & target%in%c("50","65","75")]
 rv<-setNames(rp$risk_difference_eCCA_minus_iCCA,rp$target)
 additive<-data.table(
  contrast=c("Age 65 minus age 50","Age 75 minus age 50"),
  estimate=c(rv[["65"]]-rv[["50"]],rv[["75"]]-rv[["50"]]),
  lower95=c(quantile(rw[["65"]]-rw[["50"]],.025,type=6),quantile(rw[["75"]]-rw[["50"]],.025,type=6)),
  upper95=c(quantile(rw[["65"]]-rw[["50"]],.975,type=6),quantile(rw[["75"]]-rw[["50"]],.975,type=6)),
  bootstrap_replicates=nrow(rw))
 fwrite(additive,file.path(ext_dir,"additive_scale_rd_contrasts.csv"))
}
cat("Reporting extensions complete; RMTL bootstrap:",B,"\n")
