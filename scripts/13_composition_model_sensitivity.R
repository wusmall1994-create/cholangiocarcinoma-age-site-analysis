options(stringsAsFactors=FALSE)
src<-readLines('11_horizon_anatomy_absolute_risks.R')
eval(parse(text=src[seq_len(grep('^binary_point_long <-',src)[1]-1L)]))
out_dir<-file.path(results_dir,'composition_decomposition')
z<-copy(schema_data[schema3!='iCCA'])
z[,distal:=as.integer(schema3=='Distal eCCA')]
risks<-fread(file.path(results_dir,'reporting_extensions','schema3_standardized_cif5.csv'))
r<-dcast(risks[age %in% c(50,75)],age~exposure_level,value.var='risk')
a<-r[age==50];b<-r[age==75]
alt_knots<-as.numeric(quantile(z$age_years,c(.25,.50,.75),type=2))
formulas<-list(
 original_spline=as.formula(paste('distal ~',ns_term)),
 linear_age=distal~age_years,
 quartile_spline=as.formula(paste0('distal ~ splines::ns(age_years, knots=c(',numvec(alt_knots),'), Boundary.knots=c(',numvec(age_bounds),'))'))
)
out<-rbindlist(lapply(names(formulas),function(nm){
 fit<-glm(formulas[[nm]],data=z,family=binomial())
 w<-as.numeric(predict(fit,newdata=data.frame(age_years=c(50,75)),type='response'))
 comp<-diff(w)*((a$`Distal eCCA`+b$`Distal eCCA`)-(a$`Perihilar eCCA`+b$`Perihilar eCCA`))/2
 peri<-(1-mean(w))*(b$`Perihilar eCCA`-a$`Perihilar eCCA`)
 distal<-mean(w)*(b$`Distal eCCA`-a$`Distal eCCA`)
 intra<- -(b$iCCA-a$iCCA)
 total<-comp+peri+distal+intra
 data.table(specification=nm,distal_weight50=w[1],distal_weight75=w[2],composition=comp,perihilar=peri,distal=distal,iCCA=intra,total=total,composition_fraction=comp/total,AIC=AIC(fit))
}))
fwrite(out,file.path(out_dir,'composition_model_sensitivity.csv'))
print(out)
