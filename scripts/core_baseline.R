s<-readLines(file.path(Sys.getenv('CCA_ANALYSIS_CODE_DIR','.'),'run_unified_ridge.R'));eval(parse(text=s[seq_len(grep('^mode <-',s)[1]-1L)]))
x<-dat
fwrite(rbindlist(lapply(c(cov_binary,'stage4','grade4','surgery_primary','Radiation recode','Chemotherapy recode (yes, no/unk)','Sequence number'),function(v)x[,.(n=.N),by=.(site=as.character(site_group),level=as.character(get(v)))][,variable:=v])), file.path(endpoint_out,'baseline_counts.csv'))
fwrite(x[,.(N=.N,median=median(age_years),q1=quantile(age_years,.25),q3=quantile(age_years,.75)),by=site_group],file.path(endpoint_out,'age_summary.csv'))
fwrite(x[,.(n=.N),by=.(site=as.character(site_group),level=as.character(cut(age_years,c(15,40,50,65,75,85),right=FALSE,labels=c('15–39','40–49','50–64','65–74','75–84'))))],file.path(endpoint_out,'age_counts.csv'))
