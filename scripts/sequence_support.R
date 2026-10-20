s<-readLines(file.path(Sys.getenv('CCA_ANALYSIS_CODE_DIR','.'),'run_unified_ridge.R'));eval(parse(text=s[seq_len(grep('^mode <-',s)[1]-1L)]))
dat[,sequence_detail:=ifelse(`Sequence number`=='One primary only','Single_primary',ifelse(`Sequence number`=='1st of 2 or more primaries','First_in_multiple','Later_in_multiple'))]
x<-rbindlist(lapply(defs,function(de){z<-endpoint(dat,de);z[,.(N=.N,patients=uniqueN(patient_id),endpoint_deaths=sum(cs_event),complement_deaths=sum(other_event),age_min=min(age_years),age_max=max(age_years),age40_59=sum(age_years>=40 & age_years<60),age65_84=sum(age_years>=65)),by=.(sequence_detail,site_group)][,definition:=de]}))
fwrite(x,file.path(endpoint_out,'sequence_support.csv'));print(x)
