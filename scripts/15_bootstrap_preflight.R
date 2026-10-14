s<-readLines('11_horizon_anatomy_absolute_risks.R')
eval(parse(text=s[seq_len(grep('^binary_point_long <-',s)[1]-1L)]))
for(v in cov_schema){cat('\n',v,'\n');print(schema_data[,.(N=.N,cancer=sum(cs_event),other=sum(other_event)),by=v])}
