args<-commandArgs(TRUE)
worker<-as.integer(args[1]); workers<-as.integer(args[2])
out <- args[3]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
Sys.setenv(CCA_RESULTS_DIR=out)
src <- readLines("scripts/run_unified_ridge.R")
src<-src[seq_len(grep("^mode <-",src)[1]-1L)]
src<-sub("theta=.1,scale","theta=current_penalty,scale",src,fixed=TRUE)
eval(parse(text=src))
defs<-c("SEER_original","Common_biliary")
ids<-split(seq_len(nrow(full)),full$patient_id)
for(b in c(if(worker==1) 0L, seq.int(worker,2000L,by=workers))){
 file<-file.path(out,sprintf("b%04d.csv",b))
 if(file.exists(file)) next
 z<-full
 if(b>0){set.seed(202610192L+b); z<-full[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]}
 ans<-list()
 for(current_penalty in if(b==0)c(.01,.1,1) else c(.01,1)){
 for(cohort in c("A34","A40")){
 x<-z[known==TRUE]
 if(cohort=="A40")x<-x[`Sequence number` %chin% c("One primary only","1st of 2 or more primaries")]
 warning_log<-list()
 states<-rbindlist(lapply(defs,function(label)pred(ridge_pair(endpoint(x,label),label),x,label)))
 result<-contrasts_cod(states)[state=="cancer" & definition=="Common_biliary_minus_SEER" & startsWith(contrast,"RD_")]
 stopifnot(nrow(result)==3, all(is.finite(result$estimate)))
 ans[[length(ans)+1]]<-result[,.(b,cohort,theta=current_penalty,contrast,estimate,warnings=length(warning_log))]
 }}
 fwrite(rbindlist(ans),paste0(file,".tmp"));file.rename(paste0(file,".tmp"),file)
}
