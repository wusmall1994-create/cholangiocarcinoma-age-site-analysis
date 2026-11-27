args <- commandArgs(trailingOnly=TRUE)
worker <- as.integer(args[1])
workers <- as.integer(args[2])
out <- args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
Sys.setenv(CCA_RESULTS_DIR=out)
src <- 'scripts/run_unified_ridge.R'
lines <- readLines(src)
eval(parse(text=lines[seq_len(grep('^mode <-',lines)[1]-1L)]))
eligible <- c('One primary only', '1st of 2 or more primaries')
defs <- c('SEER_original','Common_biliary')
fit_first <- function(z) {
  x <- z[known & `Sequence number` %chin% eligible]
  states <- rbindlist(lapply(defs,function(label) {
    pred(ridge_pair(endpoint(x,label),label),x,label)
  }))
  contrasts_cod(states)[state=='cancer' & contrast %chin%
    c('RD_age50','RD_age75','RD_change_75_minus_50')]
}
if (worker==1L) {
  warning_log <- list()
  point <- fit_first(full)
  stopifnot(length(warning_log)==0L)
  fwrite(point,file.path(out,'point.csv'))
  keep <- full[known & `Sequence number` %chin% eligible]
  stopifnot(nrow(keep)==12734L,uniqueN(keep$patient_id)==nrow(keep))
  counts <- rbindlist(lapply(defs,function(label) endpoint(keep,label)[,
    .(records=.N,patients=uniqueN(patient_id),endpoint_deaths=sum(cs_event),
      complement_deaths=sum(other_event)),by=.(site_group,`Sequence number`)][,definition:=label]))
  fwrite(counts,file.path(out,'counts.csv'))
}
ids <- split(seq_len(nrow(full)),full$patient_id)
for (b in seq.int(worker,2000L,by=workers)) {
  file <- file.path(out,sprintf('rep_%04d.rds',b))
  if(file.exists(file)) next
  set.seed(202610192L+b)
  z <- full[unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)]
  warning_log <- list()
  value <- fit_first(z)
  stopifnot(nrow(value)==9L,all(is.finite(value$estimate)))
  saveRDS(list(value=value,warnings=warning_log),file)
}
