library(data.table)
args <- commandArgs(TRUE)
files <- file.path(args[1], sprintf("b%04d.csv", 1:2000))
stopifnot(all(file.exists(files)))
boot <- rbindlist(lapply(files, fread))
keys <- c("cohort","theta","contrast")
stopifnot(nrow(boot)==24000L, !anyDuplicated(boot[, c("b",keys), with=FALSE]), all(is.finite(boot$estimate)))
ci <- boot[,.(lower95=quantile(estimate,.025),upper95=quantile(estimate,.975),B=.N),by=keys]
stopifnot(all(ci$B==2000L))
point <- fread(file.path(args[1],"b0000.csv"))
alt <- merge(point[theta!=.1,.(cohort,theta,contrast,estimate)],ci,by=keys)
main <- fread("results/endpoint_sensitivity/ridge_contrasts_with_CI.csv")[sequence_group=="All"]
first <- fread(args[2])
first[,B:=replicates]
base <- rbindlist(list(A34=main,A40=first),idcol="cohort",fill=TRUE)[definition=="Common_biliary_minus_SEER" & state=="cancer" & startsWith(contrast,"RD_"),.(cohort,theta=.1,contrast,estimate,lower95,upper95,B)]
check <- merge(base,point[theta==.1],by=keys)
stopifnot(nrow(check)==6L,max(abs(check$estimate.x-check$estimate.y))<1e-10)
ans <- rbindlist(list(base,alt),use.names=TRUE)
stopifnot(nrow(ans)==18L)
for (col in c("estimate","lower95","upper95")) set(ans,j=col,value=100*ans[[col]])
setorder(ans,cohort,theta,contrast)
fwrite(ans,args[3])
