# Cohort-specific reverse Kaplan-Meier potential follow-up.
# Requires the authorized derived data produced by 01_import_qc.R.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table); library(survival)})
derived_dir <- Sys.getenv('CCA_DERIVED_DIR', unset='')
results_dir <- Sys.getenv('CCA_RESULTS_DIR', unset='')
if (!nzchar(derived_dir) || !nzchar(results_dir)) {
  stop('Set CCA_DERIVED_DIR and CCA_RESULTS_DIR to local absolute paths.')
}
out <- file.path(results_dir, 'followup')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
h <- as.data.table(readRDS(file.path(derived_dir,'cca_strict_hazard_2004_2023.rds')))
cif <- as.data.table(readRDS(file.path(derived_dir,'cca_strict_cif5_2004_2018.rds')))
known <- function(z) z[!is.na(competing_event)]
exact <- function(z) known(z)[age_years <= 84]
eod <- exact(cif)[diagnosis_year >= 2010 & `EOD Schema ID Recode (2010+)` %chin%
                   c('Bile Ducts Intrahepat','Bile Ducts Perihilar','Bile Duct Distal')]
cohorts <- list(descriptive_2004_2023=h,
                cause_known_all_age_2004_2023=known(h),
                primary_hazard_2004_2023=exact(h),
                primary_risk_2004_2018=exact(cif),
                harmonized_EOD_2010_2018=eod)
expected <- c(descriptive_2004_2023=28654L, cause_known_all_age_2004_2023=28323L,
              primary_hazard_2004_2023=25400L, primary_risk_2004_2018=15578L,
              harmonized_EOD_2010_2018=11361L)
stopifnot(identical(vapply(cohorts,nrow,integer(1)),expected))
summarize_followup <- function(z,label) {
  time <- z$survival_time_months
  dead <- z$os_event
  stopifnot(all(is.finite(time)), all(time>0), all(dead %in% 0:1),
            all(dead==as.integer(z[['Vital status recode (study cutoff used)']]=='Dead')))
  # Death censors the follow-up distribution; alive/last contact is its event.
  fit <- survfit(Surv(time,1L-dead)~1)
  quartiles <- as.numeric(quantile(fit,probs=c(.25,.5,.75),conf.int=FALSE))
  # Independent product-limit reconstruction, including tied observations.
  tab <- data.table(time=time,event=1L-dead)[,.(n=.N,d=sum(event)),by=time][order(time)]
  tab[,nrisk:=rev(cumsum(rev(n)))]; tab[,s:=cumprod(1-d/nrisk)]
  check <- vapply(c(.25,.5,.75), function(p) {
    index <- which(tab$s <= 1-p)
    if(length(index)) tab$time[index[1]] else NA_real_
  }, numeric(1))
  stopifnot(isTRUE(all.equal(quartiles,check,check.attributes=FALSE)))
  observed <- as.numeric(quantile(time,c(.25,.5,.75),type=2))
  data.table(cohort=label,N=nrow(z),patients=uniqueN(z[['Patient ID']]),
             all_deaths=sum(dead),alive=sum(1L-dead),
             cancer_deaths=sum(z$competing_event==1,na.rm=TRUE),
             other_deaths=sum(z$competing_event==2,na.rm=TRUE),
             reverse_KM_q25=quartiles[1],reverse_KM_median=quartiles[2],
             reverse_KM_q75=quartiles[3],observed_q25=observed[1],
             observed_median=observed[2],observed_q75=observed[3],max_observed=max(time))
}
summary <- rbindlist(Map(summarize_followup,cohorts,names(cohorts)))
fwrite(summary,file.path(out,'followup_by_cohort.csv'))
writeLines(c(
 'Unit: tumor record, matching analysis cohorts; no patient-first-record deduplication.',
 'Full recorded survival months; zero recorded months assigned 0.5 by import.',
 'Reverse KM: all-cause deaths censored; alive/last-contact observations are events.',
 'No truncation at 60 months: follow-up availability differs from the prediction horizon.',
 'Independent product-limit calculation reproduced every survival::survfit quartile.',
 'All time summaries are in months. q25 and q75 bound the interquartile range.',
 'observed_* are empirical recorded-time quantiles (type 2), not reverse-KM estimates.',
 'N and patients are aggregate counts; cancer/other deaths omit unknown cause.',
 'No patient-level data or fitted model exported.'),file.path(out,'METHOD_AND_VERIFICATION.txt'))
capture.output(sessionInfo(),file=file.path(out,'R_SESSION_INFO.txt'))
print(summary)
