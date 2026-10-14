options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(splines)
})

source("00_config.R")

prepare_covariates <- function(x) {
  x <- copy(as.data.table(x))
  x[, patient_id := as.character(`Patient ID`)]
  x[, site_group := factor(site_group, levels = c("iCCA", "eCCA"))]
  x[, sex_factor := factor(Sex, levels = c("Female", "Male"))]
  x[, race_factor := factor(
    `Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)`,
    levels = c(
      "Non-Hispanic White", "Non-Hispanic Black",
      "Non-Hispanic Asian or Pacific Islander",
      "Non-Hispanic American Indian/Alaska Native",
      "Hispanic (All Races)", "Non-Hispanic Unknown Race"
    )
  )]
  x[, marital4 := factor(fcase(
    `Marital status at diagnosis` == "Married (including common law)", "Married",
    `Marital status at diagnosis` == "Widowed", "Widowed",
    `Marital status at diagnosis` == "Unknown", "Unknown",
    default = "Unmarried"
  ), levels = c("Married", "Unmarried", "Widowed", "Unknown"))]
  income_label <- x[["Median household income inflation adj to 2024"]]
  income_first <- suppressWarnings(as.numeric(gsub(",", "", sub(".*?([0-9][0-9,]*).*", "\\1", income_label))))
  income_first[grepl("Unknown", income_label, ignore.case = TRUE)] <- NA_real_
  x[, income3 := factor(fcase(
    is.na(income_first), "Unknown", income_first < 70000, "<70,000",
    income_first < 100000, "70,000-99,999", default = ">=100,000"
  ), levels = c("<70,000", "70,000-99,999", ">=100,000", "Unknown"))]
  rural_label <- x[["Rural-Urban Continuum Code"]]
  x[, rural3 := factor(fcase(
    grepl("^Counties in metropolitan", rural_label), "Metropolitan",
    grepl("^Nonmetropolitan", rural_label), "Nonmetropolitan",
    default = "Unknown"
  ), levels = c("Metropolitan", "Nonmetropolitan", "Unknown"))]
  x[, era4 := factor(fcase(
    diagnosis_year <= 2009L, "2004-2009",
    diagnosis_year <= 2015L, "2010-2015",
    diagnosis_year <= 2019L, "2016-2019",
    default = "2020-2023"
  ), levels = c("2004-2009", "2010-2015", "2016-2019", "2020-2023"))]
  x[, cause_known := !is.na(competing_event)]
  x[, cs_event := as.integer(competing_event == 1L)]
  x[, other_event := as.integer(competing_event == 2L)]
  x
}

d <- prepare_covariates(readRDS(file.path(derived_dir, "cca_strict_cif5_2004_2018.rds")))
d <- droplevels(d[age_topcoded_85plus == FALSE & cause_known == TRUE])
# A single record had unknown neighborhood income and no other-cause event,
# producing separation in the interval-specific competing-event model.
# Exclude that record for this model-check sensitivity analysis.
d_nonph <- droplevels(d[income3 != "Unknown"])

age_knots <- c(52, 68, 80)
age_bounds <- c(15, 84)
ns_term <- paste0(
  "splines::ns(age_years, knots=c(", paste(age_knots, collapse = ","),
  "), Boundary.knots=c(", paste(age_bounds, collapse = ","), "))"
)
cov_model1 <- c("sex_factor", "race_factor", "marital4", "income3", "rural3", "era4")
period_levels <- c("0-12", "12-36", "36-60")

make_split_data <- function(data, event_var) {
  z <- copy(data)
  z[, time60 := pmin(survival_time_months, 60)]
  z[, event60 := as.integer(get(event_var) == 1L & survival_time_months <= 60)]
  z <- survSplit(
    Surv(time60, event60) ~ ., data = as.data.frame(z),
    cut = c(12, 36), start = "tstart", end = "tstop",
    event = "event60", episode = "period_code"
  )
  z$followup_period <- factor(z$period_code, levels = 1:3, labels = period_levels)
  z
}

fit_piecewise <- function(data, event_var) {
  z <- make_split_data(data, event_var)
  rhs <- paste0(
    ns_term, " * site_group * followup_period + ",
    paste(cov_model1, collapse = " + "), " + strata(followup_period)"
  )
  fit <- coxph(
    as.formula(paste0("Surv(tstart, tstop, event60) ~ ", rhs)),
    data = z, ties = "efron", x = TRUE
  )
  list(fit = fit, basehaz = suppressWarnings(basehaz(fit, centered = FALSE)))
}

baseline_increments <- function(base_hazard) {
  z <- as.data.table(base_hazard)
  z[, period := sub("^followup_period=", "", as.character(strata))]
  setorder(z, period, time)
  z[, increment := hazard - shift(hazard, fill = 0), by = period]
  z[, .(period, time, increment)]
}

standardized_piecewise_cif <- function(fits, data, age, site) {
  nd <- copy(data)
  nd[, age_years := age]
  nd[, site_group := factor(site, levels = c("iCCA", "eCCA"))]
  survival_probability <- rep(1, nrow(nd))
  cif <- rep(0, nrow(nd))
  hc <- baseline_increments(fits$cancer$basehaz)
  ho <- baseline_increments(fits$other$basehaz)

  for (period_value in period_levels) {
    nd[, followup_period := factor(period_value, levels = period_levels)]
    rc <- exp(predict(fits$cancer$fit, newdata = as.data.frame(nd), type = "lp", reference = "zero"))
    ro <- exp(predict(fits$other$fit, newdata = as.data.frame(nd), type = "lp", reference = "zero"))
    increments <- merge(
      hc[period == period_value, .(time, dhc = increment)],
      ho[period == period_value, .(time, dho = increment)],
      by = "time", all = TRUE
    )
    increments[is.na(dhc), dhc := 0]
    increments[is.na(dho), dho := 0]
    setorder(increments, time)
    for (k in seq_len(nrow(increments))) {
      total_increment <- increments$dhc[k] * rc + increments$dho[k] * ro
      cancer_fraction <- ifelse(
        total_increment > 0,
        increments$dhc[k] * rc / total_increment,
        0
      )
      cif <- cif + survival_probability * (1 - exp(-total_increment)) * cancer_fraction
      survival_probability <- survival_probability * exp(-total_increment)
    }
  }
  mean(cif)
}

fits <- list(
  cancer = fit_piecewise(d_nonph, "cs_event"),
  other = fit_piecewise(d_nonph, "other_event")
)

ages <- c(20, 30, 40, 50, 65, 75)
piecewise_long <- rbindlist(lapply(ages, function(age) {
  rbindlist(lapply(c("iCCA", "eCCA"), function(site) {
    data.table(
      age = age,
      site_group = site,
      five_year_cif = standardized_piecewise_cif(fits, d_nonph, age, site)
    )
  }))
}))
piecewise_rd <- dcast(piecewise_long, age ~ site_group, value.var = "five_year_cif")
piecewise_rd[, piecewise_risk_difference := eCCA - iCCA]

primary <- fread(file.path(tables_dir, "diagnostic_adjusted_cif5_risk_differences.csv"))[
  target_type == "Fixed age"
]
primary[, age := as.integer(target)]
primary <- primary[age %in% ages]
comparison <- merge(
  piecewise_rd,
  primary[, .(
    age,
    proportional_hazards_risk_difference = risk_difference_eCCA_minus_iCCA,
    proportional_hazards_lower95 = rd_lower95,
    proportional_hazards_upper95 = rd_upper95
  )],
  by = "age", all.x = TRUE
)
comparison[, absolute_model_difference := piecewise_risk_difference - proportional_hazards_risk_difference]

fwrite(piecewise_long, file.path(source_dir, "nonph_piecewise_standardized_cif.csv"))
fwrite(comparison, file.path(tables_dir, "nonph_piecewise_vs_ph_risk_differences.csv"))
saveRDS(fits, file.path(models_dir, "nonph_piecewise_cause_specific_models.rds"), compress = "xz")

cat("Analysis cohort:", nrow(d_nonph), "records\n")
cat("Cancer deaths through five years:", sum(d_nonph$cs_event), "\n")
cat("Other-cause deaths through five years:", sum(d_nonph$other_event), "\n")
print(comparison)
cat("Maximum absolute model difference:", max(abs(comparison$absolute_model_difference)), "\n")
