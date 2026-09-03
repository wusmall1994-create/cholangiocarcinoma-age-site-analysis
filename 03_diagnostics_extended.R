options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(cmprsk)
  library(splines)
})

source("00_config.R")

models <- readRDS(file.path(models_dir, "fitted_models.rds"))
age_knots <- models$age_knots
age_bounds <- models$age_bounds
spline_term <- "splines::ns(age_years, knots = age_knots, Boundary.knots = age_bounds)"

preferred_levels <- list(
  sex_factor = "Female",
  race_factor = "Non-Hispanic White",
  marital4 = "Married",
  income3 = "70,000-99,999",
  rural3 = "Metropolitan",
  era4 = "2010-2015",
  stage4 = "Localized",
  grade4 = "II",
  surgery_primary = "No",
  radiation_any = "No/Unknown",
  chemo_binary = "No/Unknown",
  followup_period = "0-12 months"
)

profile_from_fit <- function(fit, ages, periods = NULL) {
  if (is.null(periods)) {
    nd <- data.frame(
      age_years = rep(ages, each = 2L),
      site_group = factor(rep(c("iCCA", "eCCA"), times = length(ages)), levels = fit$xlevels$site_group)
    )
  } else {
    grid <- CJ(age_years = ages, followup_period = periods, site = c("iCCA", "eCCA"), sorted = FALSE)
    grid[, site_order := match(site, c("iCCA", "eCCA"))]
    setorder(grid, followup_period, age_years, site_order)
    nd <- data.frame(
      age_years = grid$age_years,
      site_group = factor(grid$site, levels = fit$xlevels$site_group),
      followup_period = factor(grid$followup_period, levels = fit$xlevels$followup_period)
    )
  }
  extra_levels <- setdiff(names(fit$xlevels), names(nd))
  extra_levels <- extra_levels[!grepl("^strata\\(", extra_levels)]
  for (nm in extra_levels) {
    lv <- fit$xlevels[[nm]]
    preferred <- preferred_levels[[nm]]
    value <- if (!is.null(preferred) && preferred %in% lv) preferred else lv[1L]
    nd[[nm]] <- factor(value, levels = lv)
  }
  nd
}

cox_curve <- function(fit, ages, outcome, model) {
  nd <- profile_from_fit(fit, ages)
  mm <- model.matrix(delete.response(terms(fit)), nd, contrasts.arg = fit$contrasts)
  mm <- mm[, names(coef(fit)), drop = FALSE]
  rbindlist(lapply(seq_along(ages), function(i) {
    rows <- (2L * i - 1L):(2L * i)
    z <- mm[rows[2L], ] - mm[rows[1L], ]
    est <- sum(z * coef(fit))
    se <- sqrt(drop(z %*% vcov(fit) %*% z))
    data.table(
      outcome = outcome, model = model, age = ages[i],
      estimate = exp(est), lower95 = exp(est - 1.96 * se), upper95 = exp(est + 1.96 * se)
    )
  }))
}

fg_curve <- function(fg, ages) {
  nd <- data.frame(
    age_years = rep(ages, each = 2L),
    site_group = factor(rep(c("iCCA", "eCCA"), times = length(ages)), levels = c("iCCA", "eCCA"))
  )
  mm <- model.matrix(
    as.formula(paste0("~ ", spline_term, " * site_group")),
    nd
  )
  mm <- mm[, fg$design_columns, drop = FALSE]
  b <- fg$fit$coef
  v <- fg$fit$var
  rbindlist(lapply(seq_along(ages), function(i) {
    rows <- (2L * i - 1L):(2L * i)
    z <- mm[rows[2L], ] - mm[rows[1L], ]
    est <- sum(z * b)
    se <- sqrt(drop(z %*% v %*% z))
    data.table(
      outcome = "Cancer death subdistribution", model = "Unadjusted", age = ages[i],
      estimate = exp(est), lower95 = exp(est - 1.96 * se), upper95 = exp(est + 1.96 * se)
    )
  }))
}

cif_curve <- function(fg, ages, horizon = 60) {
  nd <- data.frame(
    age_years = rep(ages, each = 2L),
    site_group = factor(rep(c("iCCA", "eCCA"), times = length(ages)), levels = c("iCCA", "eCCA"))
  )
  mm <- model.matrix(as.formula(paste0("~ ", spline_term, " * site_group")), nd)
  mm <- mm[, fg$design_columns, drop = FALSE]
  h0 <- sum(fg$fit$bfitj[fg$fit$uftime <= horizon])
  risk <- 1 - exp(-h0 * exp(drop(mm %*% fg$fit$coef)))
  data.table(
    age = rep(ages, each = 2L),
    site_group = rep(c("iCCA", "eCCA"), times = length(ages)),
    cif5 = risk
  )
}

age_grid <- seq(15, 85, by = 1)
relative_curves <- rbindlist(list(
  cox_curve(models$cox_cs_m1$fit, age_grid, "Cancer-specific hazard", "Model 1"),
  cox_curve(models$cox_cs_m2$fit, age_grid, "Cancer-specific hazard", "Model 2"),
  cox_curve(models$cox_os_m2$fit, age_grid, "Overall mortality", "Model 2"),
  fg_curve(models$fg_unadjusted, age_grid)
))
cif_continuous <- cif_curve(models$fg_absolute, age_grid)

prepare_piecewise <- function(x) {
  x <- copy(as.data.table(x))
  x <- x[!is.na(competing_event)]
  x[, cs_event := as.integer(competing_event == 1L)]
  x[, patient_cluster := factor(`Patient ID`)]
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
    is.na(income_first), "Unknown",
    income_first < 70000, "<70,000",
    income_first < 100000, "70,000-99,999",
    default = ">=100,000"
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
  x[, stage4 := factor(stage4, levels = c("Localized", "Regional", "Distant", "Unknown"))]
  x[, grade4 := factor(grade4, levels = c("I", "II", "III", "IV", "Unknown"))]
  x[, surgery_primary := factor(surgery_primary, levels = c("No", "Yes", "Unknown"))]
  radiation_label <- x[["Radiation recode"]]
  x[, radiation_any := factor(fifelse(
    radiation_label %chin% c("None/Unknown", "Refused (1988+)", "Recommended, unknown if administered"),
    "No/Unknown", "Yes"
  ), levels = c("No/Unknown", "Yes"))]
  x[, chemo_binary := factor(`Chemotherapy recode (yes, no/unk)`, levels = c("No/Unknown", "Yes"))]
  x
}

d0 <- readRDS(file.path(derived_dir, "cca_strict_hazard_2004_2023.rds"))
dp <- prepare_piecewise(d0)
dp <- survSplit(
  Surv(survival_time_months, cs_event) ~ .,
  data = as.data.frame(dp), cut = c(12, 36),
  start = "tstart", end = "tstop", event = "cs_event", episode = "period_code"
)
dp$followup_period <- factor(
  dp$period_code,
  levels = c(1, 2, 3),
  labels = c("0-12 months", "12-36 months", ">36 months")
)

covars <- paste(
  c("sex_factor", "race_factor", "marital4", "income3", "rural3", "era4",
    "stage4", "grade4", "surgery_primary", "radiation_any", "chemo_binary"),
  collapse = " + "
)
full_formula <- as.formula(paste0(
  "Surv(tstart, tstop, cs_event) ~ (", spline_term,
  " * site_group) + (", spline_term, ":followup_period) + ",
  "site_group:followup_period + (", spline_term,
  ":site_group:followup_period) + ", covars, " + strata(followup_period)"
))
reduced_formula <- as.formula(paste0(
  "Surv(tstart, tstop, cs_event) ~ (", spline_term, " * site_group) + (",
  spline_term, ":followup_period) + site_group:followup_period + ",
  covars, " + strata(followup_period)"
))

tv_full_lrt <- coxph(full_formula, data = dp, ties = "efron", x = TRUE, model = TRUE)
tv_reduced <- coxph(reduced_formula, data = dp, ties = "efron", x = TRUE, model = TRUE)
tv_full <- coxph(
  full_formula, data = dp, ties = "efron", x = TRUE, model = TRUE,
  cluster = dp$patient_cluster, robust = TRUE
)
tv_chisq <- 2 * (tv_full_lrt$loglik[2] - tv_reduced$loglik[2])
tv_df <- length(coef(tv_full_lrt)) - length(coef(tv_reduced))
tv_p <- pchisq(tv_chisq, tv_df, lower.tail = FALSE)

piecewise_contrasts <- function(fit, ages) {
  periods <- fit$xlevels$followup_period
  nd <- profile_from_fit(fit, ages, periods)
  mm <- model.matrix(fit, data = nd)
  keep <- !is.na(coef(fit))
  beta <- coef(fit)[keep]
  variance <- vcov(fit)[keep, keep, drop = FALSE]
  mm <- mm[, names(beta), drop = FALSE]
  rbindlist(lapply(seq_along(periods), function(j) {
    rbindlist(lapply(seq_along(ages), function(i) {
      base <- (j - 1L) * length(ages) * 2L + (i - 1L) * 2L
      rows <- base + 1:2
      z <- mm[rows[2L], ] - mm[rows[1L], ]
      est <- sum(z * beta)
      se <- sqrt(drop(z %*% variance %*% z))
      data.table(
        followup_period = periods[j], age = ages[i],
        estimate = exp(est), lower95 = exp(est - 1.96 * se), upper95 = exp(est + 1.96 * se)
      )
    }))
  }))
}
piecewise_results <- piecewise_contrasts(tv_full, c(20, 30, 40, 50, 65, 75))
piecewise_test <- data.table(
  test = "Age-by-site interaction varies across follow-up periods",
  chisq = tv_chisq, df = tv_df, p = tv_p
)

dc <- prepare_piecewise(readRDS(file.path(derived_dir, "cca_strict_cif5_2004_2018.rds")))
dc[, cif_group := paste(age_group, site_group, sep = "|")]
ci <- cuminc(
  ftime = dc$survival_time_months,
  fstatus = dc$competing_event,
  group = dc$cif_group,
  cencode = 0L
)
tp <- timepoints(ci, times = 60)
event_rows <- grepl(" 1$", rownames(tp$est))
aj <- data.table(
  group = sub(" 1$", "", rownames(tp$est)[event_rows]),
  estimate = as.numeric(tp$est[event_rows, 1L]),
  variance = as.numeric(tp$var[event_rows, 1L])
)
aj[, c("age_group", "site_group") := tstrsplit(group, "\\|", fixed = FALSE)]
aj[, `:=`(
  lower95 = pmax(0, estimate - 1.96 * sqrt(variance)),
  upper95 = pmin(1, estimate + 1.96 * sqrt(variance))
)]
aj_wide <- dcast(aj, age_group ~ site_group, value.var = c("estimate", "variance"))
aj_wide[, `:=`(
  risk_difference_eCCA_minus_iCCA = estimate_eCCA - estimate_iCCA,
  rd_se = sqrt(variance_eCCA + variance_iCCA)
)]
aj_wide[, `:=`(
  rd_lower95 = risk_difference_eCCA_minus_iCCA - 1.96 * rd_se,
  rd_upper95 = risk_difference_eCCA_minus_iCCA + 1.96 * rd_se
)]

fwrite(relative_curves, file.path(source_dir, "figure2_relative_effect_curves.csv"))
fwrite(cif_continuous, file.path(source_dir, "figure3_cif_curves.csv"))
fwrite(piecewise_results, file.path(tables_dir, "piecewise_followup_effects.csv"))
fwrite(piecewise_test, file.path(tables_dir, "piecewise_interaction_test.csv"))
fwrite(aj, file.path(tables_dir, "five_year_cif_age_groups.csv"))
fwrite(aj_wide, file.path(tables_dir, "five_year_risk_differences_age_groups.csv"))
fwrite(piecewise_results, file.path(source_dir, "figure4_piecewise_effects.csv"))
saveRDS(list(full = tv_full, reduced = tv_reduced), file.path(models_dir, "piecewise_models.rds"), compress = "xz")

cat("Piecewise three-way interaction p:", format.pval(tv_p, digits = 4), "\n")
print(piecewise_results)
print(aj_wide)
