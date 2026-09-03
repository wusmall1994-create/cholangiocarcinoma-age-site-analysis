options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(cmprsk)
  library(splines)
})

source("00_config.R")

d <- as.data.table(readRDS(file.path(derived_dir, "cca_strict_hazard_2004_2023.rds")))
d_cif <- as.data.table(readRDS(file.path(derived_dir, "cca_strict_cif5_2004_2018.rds")))

prepare_covariates <- function(x) {
  x <- copy(x)
  x[, patient_cluster := factor(`Patient ID`)]
  x[, site_group := factor(site_group, levels = c("iCCA", "eCCA"))]
  x[, sex_factor := factor(Sex, levels = c("Female", "Male"))]
  x[, race_factor := factor(
    `Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)`,
    levels = c(
      "Non-Hispanic White",
      "Non-Hispanic Black",
      "Non-Hispanic Asian or Pacific Islander",
      "Non-Hispanic American Indian/Alaska Native",
      "Hispanic (All Races)",
      "Non-Hispanic Unknown Race"
    )
  )]
  x[, marital4 := factor(
    fcase(
      `Marital status at diagnosis` == "Married (including common law)", "Married",
      `Marital status at diagnosis` == "Widowed", "Widowed",
      `Marital status at diagnosis` == "Unknown", "Unknown",
      default = "Unmarried"
    ),
    levels = c("Married", "Unmarried", "Widowed", "Unknown")
  )]
  income_label <- x[["Median household income inflation adj to 2024"]]
  income_first <- suppressWarnings(as.numeric(gsub(",", "", sub(".*?([0-9][0-9,]*).*", "\\1", income_label))))
  income_first[grepl("Unknown", income_label, ignore.case = TRUE)] <- NA_real_
  x[, income3 := factor(
    fcase(
      is.na(income_first), "Unknown",
      income_first < 70000, "<70,000",
      income_first < 100000, "70,000-99,999",
      default = ">=100,000"
    ),
    levels = c("<70,000", "70,000-99,999", ">=100,000", "Unknown")
  )]
  rural_label <- x[["Rural-Urban Continuum Code"]]
  x[, rural3 := factor(
    fcase(
      grepl("^Counties in metropolitan", rural_label), "Metropolitan",
      grepl("^Nonmetropolitan", rural_label), "Nonmetropolitan",
      default = "Unknown"
    ),
    levels = c("Metropolitan", "Nonmetropolitan", "Unknown")
  )]
  x[, era4 := factor(
    fcase(
      diagnosis_year <= 2009L, "2004-2009",
      diagnosis_year <= 2015L, "2010-2015",
      diagnosis_year <= 2019L, "2016-2019",
      default = "2020-2023"
    ),
    levels = c("2004-2009", "2010-2015", "2016-2019", "2020-2023")
  )]
  x[, stage4 := factor(stage4, levels = c("Localized", "Regional", "Distant", "Unknown"))]
  x[, grade4 := factor(grade4, levels = c("I", "II", "III", "IV", "Unknown"))]
  x[, surgery_primary := factor(surgery_primary, levels = c("No", "Yes", "Unknown"))]
  radiation_label <- x[["Radiation recode"]]
  x[, radiation_any := factor(
    fifelse(
      radiation_label %chin% c(
        "None/Unknown", "Refused (1988+)",
        "Recommended, unknown if administered"
      ),
      "No/Unknown", "Yes"
    ),
    levels = c("No/Unknown", "Yes")
  )]
  x[, chemo_binary := factor(
    `Chemotherapy recode (yes, no/unk)`,
    levels = c("No/Unknown", "Yes")
  )]
  x[, positive_histology := `Diagnostic Confirmation` %chin% c(
    "Positive histology",
    "Positive microscopic confirm, method not specified"
  )]
  x[, first_primary := `Sequence number` == "One primary only" |
      startsWith(`Sequence number`, "1st of")]
  x[, cause_known := !is.na(competing_event)]
  x
}

d <- prepare_covariates(d)
d_cif <- prepare_covariates(d_cif)

age_knots <- as.numeric(quantile(
  d[age_topcoded_85plus == FALSE, age_years],
  probs = c(0.10, 0.50, 0.90), na.rm = TRUE, type = 2
))
age_bounds <- range(d$age_years, na.rm = TRUE)
spline_term <- "splines::ns(age_years, knots = age_knots, Boundary.knots = age_bounds)"

cov_model1 <- c("sex_factor", "race_factor", "marital4", "income3", "rural3", "era4")
cov_model2 <- c(cov_model1, "stage4", "grade4", "surgery_primary", "radiation_any", "chemo_binary")

rhs_text <- function(covars, interaction = TRUE) {
  age_site <- if (interaction) paste0(spline_term, " * site_group") else paste0(spline_term, " + site_group")
  paste(c(age_site, covars), collapse = " + ")
}

fit_cox_pair <- function(data, event, covars) {
  data <- droplevels(as.data.frame(data))
  f0 <- as.formula(paste0("Surv(survival_time_months, ", event, ") ~ ", rhs_text(covars, FALSE)))
  f1 <- as.formula(paste0("Surv(survival_time_months, ", event, ") ~ ", rhs_text(covars, TRUE)))
  no_int <- coxph(f0, data = data, ties = "efron", x = TRUE, model = TRUE)
  int_lrt <- coxph(f1, data = data, ties = "efron", x = TRUE, model = TRUE)
  int_robust <- coxph(
    f1, data = data, ties = "efron", x = TRUE, model = TRUE,
    cluster = data$patient_cluster, robust = TRUE
  )
  lr <- 2 * (int_lrt$loglik[2] - no_int$loglik[2])
  df_lr <- length(coef(int_lrt)) - length(coef(no_int))
  list(
    no_int = no_int,
    int_lrt = int_lrt,
    fit = int_robust,
    interaction_chisq = unname(lr),
    interaction_df = df_lr,
    interaction_p = pchisq(lr, df_lr, lower.tail = FALSE)
  )
}

profile_data <- function(data, ages) {
  valid_level <- function(x, preferred) {
    lv <- levels(x)
    if (preferred %in% lv) preferred else lv[1L]
  }
  data.frame(
    age_years = rep(ages, each = 2L),
    site_group = factor(rep(c("iCCA", "eCCA"), times = length(ages)), levels = levels(data$site_group)),
    sex_factor = factor(valid_level(data$sex_factor, "Female"), levels = levels(data$sex_factor)),
    race_factor = factor(valid_level(data$race_factor, "Non-Hispanic White"), levels = levels(data$race_factor)),
    marital4 = factor(valid_level(data$marital4, "Married"), levels = levels(data$marital4)),
    income3 = factor(valid_level(data$income3, "70,000-99,999"), levels = levels(data$income3)),
    rural3 = factor(valid_level(data$rural3, "Metropolitan"), levels = levels(data$rural3)),
    era4 = factor(valid_level(data$era4, "2010-2015"), levels = levels(data$era4)),
    stage4 = factor(valid_level(data$stage4, "Localized"), levels = levels(data$stage4)),
    grade4 = factor(valid_level(data$grade4, "II"), levels = levels(data$grade4)),
    surgery_primary = factor(valid_level(data$surgery_primary, "No"), levels = levels(data$surgery_primary)),
    radiation_any = factor(valid_level(data$radiation_any, "No/Unknown"), levels = levels(data$radiation_any)),
    chemo_binary = factor(valid_level(data$chemo_binary, "No/Unknown"), levels = levels(data$chemo_binary))
  )
}

cox_site_contrasts <- function(fit, data, ages, model_name, outcome_name) {
  nd <- profile_data(data, ages)
  mm <- model.matrix(delete.response(terms(fit)), nd, contrasts.arg = fit$contrasts)
  mm <- mm[, names(coef(fit)), drop = FALSE]
  out <- rbindlist(lapply(seq_along(ages), function(i) {
    rows <- (2L * i - 1L):(2L * i)
    z <- mm[rows[2L], ] - mm[rows[1L], ]
    est <- sum(z * coef(fit))
    se <- sqrt(drop(z %*% vcov(fit) %*% z))
    data.table(
      outcome = outcome_name,
      model = model_name,
      age = ages[i],
      contrast = "eCCA vs iCCA",
      estimate = exp(est),
      lower95 = exp(est - 1.96 * se),
      upper95 = exp(est + 1.96 * se)
    )
  }))
  out
}

fg_design <- function(data, covars, interaction = TRUE) {
  f <- as.formula(paste0("~ ", rhs_text(covars, interaction)))
  mm <- model.matrix(f, data = data)
  mm[, colnames(mm) != "(Intercept)", drop = FALSE]
}

fit_fg_pair <- function(data, covars) {
  data <- data[cause_known == TRUE]
  x0 <- fg_design(data, covars, FALSE)
  x1 <- fg_design(data, covars, TRUE)
  fit0 <- crr(
    ftime = data$survival_time_months,
    fstatus = data$competing_event,
    cov1 = x0,
    failcode = 1L, cencode = 0L, variance = FALSE
  )
  fit1 <- crr(
    ftime = data$survival_time_months,
    fstatus = data$competing_event,
    cov1 = x1,
    failcode = 1L, cencode = 0L
  )
  lr <- 2 * (fit1$loglik - fit0$loglik)
  df_lr <- length(fit1$coef) - length(fit0$coef)
  list(
    no_int = fit0,
    fit = fit1,
    design_columns = colnames(x1),
    interaction_chisq = unname(lr),
    interaction_df = df_lr,
    interaction_p = pchisq(lr, df_lr, lower.tail = FALSE)
  )
}

fg_site_contrasts <- function(fg, data, covars, ages, model_name) {
  nd <- profile_data(data, ages)
  mm <- model.matrix(as.formula(paste0("~ ", rhs_text(covars, TRUE))), nd)
  mm <- mm[, fg$design_columns, drop = FALSE]
  b <- fg$fit$coef
  v <- fg$fit$var
  rbindlist(lapply(seq_along(ages), function(i) {
    rows <- (2L * i - 1L):(2L * i)
    z <- mm[rows[2L], ] - mm[rows[1L], ]
    est <- sum(z * b)
    se <- sqrt(drop(z %*% v %*% z))
    data.table(
      outcome = "Cancer death subdistribution",
      model = model_name,
      age = ages[i],
      contrast = "eCCA vs iCCA",
      estimate = exp(est),
      lower95 = exp(est - 1.96 * se),
      upper95 = exp(est + 1.96 * se)
    )
  }))
}

analysis_ages <- c(20, 30, 40, 50, 65, 75)
d_known <- d[cause_known == TRUE]

cox_cs_m1 <- fit_cox_pair(d_known, "os_event * 0 + (competing_event == 1L)", cov_model1)
cox_cs_m2 <- fit_cox_pair(d_known, "os_event * 0 + (competing_event == 1L)", cov_model2)
cox_os_m1 <- fit_cox_pair(d, "os_event", cov_model1)
cox_os_m2 <- fit_cox_pair(d, "os_event", cov_model2)
fg_unadjusted <- fit_fg_pair(d, character(0L))

relative_results <- rbindlist(list(
  cox_site_contrasts(cox_cs_m1$fit, d_known, analysis_ages, "Model 1", "Cancer-specific hazard"),
  cox_site_contrasts(cox_cs_m2$fit, d_known, analysis_ages, "Model 2", "Cancer-specific hazard"),
  cox_site_contrasts(cox_os_m1$fit, d, analysis_ages, "Model 1", "Overall mortality"),
  cox_site_contrasts(cox_os_m2$fit, d, analysis_ages, "Model 2", "Overall mortality"),
  fg_site_contrasts(fg_unadjusted, d, character(0L), analysis_ages, "Unadjusted")
))

interaction_results <- rbindlist(list(
  data.table(outcome = "Cancer-specific hazard", model = "Model 1", chisq = cox_cs_m1$interaction_chisq, df = cox_cs_m1$interaction_df, p = cox_cs_m1$interaction_p),
  data.table(outcome = "Cancer-specific hazard", model = "Model 2", chisq = cox_cs_m2$interaction_chisq, df = cox_cs_m2$interaction_df, p = cox_cs_m2$interaction_p),
  data.table(outcome = "Overall mortality", model = "Model 1", chisq = cox_os_m1$interaction_chisq, df = cox_os_m1$interaction_df, p = cox_os_m1$interaction_p),
  data.table(outcome = "Overall mortality", model = "Model 2", chisq = cox_os_m2$interaction_chisq, df = cox_os_m2$interaction_df, p = cox_os_m2$interaction_p),
  data.table(outcome = "Cancer death subdistribution", model = "Unadjusted", chisq = fg_unadjusted$interaction_chisq, df = fg_unadjusted$interaction_df, p = fg_unadjusted$interaction_p)
))

fg_abs_data <- d_cif[cause_known == TRUE]
fg_abs <- fit_fg_pair(fg_abs_data, character(0L))
fg_absolute_risk <- function(fg, data, ages, horizon = 60) {
  nd <- profile_data(data, ages)
  mm <- model.matrix(as.formula(paste0("~ ", rhs_text(character(0L), TRUE))), nd)
  mm <- mm[, fg$design_columns, drop = FALSE]
  h0 <- sum(fg$fit$bfitj[fg$fit$uftime <= horizon])
  risk <- 1 - exp(-h0 * exp(drop(mm %*% fg$fit$coef)))
  ans <- data.table(
    age = rep(ages, each = 2L),
    site_group = rep(c("iCCA", "eCCA"), times = length(ages)),
    cif5 = risk
  )
  wide <- dcast(ans, age ~ site_group, value.var = "cif5")
  wide[, risk_difference_eCCA_minus_iCCA := eCCA - iCCA]
  wide
}
absolute_results <- fg_absolute_risk(fg_abs, fg_abs_data, analysis_ages)

zph <- cox.zph(cox_cs_m2$fit, transform = "km", terms = TRUE, singledf = TRUE)
zph_table <- data.table(
  term = rownames(zph$table),
  chisq = zph$table[, "chisq"],
  df = zph$table[, "df"],
  p = zph$table[, "p"]
)

cell_counts <- d[, .(
  n = .N,
  cancer_deaths = sum(competing_event == 1L, na.rm = TRUE),
  other_deaths = sum(competing_event == 2L, na.rm = TRUE),
  unknown_cause_deaths = sum(is.na(competing_event) & os_event == 1L),
  alive = sum(os_event == 0L)
), by = .(age_group, site_group)]

table1_rows <- list()
add_cont <- function(label, variable) {
  table1_rows[[length(table1_rows) + 1L]] <<- data.table(
    characteristic = label,
    level = "Median [IQR]",
    iCCA = sprintf("%.0f [%.0f-%.0f]", median(d[site_group == "iCCA"][[variable]], na.rm = TRUE), quantile(d[site_group == "iCCA"][[variable]], 0.25, na.rm = TRUE), quantile(d[site_group == "iCCA"][[variable]], 0.75, na.rm = TRUE)),
    eCCA = sprintf("%.0f [%.0f-%.0f]", median(d[site_group == "eCCA"][[variable]], na.rm = TRUE), quantile(d[site_group == "eCCA"][[variable]], 0.25, na.rm = TRUE), quantile(d[site_group == "eCCA"][[variable]], 0.75, na.rm = TRUE))
  )
}
add_cat <- function(label, variable) {
  levs <- levels(d[[variable]])
  if (is.null(levs)) levs <- sort(unique(as.character(d[[variable]])))
  den_i <- sum(d$site_group == "iCCA")
  den_e <- sum(d$site_group == "eCCA")
  for (lv in levs) {
    ni <- sum(d$site_group == "iCCA" & d[[variable]] == lv, na.rm = TRUE)
    ne <- sum(d$site_group == "eCCA" & d[[variable]] == lv, na.rm = TRUE)
    table1_rows[[length(table1_rows) + 1L]] <<- data.table(
      characteristic = label,
      level = lv,
      iCCA = sprintf("%d (%.1f%%)", ni, 100 * ni / den_i),
      eCCA = sprintf("%d (%.1f%%)", ne, 100 * ne / den_e)
    )
  }
}
table1_rows[[1L]] <- data.table(
  characteristic = "Cohort size", level = "N",
  iCCA = as.character(sum(d$site_group == "iCCA")),
  eCCA = as.character(sum(d$site_group == "eCCA"))
)
add_cont("Age at diagnosis, years", "age_years")
add_cat("Age group", "age_group")
add_cat("Sex", "sex_factor")
add_cat("Race/ethnicity", "race_factor")
add_cat("Marital status", "marital4")
add_cat("County household income, US$", "income3")
add_cat("County rurality", "rural3")
add_cat("Diagnosis era", "era4")
add_cat("Summary stage", "stage4")
add_cat("Grade", "grade4")
add_cat("Primary-site surgery", "surgery_primary")
add_cat("Radiotherapy", "radiation_any")
add_cat("Chemotherapy", "chemo_binary")
table1 <- rbindlist(table1_rows)

sensitivity_specs <- list(
  Main = rep(TRUE, nrow(d)),
  `Exclude C22.1 + 8162` = !d$discordant_c221_8162,
  `Positive histology only` = d$positive_histology,
  `First primary only` = d$first_primary,
  `Exclude age 85+` = !d$age_topcoded_85plus
)

fit_sensitivity <- function(label, keep) {
  ds <- droplevels(d[keep & cause_known == TRUE])
  usable_covars <- cov_model2[vapply(cov_model2, function(v) {
    z <- ds[[v]]
    length(unique(z[!is.na(z)])) > 1L
  }, logical(1))]
  pair <- fit_cox_pair(ds, "os_event * 0 + (competing_event == 1L)", usable_covars)
  hr65 <- cox_site_contrasts(pair$fit, ds, 65, "Model 2", "Cancer-specific hazard")
  data.table(
    analysis = label,
    n = nrow(ds),
    cancer_deaths = sum(ds$competing_event == 1L),
    interaction_p = pair$interaction_p,
    hr_eCCA_vs_iCCA_age65 = hr65$estimate,
    lower95 = hr65$lower95,
    upper95 = hr65$upper95
  )
}
sensitivity_results <- rbindlist(Map(fit_sensitivity, names(sensitivity_specs), sensitivity_specs))

era_sensitivity <- rbindlist(lapply(levels(d$era4), function(er) {
  ds <- droplevels(d[era4 == er & cause_known == TRUE])
  covs <- setdiff(cov_model2, "era4")
  covs <- covs[vapply(covs, function(v) {
    z <- ds[[v]]
    length(unique(z[!is.na(z)])) > 1L
  }, logical(1))]
  pair <- fit_cox_pair(ds, "os_event * 0 + (competing_event == 1L)", covs)
  hr65 <- cox_site_contrasts(pair$fit, ds, 65, "Model 2", "Cancer-specific hazard")
  data.table(
    era = er,
    n = nrow(ds),
    cancer_deaths = sum(ds$competing_event == 1L),
    interaction_p = pair$interaction_p,
    hr_eCCA_vs_iCCA_age65 = hr65$estimate,
    lower95 = hr65$lower95,
    upper95 = hr65$upper95
  )
}))

fwrite(table1, file.path(tables_dir, "table1_baseline.csv"))
fwrite(cell_counts, file.path(tables_dir, "age_site_event_cells.csv"))
fwrite(interaction_results, file.path(tables_dir, "interaction_tests.csv"))
fwrite(relative_results, file.path(tables_dir, "age_specific_relative_effects.csv"))
fwrite(absolute_results, file.path(tables_dir, "five_year_cif_age_specific.csv"))
fwrite(zph_table, file.path(tables_dir, "cox_ph_diagnostics.csv"))
fwrite(sensitivity_results, file.path(tables_dir, "sensitivity_analyses.csv"))
fwrite(era_sensitivity, file.path(tables_dir, "era_sensitivity.csv"))

fwrite(relative_results, file.path(source_dir, "figure2_relative_effects.csv"))
fwrite(absolute_results, file.path(source_dir, "figure3_five_year_cif.csv"))
fwrite(sensitivity_results, file.path(source_dir, "figure4_sensitivity.csv"))

saveRDS(
  list(
    age_knots = age_knots,
    age_bounds = age_bounds,
    cox_cs_m1 = cox_cs_m1,
    cox_cs_m2 = cox_cs_m2,
    cox_os_m1 = cox_os_m1,
    cox_os_m2 = cox_os_m2,
    fg_unadjusted = fg_unadjusted,
    fg_absolute = fg_abs
  ),
  file.path(models_dir, "fitted_models.rds"),
  compress = "xz"
)

analysis_summary <- c(
  paste0("Age spline internal knots: ", paste(age_knots, collapse = ", ")),
  paste0("Age spline boundary knots: ", paste(age_bounds, collapse = ", ")),
  paste0("Cause-specific Cox Model 2 interaction p: ", format.pval(cox_cs_m2$interaction_p, digits = 4)),
  paste0("Fine-Gray unadjusted interaction p: ", format.pval(fg_unadjusted$interaction_p, digits = 4)),
  paste0("Overall survival Cox Model 2 interaction p: ", format.pval(cox_os_m2$interaction_p, digits = 4)),
  paste0("Cox PH global p: ", format.pval(zph_table[term == "GLOBAL", p], digits = 4))
)
writeLines(analysis_summary, file.path(results_dir, "analysis_summary.txt"))
cat(paste(analysis_summary, collapse = "\n"), "\n")
