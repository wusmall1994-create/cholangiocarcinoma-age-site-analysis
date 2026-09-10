options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(splines)
})

source("00_config.R")

prepare_covariates_enhanced <- function(x) {
  x <- copy(as.data.table(x))
  x[, patient_id := as.character(`Patient ID`)]
  x[, patient_cluster := factor(patient_id)]
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
  x[, era5 := factor(
    fcase(
      diagnosis_year <= 2003L, "2000-2003",
      diagnosis_year <= 2009L, "2004-2009",
      diagnosis_year <= 2015L, "2010-2015",
      diagnosis_year <= 2019L, "2016-2019",
      default = "2020-2023"
    ),
    levels = c("2000-2003", "2004-2009", "2010-2015", "2016-2019", "2020-2023")
  )]
  x[, stage4 := factor(stage4, levels = c("Localized", "Regional", "Distant", "Unknown"))]
  x[, grade4 := factor(grade4, levels = c("I", "II", "III", "IV", "Unknown"))]
  x[, surgery_primary := factor(surgery_primary, levels = c("No", "Yes", "Unknown"))]
  radiation_label <- x[["Radiation recode"]]
  x[, radiation_any := factor(
    fifelse(
      radiation_label %chin% c("None/Unknown", "Refused (1988+)", "Recommended, unknown if administered"),
      "No/Unknown", "Yes"
    ),
    levels = c("No/Unknown", "Yes")
  )]
  x[, chemo_binary := factor(
    `Chemotherapy recode (yes, no/unk)`, levels = c("No/Unknown", "Yes")
  )]
  x[, cause_known := !is.na(competing_event)]
  x[, cs_event := as.integer(competing_event == 1L)]
  x[, other_event := as.integer(competing_event == 2L)]
  x[, age_group_exact := factor(
    fcase(age_years <= 39, "15-39", age_years <= 64, "40-64", default = "65-84"),
    levels = c("15-39", "40-64", "65-84")
  )]
  eod_schema <- x[["EOD Schema ID Recode (2010+)"]]
  x[, schema3 := factor(
    fcase(
      eod_schema == "Bile Ducts Intrahepat", "iCCA",
      eod_schema == "Bile Ducts Perihilar", "Perihilar eCCA",
      eod_schema == "Bile Duct Distal", "Distal eCCA",
      default = NA_character_
    ),
    levels = c("iCCA", "Perihilar eCCA", "Distal eCCA")
  )]
  x[, histology_code := suppressWarnings(as.integer(`Histologic Type ICD-O-3`))]
  x
}

d <- prepare_covariates_enhanced(readRDS(file.path(derived_dir, "cca_strict_hazard_2004_2023.rds")))
d_all <- prepare_covariates_enhanced(readRDS(file.path(derived_dir, "cca_strict_2000_2023.rds")))
d_cif <- prepare_covariates_enhanced(readRDS(file.path(derived_dir, "cca_strict_cif5_2004_2018.rds")))

d_exact_known <- droplevels(d[age_topcoded_85plus == FALSE & cause_known == TRUE])
d_all_exact_known <- droplevels(d_all[age_topcoded_85plus == FALSE & cause_known == TRUE])
d_cif <- droplevels(d_cif[age_topcoded_85plus == FALSE & cause_known == TRUE])

age_knots <- c(52, 68, 80)
age_bounds <- c(15, 84)
numvec <- function(x) paste(format(x, scientific = FALSE, trim = TRUE), collapse = ", ")
ns_term <- paste0(
  "splines::ns(age_years, knots = c(", numvec(age_knots),
  "), Boundary.knots = c(", numvec(age_bounds), "))"
)

cov_demographic <- c("sex_factor", "race_factor", "marital4", "income3", "rural3")
cov_model1 <- c(cov_demographic, "era4")
cov_model2_nograde <- c(cov_model1, "stage4", "surgery_primary", "radiation_any", "chemo_binary")

usable_covariates <- function(data, covars) {
  covars[vapply(covars, function(v) {
    z <- data[[v]]
    length(unique(z[!is.na(z)])) > 1L
  }, logical(1))]
}

lrt_result <- function(reduced, full) {
  chisq <- 2 * (full$loglik[2L] - reduced$loglik[2L])
  df <- sum(!is.na(coef(full))) - sum(!is.na(coef(reduced)))
  list(chisq = unname(chisq), df = df, p = pchisq(chisq, df, lower.tail = FALSE))
}

fit_interaction_tests <- function(data, exposure, covars, label) {
  data <- droplevels(as.data.frame(data))
  covars <- usable_covariates(data, covars)
  adjust <- if (length(covars)) paste0(" + ", paste(covars, collapse = " + ")) else ""
  lhs <- "Surv(survival_time_months, cs_event)"
  f_constant <- as.formula(paste0(lhs, " ~ ", ns_term, " + ", exposure, adjust))
  f_linear <- as.formula(paste0(lhs, " ~ ", ns_term, " + ", exposure, " + age_years:", exposure, adjust))
  f_spline <- as.formula(paste0(lhs, " ~ ", ns_term, " * ", exposure, adjust))
  fit_constant <- coxph(f_constant, data = data, ties = "efron")
  fit_linear <- coxph(f_linear, data = data, ties = "efron")
  fit_spline_lrt <- coxph(f_spline, data = data, ties = "efron")
  fit_spline_robust <- coxph(
    f_spline, data = data, ties = "efron", cluster = data$patient_cluster, robust = TRUE
  )
  overall <- lrt_result(fit_constant, fit_spline_lrt)
  nonlinear <- lrt_result(fit_linear, fit_spline_lrt)
  list(
    robust_fit = fit_spline_robust,
    tests = rbindlist(list(
      data.table(
        analysis = label, test = "Overall age-by-site interaction",
        n = nrow(data), cancer_deaths = sum(data$cs_event),
        chisq = overall$chisq, df = overall$df, p = overall$p
      ),
      data.table(
        analysis = label, test = "Nonlinear interaction component",
        n = nrow(data), cancer_deaths = sum(data$cs_event),
        chisq = nonlinear$chisq, df = nonlinear$df, p = nonlinear$p
      )
    ))
  )
}

reference_values <- list(
  sex_factor = "Female", race_factor = "Non-Hispanic White", marital4 = "Married",
  income3 = "70,000-99,999", rural3 = "Metropolitan", era4 = "2010-2015",
  era5 = "2010-2015", stage4 = "Localized", surgery_primary = "No",
  radiation_any = "No/Unknown", chemo_binary = "No/Unknown"
)

fill_model_factors <- function(fit, nd) {
  for (nm in setdiff(names(fit$xlevels), names(nd))) {
    lv <- fit$xlevels[[nm]]
    preferred <- reference_values[[nm]]
    value <- if (!is.null(preferred) && preferred %in% lv) preferred else lv[1L]
    nd[[nm]] <- factor(value, levels = lv)
  }
  nd
}

prediction_terms <- function(fit) delete.response(terms(fit))

site_contrasts <- function(fit, ages) {
  output <- lapply(ages, function(age) {
    nd <- data.frame(age_years = c(age, age))
    nd$site_group <- factor(c("iCCA", "eCCA"), levels = fit$xlevels$site_group)
    nd <- fill_model_factors(fit, nd)
    mm <- model.matrix(prediction_terms(fit), nd, contrasts.arg = fit$contrasts)
    keep <- !is.na(coef(fit))
    beta <- coef(fit)[keep]
    variance <- vcov(fit)[keep, keep, drop = FALSE]
    mm <- mm[, names(beta), drop = FALSE]
    z <- mm[2L, ] - mm[1L, ]
    log_hr <- sum(z * beta)
    se <- sqrt(drop(z %*% variance %*% z))
    data.table(age = age, estimate = exp(log_hr), lower95 = exp(log_hr - 1.96 * se), upper95 = exp(log_hr + 1.96 * se))
  })
  rbindlist(output)
}

# 1. Formal tests of the nonlinear component of the age interaction.
nonlinear_binary_list <- lapply(names(list(`Model 1` = cov_model1, `Model 2 without grade` = cov_model2_nograde)), function(nm) {
  specs <- list(`Model 1` = cov_model1, `Model 2 without grade` = cov_model2_nograde)
  fit <- fit_interaction_tests(d_exact_known, "site_group", specs[[nm]], "Binary anatomy, 2004-2023")
  fit$tests[, model := nm]
  fit
})
nonlinear_binary <- rbindlist(lapply(nonlinear_binary_list, `[[`, "tests"), fill = TRUE)

schema_data <- droplevels(d_exact_known[diagnosis_year >= 2010L & !is.na(schema3)])
schema_data[, era_schema := factor(
  fcase(diagnosis_year <= 2015L, "2010-2015", diagnosis_year <= 2019L, "2016-2019", default = "2020-2023"),
  levels = c("2010-2015", "2016-2019", "2020-2023")
)]
schema_tests_list <- lapply(
  list(
    `Model 1` = c(cov_demographic, "era_schema"),
    `Model 2 without grade` = c(cov_demographic, "era_schema", "stage4", "surgery_primary", "radiation_any", "chemo_binary")
  ),
  function(covars) fit_interaction_tests(schema_data, "schema3", covars, "Three-category EOD anatomy, 2010-2023")
)
for (i in seq_along(schema_tests_list)) schema_tests_list[[i]]$tests[, model := names(schema_tests_list)[i]]
nonlinear_schema <- rbindlist(lapply(schema_tests_list, `[[`, "tests"), fill = TRUE)

# 2. Full-period sensitivity analysis for relative effects.
full_covars <- c(cov_demographic, "era5")
full_period <- fit_interaction_tests(
  d_all_exact_known, "site_group", full_covars, "Binary anatomy, 2000-2023"
)
full_period$tests[, model := "Model 1"]
full_period_effects <- site_contrasts(full_period$robust_fit, c(20, 30, 40, 50, 65, 75))
full_period_effects[, `:=`(analysis = "Binary anatomy, 2000-2023", model = "Model 1")]

# 3. Descriptive sequential adjustment at fixed ages.
sequential_specs <- list(
  `Demographic and era` = cov_model1,
  `Plus stage` = c(cov_model1, "stage4"),
  `Plus stage and surgery` = c(cov_model1, "stage4", "surgery_primary"),
  `Plus stage and recorded treatment` = cov_model2_nograde
)
sequential_results <- rbindlist(lapply(names(sequential_specs), function(nm) {
  fit <- fit_interaction_tests(d_exact_known, "site_group", sequential_specs[[nm]], paste0("Sequential adjustment: ", nm))
  ans <- site_contrasts(fit$robust_fit, c(40, 50, 65, 75))
  ans[, model := nm]
  ans
}))

# 4. Distribution of histology 8162 across topography and EOD schema.
coding_data <- copy(d_all[diagnosis_year >= 2010L])
coding_data[, topography := factor(
  fifelse(site_group == "iCCA", "C22.1 (intrahepatic bile duct)", "C24.0 (extrahepatic bile duct)")
)]
coding_data[, histology := factor(paste0(histology_code, ": ", fifelse(histology_code == 8162L, "Klatskin tumor", "cholangiocarcinoma")))]
coding_data[, eod_anatomy := fcase(
  `EOD Schema ID Recode (2010+)` == "Bile Ducts Intrahepat", "Intrahepatic",
  `EOD Schema ID Recode (2010+)` == "Bile Ducts Perihilar", "Perihilar",
  `EOD Schema ID Recode (2010+)` == "Bile Duct Distal", "Distal",
  `EOD Schema ID Recode (2010+)` == "Cystic Duct", "Cystic duct",
  default = "Other or unavailable"
)]
histology_schema_distribution <- coding_data[, .N, by = .(topography, histology, eod_anatomy)]
histology_schema_distribution[, percent_within_topography_histology := 100 * N / sum(N), by = .(topography, histology)]
setorder(histology_schema_distribution, topography, histology, -N)
histology_topography_summary <- d_all[, .N, by = .(
  period = fifelse(diagnosis_year >= 2010L, "2010-2023", "2000-2009"),
  topography = fifelse(site_group == "iCCA", "C22.1", "C24.0"),
  histology_code
)]
histology_topography_summary[, percent_within_topography := 100 * N / sum(N), by = .(period, topography)]
setorder(histology_topography_summary, period, topography, histology_code)

full_period_age_counts <- rbindlist(list(
  d_exact_known[, .(N = .N, cancer_deaths = sum(cs_event)), by = .(age_group_exact, site_group)][, period := "2004-2023"],
  d_all_exact_known[, .(N = .N, cancer_deaths = sum(cs_event)), by = .(age_group_exact, site_group)][, period := "2000-2023"]
))
setorder(full_period_age_counts, period, age_group_exact, site_group)

# 5. Quantify the absolute-risk distortion from imposing a constant site effect.
fit_cif_models <- function(data, interaction = TRUE) {
  data <- droplevels(as.data.frame(data))
  age_site <- if (interaction) paste0(ns_term, " * site_group") else paste0(ns_term, " + site_group")
  rhs <- paste(c(age_site, cov_model1), collapse = " + ")
  fit_cancer <- coxph(as.formula(paste0("Surv(survival_time_months, cs_event) ~ ", rhs)), data = data, ties = "efron")
  fit_other <- coxph(as.formula(paste0("Surv(survival_time_months, other_event) ~ ", rhs)), data = data, ties = "efron")
  list(
    cancer = fit_cancer, other = fit_other,
    bh_cancer = suppressWarnings(basehaz(fit_cancer, centered = TRUE)),
    bh_other = suppressWarnings(basehaz(fit_other, centered = TRUE))
  )
}

cumhaz_at <- function(base_hazard, times) {
  idx <- findInterval(times, base_hazard$time)
  c(0, base_hazard$hazard)[idx + 1L]
}

standardized_cif <- function(fits, newdata, horizon = 60) {
  event_times <- sort(unique(c(
    fits$bh_cancer$time[fits$bh_cancer$time <= horizon],
    fits$bh_other$time[fits$bh_other$time <= horizon], horizon
  )))
  hc <- cumhaz_at(fits$bh_cancer, event_times)
  ho <- cumhaz_at(fits$bh_other, event_times)
  dhc <- diff(c(0, hc))
  dho <- diff(c(0, ho))
  rc <- exp(predict(fits$cancer, newdata = newdata, type = "lp", reference = "sample"))
  ro <- exp(predict(fits$other, newdata = newdata, type = "lp", reference = "sample"))
  survival_probability <- rep(1, nrow(newdata))
  cif <- rep(0, nrow(newdata))
  for (k in seq_along(event_times)) {
    total_increment <- dhc[k] * rc + dho[k] * ro
    cancer_fraction <- ifelse(total_increment > 0, dhc[k] * rc / total_increment, 0)
    cif <- cif + survival_probability * (1 - exp(-total_increment)) * cancer_fraction
    survival_probability <- survival_probability * exp(-total_increment)
  }
  mean(cif)
}

cif_fixed_ages <- function(data, fits, model_type, ages = c(40, 50, 65, 75)) {
  rbindlist(lapply(ages, function(age) {
    risks <- sapply(levels(data$site_group), function(site) {
      nd <- copy(data)
      nd[, age_years := age]
      nd[, site_group := factor(site, levels = levels(data$site_group))]
      standardized_cif(fits, as.data.frame(nd), 60)
    })
    data.table(
      age = age, model_type = model_type,
      iCCA = unname(risks["iCCA"]), eCCA = unname(risks["eCCA"]),
      risk_difference_eCCA_minus_iCCA = unname(risks["eCCA"] - risks["iCCA"])
    )
  }))
}

flexible_fits <- fit_cif_models(d_cif, interaction = TRUE)
constant_fits <- fit_cif_models(d_cif, interaction = FALSE)
cif_model_comparison_point <- rbindlist(list(
  cif_fixed_ages(d_cif, flexible_fits, "Age-specific spline interaction"),
  cif_fixed_ages(d_cif, constant_fits, "Constant site effect")
))
cif_distortion_point <- dcast(
  cif_model_comparison_point, age ~ model_type,
  value.var = "risk_difference_eCCA_minus_iCCA"
)
setnames(
  cif_distortion_point,
  c("Age-specific spline interaction", "Constant site effect"),
  c("rd_flexible", "rd_constant")
)
cif_distortion_point[, constant_minus_flexible := rd_constant - rd_flexible]

bootstrap_replicates <- as.integer(Sys.getenv("CCA_ENHANCEMENT_BOOTSTRAP_B", unset = "500"))
bootstrap_seed <- 20260910L
cluster_rows <- split(seq_len(nrow(d_cif)), d_cif$patient_id)

bootstrap_one <- function(iteration) {
  sampled_clusters <- sample.int(length(cluster_rows), length(cluster_rows), replace = TRUE)
  sampled_rows <- unlist(cluster_rows[sampled_clusters], use.names = FALSE)
  db <- droplevels(d_cif[sampled_rows])
  tryCatch({
    flexible <- cif_fixed_ages(db, fit_cif_models(db, TRUE), "Age-specific spline interaction")
    constant <- cif_fixed_ages(db, fit_cif_models(db, FALSE), "Constant site effect")
    point <- rbindlist(list(flexible, constant))
    wide <- dcast(point, age ~ model_type, value.var = "risk_difference_eCCA_minus_iCCA")
    setnames(wide, c("Age-specific spline interaction", "Constant site effect"), c("rd_flexible", "rd_constant"))
    wide[, `:=`(constant_minus_flexible = rd_constant - rd_flexible, replicate = iteration)]
    list(comparison = point[, replicate := iteration], distortion = wide)
  }, error = function(e) NULL)
}

if (bootstrap_replicates > 0L) {
  worker_count <- min(3L, max(1L, parallel::detectCores(logical = FALSE) - 1L))
  cl <- parallel::makeCluster(worker_count)
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(data.table)
      library(survival)
      library(splines)
    })
    NULL
  })
  parallel::clusterExport(
    cl,
    c(
      "d_cif", "cluster_rows", "fit_cif_models", "cif_fixed_ages",
      "standardized_cif", "cumhaz_at", "ns_term", "cov_model1"
    ),
    envir = environment()
  )
  parallel::clusterSetRNGStream(cl, iseed = bootstrap_seed)
  bootstrap_results <- tryCatch(
    parallel::parLapplyLB(cl, seq_len(bootstrap_replicates), bootstrap_one),
    finally = parallel::stopCluster(cl)
  )
  bootstrap_results <- Filter(Negate(is.null), bootstrap_results)
} else {
  bootstrap_results <- list()
}

if (length(bootstrap_results)) {
  boot_comparison <- rbindlist(lapply(bootstrap_results, `[[`, "comparison"))
  boot_distortion <- rbindlist(lapply(bootstrap_results, `[[`, "distortion"))
  comparison_ci <- boot_comparison[, .(
    rd_lower95 = quantile(risk_difference_eCCA_minus_iCCA, 0.025, na.rm = TRUE, type = 6),
    rd_upper95 = quantile(risk_difference_eCCA_minus_iCCA, 0.975, na.rm = TRUE, type = 6)
  ), by = .(age, model_type)]
  distortion_ci <- boot_distortion[, .(
    distortion_lower95 = quantile(constant_minus_flexible, 0.025, na.rm = TRUE, type = 6),
    distortion_upper95 = quantile(constant_minus_flexible, 0.975, na.rm = TRUE, type = 6)
  ), by = age]
  cif_model_comparison <- merge(cif_model_comparison_point, comparison_ci, by = c("age", "model_type"), all.x = TRUE)
  cif_distortion <- merge(cif_distortion_point, distortion_ci, by = "age", all.x = TRUE)
} else {
  boot_comparison <- data.table()
  boot_distortion <- data.table()
  cif_model_comparison <- copy(cif_model_comparison_point)[, `:=`(rd_lower95 = NA_real_, rd_upper95 = NA_real_)]
  cif_distortion <- copy(cif_distortion_point)[, `:=`(distortion_lower95 = NA_real_, distortion_upper95 = NA_real_)]
}

fwrite(nonlinear_binary, file.path(tables_dir, "enhancement_nonlinearity_binary.csv"))
fwrite(nonlinear_schema, file.path(tables_dir, "enhancement_nonlinearity_schema3.csv"))
fwrite(full_period$tests, file.path(tables_dir, "enhancement_fullperiod_tests.csv"))
fwrite(full_period_effects, file.path(tables_dir, "enhancement_fullperiod_effects.csv"))
fwrite(sequential_results, file.path(tables_dir, "enhancement_sequential_adjustment.csv"))
fwrite(histology_schema_distribution, file.path(tables_dir, "enhancement_histology_schema_distribution.csv"))
fwrite(histology_topography_summary, file.path(tables_dir, "enhancement_histology_topography_summary.csv"))
fwrite(full_period_age_counts, file.path(tables_dir, "enhancement_fullperiod_age_counts.csv"))
fwrite(cif_model_comparison, file.path(tables_dir, "enhancement_cif_constant_vs_flexible.csv"))
fwrite(cif_distortion, file.path(tables_dir, "enhancement_cif_distortion.csv"))
if (nrow(boot_comparison)) fwrite(boot_comparison, file.path(source_dir, "enhancement_cif_model_bootstrap.csv.gz"))
if (nrow(boot_distortion)) fwrite(boot_distortion, file.path(source_dir, "enhancement_cif_distortion_bootstrap.csv.gz"))

saveRDS(
  list(
    nonlinear_binary = nonlinear_binary_list,
    nonlinear_schema = schema_tests_list,
    full_period = full_period,
    flexible_cif_fits = flexible_fits,
    constant_cif_fits = constant_fits,
    bootstrap_requested = bootstrap_replicates,
    bootstrap_successful = length(bootstrap_results),
    bootstrap_seed = bootstrap_seed
  ),
  file.path(models_dir, "manuscript_enhancement_models.rds")
)

qc <- data.table(
  check = c(
    "Primary nonlinear test available", "Schema nonlinear test available",
    "Full-period sensitivity N", "C22.1 plus 8162 count, 2010-2023",
    "CIF bootstrap requested", "CIF bootstrap successful"
  ),
  value = c(
    nrow(nonlinear_binary[test == "Nonlinear interaction component" & model == "Model 1"]),
    nrow(nonlinear_schema[test == "Nonlinear interaction component" & model == "Model 1"]),
    nrow(d_all_exact_known),
    coding_data[site_group == "iCCA" & histology_code == 8162L, .N],
    bootstrap_replicates, length(bootstrap_results)
  )
)
fwrite(qc, file.path(results_dir, "manuscript_enhancement_qc.csv"))

cat("Binary nonlinearity tests:\n")
print(nonlinear_binary)
cat("Three-category nonlinearity tests:\n")
print(nonlinear_schema)
cat("Full-period sensitivity:\n")
print(full_period$tests)
print(full_period_effects)
cat("Constant-effect distortion of standardized five-year risk differences:\n")
print(cif_distortion)
cat("Bootstrap successful:", length(bootstrap_results), "of", bootstrap_replicates, "\n")
