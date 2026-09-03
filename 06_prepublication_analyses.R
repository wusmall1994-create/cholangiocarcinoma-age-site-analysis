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
  x[, patient_cluster := factor(patient_id)]
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
  x[, cause_known := !is.na(competing_event)]
  x[, cs_event := as.integer(competing_event == 1L)]
  x[, other_event := as.integer(competing_event == 2L)]
  x[, age85_flag := factor(
    fifelse(age_topcoded_85plus, "85+", "15-84"),
    levels = c("15-84", "85+")
  )]
  x[, age_group_exact := factor(
    fcase(
      age_years <= 39, "15-39",
      age_years <= 64, "40-64",
      default = "65-84"
    ),
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

d <- prepare_covariates(readRDS(file.path(derived_dir, "cca_strict_hazard_2004_2023.rds")))
d_cif <- prepare_covariates(readRDS(file.path(derived_dir, "cca_strict_cif5_2004_2018.rds")))

d_exact <- d[age_topcoded_85plus == FALSE]
d_exact_known <- d_exact[cause_known == TRUE]
age_knots <- as.numeric(quantile(
  d_exact$age_years,
  probs = c(0.10, 0.50, 0.90), na.rm = TRUE, type = 2
))
age_bounds <- range(d_exact$age_years, na.rm = TRUE)
numvec <- function(x) paste(format(x, scientific = FALSE, trim = TRUE), collapse = ", ")
ns_term <- paste0(
  "splines::ns(age_years, knots = c(", numvec(age_knots),
  "), Boundary.knots = c(", numvec(age_bounds), "))"
)
ns_all_term <- paste0(
  "splines::ns(age_years, knots = c(", numvec(age_knots),
  "), Boundary.knots = c(15, 85))"
)

cov_demographic <- c("sex_factor", "race_factor", "marital4", "income3", "rural3")
cov_model1 <- c(cov_demographic, "era4")
cov_model2 <- c(cov_model1, "stage4", "grade4", "surgery_primary", "radiation_any", "chemo_binary")
cov_model2_nograde <- setdiff(cov_model2, "grade4")

usable_covariates <- function(data, covars) {
  covars[vapply(covars, function(v) {
    z <- data[[v]]
    length(unique(z[!is.na(z)])) > 1L
  }, logical(1))]
}

fit_age_exposure <- function(data, exposure, covars, time_var = "survival_time_months", event_var = "cs_event") {
  data <- droplevels(as.data.frame(data))
  covars <- usable_covariates(data, covars)
  adjust <- if (length(covars)) paste0(" + ", paste(covars, collapse = " + ")) else ""
  survival_lhs <- paste0("Surv(", time_var, ", ", event_var, ")")
  reduced_formula <- as.formula(paste0(survival_lhs, " ~ ", ns_term, " + ", exposure, adjust))
  full_formula <- as.formula(paste0(survival_lhs, " ~ ", ns_term, " * ", exposure, adjust))
  reduced <- coxph(reduced_formula, data = data, ties = "efron")
  full_lrt <- coxph(full_formula, data = data, ties = "efron")
  full <- coxph(
    full_formula, data = data, ties = "efron",
    cluster = data$patient_cluster, robust = TRUE
  )
  chisq <- 2 * (full_lrt$loglik[2L] - reduced$loglik[2L])
  test_df <- sum(!is.na(coef(full_lrt))) - sum(!is.na(coef(reduced)))
  list(
    fit = full,
    lrt_fit = full_lrt,
    reduced = reduced,
    chisq = unname(chisq),
    df = test_df,
    p = pchisq(chisq, test_df, lower.tail = FALSE),
    n = nrow(data),
    events = sum(data[[event_var]] == 1L)
  )
}

reference_values <- list(
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
  chemo_binary = "No/Unknown"
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

prediction_terms <- function(fit) {
  tt <- delete.response(terms(fit))
  strata_terms <- grep("^strata\\(", attr(tt, "term.labels"))
  if (length(strata_terms)) {
    tt <- drop.terms(tt, strata_terms, keep.response = FALSE)
  }
  tt
}

exposure_contrasts <- function(fit, exposure, ages, modifier = NULL, modifier_values = NULL) {
  exposure_levels <- fit$xlevels[[exposure]]
  ref <- exposure_levels[1L]
  if (is.null(modifier)) modifier_values <- NA_character_
  output <- list()
  counter <- 0L
  for (mod_value in modifier_values) {
    for (age in ages) {
      for (contrast_level in exposure_levels[-1L]) {
        nd <- data.frame(age_years = c(age, age))
        nd[[exposure]] <- factor(c(ref, contrast_level), levels = exposure_levels)
        if (!is.null(modifier)) {
          nd[[modifier]] <- factor(rep(mod_value, 2L), levels = fit$xlevels[[modifier]])
        }
        nd <- fill_model_factors(fit, nd)
        mm <- model.matrix(prediction_terms(fit), nd, contrasts.arg = fit$contrasts)
        keep <- !is.na(coef(fit))
        beta <- coef(fit)[keep]
        variance <- vcov(fit)[keep, keep, drop = FALSE]
        mm <- mm[, names(beta), drop = FALSE]
        z <- mm[2L, ] - mm[1L, ]
        estimate <- sum(z * beta)
        se <- sqrt(drop(z %*% variance %*% z))
        counter <- counter + 1L
        output[[counter]] <- data.table(
          age = age,
          reference = ref,
          contrast = contrast_level,
          modifier = if (is.null(modifier)) NA_character_ else mod_value,
          estimate = exp(estimate),
          lower95 = exp(estimate - 1.96 * se),
          upper95 = exp(estimate + 1.96 * se)
        )
      }
    }
  }
  rbindlist(output)
}

# 1. Primary exact-age models and an explicit 85+ indicator sensitivity model.
exact_specs <- list(
  `Model 1` = cov_model1,
  `Model 2` = cov_model2,
  `Model 2 without grade` = cov_model2_nograde
)
exact_fits <- lapply(exact_specs, function(covars) {
  fit_age_exposure(d_exact_known, "site_group", covars)
})
exact_tests <- rbindlist(lapply(names(exact_fits), function(nm) {
  fit <- exact_fits[[nm]]
  data.table(
    analysis = "Exact ages 15-84",
    model = nm,
    n = fit$n,
    cancer_deaths = fit$events,
    chisq = fit$chisq,
    df = fit$df,
    p = fit$p
  )
}))
exact_effects <- rbindlist(lapply(names(exact_fits), function(nm) {
  ans <- exposure_contrasts(exact_fits[[nm]]$fit, "site_group", c(20, 30, 40, 50, 65, 75))
  ans[, `:=`(analysis = "Exact ages 15-84", model = nm)]
  ans
}))

fit_age85_indicator <- function(data, covars) {
  data <- droplevels(as.data.frame(data[cause_known == TRUE]))
  covars <- usable_covariates(data, covars)
  adjust <- paste0(" + ", paste(covars, collapse = " + "))
  reduced_formula <- as.formula(paste0(
    "Surv(survival_time_months, cs_event) ~ ", ns_all_term,
    " + site_group + age85_flag", adjust
  ))
  full_formula <- as.formula(paste0(
    "Surv(survival_time_months, cs_event) ~ ", ns_all_term,
    " * site_group + age85_flag * site_group", adjust
  ))
  reduced <- coxph(reduced_formula, data = data, ties = "efron")
  full_lrt <- coxph(full_formula, data = data, ties = "efron")
  full <- coxph(full_formula, data = data, ties = "efron", cluster = data$patient_cluster, robust = TRUE)
  chisq <- 2 * (full_lrt$loglik[2L] - reduced$loglik[2L])
  test_df <- sum(!is.na(coef(full_lrt))) - sum(!is.na(coef(reduced)))
  list(
    fit = full,
    chisq = unname(chisq), df = test_df,
    p = pchisq(chisq, test_df, lower.tail = FALSE),
    n = nrow(data), events = sum(data$cs_event)
  )
}

age85_fits <- list(
  `Model 1` = fit_age85_indicator(d, cov_model1),
  `Model 2 without grade` = fit_age85_indicator(d, cov_model2_nograde)
)
age85_tests <- rbindlist(lapply(names(age85_fits), function(nm) {
  fit <- age85_fits[[nm]]
  data.table(
    analysis = "All ages; 85+ modeled as an indicator",
    model = nm, n = fit$n, cancer_deaths = fit$events,
    chisq = fit$chisq, df = fit$df, p = fit$p
  )
}))
age85_counts <- d[, .(
  n = .N,
  cancer_deaths = sum(competing_event == 1L, na.rm = TRUE),
  other_deaths = sum(competing_event == 2L, na.rm = TRUE),
  unknown_cause_deaths = sum(os_event == 1L & is.na(competing_event))
), by = .(age85_flag, site_group)]

# 2. Formal age-by-site-by-era interaction with common follow-up horizons.
fit_era_threeway <- function(data, covars, horizon, last_year) {
  ds <- copy(data[diagnosis_year <= last_year])
  last_label <- if (horizon == 24) "2020-2021" else "2020"
  ds[, era_common := factor(
    fcase(
      diagnosis_year <= 2009L, "2004-2009",
      diagnosis_year <= 2015L, "2010-2015",
      diagnosis_year <= 2019L, "2016-2019",
      default = last_label
    ),
    levels = c("2004-2009", "2010-2015", "2016-2019", last_label)
  )]
  ds[, common_time := pmin(survival_time_months, horizon)]
  ds[, common_event := as.integer(cs_event == 1L & survival_time_months <= horizon)]
  ds <- droplevels(as.data.frame(ds))
  covars <- setdiff(usable_covariates(ds, covars), "era4")
  adjust <- if (length(covars)) paste0(" + ", paste(covars, collapse = " + ")) else ""
  reduced_formula <- as.formula(paste0(
    "Surv(common_time, common_event) ~ ", ns_term, " * site_group + ",
    ns_term, " * era_common + site_group * era_common", adjust
  ))
  full_formula <- as.formula(paste0(
    "Surv(common_time, common_event) ~ ", ns_term,
    " * site_group * era_common", adjust
  ))
  reduced <- coxph(reduced_formula, data = ds, ties = "efron")
  full_lrt <- coxph(full_formula, data = ds, ties = "efron")
  full <- coxph(full_formula, data = ds, ties = "efron", cluster = ds$patient_cluster, robust = TRUE)
  chisq <- 2 * (full_lrt$loglik[2L] - reduced$loglik[2L])
  test_df <- sum(!is.na(coef(full_lrt))) - sum(!is.na(coef(reduced)))
  list(
    fit = full, lrt_fit = full_lrt, reduced = reduced,
    chisq = unname(chisq), df = test_df,
    p = pchisq(chisq, test_df, lower.tail = FALSE),
    n = nrow(ds), events = sum(ds$common_event), data = ds,
    horizon = horizon, last_year = last_year
  )
}

era_specs <- list(
  `Model 1` = cov_model1,
  `Model 2 without grade` = cov_model2_nograde
)
era_fits_24 <- lapply(era_specs, function(covars) {
  fit_era_threeway(d_exact_known, covars, horizon = 24, last_year = 2021)
})
era_fits_36 <- lapply(era_specs, function(covars) {
  fit_era_threeway(d_exact_known, covars, horizon = 36, last_year = 2020)
})
era_tests <- rbindlist(c(
  lapply(names(era_fits_24), function(nm) {
    fit <- era_fits_24[[nm]]
    data.table(
      horizon_months = 24L, through_year = 2021L, model = nm,
      n = fit$n, cancer_deaths = fit$events,
      chisq = fit$chisq, df = fit$df, p = fit$p
    )
  }),
  lapply(names(era_fits_36), function(nm) {
    fit <- era_fits_36[[nm]]
    data.table(
      horizon_months = 36L, through_year = 2020L, model = nm,
      n = fit$n, cancer_deaths = fit$events,
      chisq = fit$chisq, df = fit$df, p = fit$p
    )
  })
))
era_effect_grid <- seq(20, 84, by = 1)
era_effects_24 <- rbindlist(lapply(names(era_fits_24), function(nm) {
  fit <- era_fits_24[[nm]]
  ans <- exposure_contrasts(
    fit$fit, "site_group", era_effect_grid,
    modifier = "era_common", modifier_values = fit$fit$xlevels$era_common
  )
  ans[, `:=`(horizon_months = 24L, model = nm)]
  ans
}))
era_counts_24 <- as.data.table(era_fits_24[[1L]]$data)[, .(
  n = .N, cancer_deaths_within_horizon = sum(common_event)
), by = .(era_common, site_group, age_group_exact)]

# 3. Three-category anatomic schema analysis, restricted to reliable EOD schema years.
schema_data <- droplevels(d_exact_known[
  diagnosis_year >= 2010L & !is.na(schema3)
])
schema_data[, era_schema := factor(
  fcase(
    diagnosis_year <= 2015L, "2010-2015",
    diagnosis_year <= 2019L, "2016-2019",
    default = "2020-2023"
  ),
  levels = c("2010-2015", "2016-2019", "2020-2023")
)]
schema_cov_m1 <- c(cov_demographic, "era_schema")
schema_cov_m2_nograde <- c(schema_cov_m1, "stage4", "surgery_primary", "radiation_any", "chemo_binary")
schema_fits <- list(
  `Model 1` = fit_age_exposure(schema_data, "schema3", schema_cov_m1),
  `Model 2 without grade` = fit_age_exposure(schema_data, "schema3", schema_cov_m2_nograde)
)
schema_tests <- rbindlist(lapply(names(schema_fits), function(nm) {
  fit <- schema_fits[[nm]]
  data.table(
    analysis = "EOD schema 2010-2023; cystic duct excluded",
    model = nm, n = fit$n, cancer_deaths = fit$events,
    chisq = fit$chisq, df = fit$df, p = fit$p
  )
}))
schema_grid <- seq(20, 84, by = 1)
schema_effects <- rbindlist(lapply(names(schema_fits), function(nm) {
  ans <- exposure_contrasts(schema_fits[[nm]]$fit, "schema3", schema_grid)
  ans[, model := nm]
  ans
}))
schema_counts <- schema_data[, .(
  n = .N,
  cancer_deaths = sum(cs_event),
  other_deaths = sum(other_event)
), by = .(schema3, age_group_exact)]

schema_crosswalk <- d[diagnosis_year >= 2010L, .N, by = .(
  site_group,
  eod_schema = `EOD Schema ID Recode (2010+)`,
  tnm_schema = `TNM 7/CS v0204+ Schema recode`
)]
setorder(schema_crosswalk, site_group, eod_schema, tnm_schema)

hist8160_data <- droplevels(d_exact_known[histology_code == 8160L])
hist8160_fits <- list(
  `Model 1` = fit_age_exposure(hist8160_data, "site_group", cov_model1),
  `Model 2 without grade` = fit_age_exposure(hist8160_data, "site_group", cov_model2_nograde)
)
hist8160_tests <- rbindlist(lapply(names(hist8160_fits), function(nm) {
  fit <- hist8160_fits[[nm]]
  data.table(
    analysis = "Histology 8160 only; exact ages 15-84",
    model = nm, n = fit$n, cancer_deaths = fit$events,
    chisq = fit$chisq, df = fit$df, p = fit$p
  )
}))

# 4. Adjusted five-year CIF using cause-specific Cox g-formula standardization.
cif_data <- droplevels(d_cif[age_topcoded_85plus == FALSE & cause_known == TRUE])

fit_cif_models <- function(data) {
  data <- droplevels(as.data.frame(data))
  rhs <- paste(c(paste0(ns_term, " * site_group"), cov_model1), collapse = " + ")
  fit_cancer <- coxph(
    as.formula(paste0("Surv(survival_time_months, cs_event) ~ ", rhs)),
    data = data, ties = "efron"
  )
  fit_other <- coxph(
    as.formula(paste0("Surv(survival_time_months, other_event) ~ ", rhs)),
    data = data, ties = "efron"
  )
  list(
    cancer = fit_cancer,
    other = fit_other,
    bh_cancer = suppressWarnings(basehaz(fit_cancer, centered = TRUE)),
    bh_other = suppressWarnings(basehaz(fit_other, centered = TRUE))
  )
}

cumhaz_at <- function(base_hazard, times) {
  idx <- findInterval(times, base_hazard$time)
  c(0, base_hazard$hazard)[idx + 1L]
}

standardized_cif <- function(fits, newdata, horizon = 60) {
  bh_cancer <- fits$bh_cancer
  bh_other <- fits$bh_other
  event_times <- sort(unique(c(
    bh_cancer$time[bh_cancer$time <= horizon],
    bh_other$time[bh_other$time <= horizon],
    horizon
  )))
  hc <- cumhaz_at(bh_cancer, event_times)
  ho <- cumhaz_at(bh_other, event_times)
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

cif_targets <- function(data, fits, fixed_ages = c(20, 30, 40, 50, 65, 75)) {
  output <- list()
  counter <- 0L
  for (age in fixed_ages) {
    for (site in levels(data$site_group)) {
      nd <- copy(data)
      nd[, age_years := age]
      nd[, site_group := factor(site, levels = levels(data$site_group))]
      counter <- counter + 1L
      output[[counter]] <- data.table(
        target_type = "Fixed age",
        target = as.character(age),
        site_group = site,
        estimate = standardized_cif(fits, as.data.frame(nd), 60)
      )
    }
  }
  for (group_label in levels(data$age_group_exact)) {
    subgroup <- data[age_group_exact == group_label]
    for (site in levels(data$site_group)) {
      nd <- copy(subgroup)
      nd[, site_group := factor(site, levels = levels(data$site_group))]
      counter <- counter + 1L
      output[[counter]] <- data.table(
        target_type = "Age group",
        target = group_label,
        site_group = site,
        estimate = standardized_cif(fits, as.data.frame(nd), 60)
      )
    }
  }
  rbindlist(output)
}

add_risk_differences <- function(long_data, value_col = "estimate") {
  wide <- dcast(long_data, target_type + target ~ site_group, value.var = value_col)
  wide[, risk_difference_eCCA_minus_iCCA := eCCA - iCCA]
  wide
}

cif_fits <- fit_cif_models(cif_data)
cif_point <- cif_targets(cif_data, cif_fits)
cif_point_rd <- add_risk_differences(cif_point)

bootstrap_replicates <- as.integer(Sys.getenv("CCA_BOOTSTRAP_B", unset = "500"))
bootstrap_seed <- 20260902L
cluster_rows <- split(seq_len(nrow(cif_data)), cif_data$patient_id)

bootstrap_one <- function(iteration) {
  sampled_clusters <- sample.int(length(cluster_rows), length(cluster_rows), replace = TRUE)
  sampled_rows <- unlist(cluster_rows[sampled_clusters], use.names = FALSE)
  db <- droplevels(cif_data[sampled_rows])
  result <- tryCatch({
    fits <- fit_cif_models(db)
    long <- cif_targets(db, fits)
    wide <- add_risk_differences(long)
    list(long = long, wide = wide)
  }, error = function(e) NULL)
  result
}

worker_count <- min(3L, max(1L, parallel::detectCores(logical = FALSE) - 1L))
if (bootstrap_replicates > 0L) {
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
      "cif_data", "cluster_rows", "fit_cif_models", "cif_targets",
      "standardized_cif", "cumhaz_at", "add_risk_differences",
      "ns_term", "cov_model1"
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
  cif_boot_long <- rbindlist(lapply(seq_along(bootstrap_results), function(i) {
    ans <- copy(bootstrap_results[[i]]$long)
    ans[, replicate := i]
    ans
  }))
  cif_boot_rd <- rbindlist(lapply(seq_along(bootstrap_results), function(i) {
    ans <- copy(bootstrap_results[[i]]$wide)
    ans[, replicate := i]
    ans
  }))
  cif_ci <- cif_boot_long[, .(
    lower95 = quantile(estimate, 0.025, na.rm = TRUE, type = 6),
    upper95 = quantile(estimate, 0.975, na.rm = TRUE, type = 6)
  ), by = .(target_type, target, site_group)]
  cif_rd_ci <- cif_boot_rd[, .(
    rd_lower95 = quantile(risk_difference_eCCA_minus_iCCA, 0.025, na.rm = TRUE, type = 6),
    rd_upper95 = quantile(risk_difference_eCCA_minus_iCCA, 0.975, na.rm = TRUE, type = 6)
  ), by = .(target_type, target)]
  cif_adjusted <- merge(cif_point, cif_ci, by = c("target_type", "target", "site_group"), all.x = TRUE)
  cif_adjusted_rd <- merge(cif_point_rd, cif_rd_ci, by = c("target_type", "target"), all.x = TRUE)
} else {
  cif_boot_long <- data.table()
  cif_boot_rd <- data.table()
  cif_adjusted <- copy(cif_point)[, `:=`(lower95 = NA_real_, upper95 = NA_real_)]
  cif_adjusted_rd <- copy(cif_point_rd)[, `:=`(rd_lower95 = NA_real_, rd_upper95 = NA_real_)]
}

# 5. Six-interval time-varying age-by-site effects.
fit_time_varying <- function(data, covars) {
  dp <- survSplit(
    Surv(survival_time_months, cs_event) ~ .,
    data = as.data.frame(data), cut = c(6, 12, 24, 36, 60),
    start = "tstart", end = "tstop", event = "cs_event", episode = "period_code"
  )
  dp$followup_period <- factor(
    dp$period_code,
    levels = 1:6,
    labels = c("0-6", "6-12", "12-24", "24-36", "36-60", ">60")
  )
  covars <- usable_covariates(dp, covars)
  adjust <- paste0(" + ", paste(covars, collapse = " + "))
  reduced_formula <- as.formula(paste0(
    "Surv(tstart, tstop, cs_event) ~ ", ns_term, " * site_group + ",
    ns_term, " * followup_period + site_group * followup_period",
    adjust, " + strata(followup_period)"
  ))
  full_formula <- as.formula(paste0(
    "Surv(tstart, tstop, cs_event) ~ ", ns_term,
    " * site_group * followup_period", adjust,
    " + strata(followup_period)"
  ))
  reduced <- coxph(reduced_formula, data = dp, ties = "efron")
  full_lrt <- coxph(full_formula, data = dp, ties = "efron")
  full <- coxph(
    full_formula, data = dp, ties = "efron",
    cluster = dp$patient_cluster, robust = TRUE
  )
  chisq <- 2 * (full_lrt$loglik[2L] - reduced$loglik[2L])
  test_df <- sum(!is.na(coef(full_lrt))) - sum(!is.na(coef(reduced)))
  list(
    fit = full, lrt_fit = full_lrt, reduced = reduced,
    chisq = unname(chisq), df = test_df,
    p = pchisq(chisq, test_df, lower.tail = FALSE),
    n = nrow(data), events = sum(data$cs_event), split_rows = nrow(dp)
  )
}

time_fits <- list(
  `Model 1` = fit_time_varying(d_exact_known, cov_model1),
  `Model 2 without grade` = fit_time_varying(d_exact_known, cov_model2_nograde)
)
time_tests <- rbindlist(lapply(names(time_fits), function(nm) {
  fit <- time_fits[[nm]]
  data.table(
    analysis = "Six-interval time-varying age-by-site effect",
    model = nm, n = fit$n, cancer_deaths = fit$events,
    split_rows = fit$split_rows,
    chisq = fit$chisq, df = fit$df, p = fit$p
  )
}))
time_effects <- rbindlist(lapply(names(time_fits), function(nm) {
  fit <- time_fits[[nm]]
  ans <- exposure_contrasts(
    fit$fit, "site_group", c(20, 30, 40, 50, 65, 75),
    modifier = "followup_period",
    modifier_values = fit$fit$xlevels$followup_period
  )
  ans[, model := nm]
  ans
}))
setnames(time_effects, "modifier", "followup_period")

# Outputs.
fwrite(rbindlist(list(exact_tests, age85_tests), fill = TRUE), file.path(tables_dir, "prepub_age85_interaction_tests.csv"))
fwrite(exact_effects, file.path(tables_dir, "prepub_exact_age_effects.csv"))
fwrite(age85_counts, file.path(tables_dir, "prepub_age85_counts.csv"))
fwrite(era_tests, file.path(tables_dir, "prepub_era_threeway_tests.csv"))
fwrite(era_counts_24, file.path(tables_dir, "prepub_era_common_horizon_counts.csv"))
fwrite(era_effects_24, file.path(tables_dir, "prepub_era_24m_effect_curves.csv"))
fwrite(schema_tests, file.path(tables_dir, "prepub_schema3_interaction_tests.csv"))
fwrite(schema_counts, file.path(tables_dir, "prepub_schema3_counts.csv"))
fwrite(schema_effects, file.path(tables_dir, "prepub_schema3_effect_curves.csv"))
fwrite(schema_crosswalk, file.path(tables_dir, "prepub_schema_crosswalk.csv"))
fwrite(hist8160_tests, file.path(tables_dir, "prepub_histology8160_tests.csv"))
fwrite(cif_adjusted, file.path(tables_dir, "prepub_adjusted_cif5.csv"))
fwrite(cif_adjusted_rd, file.path(tables_dir, "prepub_adjusted_cif5_risk_differences.csv"))
fwrite(time_tests, file.path(tables_dir, "prepub_time_varying_tests.csv"))
fwrite(time_effects, file.path(tables_dir, "prepub_time_varying_effects.csv"))

fwrite(era_effects_24[model == "Model 1"], file.path(source_dir, "figure5_era_24m_curves.csv"))
fwrite(schema_effects[model == "Model 1"], file.path(source_dir, "figure5_schema3_curves.csv"))
fwrite(cif_adjusted[target_type == "Fixed age"], file.path(source_dir, "figure5_adjusted_cif.csv"))
fwrite(cif_adjusted_rd[target_type == "Fixed age"], file.path(source_dir, "figure5_adjusted_cif_rd.csv"))
fwrite(time_effects[model == "Model 1"], file.path(source_dir, "figure5_time_varying.csv"))
if (nrow(cif_boot_long)) {
  fwrite(cif_boot_long, file.path(source_dir, "adjusted_cif_bootstrap_replicates.csv.gz"))
  fwrite(cif_boot_rd, file.path(source_dir, "adjusted_cif_rd_bootstrap_replicates.csv.gz"))
}

saveRDS(
  list(
    age_knots = age_knots,
    age_bounds = age_bounds,
    exact_fits = exact_fits,
    age85_fits = age85_fits,
    era_fits_24 = era_fits_24,
    era_fits_36 = era_fits_36,
    schema_fits = schema_fits,
    hist8160_fits = hist8160_fits,
    cif_fits = cif_fits,
    time_fits = time_fits,
    bootstrap_requested = bootstrap_replicates,
    bootstrap_successful = length(bootstrap_results),
    bootstrap_seed = bootstrap_seed
  ),
  file.path(models_dir, "prepublication_models.rds"),
  compress = "xz"
)

qc <- data.table(
  item = c(
    "Exact-age primary cohort N", "Exact-age primary cancer deaths",
    "Top-coded 85+ records", "24-month era cohort N",
    "Schema-3 cohort N", "Adjusted CIF cohort N",
    "Bootstrap requested", "Bootstrap successful",
    "Time-varying split rows"
  ),
  value = c(
    nrow(d_exact_known), sum(d_exact_known$cs_event),
    sum(d$age_topcoded_85plus), era_fits_24[[1L]]$n,
    nrow(schema_data), nrow(cif_data),
    bootstrap_replicates, length(bootstrap_results),
    time_fits[[1L]]$split_rows
  )
)
fwrite(qc, file.path(results_dir, "prepublication_qc.csv"))

cat("Age knots:", paste(age_knots, collapse = ", "), "\n")
cat("Exact-age interaction tests:\n")
print(exact_tests)
cat("Era three-way tests:\n")
print(era_tests)
cat("Schema-3 interaction tests:\n")
print(schema_tests)
cat("Histology 8160-only tests:\n")
print(hist8160_tests)
cat("Time-varying tests:\n")
print(time_tests)
cat("Adjusted CIF risk differences:\n")
print(cif_adjusted_rd)
cat("Bootstrap successful:", length(bootstrap_results), "of", bootstrap_replicates, "\n")
