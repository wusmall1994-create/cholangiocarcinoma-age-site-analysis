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
  x[, era_schema_risk := factor(
    fifelse(diagnosis_year <= 2015L, "2010-2015", "2016-2018"),
    levels = c("2010-2015", "2016-2018")
  )]
  x[, schema3 := factor(fcase(
    `EOD Schema ID Recode (2010+)` == "Bile Ducts Intrahepat", "iCCA",
    `EOD Schema ID Recode (2010+)` == "Bile Ducts Perihilar", "Perihilar eCCA",
    `EOD Schema ID Recode (2010+)` == "Bile Duct Distal", "Distal eCCA",
    default = NA_character_
  ), levels = c("iCCA", "Perihilar eCCA", "Distal eCCA"))]
  x[, cause_known := !is.na(competing_event)]
  x[, cs_event := as.integer(competing_event == 1L)]
  x[, other_event := as.integer(competing_event == 2L)]
  x
}

d_hazard <- prepare_covariates(readRDS(file.path(derived_dir, "cca_strict_hazard_2004_2023.rds")))
d_cif <- prepare_covariates(readRDS(file.path(derived_dir, "cca_strict_cif5_2004_2018.rds")))

age_reference <- d_hazard[age_topcoded_85plus == FALSE]
age_knots <- as.numeric(quantile(age_reference$age_years, c(0.10, 0.50, 0.90), na.rm = TRUE, type = 2))
age_bounds <- range(age_reference$age_years, na.rm = TRUE)
numvec <- function(x) paste(format(x, scientific = FALSE, trim = TRUE), collapse = ", ")
ns_term <- paste0(
  "splines::ns(age_years, knots = c(", numvec(age_knots),
  "), Boundary.knots = c(", numvec(age_bounds), "))"
)

cov_binary <- c("sex_factor", "race_factor", "marital4", "income3", "rural3", "era4")
cov_schema <- c("sex_factor", "race_factor", "marital4", "income3", "rural3", "era_schema_risk")
ages <- c(40, 50, 65, 75)
horizons <- c(12, 36, 60)

fit_competing_models <- function(data, exposure, covars) {
  rhs <- paste(c(paste0(ns_term, " * ", exposure), covars), collapse = " + ")
  cancer <- coxph(as.formula(paste0("Surv(survival_time_months, cs_event) ~ ", rhs)),
                   data = data, ties = "efron")
  other <- coxph(as.formula(paste0("Surv(survival_time_months, other_event) ~ ", rhs)),
                  data = data, ties = "efron")
  list(
    cancer = cancer,
    other = other,
    bh_cancer = suppressWarnings(basehaz(cancer, centered = TRUE)),
    bh_other = suppressWarnings(basehaz(other, centered = TRUE))
  )
}

cumhaz_at <- function(base_hazard, times) {
  idx <- findInterval(times, base_hazard$time)
  c(0, base_hazard$hazard)[idx + 1L]
}

standardized_cif <- function(fits, newdata, horizon) {
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
  surv <- rep(1, nrow(newdata))
  cif <- rep(0, nrow(newdata))
  for (k in seq_along(event_times)) {
    total <- dhc[k] * rc + dho[k] * ro
    fraction <- ifelse(total > 0, dhc[k] * rc / total, 0)
    cif <- cif + surv * (1 - exp(-total)) * fraction
    surv <- surv * exp(-total)
  }
  mean(cif)
}

predict_grid <- function(data, fits, exposure, exposure_levels, horizons, ages) {
  rbindlist(lapply(horizons, function(horizon) {
    rbindlist(lapply(ages, function(age) {
      rbindlist(lapply(exposure_levels, function(level) {
        nd <- copy(data)
        nd[, age_years := age]
        set(nd, j = exposure, value = factor(level, levels = exposure_levels))
        data.table(
          horizon_months = horizon,
          age = age,
          exposure_level = level,
          risk = standardized_cif(fits, as.data.frame(nd), horizon)
        )
      }))
    }))
  }))
}

binary_differences <- function(long) {
  wide <- dcast(long, horizon_months + age ~ exposure_level, value.var = "risk")
  wide[, risk_difference_eCCA_minus_iCCA := eCCA - iCCA]
  wide
}

schema_differences <- function(long) {
  wide <- dcast(long, horizon_months + age ~ exposure_level, value.var = "risk")
  wide[, `:=`(
    rd_perihilar_minus_iCCA = `Perihilar eCCA` - iCCA,
    rd_distal_minus_iCCA = `Distal eCCA` - iCCA,
    contrast_distal_minus_perihilar = `Distal eCCA` - `Perihilar eCCA`
  )]
  wide
}

binary_data <- droplevels(d_cif[age_topcoded_85plus == FALSE & cause_known == TRUE])
schema_data <- droplevels(binary_data[diagnosis_year >= 2010L & !is.na(schema3)])

binary_point_long <- predict_grid(
  binary_data, fit_competing_models(binary_data, "site_group", cov_binary),
  "site_group", levels(binary_data$site_group), horizons, ages
)
binary_point_rd <- binary_differences(binary_point_long)
schema_point_long <- predict_grid(
  schema_data, fit_competing_models(schema_data, "schema3", cov_schema),
  "schema3", levels(schema_data$schema3), 60, ages
)
schema_point_rd <- schema_differences(schema_point_long)

bootstrap_replicates <- as.integer(Sys.getenv("CCA_HORIZON_BOOTSTRAP_B", unset = "1000"))
bootstrap_seed <- 20261010L
binary_clusters <- split(seq_len(nrow(binary_data)), binary_data$patient_id)
schema_clusters <- split(seq_len(nrow(schema_data)), schema_data$patient_id)

bootstrap_one <- function(iteration) {
  binary_ids <- sample.int(length(binary_clusters), length(binary_clusters), replace = TRUE)
  schema_ids <- sample.int(length(schema_clusters), length(schema_clusters), replace = TRUE)
  db <- droplevels(binary_data[unlist(binary_clusters[binary_ids], use.names = FALSE)])
  ds <- droplevels(schema_data[unlist(schema_clusters[schema_ids], use.names = FALSE)])
  tryCatch({
    binary_long <- predict_grid(
      db, fit_competing_models(db, "site_group", cov_binary),
      "site_group", levels(binary_data$site_group), horizons, ages
    )
    schema_long <- predict_grid(
      ds, fit_competing_models(ds, "schema3", cov_schema),
      "schema3", levels(schema_data$schema3), 60, ages
    )
    list(
      binary_long = binary_long,
      binary_rd = binary_differences(binary_long),
      schema_long = schema_long,
      schema_rd = schema_differences(schema_long)
    )
  }, error = function(e) NULL)
}

if (bootstrap_replicates > 0L) {
  worker_count <- min(10L, max(1L, parallel::detectCores(logical = TRUE) - 2L))
  cl <- parallel::makeCluster(worker_count)
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({library(data.table); library(survival); library(splines)})
    NULL
  })
  parallel::clusterExport(cl, c(
    "binary_data", "schema_data", "binary_clusters", "schema_clusters",
    "fit_competing_models", "predict_grid", "standardized_cif", "cumhaz_at",
    "binary_differences", "schema_differences", "cov_binary", "cov_schema",
    "ns_term", "horizons", "ages"
  ), envir = environment())
  parallel::clusterSetRNGStream(cl, bootstrap_seed)
  boot <- tryCatch(
    parallel::parLapplyLB(cl, seq_len(bootstrap_replicates), bootstrap_one),
    finally = parallel::stopCluster(cl)
  )
  boot <- Filter(Negate(is.null), boot)
} else {
  boot <- list()
}

add_rep <- function(x, i) copy(x)[, replicate := i]
if (length(boot)) {
  boot_binary_long <- rbindlist(lapply(seq_along(boot), function(i) add_rep(boot[[i]]$binary_long, i)))
  boot_binary_rd <- rbindlist(lapply(seq_along(boot), function(i) add_rep(boot[[i]]$binary_rd, i)))
  boot_schema_long <- rbindlist(lapply(seq_along(boot), function(i) add_rep(boot[[i]]$schema_long, i)))
  boot_schema_rd <- rbindlist(lapply(seq_along(boot), function(i) add_rep(boot[[i]]$schema_rd, i)))

  binary_ci <- boot_binary_long[, .(
    lower95 = quantile(risk, 0.025, na.rm = TRUE, type = 6),
    upper95 = quantile(risk, 0.975, na.rm = TRUE, type = 6)
  ), by = .(horizon_months, age, exposure_level)]
  binary_rd_ci <- boot_binary_rd[, .(
    rd_lower95 = quantile(risk_difference_eCCA_minus_iCCA, 0.025, na.rm = TRUE, type = 6),
    rd_upper95 = quantile(risk_difference_eCCA_minus_iCCA, 0.975, na.rm = TRUE, type = 6)
  ), by = .(horizon_months, age)]
  schema_ci <- boot_schema_long[, .(
    lower95 = quantile(risk, 0.025, na.rm = TRUE, type = 6),
    upper95 = quantile(risk, 0.975, na.rm = TRUE, type = 6)
  ), by = .(horizon_months, age, exposure_level)]
  schema_rd_ci <- boot_schema_rd[, .(
    perihilar_lower95 = quantile(rd_perihilar_minus_iCCA, 0.025, na.rm = TRUE, type = 6),
    perihilar_upper95 = quantile(rd_perihilar_minus_iCCA, 0.975, na.rm = TRUE, type = 6),
    distal_lower95 = quantile(rd_distal_minus_iCCA, 0.025, na.rm = TRUE, type = 6),
    distal_upper95 = quantile(rd_distal_minus_iCCA, 0.975, na.rm = TRUE, type = 6),
    contrast_lower95 = quantile(contrast_distal_minus_perihilar, 0.025, na.rm = TRUE, type = 6),
    contrast_upper95 = quantile(contrast_distal_minus_perihilar, 0.975, na.rm = TRUE, type = 6),
    contrast_p = {
      x <- contrast_distal_minus_perihilar
      min(1, 2 * min((sum(x <= 0) + 1) / (length(x) + 1), (sum(x >= 0) + 1) / (length(x) + 1)))
    }
  ), by = .(horizon_months, age)]
} else {
  binary_ci <- binary_point_long[, .(horizon_months, age, exposure_level, lower95 = NA_real_, upper95 = NA_real_)]
  binary_rd_ci <- binary_point_rd[, .(horizon_months, age, rd_lower95 = NA_real_, rd_upper95 = NA_real_)]
  schema_ci <- schema_point_long[, .(horizon_months, age, exposure_level, lower95 = NA_real_, upper95 = NA_real_)]
  schema_rd_ci <- schema_point_rd[, .(
    horizon_months, age, perihilar_lower95 = NA_real_, perihilar_upper95 = NA_real_,
    distal_lower95 = NA_real_, distal_upper95 = NA_real_, contrast_lower95 = NA_real_,
    contrast_upper95 = NA_real_, contrast_p = NA_real_
  )]
}

binary_out <- merge(binary_point_long, binary_ci, by = c("horizon_months", "age", "exposure_level"), all.x = TRUE)
binary_rd_out <- merge(binary_point_rd, binary_rd_ci, by = c("horizon_months", "age"), all.x = TRUE)
schema_out <- merge(schema_point_long, schema_ci, by = c("horizon_months", "age", "exposure_level"), all.x = TRUE)
schema_rd_out <- merge(schema_point_rd, schema_rd_ci, by = c("horizon_months", "age"), all.x = TRUE)

# Preserve the prespecified primary five-year bootstrap intervals so the
# horizon extension reproduces the published primary analysis exactly.
primary_five_year_file <- file.path(tables_dir, "diagnostic_adjusted_cif5_risk_differences.csv")
if (file.exists(primary_five_year_file)) {
  primary_five_year <- fread(primary_five_year_file)[
    target_type == "Fixed age" & as.integer(target) %in% ages,
    .(age = as.integer(target), primary_lower95 = rd_lower95, primary_upper95 = rd_upper95)
  ]
  binary_rd_out[primary_five_year, on = .(age),
                `:=`(rd_lower95 = fifelse(horizon_months == 60, i.primary_lower95, rd_lower95),
                     rd_upper95 = fifelse(horizon_months == 60, i.primary_upper95, rd_upper95))]
}

ext_dir <- file.path(results_dir, "reporting_extensions")
dir.create(ext_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(binary_out, file.path(ext_dir, "horizon_standardized_cif.csv"))
fwrite(binary_rd_out, file.path(ext_dir, "horizon_standardized_risk_differences.csv"))
fwrite(schema_out, file.path(ext_dir, "schema3_standardized_cif5.csv"))
fwrite(schema_rd_out, file.path(ext_dir, "schema3_standardized_risk_differences.csv"))
if (length(boot)) {
  fwrite(boot_binary_rd, file.path(source_dir, "horizon_standardized_rd_bootstrap.csv.gz"))
  fwrite(boot_schema_rd, file.path(source_dir, "schema3_standardized_rd_bootstrap.csv.gz"))
}

qc <- data.table(
  item = c("Binary cohort N", "Schema cohort N", "Bootstrap requested", "Bootstrap successful"),
  value = c(nrow(binary_data), nrow(schema_data), bootstrap_replicates, length(boot))
)
fwrite(qc, file.path(ext_dir, "horizon_anatomy_qc.csv"))

cat("Horizon and anatomy absolute-risk extensions complete.\n")
print(qc)
print(binary_rd_out)
print(schema_rd_out)
