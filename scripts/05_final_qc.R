options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(data.table))

source("00_config.R")

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name, status = if (isTRUE(passed)) "PASS" else "FAIL", detail = as.character(detail)
  )
}

raw <- input_csv
add_check("Raw SEER export exists", file.exists(raw) && file.info(raw)$size > 0, if (file.exists(raw)) file.info(raw)$size else "missing")

d <- readRDS(file.path(derived, "cca_strict_hazard_2004_2023.rds"))
dc <- readRDS(file.path(derived, "cca_strict_cif5_2004_2018.rds"))
add_check("Primary cohort N", nrow(d) == 28654, nrow(d))
add_check("Five-year cohort N", nrow(dc) == 17751, nrow(dc))
add_check("Cause-known cohort N", sum(!is.na(d$competing_event)) == 28323, sum(!is.na(d$competing_event)))

interactions <- fread(file.path(results, "tables", "interaction_tests.csv"))
p_legacy <- interactions[outcome == "Cancer-specific hazard" & model == "Model 2", p]
add_check("Legacy full-age interaction available", length(p_legacy) == 1L && is.finite(p_legacy), format(p_legacy, scientific = TRUE))
add_check("Legacy full-age interaction reproduces", abs(p_legacy - 1.17137293007335e-08) < 1e-14, format(p_legacy, digits = 16))

prepub <- fread(file.path(results, "tables", "diagnostic_age85_interaction_tests.csv"))
p_main <- prepub[analysis == "Exact ages 15-84" & model == "Model 1", p]
add_check("Publication-primary exact-age interaction", length(p_main) == 1L && is.finite(p_main) && p_main < 0.05, format(p_main, scientific = TRUE))

era <- fread(file.path(results, "tables", "diagnostic_era_threeway_tests.csv"))
add_check("Formal era tests complete", nrow(era) == 4L && !anyNA(era[, .(chisq, df, p)]), paste(nrow(era), "rows"))

schema <- fread(file.path(results, "tables", "diagnostic_schema3_interaction_tests.csv"))
add_check("Schema-3 tests complete", nrow(schema) == 2L && !anyNA(schema[, .(chisq, df, p)]), paste(nrow(schema), "rows"))

nonlinearity <- fread(file.path(results, "tables", "sensitivity_nonlinearity_binary.csv"))
p_nonlinear <- nonlinearity[test == "Nonlinear interaction component" & model == "Model 1", p]
add_check("Primary nonlinear interaction component", length(p_nonlinear) == 1L && is.finite(p_nonlinear) && p_nonlinear < 0.05,
          format(p_nonlinear, scientific = TRUE))

full_period <- fread(file.path(results, "tables", "sensitivity_fullperiod_tests.csv"))
p_full_nonlinear <- full_period[test == "Nonlinear interaction component" & model == "Model 1", p]
add_check("Expanded-period nonlinear interaction", length(p_full_nonlinear) == 1L && is.finite(p_full_nonlinear) && p_full_nonlinear < 0.05,
          format(p_full_nonlinear, scientific = TRUE))

time_varying <- fread(file.path(results, "tables", "diagnostic_time_varying_tests.csv"))
add_check("Time-varying tests complete", nrow(time_varying) == 2L && !anyNA(time_varying[, .(chisq, df, p)]), paste(nrow(time_varying), "rows"))

diagnostic_qc <- fread(file.path(results, "absolute_risk_qc.csv"))
boot_requested <- diagnostic_qc[item == "Bootstrap requested", value]
boot_successful <- diagnostic_qc[item == "Bootstrap successful", value]
expected_bootstrap <- as.integer(Sys.getenv("CCA_BOOTSTRAP_B", unset = "1000"))
add_check("Adjusted-CIF bootstrap complete", boot_requested == expected_bootstrap && boot_successful == expected_bootstrap,
          paste(boot_successful, "of", boot_requested))

sensitivity_qc <- fread(file.path(results, "sensitivity_qc.csv"))
sensitivity_requested <- sensitivity_qc[check == "CIF bootstrap requested", as.integer(value)]
sensitivity_successful <- sensitivity_qc[check == "CIF bootstrap successful", as.integer(value)]
expected_sensitivity_bootstrap <- as.integer(Sys.getenv("CCA_ENHANCEMENT_BOOTSTRAP_B", unset = "1000"))
add_check("Constant-effect comparison bootstrap complete",
          sensitivity_requested == expected_sensitivity_bootstrap && sensitivity_successful == expected_sensitivity_bootstrap,
          paste(sensitivity_successful, "of", sensitivity_requested))

cif_adjusted <- fread(file.path(results, "tables", "diagnostic_adjusted_cif5_risk_differences.csv"))
add_check("Adjusted CIF intervals complete", nrow(cif_adjusted) == 9L &&
          !anyNA(cif_adjusted[, .(risk_difference_eCCA_minus_iCCA, rd_lower95, rd_upper95)]),
          paste(nrow(cif_adjusted), "rows"))

horizon_qc <- fread(file.path(results, "reporting_extensions", "horizon_anatomy_qc.csv"))
horizon_requested <- horizon_qc[item == "Bootstrap requested", as.integer(value)]
horizon_successful <- horizon_qc[item == "Bootstrap successful", as.integer(value)]
expected_horizon_bootstrap <- as.integer(Sys.getenv("CCA_HORIZON_BOOTSTRAP_B", unset = "1000"))
add_check("Horizon and anatomy bootstrap complete",
          horizon_requested == expected_horizon_bootstrap && horizon_successful == expected_horizon_bootstrap,
          paste(horizon_successful, "of", horizon_requested))

horizon_rd <- fread(file.path(results, "reporting_extensions", "horizon_standardized_risk_differences.csv"))
schema_rd <- fread(file.path(results, "reporting_extensions", "schema3_standardized_risk_differences.csv"))
add_check("One-, three-, and five-year risk differences complete",
          nrow(horizon_rd) == 12L && !anyNA(horizon_rd[, .(risk_difference_eCCA_minus_iCCA, rd_lower95, rd_upper95)]),
          paste(nrow(horizon_rd), "rows"))
add_check("Detailed-anatomy absolute-risk contrasts complete",
          nrow(schema_rd) == 4L && !anyNA(schema_rd[, .(rd_perihilar_minus_iCCA, rd_distal_minus_iCCA,
                                                        contrast_distal_minus_perihilar, contrast_lower95,
                                                        contrast_upper95, contrast_p)]),
          paste(nrow(schema_rd), "rows"))

pw <- fread(file.path(results, "tables", "piecewise_followup_effects.csv"))
add_check("Piecewise estimates complete", nrow(pw) == 18L && !anyNA(pw[, .(estimate, lower95, upper95)]), paste(nrow(pw), "rows"))

source_files <- list.files(file.path(results, "source_data"), pattern = "\\.csv$", full.names = TRUE)
source_ok <- vapply(source_files, function(f) {
  x <- tryCatch(fread(f), error = function(e) NULL)
  !is.null(x) && nrow(x) > 0L && ncol(x) > 0L
}, logical(1))
add_check("All figure source-data files nonempty", all(source_ok), paste(sum(source_ok), "of", length(source_ok)))

fig_stems <- c(
  "Figure1_cohort_and_events", "Figure2_age_relative_effects",
  "Figure3_five_year_absolute_risk", "Figure4_sensitivity_analyses",
  "Figure5_absolute_risk_validity_checks"
)
fig_ext <- c("tiff")
figure_paths <- as.vector(outer(fig_stems, fig_ext, function(s, e) file.path(results, "figures", paste0(s, ".", e))))
figure_ok <- file.exists(figure_paths) & file.info(figure_paths)$size > 0
add_check("All figure formats exist", all(figure_ok), paste(sum(figure_ok), "of", length(figure_ok)))

scripts <- file.path(workspace, c(
  "01_import_qc.R", "02_models.R", "03_diagnostics_extended.R",
  "04_figures.R", "06_absolute_risk_analysis.R",
  "08_sensitivity_analyses.R", "09_competing_risk_summaries.R",
  "10_nonproportional_cif_sensitivity.R", "11_horizon_anatomy_absolute_risks.R",
  "07_absolute_risk_figures.R", "05_final_qc.R"
))
parse_ok <- vapply(scripts, function(f) {
  tryCatch({ parse(file = f); TRUE }, error = function(e) FALSE)
}, logical(1))
add_check("All R scripts parse", all(parse_ok), paste(sum(parse_ok), "of", length(parse_ok)))

out <- rbindlist(checks)
fwrite(out, file.path(results, "final_qc_report.csv"))
writeLines(c(
  "FINAL QUALITY-CONTROL REPORT",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  paste(out$status, out$check, "-", out$detail),
  "",
  paste("Overall:", if (all(out$status == "PASS")) "PASS" else "FAIL")
), file.path(results, "final_qc_report.txt"))

print(out)
if (!all(out$status == "PASS")) quit(status = 1L)
