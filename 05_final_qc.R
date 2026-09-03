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
add_check("Legacy full-age interaction reproduces", abs(p_legacy - 1.59661109144431e-08) < 1e-14, format(p_legacy, digits = 16))

prepub <- fread(file.path(results, "tables", "prepub_age85_interaction_tests.csv"))
p_main <- prepub[analysis == "Exact ages 15-84" & model == "Model 1", p]
add_check("Publication-primary exact-age interaction", length(p_main) == 1L && is.finite(p_main) && p_main < 0.05, format(p_main, scientific = TRUE))

era <- fread(file.path(results, "tables", "prepub_era_threeway_tests.csv"))
add_check("Formal era tests complete", nrow(era) == 4L && !anyNA(era[, .(chisq, df, p)]), paste(nrow(era), "rows"))

schema <- fread(file.path(results, "tables", "prepub_schema3_interaction_tests.csv"))
add_check("Schema-3 tests complete", nrow(schema) == 2L && !anyNA(schema[, .(chisq, df, p)]), paste(nrow(schema), "rows"))

time_varying <- fread(file.path(results, "tables", "prepub_time_varying_tests.csv"))
add_check("Time-varying tests complete", nrow(time_varying) == 2L && !anyNA(time_varying[, .(chisq, df, p)]), paste(nrow(time_varying), "rows"))

prepub_qc <- fread(file.path(results, "prepublication_qc.csv"))
boot_requested <- prepub_qc[item == "Bootstrap requested", value]
boot_successful <- prepub_qc[item == "Bootstrap successful", value]
expected_bootstrap <- as.integer(Sys.getenv("CCA_BOOTSTRAP_B", unset = "500"))
add_check("Adjusted-CIF bootstrap complete", boot_requested == expected_bootstrap && boot_successful == expected_bootstrap,
          paste(boot_successful, "of", boot_requested))

cif_adjusted <- fread(file.path(results, "tables", "prepub_adjusted_cif5_risk_differences.csv"))
add_check("Adjusted CIF intervals complete", nrow(cif_adjusted) == 9L &&
          !anyNA(cif_adjusted[, .(risk_difference_eCCA_minus_iCCA, rd_lower95, rd_upper95)]),
          paste(nrow(cif_adjusted), "rows"))

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
  "Figure5_prepublication_validity_checks"
)
fig_ext <- c("svg", "pdf", "tiff", "png")
figure_paths <- as.vector(outer(fig_stems, fig_ext, function(s, e) file.path(results, "figures", paste0(s, ".", e))))
figure_ok <- file.exists(figure_paths) & file.info(figure_paths)$size > 0
add_check("All figure formats exist", all(figure_ok), paste(sum(figure_ok), "of", length(figure_ok)))

scripts <- file.path(workspace, c(
  "01_import_qc.R", "02_models.R", "03_diagnostics_extended.R",
  "04_figures.R", "06_prepublication_analyses.R",
  "07_prepublication_figures.R", "05_final_qc.R"
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
