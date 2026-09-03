options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

source("00_config.R")

required_columns <- c(
  "Patient ID",
  "Record number recode",
  "Year of diagnosis",
  "Age recode with single ages and 85+",
  "Primary Site - labeled",
  "Histologic Type ICD-O-3",
  "Type of Reporting Source",
  "Survival months",
  "Vital status recode (study cutoff used)",
  "SEER cause-specific death classification",
  "SEER other cause of death classification",
  "Combined Summary Stage with Expanded Regional Codes (2004+)",
  "Grade Recode (thru 2017)",
  "Derived Summary Grade 2018 (2018+)",
  "RX Summ--Surg Prim Site (1998-2022)",
  "RX Summ--Surg Prim Site 2023 (2023+)"
)

raw <- fread(
  input_csv,
  colClasses = "character",
  na.strings = "NA",
  encoding = "UTF-8",
  quote = "\"",
  check.names = FALSE
)

missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

stopifnot(ncol(raw) == 42L)
stopifnot(nrow(raw) == 45243L)

raw[, diagnosis_year := as.integer(`Year of diagnosis`)]
raw[, age_years := as.integer(sub("^([0-9]+).*", "\\1", `Age recode with single ages and 85+`))]
raw[, age_topcoded_85plus := grepl("85\\+", `Age recode with single ages and 85+`)]
raw[, age_group := factor(
  fcase(
    age_years <= 39L, "15-39",
    age_years <= 64L, "40-64",
    default = "65+"
  ),
  levels = c("15-39", "40-64", "65+")
)]

raw[, site_group := factor(
  fcase(
    startsWith(`Primary Site - labeled`, "C22.1"), "iCCA",
    startsWith(`Primary Site - labeled`, "C24.0"), "eCCA",
    default = NA_character_
  ),
  levels = c("iCCA", "eCCA")
)]

raw[, strict_histology := `Histologic Type ICD-O-3` %chin% c("8160", "8162")]
raw[, discordant_c221_8162 := site_group == "iCCA" & `Histologic Type ICD-O-3` == "8162"]

raw[, survival_months_num := as.integer(`Survival months`)]
raw[, survival_time_months := fifelse(survival_months_num == 0L, 0.5, as.numeric(survival_months_num))]
raw[, os_event := as.integer(`Vital status recode (study cutoff used)` == "Dead")]
raw[, competing_event := fcase(
  `SEER cause-specific death classification` == "Dead (attributable to this cancer dx)", 1L,
  `SEER other cause of death classification` == "Dead (attributable to causes other than this cancer dx)", 2L,
  `Vital status recode (study cutoff used)` == "Alive", 0L,
  default = NA_integer_
)]

raw[, stage4 := factor(
  fcase(
    `Combined Summary Stage with Expanded Regional Codes (2004+)` == "Localized only", "Localized",
    grepl("^Regional", `Combined Summary Stage with Expanded Regional Codes (2004+)`), "Regional",
    `Combined Summary Stage with Expanded Regional Codes (2004+)` == "Distant site(s)/node(s) involved", "Distant",
    default = "Unknown"
  ),
  levels = c("Localized", "Regional", "Distant", "Unknown")
)]

raw[, grade_source := fifelse(
  diagnosis_year <= 2017L,
  `Grade Recode (thru 2017)`,
  `Derived Summary Grade 2018 (2018+)`
)]
raw[, grade4 := factor(
  fcase(
    grepl("Well differentiated|Grade I([^I]|$)", grade_source, ignore.case = TRUE), "I",
    grepl("Moderately differentiated|Grade II([^I]|$)", grade_source, ignore.case = TRUE), "II",
    grepl("Poorly differentiated|Grade III([^I]|$)", grade_source, ignore.case = TRUE), "III",
    grepl("Undifferentiated|anaplastic|Grade IV", grade_source, ignore.case = TRUE), "IV",
    default = "Unknown"
  ),
  levels = c("I", "II", "III", "IV", "Unknown")
)]

raw[, surgery_primary_raw := fifelse(
  diagnosis_year <= 2022L,
  `RX Summ--Surg Prim Site (1998-2022)`,
  `RX Summ--Surg Prim Site 2023 (2023+)`
)]
raw[, surgery_primary := factor(
  fcase(
    surgery_primary_raw %chin% c("00", "A000"), "No",
    is.na(surgery_primary_raw) | surgery_primary_raw %chin% c("99", "A990"), "Unknown",
    default = "Yes"
  ),
  levels = c("No", "Yes", "Unknown")
)]

if (anyNA(raw$site_group)) stop("Unexpected primary site outside C22.1/C24.0")
if (min(raw$age_years, na.rm = TRUE) < 15L) stop("Age selection check failed")
if (range(raw$diagnosis_year, na.rm = TRUE)[1L] != 2000L ||
    range(raw$diagnosis_year, na.rm = TRUE)[2L] != 2023L) {
  stop("Diagnosis-year selection check failed")
}
if (any(grepl("Autopsy|Death certificate", raw$`Type of Reporting Source`, ignore.case = TRUE))) {
  stop("DCO/autopsy exclusion check failed")
}

strict_2000_2023 <- raw[strict_histology == TRUE]
strict_2004_2023 <- raw[strict_histology & diagnosis_year >= 2004L]
strict_2004_2018 <- raw[strict_histology & diagnosis_year >= 2004L & diagnosis_year <= 2018L]

saveRDS(raw, file.path(derived_dir, "cca_broad_2000_2023.rds"), compress = "xz")
saveRDS(strict_2000_2023, file.path(derived_dir, "cca_strict_2000_2023.rds"), compress = "xz")
saveRDS(strict_2004_2023, file.path(derived_dir, "cca_strict_hazard_2004_2023.rds"), compress = "xz")
saveRDS(strict_2004_2018, file.path(derived_dir, "cca_strict_cif5_2004_2018.rds"), compress = "xz")

qc_lines <- c(
  "SEER cholangiocarcinoma import/QC",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Input: ", input_csv),
  "",
  paste0("Broad cohort 2000-2023: ", nrow(raw)),
  paste0("Unique patients: ", uniqueN(raw$`Patient ID`)),
  paste0("Strict histology 8160/8162, 2000-2023: ", nrow(strict_2000_2023)),
  paste0("Strict hazard cohort, 2004-2023: ", nrow(strict_2004_2023)),
  paste0("Strict 5-year CIF cohort, 2004-2018: ", nrow(strict_2004_2018)),
  paste0("Discordant C22.1 + 8162 records: ", sum(raw$discordant_c221_8162)),
  paste0("Unknown cause of death in strict hazard cohort: ", sum(is.na(strict_2004_2023$competing_event))),
  paste0("Zero survival months in strict hazard cohort: ", sum(strict_2004_2023$survival_months_num == 0L)),
  "",
  "Primary exposure follows topography: C22.1=iCCA, C24.0=eCCA.",
  "C22.1 + 8162 records are retained and flagged; exclude them in a prespecified sensitivity analysis.",
  "For survival modeling, zero-month survival is represented as 0.5 months.",
  "Deaths with unknown cause are retained for OS and excluded from cause-specific/competing-risk analyses."
)
writeLines(qc_lines, file.path(derived_dir, "qc_summary.txt"), useBytes = TRUE)

cat(paste(qc_lines, collapse = "\n"), "\n")
