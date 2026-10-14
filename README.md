# Cholangiocarcinoma age and anatomical-site prognosis

Analysis code and non-identifiable aggregate results for **Age-related anatomical-site differences in cholangiocarcinoma prognosis: cancer-death, competing mortality, and all-cause risk**.

The repository presents the analysis specifications used in the manuscript. It contains no raw or derived patient records, identifiers, fitted patient-level model objects, or manuscript drafts. SEER data require registration and a data-use agreement; obtain them from https://seer.cancer.gov/data/.

## Analysis structure

| Analysis | Diagnosis window | Target and estimator | Principal outputs |
| --- | --- | --- | --- |
| Primary age-by-site hypothesis | 2004–2023 | 25,400 records; unpenalized cause-specific Cox model | Relative effects and interaction tests |
| Primary standardized outcomes | 2004–2018 | 15,578 records; own empirical target; joint unpenalized cause-specific models | Cancer death, other death, survival and all-cause death; Table 2, Figure 2, Tables S11, S13 and S24 |
| Harmonized definitions | 2010–2018 | Six fitting definitions; common 11,361-record EOD target; ridge theta=0.1 | Tables 3–4, Figure 3, Tables S19–S21 |
| Paired model and stage checks | 2010–2018 | PH versus interval-specific models; separate targets within each stage | Tables S25–S27 and Figure S4 |

Risk differences are eCCA minus iCCA. All-cause death equals cancer death plus other death within the cause-known cohort. Synthetic mixtures and their decomposition are descriptive and are not causal mediation or an exact decomposition of a separately fitted binary model. Point-only checks are not interval estimates.

## Reproduction

Use an existing R installation (analysis environment: R 4.6.1) and the packages recorded in `results/clinical_extension/R_SESSION_INFO.txt`. Python with matplotlib is needed for publication figures. No dependency installation is run automatically.

From `scripts/`, set absolute paths for `CCA_PROJECT_DIR`, `CCA_INPUT_CSV`, `CCA_DERIVED_DIR` and `CCA_RESULTS_DIR`. The input is the authorized SEER 17 November 2025 submission (2000–2023) export; required fields and transformations are explicit in `01_import_qc.R`. Use a separate, empty local results directory, never the distributed aggregate-results directory.

1. Run `01_import_qc.R`, `02_models.R`, `03_diagnostics_extended.R`, `06_absolute_risk_analysis.R`, `08_sensitivity_analyses.R`, `09_competing_risk_summaries.R`, `10_nonproportional_cif_sensitivity.R`, and `11_horizon_anatomy_absolute_risks.R` for the primary and specified sensitivity analyses. Inspect the scripts' environment options before long bootstrap jobs.
2. For the six-definition analysis, run `16_paired_bootstrap.R` with `CCA_BOOT_MODE=preflight`, then `CCA_BOOT_MODE=run`, `CCA_PAIRED_B=1000`, `CCA_WORKERS=1`, and `CCA_WORKER=1`. Run `18_audit_paired_bootstrap.R` with `CCA_AUDIT_MODE=final` and `19_validate_ridge_prediction.R`.
3. Set `CCA_CLINICAL_SCRIPTS_DIR` to the absolute scripts directory, `CCA_CLINICAL_OUT` to a separate empty output directory, and `CCA_AGGREGATE_RESULTS` to the distributed results directory. Run `competing_outcomes.R` with `CCA_CLINICAL_MODE=point`, then `descriptive`, then `bootstrap` with `CCA_CLINICAL_B=1000`. Run `clinical_auxiliary_bootstrap.R`, `original_ridge_check.R`, `conditional_composition_check.R`, `summarize_competing_outcomes.R`, `summarize_original_ridge.R`, and `age_risk_contrasts.R`.
4. Run `audit_clinical_points.R`, `diagnose_clinical_boot_warnings.R`, and `final_clinical_audit.R`. Numerical warnings are retained and audited, not silently excluded. Scripts 12–15, 17 and 20 support decomposition specification and numerical diagnostics; they do not replace the final paired intervals.

Resume requires identical data, model, code, seed and worker settings. Multiple workers require distinct indices and the same worker count for a stream. Do not redistribute locally generated RDS files.

## Figures from distributed aggregates

From the repository root:

```text
python scripts/figures_cohort_and_definitions.py --results results --out figures
python scripts/figures_competing_outcomes.py --clinical results/clinical_extension --results results --out figures
```

The first command generates Figures 1 and 3; the second generates Figure 2 and Supplementary Figures S3 and S4. Outputs are 600-dpi TIFF. The primary R figure scripts support the other supplementary panels. Lines connect displayed ages, and intervals are pointwise rather than simultaneous.

## Aggregate sources and uncertainty

`results/clinical_extension/original_*_with_CI.csv` supplies the primary competing outcomes; `original_age65_age75_contrasts_with_CI.csv` supplies paired age contrasts. `results/composition_decomposition/paired_1000/estimates_with_ci.csv`, `paired_differences_with_ci.csv` and `noncomposition_ci.csv` supply the harmonized definition analysis. The PH/interval paired comparison uses `interval_minus_PH_paired_CI.csv`; do not subtract marginal CI endpoints.

All reported probabilities and differences in these files are in probability units unless a column explicitly states percentages or percentage points. Multiply by 100 for manuscript risk percentages and percentage-point RDs. `lower95`/`upper95` denote percentile 95% limits; replicate-level files contain aggregate estimates, not patient records.

The six-definition stream uses seed 202610110+b, the PH/interval and stage stream 202610140+b, and the primary competing-outcome/ridge stream 202610150+b, with 1,000 patient-cluster samples each. Different streams can yield slightly different marginal intervals for identical point specifications. All-cause and age contrasts use within-replicate sums/differences. The original unpenalized outcome stream retained 677 warning-bearing samples among 1,000, predominantly sparse zero-event other-death categories. The paired ridge comparison changed evaluated RDs by at most 0.86 percentage points. See `FINAL_NUMERICAL_AUDIT.csv` and the diagnostic aggregates for limits of this check.

## Access and citation

Code uses the existing MIT license. This license does not grant rights to SEER patient data. Cite the repository URL and exact commit used. `MANIFEST.json` provides SHA-256 checksums of the distributed files. The repository and accompanying Supplementary Data 1 contain the same final analysis files. A GitHub update is not a new Zenodo deposit; no DOI for this exact code snapshot is claimed.
