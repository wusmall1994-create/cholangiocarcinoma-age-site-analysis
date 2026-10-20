# Cholangiocarcinoma prognosis across death endpoints and primary sequence

Analysis code and aggregate results for **Age-specific prognosis comparisons between intrahepatic and extrahepatic cholangiocarcinoma depend on death endpoint definitions and primary sequence eligibility**.

The study compares age-specific extrahepatic versus intrahepatic cholangiocarcinoma prognosis under SEER attribution, common biliary and common hepatobiliary death definitions, registry-recorded primary-sequence restrictions, and direct all-cause mortality. Individual-level SEER records, identifiers and fitted patient-level models are not distributed. Data access requires an authorized SEER export and data-use agreement: https://seer.cancer.gov/data/.

## Analysis populations and estimands

| Analysis | Diagnoses | Records / patients | Specification |
| --- | --- | --- | --- |
| Primary age-by-site hazard hypothesis | 2004–2023 | 25,400 / 25,353 | Unpenalized cause-specific Cox interaction test |
| Endpoint and primary-sequence sensitivity | 2004–2018, ages 15–84, known cause classification | 15,578 / 15,561 | Fixed-ridge competing-risk models; shared target within each population |
| Registry-recorded single primary | Same risk eligibility | 12,222 / 12,222 | Stratum-specific empirical target |
| Multiple-primary sequence | Same risk eligibility | 3,356 / 3,339 | Stratum-specific empirical target |
| Direct all-cause mortality | Same risk eligibility including unknown-cause deaths | 15,780 / 15,763 | Direct fixed-ridge all-cause Cox model |
| Anatomical supplementary analysis | 2010–2018, EOD classifiable | 11,361 records | Common EOD target and paired anatomical definitions |

Risk differences (RDs) are eCCA minus iCCA; the age contrast is RD75 minus RD50. Common biliary death comprises COD-to-site categories 21072 and 21090; the hepatobiliary set adds 21071 and 21080. These sets are applied identically across anatomical sites and change the estimand. They are not adjudicated causes of death. Complementary deaths include cancers outside the selected endpoint.

The primary hypothesis is the 2004–2023 age-by-site hazard interaction. Endpoint, sequence, directed-reclassification and anatomical analyses are exploratory. Sequence restrictions change the population and target. The 621-record scenario changes endpoint membership while retaining records, follow-up, covariates and target; it is a classification sensitivity analysis, not bias correction. Direct all-cause mortality avoids cause allocation; intervals spanning zero do not establish equivalence.

## Software and input preparation

Use an existing R installation with `data.table`, `survival`, `splines`, `cmprsk` and the plotting packages referenced by the supplementary scripts. The analysis environment is recorded in `results/endpoint_sensitivity/sessionInfo.txt`. Python figure generation requires matplotlib and Pillow. No installation command is run automatically.

From `scripts/`, set these environment variables to absolute paths:

| Variable | Content |
| --- | --- |
| `CCA_PROJECT_DIR` | Local project directory |
| `CCA_INPUT_CSV` | Authorized SEER 17 November 2025 export, 2000–2023; required columns are explicit in `01_import_qc.R` |
| `CCA_DERIVED_DIR` | Local directory for derived patient-level objects |
| `CCA_RESULTS_DIR` | Writable local output directory, separate from the distributed `results/` directory |
| `CCA_COD_CSV` | Matched export with `COD to site recode` and `COD to site recode ICD-O-3 2023 Revision Expanded (1999+)` |
| `CCA_ANALYSIS_CODE_DIR` | Absolute path to `scripts/` |

Run `Rscript 01_import_qc.R`. This produces `cca_strict_hazard_2004_2023.rds` and `cca_strict_cif5_2004_2018.rds`. Keep the covariates from this input fixed when matching the COD export by patient identifier and record number. The COD export adds endpoint fields; it does not replace the covariate source. Both input exports must have compatible record keys. All patient-level data and generated RDS files remain local.

## Endpoint and sequence analysis

1. Set `CCA_ENDPOINT_MODE=point` and run `Rscript run_unified_ridge.R`.
2. Run `Rscript core_baseline.R` and `python format_core_baseline.py` for the 15,578-record baseline summaries and formatted Table 1. Run `Rscript sequence_support.R` for descriptive first/later-primary event support (Table S39).
3. Set `CCA_ENDPOINT_MODE=bootstrap`, `CCA_ENDPOINT_B=2000`, `CCA_ENDPOINT_WORKERS=1`, and `CCA_ENDPOINT_WORKER=1`, then run `Rscript run_unified_ridge.R`.
4. Run `Rscript summarize_unified_ridge.R`.
5. Run `python audit_bootstrap.py` and `Rscript check_predictions_and_penalties.R`.

For parallel processing, launch N processes with the same `CCA_ENDPOINT_WORKERS=N` and distinct `CCA_ENDPOINT_WORKER=1,…,N`. Each sample uses seed `202610192+b`. Do not change data, code, seeds or worker settings when resuming completed samples. Existing replicate files are reused; select a fresh output directory for a different specification.

Patient clusters are sampled from the expanded 15,780-record cohort. Cause-known eligibility and sequence restrictions are reapplied within each sample. Every endpoint, sequence restriction, directed scenario and direct all-cause analysis shares the patient draw. Endpoint definitions within a population share the same resampled target; sequence strata have separate targets. Joint cause-specific hazard integration produces cancer death, complementary death and survival probabilities.

Ridge theta is fixed at 0.1 on all unscaled design columns, including spline/site interaction columns (`scale=FALSE`, Efron ties, `eps=1e-11`, `iter.max=100`). There is no cross-validation or outcome-based tuning. Alternative theta values 0.01 and 1 are point checks only. The bootstrap refits the fixed specification and computes percentile pointwise 95% intervals (R type-7 quantiles). Paired endpoint and age differences are calculated within each sample; marginal interval endpoints must not be subtracted. All 2,000 distributed samples completed without fitting warnings or failures.

## Distributed results and manuscript mapping

`results/endpoint_sensitivity/` contains the final core estimates and reproducibility summaries:

| Manuscript element | Aggregate source |
| --- | --- |
| Table 1, core risk-cohort baseline | `Table1_source.json`, `baseline_counts.csv`, `age_summary.csv`, `age_counts.csv` |
| Tables 2–3; endpoint and sequence estimates; Tables S34–S36 | `ridge_contrasts_with_CI.csv` |
| Figure 3 | `ridge_contrasts_with_CI.csv`; `scripts/draw_figures.py` |
| Bidirectional classification crosswalk and age composition | `bidirectional_transitions_by_age.csv` |
| 621-record directed scenario | `scenario621_audit.csv` and the Scenario621 rows in the core contrast table |
| Numerical diagnostics | `bootstrap_diagnostics.csv`, `independent_prediction_QA.csv`, `alternative_penalty_point_checks.csv`, `quantile_MCSE.csv`, `half_sample_interval_stability.csv` |
| Complete paired replicate estimates | `aggregate_bootstrap_contrasts.csv.gz` |
| Table S38, broad-cohort baseline | `results/tables/broad_cohort_baseline.csv` |

The replicate matrix has 2,000 samples × 447 contrasts and contains model estimates only. It is gzip-compressed; `audit_bootstrap.py` reads it directly. Local summarization produces an uncompressed CSV. Estimates are probabilities unless the field explicitly says percentage points; multiply by 100 for manuscript risk percentages and RDs. Figure 3 labels the paired common-biliary-minus-SEER age contrast and its own bootstrap interval, as well as the single-primary estimates. `sequence_support.csv` and `S39_source.json` describe first versus later primary status within multiple-primary records. First-primary eCCA has 188 records and 21 SEER-attributed deaths; these support counts do not establish subgroup endpoint sensitivity. No additional subgroup bootstrap analysis is claimed.

To reproduce the main figures from distributed aggregates, set `CCA_RESULTS_DIR` to the absolute `results/endpoint_sensitivity` directory and `CCA_FIGURE_OUTPUT_DIR` to a separate writable output directory, then run:

```text
python scripts/draw_figures.py
```

The three outputs are 600-dpi LZW-compressed TIFF files. Figure 1 cohort counts and Figure 2 classification rules are fixed study annotations; Figure 3 reads the distributed contrast table. Arrows in Figure 2 represent classification paths, not causal mechanisms.

## Supplementary analyses

The numbered scripts implement import, the primary hazard test, spline and model diagnostics, competing-risk summaries, horizons, EOD case definitions, descriptive anatomical decomposition and cohort-specific potential follow-up. From `scripts/`, use `02_models.R`, `03_diagnostics_extended.R`, `06_absolute_risk_analysis.R`, `08_sensitivity_analyses.R`, `09_competing_risk_summaries.R`, `10_nonproportional_cif_sensitivity.R` and `11_horizon_anatomy_absolute_risks.R` for these analyses. `04_figures.R` and `07_absolute_risk_figures.R` support supplementary plotting.

For paired anatomical definitions, use `16_paired_bootstrap.R` (`CCA_BOOT_MODE=preflight`, then `run`; `CCA_PAIRED_B=1000`, `CCA_WORKERS=1`, `CCA_WORKER=1`), followed by `18_audit_paired_bootstrap.R` (`CCA_AUDIT_MODE=final`) and `19_validate_ridge_prediction.R`. The composition and diagnostics scripts 12–15, 17 and 20 support this specification. `results/composition_decomposition/paired_1000/` contains its paired estimates and audits.

For supplementary unpenalized competing outcomes and PH/interval/stage comparisons, set `CCA_CLINICAL_SCRIPTS_DIR` to `scripts/`, `CCA_CLINICAL_OUT` to a separate writable directory, and `CCA_AGGREGATE_RESULTS` to the distributed `results/`. Run `competing_outcomes.R` with `CCA_CLINICAL_MODE=point`, `descriptive`, then `bootstrap` (`CCA_CLINICAL_B=1000`); run `clinical_auxiliary_bootstrap.R`, `original_ridge_check.R`, `conditional_composition_check.R`, `summarize_competing_outcomes.R`, `summarize_original_ridge.R`, `age_risk_contrasts.R`, and the clinical audit scripts. Their outputs are in `results/clinical_extension/`. Filenames containing `original` identify the unpenalized SEER-endpoint comparator, not the main fixed-ridge estimates. Numerical warnings and point-only sensitivity checks are retained as scientific diagnostics.

The separate all-cause eligibility comparison in Table S36 uses a different model/target question and bootstrap stream. Set `CCA_ALLCAUSE_RESULTS_DIR` to a separate output directory; run `allcause_eligibility.R` with `CCA_ALLCAUSE_MODE=point`, then `bootstrap` (`CCA_ALLCAUSE_B=1000`, `CCA_WORKERS=1`, `CCA_WORKER=1`), followed by `summarize_allcause_eligibility.R`. Its distributed contrasts are in `results/allcause_eligibility/`. Do not pool these estimates with the main ridge stream.

Run `21_followup_summary.R` for reverse Kaplan–Meier potential follow-up using full recorded times. Deaths are censored and living/last-contact observations treated as events for this calculation; the five-year prediction horizon does not truncate it. Outputs and method details are in `results/followup/`.

## Integrity and citation

`MANIFEST.json` records SHA-256 checksums. The MIT license applies to code, not to SEER patient data. Cite software release [v1.5.1](https://github.com/wusmall1994-create/cholangiocarcinoma-age-site-analysis/releases/tag/v1.5.1) and its Zenodo archive record. Citation metadata are in `CITATION.cff`.
