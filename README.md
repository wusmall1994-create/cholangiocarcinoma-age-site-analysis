# Age-specific anatomical-site contrasts in cholangiocarcinoma

This repository contains the R code used for a population-based analysis of whether age at diagnosis modifies cancer-specific mortality contrasts between intrahepatic and extrahepatic cholangiocarcinoma.

## Repository scope

Only statistical analysis code is included. The repository does not contain SEER records, derived patient-level data, model objects, numerical results, figures, tables, manuscripts, submission files, credentials, or personal information.

## Data access

The analysis uses the Incidence - SEER Research Data, 17 Registries, November 2025 submission (2000-2023), linked to time-dependent county income and rurality attributes. SEER data are subject to registration and a data-use agreement and cannot be redistributed by this repository. Researchers must obtain the corresponding database directly from the US National Cancer Institute SEER Program and create a case-listing export with the variables required by `01_import_qc.R`.

The default input filename is `CCA_SEER17_2000_2023_raw.csv`. Store it outside version control and provide its path with the `CCA_INPUT_CSV` environment variable.

## Software

- R 4.6.1
- data.table
- survival 3.8-6
- cmprsk 2.2-12
- ggplot2
- patchwork
- scales
- ragg
- svglite

Install missing packages with:

```r
Rscript install_dependencies.R
```

## Configuration

Run the code from the repository root. Paths can be set through environment variables:

- `CCA_INPUT_CSV`: SEER case-listing CSV
- `CCA_DERIVED_DIR`: private directory for derived patient-level files
- `CCA_RESULTS_DIR`: private directory for tables, models, and figures
- `CCA_BOOTSTRAP_B`: number of bootstrap replicates; default `500`

Example in PowerShell:

```powershell
$env:CCA_INPUT_CSV = "D:/restricted-data/CCA_SEER17_2000_2023_raw.csv"
$env:CCA_DERIVED_DIR = "D:/restricted-data/derived"
$env:CCA_RESULTS_DIR = "D:/restricted-data/results"
Rscript run_all.R
```

## Analysis sequence

1. `01_import_qc.R`: imports the SEER export, harmonizes variables, defines cohorts, and performs selection checks.
2. `02_models.R`: fits primary cause-specific, overall-survival, and competing-risk models.
3. `03_diagnostics_extended.R`: produces continuous curves, proportional-hazards diagnostics, and complementary absolute-risk estimates.
4. `04_figures.R`: generates the main figures from analysis outputs.
5. `06_prepublication_analyses.R`: evaluates age top-coding, diagnosis-era interaction, detailed anatomy, adjusted cumulative incidence, and time-varying effects.
6. `07_prepublication_figures.R`: generates figures for the extended analyses.
7. `05_final_qc.R`: checks expected cohort counts, statistical outputs, figures, and script syntax.

`run_all.R` executes these scripts in dependency order. All generated files remain in the private directories specified above and are excluded by `.gitignore`.

## Reproducibility note

The default analysis uses 500 patient-cluster bootstrap replicates with seed 20260902. A small value of `CCA_BOOTSTRAP_B` may be used only for a local smoke test; publication estimates require the default value.

## Licence

The analysis code is released under the MIT License. The licence does not apply to SEER data or override SEER data-use conditions.
