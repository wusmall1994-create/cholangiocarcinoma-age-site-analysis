options(stringsAsFactors = FALSE)

resolve_path <- function(variable, default) {
  value <- Sys.getenv(variable, unset = "")
  if (!nzchar(value)) value <- default
  normalizePath(value, winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_path("CCA_PROJECT_DIR", getwd())
data_dir <- resolve_path("CCA_DATA_DIR", file.path(project_dir, "data"))
derived_dir <- resolve_path("CCA_DERIVED_DIR", file.path(project_dir, "derived"))
results_dir <- resolve_path("CCA_RESULTS_DIR", file.path(project_dir, "results"))
input_csv <- resolve_path(
  "CCA_INPUT_CSV",
  file.path(data_dir, "CCA_SEER17_2000_2023_raw.csv")
)

tables_dir <- file.path(results_dir, "tables")
models_dir <- file.path(results_dir, "models")
source_dir <- file.path(results_dir, "source_data")
figures_dir <- file.path(results_dir, "figures")
workspace <- project_dir
root <- project_dir
results <- results_dir
derived <- derived_dir

for (path in c(derived_dir, results_dir, tables_dir, models_dir, source_dir, figures_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
