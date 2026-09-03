options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
project_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}
setwd(project_dir)
Sys.setenv(CCA_PROJECT_DIR = project_dir)

scripts <- c(
  "01_import_qc.R",
  "02_models.R",
  "03_diagnostics_extended.R",
  "04_figures.R",
  "06_prepublication_analyses.R",
  "07_prepublication_figures.R",
  "05_final_qc.R"
)

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)
for (script in scripts) {
  message("Running ", script)
  status <- system2(rscript, script)
  if (!identical(status, 0L)) stop(script, " failed with exit status ", status)
}

message("Analysis pipeline completed.")
