required <- c(
  "data.table", "survival", "cmprsk", "ggplot2", "patchwork",
  "scales", "ragg", "svglite"
)

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

message("Required R packages are available.")
