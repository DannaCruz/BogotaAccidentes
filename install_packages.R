required <- c("sf", "sfnetworks", "igraph", "tidygraph", "dplyr", "ggplot2",
              "units", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages are installed.")
}
