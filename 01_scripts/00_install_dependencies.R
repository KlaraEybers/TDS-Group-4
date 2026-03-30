# ==============================================================================
# Install Required R Packages
# Run once: Rscript install_dependencies.R
# ==============================================================================

packages <- c(
  # Core tidyverse
  "dplyr", "ggplot2", "tidyr", "purrr", "stringr", "tibble", "tidyverse",
  "readr",
  # Project management
  "here",
  # Data import
  "readxl", "data.table",
  # Imputation
  "miceRanger", "future",
  # Plotting
  "patchwork", "corrplot",
  # Statistics
  "e1071", "FactoMineR",
  # Variable selection
  "sharp",
  # Analysis
  "boot", "pROC", "broom", "scales"
)

# Only install missing packages
missing <- packages[!packages %in% installed.packages()[, "Package"]]

if (length(missing) > 0) {
  cat("Installing", length(missing), "missing packages...\n")
  install.packages(missing, repos = "https://cloud.r-project.org", Ncpus = 4)
} else {
  cat("All packages already installed.\n")
}