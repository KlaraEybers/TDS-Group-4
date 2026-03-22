# ------------------------------------------------------------------------------
# Script: run_processing_pipeline.R
# Purpose: Run the full data processing pipeline in the correct order.
# ------------------------------------------------------------------------------

rm(list = ls())

library(here)

scripts <- c(
  "00_initial_clean.R",
  "01_socio_demographics.R",
  "02_physical_measures.R",
  "03_diet.R",
  "04_biological_samples.R",
  "05_health_and_medical_history.R",
  "06_early_life_factors.R",
  "07_environment.R",
  "08_lifestyle.R",
  "09_mental_health.R",
  "10_final_clean.R",
  "11_remove_missing.R",
  "12_imputation.R"
)

cat("\nStarting processing pipeline...\n")

for (i in seq_along(scripts)) {
  script_name <- scripts[i]
  script_path <- here("01_scripts", "01_processing", script_name)
  
  cat("\n[", i, "/", length(scripts), "] Running: ", script_name, "\n", sep = "")
  
  tryCatch(
    {
      source(script_path, echo = FALSE)
      cat("Completed: ", script_name, "\n", sep = "")
    },
    error = function(e) {
      cat("Error in: ", script_name, "\n", sep = "")
      stop(e)
    }
  )
}

cat("\nProcessing pipeline completed successfully.\n")