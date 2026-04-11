# ==============================================================================
# Script: Extract Variable Names for 3 Scenarios
# ==============================================================================
rm(list = ls())

library(here)

cat("=======================================================\n")
cat("Extracting variable names for 3 scenarios...\n")
cat("=======================================================\n")

# 1. Define paths 
main_path   <- here("00_data", "01_processed", "03_imputation", "df_test_imputed.rds")
anthro_path <- here("00_data", "01_processed", "03_imputation", "df_test_imputed.rds") 
cc_path     <- here("00_data", "01_processed", "04_complete_case", "cc_test.rds")

# Read datasets
main_df   <- readRDS(main_path)
anthro_df <- readRDS(anthro_path)
cc_df     <- readRDS(cc_path)

# 2. Extract variable names (column names)
main_vars   <- names(main_df)
anthro_vars <- names(anthro_df)
cc_vars     <- names(cc_df)

# 3. Standardize lengths: Pad with NAs so they can be merged into a single dataframe
max_len <- max(length(main_vars), length(anthro_vars), length(cc_vars))

length(main_vars)   <- max_len
length(anthro_vars) <- max_len
length(cc_vars)     <- max_len

# 4. Combine into a single data frame with 3 columns
scenario_vars_df <- data.frame(
  Main_Scenario = main_vars,
  Anthropometric_Scenario = anthro_vars,
  Complete_Case_Scenario = cc_vars
)

# 5. Save to the specific results folder
results_dir <- here("02_resultsv2")

output_file <- file.path(results_dir, "scenario_variable_names.csv")

# Export to CSV. 'na = ""' ensures padded cells appear clean and empty in Excel/CSV
write.csv(scenario_vars_df, output_file, row.names = FALSE, na = "")

cat("SUCCESS! CSV file has been saved to:\n")
cat("->", output_file, "\n")
cat("=======================================================\n")