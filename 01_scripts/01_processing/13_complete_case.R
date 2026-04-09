# ------------------------------------------------------------------------------
# Script: 13_complete_case.R
# Purpose: Create complete case dataset (no imputation) for sensitivity analysis.
#          Drop all participants with any missing predictor values, then split
#          into selection (50%), train (37.5%), and test (12.5%) sets.
#
# Inputs:  00_data/01_processed/ukb_final_merged.rds
#
# Outputs:
#   - 04_complete_case/cc_selection.rds  (50% — for LASSO stability selection)
#   - 04_complete_case/cc_train.rds      (37.5% — for model fitting)
#   - 04_complete_case/cc_test.rds       (12.5% — for evaluation)
#
# Upstream:  10_final_clean.R
# Downstream: LASSO selection, logistic regression, ML/AI models
#
# Notes:
#   - Admin/outcome columns (eid, recruitment_date, death_date, first_cvd_date,
#     has_cvd) are protected from the complete case filter
#   - All splits are stratified by has_cvd to preserve outcome prevalence
# ------------------------------------------------------------------------------
rm(list = ls())

# Import packages ---------------------------------------------------------

library(dplyr)
library(here)

set.seed(42)

SCENARIO <- "complete_case"
source(here("01_scripts", "00_variable_config.R"))

# Paths -------------------------------------------------------------------

output_dir <- data_input_dir

# Import data -------------------------------------------------------------

df <- readRDS(here("00_data", "01_processed", "ukb_final_cleaned.rds"))

cat("Starting dimensions:", nrow(df), "x", ncol(df), "\n")


# Remove anthropometrics --------------------------------------------------

# Remove excluded variables (anthropometrics) before complete case filter
if (length(exclude_vars) > 0) {
  exclude_in_data <- intersect(exclude_vars, names(df))
  if (length(exclude_in_data) > 0) {
    cat("Removing excluded variables before NA filter:", 
        paste(exclude_in_data, collapse = ", "), "\n")
    df <- df %>% select(-any_of(exclude_in_data))
  }
}

# Complete case filtering -------------------------------------------------

# Admin/outcome columns — protect from NA check
protected_cols <- admin_vars

# Predictor columns — these must be complete
predictor_cols <- setdiff(names(df), protected_cols)

# Count NAs per row across predictor columns only
row_na_count <- rowSums(is.na(df[, predictor_cols]))

cat("\nMissingness distribution across predictor columns:\n")
cat("  Participants with 0 NAs:", sum(row_na_count == 0), "\n")
cat("  Participants with 1+ NAs:", sum(row_na_count > 0), "\n")
cat("  Median NAs per participant:", median(row_na_count), "\n")
cat("  Max NAs per participant:", max(row_na_count), "\n")

# Keep only complete cases
df_cc <- df[row_na_count == 0, ]

cat("\nAfter complete case filter:", nrow(df_cc), "of", nrow(df),
    "participants retained (",
    round(nrow(df_cc) / nrow(df) * 100, 1), "%)\n")

cat("CVD prevalence:", round(mean(df_cc$has_cvd) * 100, 2), "%\n")

# Stratified splitting function -------------------------------------------

stratified_split <- function(data, prop, strat_var = "has_cvd") {
  idx_cases <- which(data[[strat_var]] == 1)
  idx_controls <- which(data[[strat_var]] == 0)
  
  n_cases <- round(length(idx_cases) * prop)
  n_controls <- round(length(idx_controls) * prop)
  
  sel_cases <- sample(idx_cases, n_cases)
  sel_controls <- sample(idx_controls, n_controls)
  
  split_idx <- sort(c(sel_cases, sel_controls))
  
  list(
    set1 = data[split_idx, ],
    set2 = data[-split_idx, ]
  )
}

# 50/50 split: selection vs inference -------------------------------------

split1 <- stratified_split(df_cc, prop = 0.5)
cc_selection <- split1$set1
cc_inference <- split1$set2

cat("\n--- 50/50 Split ---\n")
cat("Selection set:", nrow(cc_selection), "| CVD prevalence:",
    round(mean(cc_selection$has_cvd) * 100, 2), "%\n")
cat("Inference set:", nrow(cc_inference), "| CVD prevalence:",
    round(mean(cc_inference$has_cvd) * 100, 2), "%\n")

# 75/25 split on inference: train vs test ---------------------------------

split2 <- stratified_split(cc_inference, prop = 0.75)
cc_train <- split2$set1
cc_test  <- split2$set2

cat("\n--- 75/25 Split (Inference) ---\n")
cat("Train set:", nrow(cc_train), "| CVD prevalence:",
    round(mean(cc_train$has_cvd) * 100, 2), "%\n")
cat("Test set:", nrow(cc_test), "| CVD prevalence:",
    round(mean(cc_test$has_cvd) * 100, 2), "%\n")

# Export ------------------------------------------------------------------

saveRDS(cc_selection, file.path(output_dir, selection_file))
saveRDS(cc_train, file.path(output_dir, train_file))
saveRDS(cc_test, file.path(output_dir, test_file))

cat("\nFiles saved to:", output_dir, "\n")
