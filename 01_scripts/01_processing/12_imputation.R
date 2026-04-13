# ------------------------------------------------------------------------------
# Script: 12_imputation.R
# Purpose: Run imputation pipelines for full, selection, training, and test datasets.
#          Scenario-aware: controlled by SCENARIO variable via 00_variable_config.R
#
# Inputs:
#   - 00_data/01_processed/ukb_final_cleaned.rds
#   - 01_scripts/00_variable_config.R
#   - 01_scripts/01_processing/12_imputation_functions.R
#   - 01_scripts/00_variable_dictionary.R
#
# Outputs (paths set by config, vary by scenario):
#   Main:
#     - 00_data/01_processed/03_imputation/df_full_imputed.rds
#     - 00_data/01_processed/03_imputation/df_selection_imputed.rds
#     - 00_data/01_processed/03_imputation/df_train_imputed.rds
#     - 00_data/01_processed/03_imputation/df_test_imputed.rds
#     - 02_results/00_figures/01_imputation/<numeric|categorical>_diagnostics_<label>.pdf
#   Sensitivity (anthro as confounders):
#     - 00_data/01_processed/03_imputation/sensitivity_anthro_confounder/df_*_imputed.rds
#     - 02_results/00_figures/08_sensitivity_anthro_confounder/imputation/*_diagnostics_<label>.pdf
#
# Upstream:
#   - 11_remove_missing.R
#   - 12_imputation_functions.R
#
# Downstream:
#   - 02_analysis/01_descriptive/
#   - 02_analysis/02_univariate/
#   - 02_analysis/03_stability_selection/
#   - 02_analysis/04_multivariate/
#   - 02_analysis/05_neural_network/
#   - 02_analysis/06_decision_trees/
#
# Key steps:
#   - Load variable classifications from 00_variable_config.R (scenario-aware)
#   - Validate config variables against actual data columns
#   - Remove participants with missing confounders (age, sex, ethnicity in main;
#     also anthropometrics in sensitivity)
#   - Exclude variables based on scenario (anthropometrics dropped in main)
#   - Identify skewed biomarkers and specify variables for log transformation
#   - Split data into selection, training, and test sets
#   - Run full-dataset imputation for descriptive and univariate analyses
#   - Run DAG-respecting staged imputation for selection, train, and test datasets
#   - Save imputed datasets and diagnostic plots to scenario-specific directories
#
# Notes:
#   - Scenario set via: SCENARIO env var, or hardcode before source()
#   - TEST_MODE allows reduced-size runs for debugging
#   - Confounders are never imputed — rows with missing confounders are dropped
#   - Log-transformed variables (crp, ggt, bilirubin, lpa) are back-transformed
#     after imputation
#   - Test-set imputation uses models trained on the training set
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

rm(list = ls())

# Import packages ---------------------------------------------------------

library(dplyr)
library(here)
library(miceRanger)
library(future)
library(ggplot2)
library(tidyr)
library(patchwork)
library(e1071)

# Import data and functions  ----------------------------------------------

df <- readRDS(here("00_data", "01_processed", "ukb_final_cleaned.rds"))
source(here("01_scripts", "01_processing", "12_imputation_functions.R"))
source(here("01_scripts", "00_variable_dictionary.R"))

# Toggle for test run or job script run -----------------------------------

TEST_MODE <- TRUE  # Set to FALSE for full run
TEST_N <- 10000


# Classify variables ------------------------------------------------------
# Classifications: admin/outcomes, confounders, exposures, biological

source(here("01_scripts", "00_variable_config.R"))

# Removed: "diabetes_doctor_diagnosed", "serious_dx",
#"age_high_bp_dx", "age_diabetes_dx", "angina_dx_age", "mi_dx_age",
# "dvt_dx_age", "pe_dx_age", "stroke_dx_age",
# "heart_problems_collapsed", "meds_chol_bp_diabetes_hormone_collapsed",
# "meds_chol_bp_diabetes_collapsed", "long_ill", "n_previous_pregnancies_collapsed"


# Validate variable lists against actual data columns ---------------------
# To deal with columns removed due to high missingness

available_cols <- names(df)

exposure_vars_full <- exposure_vars
bio_vars_full      <- bio_vars

exposure_vars <- intersect(exposure_vars, available_cols)
bio_vars      <- intersect(bio_vars, available_cols)

# Report any variables defined in config but missing from data
dropped_exp <- setdiff(exposure_vars_full, exposure_vars)
dropped_bio <- setdiff(bio_vars_full, bio_vars)

if (length(dropped_exp) > 0) {
  cat("WARNING: Config exposures not in data (removed by preprocessing?):\n")
  cat(" ", paste(dropped_exp, collapse = ", "), "\n")
}
if (length(dropped_bio) > 0) {
  cat("WARNING: Config biomarkers not in data (removed by preprocessing?):\n")
  cat(" ", paste(dropped_bio, collapse = ", "), "\n")
}

# Drop excluded variables (anthro in main analysis) -----------------------

if (length(exclude_vars) > 0) {
  cat("Dropping excluded variables:", paste(exclude_vars, collapse = ", "), "\n")
  df <- df %>% select(-any_of(exclude_vars))
}

# Remove rows with missing confounders ------------------------------------
# Confounders shouldn't be imputed — drop rows with any NA

cat("NAs per confounder before drop:\n")
# print(colSums(is.na(df[, confounder_vars])))

n_before <- nrow(df)
df <- df %>% drop_na(any_of(confounder_vars))
n_dropped <- n_before - nrow(df)

cat("Dropped", n_dropped, "rows with missing confounders\n")
cat("N:", nrow(df), "\n")
cat("CVD events:", sum(df$has_cvd == 1), "\n")
cat("No CVD events:", sum(df$has_cvd == 0), "\n")


# Investigate  skewed variables -------------------------------------------

# Skewness test
df %>%
  select(any_of(bio_vars)) %>%
  select(where(is.numeric)) %>%
  sapply(function(x) skewness(x, na.rm = TRUE)) %>%
  sort(decreasing = TRUE)

# Log transform significantly skew variables
log_vars <- c("crp", "ggt", "bilirubin", "lpa")

# Create data splits ------------------------------------------------------
set.seed(123)

stratified_split <- function(data, prop, strat_var = "has_cvd") {
  idx_cases <- which(data[[strat_var]] == 1)
  idx_controls <- which(data[[strat_var]] == 0)
  
  n_cases <- round(length(idx_cases) * prop)
  n_controls <- round(length(idx_controls) * prop)
  
  sel_cases <- sample(idx_cases, n_cases)
  sel_controls <- sample(idx_controls, n_controls)
  
  split_idx <- sort(c(sel_cases, sel_controls))
  
  list(set1 = data[split_idx, ], set2 = data[-split_idx, ])
}

# 50/50 split: selection vs inference
split1 <- stratified_split(df, prop = 0.5)

df_selection <- split1$set1
df_inference <- split1$set2

# 75/25 split on inference: train vs test
split2 <- stratified_split(df_inference, prop = 0.75)
df_train <- split2$set1
df_test  <- split2$set2

cat("Selection:", nrow(df_selection), "| CVD:", round(mean(df_selection$has_cvd) * 100, 2), "%\n")
cat("Train:", nrow(df_train), "| CVD:", round(mean(df_train$has_cvd) * 100, 2), "%\n")
cat("Test:", nrow(df_test), "| CVD:", round(mean(df_test$has_cvd) * 100, 2), "%\n")

############################################################################
# RUN IMPUTATION PIPELINES
############################################################################

# Imputation 1 Pipeline ---------------------------------------------------

# Run imputation on the full dataset for univariate and descriptive analysis
cat("\n========== IMPUTATION 1: Full dataset ==========\n")

all_vars <- c(confounder_vars, exposure_vars, bio_vars)
df_full_impute <- prep_for_imputation(df, vars_to_impute = all_vars, predictor_vars = c())

num_vars <- names(df_full_impute)[sapply(df_full_impute, is.numeric) & colSums(is.na(df_full_impute)) > 0]



# Subset if testing
if (TEST_MODE) {
  df_full_impute <- df_full_impute[1:TEST_N, ] %>% droplevels()
}

# Run
cat("[1/4] Running imputation...\n")
mr_full <- run_imputation(df_full_impute, maxiter = ifelse(TEST_MODE, 1, 5))

# Extract and back-transform
cat("[2/4] Extracting completed data...\n")
df_full_filled <- completeData(mr_full)[[1]] %>% back_transform(log_vars)


# Diagnostics
cat("[3/4] Saving diagnostics...\n")
save_diagnostics(
  back_transform(df_full_impute, log_vars),
  df_full_filled,
  "full"
)

# Add back admin columns
cat("[4/4] Adding columns...\n")
df_full_out <- add_back_columns(
  if (TEST_MODE) df[1:TEST_N, ] else df,
  df_full_filled
)

# Save dataset
saveRDS(df_full_out, file.path(imputation_dir, "df_full_imputed.rds"))
cat("Saved: df_full_imputed.rds\n")


# Imputation 2 ------------------------------------------------------------
# LASSO variable selection set (DAG-respecting)

# Use confounders and exposures to impute exposures
# Use everything to impute biological variables and anthropometrics

cat("\n========== IMPUTATION 2: Selection set ==========\n")
df_sel <- if (TEST_MODE) df_selection[1:TEST_N, ] else df_selection

result_sel <- run_dag_imputation(df_sel, "Selection")

save_diagnostics(result_sel$original, result_sel$filled, "selection")
df_sel_out <- add_back_columns(df_sel, result_sel$filled)


# Save dataset
saveRDS(df_sel_out, file.path(imputation_dir, "df_selection_imputed.rds"))
cat("Saved: df_selection_imputed.rds\n")

# Imputation 3 ------------------------------------------------------------
# Prediction set (DAG-respecting) 

# Run for training set
cat("\n========== IMPUTATION 3: Train set ==========\n")
df_tr <- if (TEST_MODE) df_train[1:TEST_N, ] else df_train

result_train <- run_dag_imputation(df_tr, "Train", return_models = TRUE)

save_diagnostics(result_train$original, result_train$filled, "train")
df_train_out <- add_back_columns(df_tr, result_train$filled)

# Save dataset
saveRDS(df_train_out, file.path(imputation_dir, "df_train_imputed.rds"))
cat("Saved: df_train_imputed.rds\n")

# Run for test set using train model
cat("\n========== IMPUTATION 3: Test set ==========\n")
df_te <- if (TEST_MODE) df_test[1:TEST_N, ] else df_test

result_test <- apply_dag_imputation(df_te, result_train, "Test")

save_diagnostics(result_test$original, result_test$filled, "test")
df_test_out <- add_back_columns(df_te, result_test$filled)

# Save dataset
saveRDS(df_test_out, file.path(imputation_dir, "df_test_imputed.rds"))
cat("Saved: df_test_imputed.rds\n")


# Export imputed datasets -----------------------------------------------------

# saveRDS(df_full_out, file.path(imputation_dir, "df_full_imputed.rds"))
# saveRDS(df_sel_out, file.path(imputation_dir, "df_selection_imputed.rds"))
# saveRDS(df_train_out, file.path(imputation_dir, "df_train_imputed.rds"))
# saveRDS(df_test_out, file.path(imputation_dir, "df_test_imputed.rds"))

