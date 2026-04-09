# ------------------------------------------------------------------------------
# Script: 12_imputation.R
# Purpose: Run imputation pipelines for full, selection, training, and test datasets.
#
# Inputs:
#   - 00_data/01_processed/ukb_final_cleaned.rds
#   - 01_scripts/01_processing/12_imputation_functions.R
#
# Outputs:
#   - 00_data/01_processed/03_imputation/df_selection_imputed.rds
#   - 00_data/01_processed/03_imputation/df_train_imputed.rds
#   - 00_data/01_processed/03_imputation/df_test_imputed.rds
#   - 02_results/00_figures/01_imputation/numeric_diagnostics_<label>.pdf
#   - 02_results/00_figures/01_imputation/categorical_diagnostics_<label>.pdf
#
# Upstream:
#   - 11_remove_missing.R
#   - 12_imputation_functions.R
#
# Downstream:
#   - 02_analysis/01_imputation/
#   - 02_analysis/02_univariate/
#   - 02_analysis/04_multivariate/
#   - 02_analysis/05_neural_network/
#   - 02_analysis/06_decision_trees/
#
# Key steps:
#   - Import cleaned data and imputation helper functions
#   - Define administrative, confounder, exposure, and biological variable sets
#   - Remove participants with missing ethnicity
#   - Identify skewed biomarkers and specify variables for log transformation
#   - Split data into selection, training, and test sets
#   - Run full-dataset imputation for descriptive and univariate analyses
#   - Run DAG-respecting staged imputation for selection, train, and test datasets
#   - Save imputed datasets and diagnostic plots
#
# Notes:
#   - TEST_MODE allows reduced-size runs for debugging
#   - Ethnicity is not imputed and rows with missing ethnicity are removed
#   - Log-transformed variables are back-transformed after imputation
#   - Test-set imputation uses models trained on the training set
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

# Toggle for test run or job script run -----------------------------------

TEST_MODE <- F  # Set to FALSE for full run
TEST_N <- 40000


# Classify variables ------------------------------------------------------
# Classifications: admin/outcomes, confounders, exposures, biological

admin_vars <- c("eid", "recruitment_date", "death_date", "first_cvd_date", "has_cvd")

confounder_vars <- c("age", "sex", "ethnicity")

exposure_vars <- c(
  "household_income", "highest_qualification", "employment", "housing", "household_size",
  "diet_score", "alcohol_consumption", "salt_added",
  "sleep_hrs", "sleep_cat", "insomnia",
  "smoking_status_refined", "pack_yrs", "passive_total_hrs", "passive_cat",
  "met_min_week_all",
  "mh_sought_help", "work_job_satisfaction_num", "neuroticism_score",
  "breastfed", "comparative_body_size_10", "maternal_smoking_birth", "birth_weight",
  "no2_2010", "pm2.5_2010", "pm10_2010", "urban_rural_3lvl"
)

bio_vars <- c(
  "sbp", "dbp", "grip_strength", "whr",
  "wbc", "platelets", "albumin", "alp", "alt", "apoa1", "apob", "ast",
  "bilirubin", "total_cholesterol", "creatinine", "crp", "cystatin_c",
  "ggt", "hba1c", "hdl_cholesterol", "ldl_direct", "lpa", "triglycerides", 
  "neutrophill_count", "lymphocyte_count", "monocyte_count"
  )

# Removed: "diabetes_doctor_diagnosed", "serious_dx",
#"age_high_bp_dx", "age_diabetes_dx", "angina_dx_age", "mi_dx_age",
# "dvt_dx_age", "pe_dx_age", "stroke_dx_age",
# "heart_problems_collapsed", "meds_chol_bp_diabetes_hormone_collapsed",
# "meds_chol_bp_diabetes_collapsed", "long_ill", "n_previous_pregnancies_collapsed"
# "heel_bone_density",  "body_fat_mass",


# Remove NA ethnicity -----------------------------------------------------
# Remove missing ethnicity (avoid imputing confounders)

# Confirm how many NAs for confounders

colSums(is.na(df[, confounder_vars]))

# Real UK Biobank has only around 2 000 NAs for ethnicity
df <- df %>% drop_na(ethnicity)

# Output numbers for flow diagram
cat("N:", nrow(df), "\n")
cat("CVD events:", sum(df$has_cvd == 1), "\n")
cat("No CVD events:", sum(df$has_cvd == 0), "\n")

# Check no missing confounder variables 
colSums(is.na(df[, confounder_vars]))


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
idx_selection <- sample(nrow(df), size = 0.5 * nrow(df))

df_selection <- df[idx_selection, ]     # 50% for LASSO
df_inference <- df[-idx_selection, ]    # 50% for prediction models

# Split prediction set into train/test set
idx_train <- sample(nrow(df_inference), size = 0.75 * nrow(df_inference))
df_train <- df_inference[idx_train, ]   # refit
df_test  <- df_inference[-idx_train, ]  # evaluation


############################################################################
# RUN IMPUTATION PIPELINES
############################################################################

# Imputation 1 Pipeline ---------------------------------------------------

# Run imputation on the full dataset for univariate and descriptive analysis
cat("\n========== IMPUTATION 1: Full dataset ==========\n")

all_vars <- c(confounder_vars, exposure_vars, bio_vars)
df_full_impute <- prep_for_imputation(df, vars_to_impute = all_vars, predictor_vars = c())

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


# Imputation 2 ------------------------------------------------------------
# LASSO variable selection set (DAG-respecting)

# Use confounders and exposures to impute exposures
# Use everything to impute biological variables and anthropometrics

cat("\n========== IMPUTATION 2: Selection set ==========\n")
df_sel <- if (TEST_MODE) df_selection[1:TEST_N, ] else df_selection

result_sel <- run_dag_imputation(df_sel, "Selection")

save_diagnostics(result_sel$original, result_sel$filled, "selection")
df_sel_out <- add_back_columns(df_sel, result_sel$filled)


sum(is.na(df_selection[1:TEST_N, ][["heart_problems_collapsed"]]))
class(result_sel$original[["heart_problems_collapsed"]])
class(result_sel$filled[["heart_problems_collapsed"]])

# Imputation 3 ------------------------------------------------------------
# Prediction set (DAG-respecting) 

# Run for training set
cat("\n========== IMPUTATION 3: Train set ==========\n")
df_tr <- if (TEST_MODE) df_train[1:TEST_N, ] else df_train

result_train <- run_dag_imputation(df_tr, "Train", return_models = TRUE)

save_diagnostics(result_train$original, result_train$filled, "train")
df_train_out <- add_back_columns(df_tr, result_train$filled)

# Run for test set using train model
cat("\n========== IMPUTATION 3: Test set ==========\n")
df_te <- if (TEST_MODE) df_test[1:TEST_N, ] else df_test

result_test <- apply_dag_imputation(df_te, result_train, "Test")

save_diagnostics(result_test$original, result_test$filled, "test")
df_test_out <- add_back_columns(df_te, result_test$filled)


# Export imputed datasets -----------------------------------------------------

# saveRDS(df_full_out, here("00_data", "01_processed", "03_imputation", "df_full_imputed.rds"))
saveRDS(df_sel_out, here("00_data", "01_processed", "03_imputation", "df_selection_imputed.rds"))
saveRDS(df_train_out, here("00_data", "01_processed", "03_imputation", "df_train_imputed.rds"))
saveRDS(df_test_out, here("00_data", "01_processed", "03_imputation", "df_test_imputed.rds"))


