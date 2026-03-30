## =============================================================================
# 00_variable_config.R
# Central configuration for variable classifications and output paths
# Source this file at the top of every analysis script instead of
# hardcoding variable lists.
#
# Usage:
#   source(here("01_scripts", "00_variable_config.R"))
#
# To switch scenarios, either:
#   1. Set SCENARIO before sourcing:  SCENARIO <- "sensitivity_anthro_confounder"
#   2. Pass via PBS/environment:      export SCENARIO=sensitivity_anthro_confounder
#   3. Default is "main" if not set
#
# Scenarios:
#   "main"                         — Anthropometrics excluded entirely. Imputed data.
#   "sensitivity_anthro_confounder" — Anthropometrics treated as confounders
#                                     (unpenalised in LASSO). Separate imputation run.
#   "complete_case"                — Anthropometrics excluded (same as main).
#                                     No imputation; only participants with zero
#                                     missing predictor values are retained.
#
# Exports:
#   Variable lists:
#     confounder_vars, exposure_vars, bio_vars, anthro_vars, admin_vars,
#     exclude_vars, confounders (alias for confounder_vars)
#
#   Data input:
#     data_input_dir  — directory containing selection/train/test datasets
#     selection_file, train_file, test_file — dataset filenames
#     imputation_dir, imputation_fig_dir — imputation-specific (NULL for complete_case)
#
#   Output paths:
#     lasso_outcome_dir, lasso_mediator_dir — LASSO data outputs
#     lasso_fig_dir, lasso_tbl_dir — LASSO figures and tables
#     incremental_fig_dir, incremental_tbl_dir — incremental AUC analysis
#     figures_multi_dir, tables_multi_dir — prediction and mediation results
#     sex_fig_dir, sex_tbl_dir — sex-stratified sensitivity analysis
#
# Notes:
#   - sbp is always removed from modelling in individual scripts (not here)
#     because it dominates variable selection regardless of scenario
#   - rm(list = ls()) is safe before sourcing because SCENARIO falls back
#     to Sys.getenv("SCENARIO") with default "main"
#   - All output directories are auto-created at the end of this file
# =============================================================================
library(here)

# ── Scenario switch ──────────────────────────────────────────────────────────
# Only set SCENARIO if it hasn't been set already (allows override)
if (!exists("SCENARIO")) {
  SCENARIO <- Sys.getenv("SCENARIO", unset = "main")
}

cat("=== CONFIG: SCENARIO =", SCENARIO, "===\n")

# ── Anthropometric variables (the contested group) ───────────────────────────
anthro_vars <- c("sbp", "dbp", "grip_strength", "whr",
                 "heel_bone_density", "body_fat_mass")

# ── Admin / outcome variables (always the same) ─────────────────────────────
admin_vars <- c("eid", "recruitment_date", "death_date",
                "first_cvd_date", "has_cvd")

# ── Base variable definitions (without anthropometrics) ──────────────────────
base_confounder_vars <- c("age", "sex", "ethnicity")

base_exposure_vars <- c(
  "household_income", "highest_qualification", "employment",
  "housing", "household_size",
  "diet_score", "alcohol_consumption", "salt_added",
  "sleep_hrs", "sleep_cat", "insomnia",
  "smoking_status_refined", "pack_yrs", "passive_total_hrs", "passive_cat",
  "met_min_week_all",
  "mh_sought_help", "satisfaction_mean", "neuroticism_score",
  "breastfed", "comparative_body_size_10", "maternal_smoking_birth", "birth_weight",
  "no2_2010", "pm2.5_2010", "pm10_2010", "urban_rural_3lvl",
  "snoring", "has_vehicle"
)

base_bio_vars <- c(
  "wbc", "platelets", "albumin", "alp", "alt", "apoa1", "apob", "ast",
  "bilirubin", "total_cholesterol", "creatinine", "crp", "cystatin_c",
  "ggt", "hba1c", "hdl_cholesterol", "ldl_direct", "lpa", "triglycerides",
  "neutrophill_count", "lymphocyte_count", "monocyte_count"
)

# ── Scenario-specific assembly ───────────────────────────────────────────────
if (SCENARIO == "main") {
  # Main analysis: anthropometrics excluded entirely
  confounder_vars <- base_confounder_vars
  exposure_vars   <- base_exposure_vars
  bio_vars        <- base_bio_vars
  exclude_vars    <- anthro_vars
  
} else if (SCENARIO == "sensitivity_anthro_confounder") {
  # Sensitivity: anthropometrics treated as confounders (unpenalised)
  confounder_vars <- c(base_confounder_vars, anthro_vars)
  exposure_vars   <- base_exposure_vars
  bio_vars        <- base_bio_vars
  exclude_vars    <- character(0)
  
} else if (SCENARIO == "complete_case") {
  # Sensitivity: Complete case
  confounder_vars <- base_confounder_vars
  exposure_vars   <- base_exposure_vars
  bio_vars        <- base_bio_vars
  exclude_vars    <- anthro_vars
  
} else {
  stop("Unknown SCENARIO: ", SCENARIO,
       "\n  Valid options: 'main', 'sensitivity_anthro_confounder'")
}

# Convenience alias (some scripts use "confounders" without _vars suffix)
confounders <- confounder_vars

# ── Output paths ─────────────────────────────────────────────────────────────
if (SCENARIO == "main") {
  
  # Imputation
  imputation_dir     <- here("00_data", "01_processed", "03_imputation")
  imputation_fig_dir <- here("02_results", "00_figures", "01_imputation")
  data_input_dir     <- imputation_dir
  full_imputed_file  <- "df_full_imputed.rds"
  selection_file     <- "df_selection_imputed.rds"
  train_file         <- "df_train_imputed.rds"
  test_file          <- "df_test_imputed.rds"
  
  # Stability selection
  lasso_outcome_dir  <- here("00_data", "02_analysis", "stability_selection",
                             "stability_selection_runs_outcome")
  lasso_mediator_dir <- here("00_data", "02_analysis", "stability_selection",
                             "stability_selection_runs_mediator")
  lasso_fig_dir <- here("02_results", "00_figures", "04_exploratory",
                        "stability_selection_lasso")
  lasso_tbl_dir <- here("02_results", "01_tables", "08_lasso", "stability_selection_lasso")
  
  # Incremental analysis
  incremental_fig_dir <- here("02_results", "00_figures", "04_exploratory",
                              "lasso_incremental", "incremental_analysis")
  incremental_tbl_dir <- here("02_results", "01_tables", "08_lasso",
                              "incremental_analysis")
  
  # Multivariate results
  figures_multi_dir  <- here("02_results", "00_figures", "03_multivariate")
  tables_multi_dir   <- here("02_results", "01_tables", "01_multivariate")
  
  # Sex-stratified sensitivity
  sex_fig_dir <- here("02_results", "00_figures", "05_sensitivity_sex")
  sex_tbl_dir <- here("02_results", "01_tables", "05_sensitivity_sex")
  
} else if (SCENARIO == "sensitivity_anthro_confounder") {
  
  # Imputation
  imputation_dir     <- here("00_data", "01_processed", "03_imputation",
                             "sensitivity_anthro_confounder")
  imputation_fig_dir <- here("02_results", "00_figures",
                             "08_sensitivity_anthro_confounder", "imputation")
  data_input_dir     <- imputation_dir
  full_imputed_file  <- "df_full_imputed.rds"
  selection_file     <- "df_selection_imputed.rds"
  train_file         <- "df_train_imputed.rds"
  test_file          <- "df_test_imputed.rds"
  
  # FAMD
  famd_fig_dir <- here("02_results", "00_figures", "04_exploratory", "famd_correlation")
  famd_tbl_dir <- here("02_results", "01_tables",  "07_famd_correlation")
  
  # Stability selection
  lasso_outcome_dir  <- here("00_data", "02_analysis", "sensitivity_analysis",
                             "anthro_as_confounders", "lasso1_outcome")
  lasso_mediator_dir <- here("00_data", "02_analysis", "sensitivity_analysis",
                             "anthro_as_confounders", "lasso2_mediator")
  lasso_fig_dir <- here("02_results", "00_figures",
                        "08_sensitivity_anthro_confounder", "stability_selection_lasso")
  lasso_tbl_dir <- here("02_results", "01_tables", "06_sensitivity_anthro_confounder")
  
  # Incremental analysis
  incremental_fig_dir <- here("02_results", "00_figures",
                              "08_sensitivity_anthro_confounder",
                              "incremental_analysis")
  incremental_tbl_dir <- here("02_results", "01_tables", "06_sensitivity_anthro_confounder",
                              "incremental_analysis")
  
  
  # Multivariate results
  figures_multi_dir  <- here("02_results", "00_figures",
                             "08_sensitivity_anthro_confounder")
  tables_multi_dir   <- here("02_results", "01_tables",
                             "06_sensitivity_anthro_confounder")
  
  # Sex-stratified sensitivity
  sex_fig_dir <- here("02_results", "00_figures",
                      "08_sensitivity_anthro_confounder", "sensitivity_sex")
  sex_tbl_dir <- here("02_results", "01_tables",
                      "08_sensitivity_anthro_confounder", "sensitivity_sex")
  
} else if (SCENARIO == "complete_case") {
  
  # Imputation (not used, but kept for consistency)
  imputation_dir     <- NULL
  imputation_fig_dir <- NULL
  data_input_dir     <- here("00_data", "01_processed", "04_complete_case")
  full_imputed_file <- NULL
  selection_file     <- "cc_selection.rds"
  train_file         <- "cc_train.rds"
  test_file          <- "cc_test.rds"
  
  # Stability selection
  lasso_outcome_dir  <- here("00_data", "02_analysis", "sensitivity_analysis",
                             "complete_case", "lasso1_outcome")
  lasso_mediator_dir <- here("00_data", "02_analysis", "sensitivity_analysis",
                             "complete_case", "lasso2_mediator")
  
  # Figures and tables
  lasso_fig_dir      <- here("02_results", "00_figures",
                             "06_sensitivity_complete_cases", "stability_selection_lasso")
  lasso_tbl_dir <- here("02_results", "01_tables", "03_sensitivity_complete_case", "stability_selection_lasso")
  incremental_fig_dir <- here("02_results", "00_figures",
                              "06_sensitivity_complete_cases", "incremental_analysis")
  incremental_tbl_dir <- here("02_results", "01_tables", "03_sensitivity_complete_case",
                              "incremental_analysis")
  figures_multi_dir  <- here("02_results", "00_figures",
                             "06_sensitivity_complete_cases")
  tables_multi_dir   <- here("02_results", "01_tables",
                             "03_sensitivity_complete_case")
  
  # Sex-stratified
  sex_fig_dir <- here("02_results", "00_figures",
                      "06_sensitivity_complete_cases", "sensitivity_sex")
  sex_tbl_dir <- here("02_results", "01_tables",
                      "03_sensitivity_complete_case", "sensitivity_sex")
}
# Create directories if they don't exist
for (d in c(imputation_dir, imputation_fig_dir, lasso_outcome_dir, lasso_mediator_dir,
            incremental_fig_dir, incremental_tbl_dir, figures_multi_dir, tables_multi_dir, lasso_fig_dir,
            lasso_tbl_dir, sex_fig_dir, sex_tbl_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── Summary ──────────────────────────────────────────────────────────────────
cat("  Confounders:", length(confounder_vars), "-", paste(confounder_vars, collapse = ", "), "\n")
cat("  Exposures:  ", length(exposure_vars), "\n")
cat("  Biomarkers: ", length(bio_vars), "\n")
cat("  Excluded:   ", length(exclude_vars),
    if (length(exclude_vars) > 0) paste("-", paste(exclude_vars, collapse = ", ")) else "", "\n")
cat("  Imputation dir: ", imputation_dir, "\n")
cat("  LASSO outcome:  ", lasso_outcome_dir, "\n")
cat("  LASSO mediator: ", lasso_mediator_dir, "\n")