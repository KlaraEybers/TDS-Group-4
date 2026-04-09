# ------------------------------------------------------------------------------
# Script: 10_final_clean.R
# Purpose: Merge cleaned domain datasets, attach CVD events, and create the final analysis dataset.
#
# Inputs:
#   - 00_data/01_processed/bio_samples_cleaned.rds
#   - 00_data/01_processed/diet_cleaned.rds
#   - 00_data/01_processed/early_life_cleaned.rds
#   - 00_data/01_processed/environment_cleaned.rds
#   - 00_data/01_processed/lifestyle_cleaned.rds
#   - 00_data/01_processed/medical_history_cleaned.rds
#   - 00_data/01_processed/mental_health_cleaned.rds
#   - 00_data/01_processed/physical_cleaned.rds
#   - 00_data/01_processed/sociodemographic_cleaned.rds
#   - 00_data/00_extracted/cvd_events.rds
#
# Outputs: 00_data/01_processed/ukb_final_merged.rds
#
# Upstream:
#   - 01_socio_demographics.R
#   - 02_physical_measures.R
#   - 03_diet.R
#   - 04_biological_samples.R
#   - 05_health_and_medical_history.R
#   - 06_early_life_factors.R
#   - 07_environment.R
#   - 08_lifestyle.R
#   - 09_mental_health.R
#
# Downstream:
#   - 11_remove_missing.R
#   - 12_imputation.R
#
# Key steps:
#   - Import all cleaned domain datasets and CVD event data
#   - Derive earliest CVD event date for each participant
#   - Merge all datasets by eid
#   - Create binary CVD event indicator
#   - Exclude participants with CVD before or at recruitment
#   - Remove .0.0 suffixes from remaining variable names
#   - Summarise missingness across variables and participants
#
# Notes:
#   - Earliest recorded CVD event is retained per participant
#   - Participants with first CVD event before or on recruitment date are excluded
#   - Missingness summaries are used for QC but not saved separately in this script
# ------------------------------------------------------------------------------
rm(list = ls())

# Import packages ---------------------------------------------------------

library(stringr)
library(dplyr)
library(here)

# Import data -------------------------------------------------------------

# Import cleaned domains
bio_samples_cleaned <- readRDS(here("00_data", "01_processed", "bio_samples_cleaned.rds"))
diet_cleaned <- readRDS(here("00_data", "01_processed", "diet_cleaned.rds"))
early_life_cleaned <- readRDS(here("00_data", "01_processed", "early_life_cleaned.rds"))
environment_cleaned <- readRDS(here("00_data", "01_processed", "environment_cleaned.rds"))
lifestyle_cleaned <- readRDS(here("00_data", "01_processed", "lifestyle_cleaned.rds"))
medical_history_cleaned <- readRDS(here("00_data", "01_processed", "medical_history_cleaned.rds"))
mental_health_cleaned <- readRDS(here("00_data", "01_processed", "mental_health_cleaned.rds"))
physical_cleaned <- readRDS(here("00_data", "01_processed", "physical_cleaned.rds"))
sociodemographic_cleaned <- readRDS(here("00_data", "01_processed", "sociodemographic_cleaned.rds"))

# Import cvd events

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  cvd_events_path <- args[1]
} else {
  cvd_events_path <- here("00_data", "00_extracted", "cvd_events.rds")
}

cvd_events <- readRDS(cvd_events_path)


# Clean CVD events --------------------------------------------------------

cvd_cleaned <- cvd_events %>%
  mutate(date = as.Date(date)) %>%
  group_by(eid) %>%
  arrange(date) %>%
  slice(1) %>%                      # keep earliest CVD event per participant
  ungroup() %>%
  rename(first_cvd_date = date) %>%
  select(eid, first_cvd_date)


# Merge datasets on eid ---------------------------------------------------

domains_merged <- sociodemographic_cleaned %>%
  left_join(physical_cleaned, by = "eid") %>%
  left_join(diet_cleaned, by = "eid") %>%
  left_join(lifestyle_cleaned, by = "eid") %>%
  left_join(medical_history_cleaned, by = "eid") %>%
  left_join(mental_health_cleaned, by = "eid") %>%
  left_join(early_life_cleaned, by = "eid") %>%
  left_join(environment_cleaned, by = "eid") %>%
  left_join(bio_samples_cleaned, by = "eid")

# Merge cvd events to ukb -------------------------------------------------

final_merged <- domains_merged %>%
  left_join(cvd_cleaned, by = "eid") %>%
  mutate(has_cvd = ifelse(!is.na(first_cvd_date), 1, 0))



# Remove CVD before recruitment and pregnant women ------------------------


# Remove participants with first CVD event before or at recruitment
# And pregnant women

final_filtered <- final_merged %>%
  filter(is.na(first_cvd_date) | first_cvd_date >= recruitment_date) %>%
  filter(is.na(`pregnancy.0.0`) | `pregnancy.0.0` != "Yes") %>%
  select(-`pregnancy.0.0`)


# Remove .0.0 suffixes  ----------------------------------------------------

final_filtered <- final_filtered %>%
  rename_with(~str_replace(., "\\.0\\.0$", ""))


# Drop unneccesary columns ------------------------------------------------

cols_to_drop <- c("diabetes_doctor_diagnosed", "serious_dx",
                  "age_high_bp_dx", "age_diabetes_dx", "angina_dx_age",
                  "mi_dx_age", "dvt_dx_age", "pe_dx_age",
                  "stroke_dx_age", "heart_problems_collapsed",
                  "meds_chol_bp_diabetes_hormone_collapsed",
                  "meds_chol_bp_diabetes_collapsed", "long_ill",
                  "n_previous_pregnancies_collapsed")

ncol_before <- ncol(final_filtered)
final_filtered <- final_filtered %>% select(-any_of(cols_to_drop))
ncol_after <- ncol(final_filtered)

cat("Columns before:", ncol_before, "\n")
cat("Columns after:", ncol_after, "\n")
cat("Dropped:", ncol_before - ncol_after, "| Expected:", length(cols_to_drop), "\n")


# Export dataset ----------------------------------------------------------

saveRDS(final_filtered, here("00_data", "01_processed", "ukb_final_merged.rds"))

