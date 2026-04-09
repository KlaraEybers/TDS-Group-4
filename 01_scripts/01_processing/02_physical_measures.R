# ------------------------------------------------------------------------------
# Script: 02_physical_measures.R
# Purpose: Extract, clean, and derive physical measure variables.
#
# Inputs: 00_data/01_processed/ukb_initial_cleaned.rds
#
# Outputs: 00_data/01_processed/physical_cleaned.rds
#
# Upstream: 00_initial_clean.R
#
# Downstream: 10_final_clean.R
#
# Key steps:
#   - Select physical measure variables from UKB dataset
#   - Clean and average blood pressure readings (SBP, DBP)
#   - Derive average grip strength from left and right hand
#   - Calculate waist-to-hip ratio
#   - Clean and recode body fat mass and heel bone density
#
# Notes:
#   - Blood pressure values outside physiological ranges set to NA
#   - Grip strength missing code (-1) recoded to NA
#   - BMI excluded in favour of waist-to-hip ratio
# ------------------------------------------------------------------------------

rm(list = ls())

# Import packages ---------------------------------------------------------

library(here)
library(dplyr)


# Import data -------------------------------------------------------------

ukb_clean <- readRDS(here("00_data", "01_processed", "ukb_initial_cleaned.rds"))

# Extract relevant columns ------------------------------------------------

# Pattern to match physical measures variables
pattern <- paste(c(
  "hip_circumference",
  "waist_circumference", 
  "bmi",
  "whole_b_f_mass",
  "heel_bdensity_tscore",
  "dbp_automated",
  "sbp_automated",
  "grip_str_l",
  "grip_str_r",
  "sbp_manual",
  "dbp_manual"
), collapse = "|")


# Extract all columns that match any of these patterns
matching_cols <- grep(pattern, colnames(ukb_clean), value = TRUE)

# Select relevant columnns
selected_data <- ukb_clean %>% select(eid, any_of(matching_cols))

#--------------------------------------------------------------------------
# Data pre-processing -----------------------------------------------------
#--------------------------------------------------------------------------

# Blood pressure ----------------------------------------------------------
# Remove outliers for SBP (70-250) and for DBP (40-150) 
# Average over all 4 columns for SBP and DBP respectively

physical_clean <- selected_data %>%
  mutate(
    # Set outliers to NA for SBP (valid range: 70-250 mmHg)
    sbp_manual.0.0_clean = ifelse(sbp_manual.0.0 >= 70 & sbp_manual.0.0 <= 250, 
                                  sbp_manual.0.0, NA),
    sbp_manual.0.1_clean = ifelse(sbp_manual.0.1 >= 70 & sbp_manual.0.1 <= 250, 
                                  sbp_manual.0.1, NA),
    sbp_automated.0.0_clean = ifelse(sbp_automated.0.0 >= 70 & sbp_automated.0.0 <= 250, 
                                     sbp_automated.0.0, NA),
    sbp_automated.0.1_clean = ifelse(sbp_automated.0.1 >= 70 & sbp_automated.0.1 <= 250, 
                                     sbp_automated.0.1, NA),
    
    # Set outliers to NA for DBP (valid range: 40-150 mmHg)
    dbp_manual.0.0_clean = ifelse(dbp_manual.0.0 >= 40 & dbp_manual.0.0 <= 150, 
                                  dbp_manual.0.0, NA),
    dbp_manual.0.1_clean = ifelse(dbp_manual.0.1 >= 40 & dbp_manual.0.1 <= 150, 
                                  dbp_manual.0.1, NA),
    dbp_automated.0.0_clean = ifelse(dbp_automated.0.0 >= 40 & dbp_automated.0.0 <= 150, 
                                     dbp_automated.0.0, NA),
    dbp_automated.0.1_clean = ifelse(dbp_automated.0.1 >= 40 & dbp_automated.0.1 <= 150, 
                                     dbp_automated.0.1, NA),

    # Average all valid SBP readings
    sbp = rowMeans(cbind(sbp_manual.0.0_clean, sbp_manual.0.1_clean, 
                          sbp_automated.0.0_clean, sbp_automated.0.1_clean), 
                   na.rm = TRUE),
    
    # Average all valid DBP readings
    dbp = rowMeans(cbind(dbp_manual.0.0_clean, dbp_manual.0.1_clean,
                          dbp_automated.0.0_clean, dbp_automated.0.1_clean),
                   na.rm = TRUE),

    # Set to NA if all readings were missing
    sbp = ifelse(is.nan(sbp), NA, sbp),
    dbp = ifelse(is.nan(dbp), NA, dbp)
  )

# Drop intermediate cleaning columns and original BP columns
physical_clean <- physical_clean %>%
  select(-contains("_clean"),
         -sbp_manual.0.0, -sbp_manual.0.1, 
         -sbp_automated.0.0, -sbp_automated.0.1,
         -dbp_manual.0.0, -dbp_manual.0.1,
         -dbp_automated.0.0, -dbp_automated.0.1)


# Hand grip strength ------------------------------------------------------

# Average left and right, set -1 to NA
physical_clean <- physical_clean %>%
  mutate(
    # Set -1 (missing code) to NA
    grip_str_l_clean = ifelse(grip_str_l.0.0 == -1, NA, grip_str_l.0.0),
    grip_str_r_clean = ifelse(grip_str_r.0.0 == -1, NA, grip_str_r.0.0),
    
    # Average left and right grip strength
    grip_strength = rowMeans(cbind(grip_str_l_clean, grip_str_r_clean), na.rm = TRUE),
    
    # Set to NA if both hands were missing
    grip_strength = ifelse(is.nan(grip_strength), NA, grip_strength)
  )

# Drop intermediate columns and originals
physical_clean <- physical_clean %>%
  select(-grip_str_l_clean, -grip_str_r_clean,
         -grip_str_l.0.0, -grip_str_r.0.0)


# Waist-to-Hip Ratio ------------------------------------------------------

# Calculate ratio and removed BMI for now
# No cutoffs chosen for outliers -> unclear physiologically impossible thresholds

# Create waist-to-hip ratio
physical_clean <- physical_clean %>%
  mutate(
    whr = waist_circumference.0.0 / hip_circumference.0.0
  ) %>%
  select(-bmi.0.0, -waist_circumference.0.0, -hip_circumference.0.0)


# Whole body fat mass and heel bone density -------------------------------

# Clean heel bone density and body fat mass
physical_clean <- physical_clean %>%
  mutate(
    # Convert heel bone density from character to numeric
    heel_bone_density = as.numeric(heel_bdensity_tscore.0.0),
    
    # Body fat mass: set negative values to NA (physiologically impossible)
    body_fat_mass = ifelse(whole_b_f_mass.0.0 < 0, NA, whole_b_f_mass.0.0)
  )

# Drop original columns
physical_clean <- physical_clean %>%
  select(-heel_bdensity_tscore.0.0, -whole_b_f_mass.0.0)


# Save RDS ----------------------------------------------------------------

saveRDS(physical_clean, here("00_data", "01_processed", "physical_cleaned.rds"))

