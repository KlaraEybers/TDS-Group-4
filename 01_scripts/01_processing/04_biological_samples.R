rm(list = ls())

# Import packages

library(here)
library(dplyr)
library(ggplot2)
library(purrr)

# Load data

ukb_clean <- readRDS(here("00_data", "01_processed", "ukb_initial_cleaned.rds"))

# Extract relevant columns

biological_samples_cols <- c(
  "wbc.0.0","platelets.0.0","albumin.0.0","alp.0.0","alt.0.0",
  "apoa1.0.0","apob.0.0","ast.0.0","bilirubin.0.0","total_cholesterol.0.0",
  "creatinine.0.0","crp.0.0","cystatin_c.0.0","ggt.0.0","hba1c.0.0",
  "hdl_cholesterol.0.0","ldl_direct.0.0","lpa.0.0","triglycerides.0.0",
  "neutrophill_count.0.0", "lymphocyte_count.0.0", "monocyte_count.0.0"
)

ukb_clean <- ukb_clean %>% select(eid, any_of(biological_samples_cols))

# Data cleaning 

# 1) Check column classes
sapply(biological_samples_cols, function(v) class(ukb_clean[[v]]))

# 2) Identify non-numeric values in character columns
non_numeric_report <- lapply(biological_samples_cols, function(v) {
  
  col <- ukb_clean[[v]]
  
  if (!is.character(col)) return(NULL)
  
  bad_vals <- unique(col[!is.na(col) & is.na(suppressWarnings(as.numeric(col)))])
  if (length(bad_vals) == 0) return(NULL)
  
  data.frame(variable = v, bad_value = bad_vals)
})

non_numeric_report <- bind_rows(non_numeric_report)
non_numeric_report

# 3) Change character columns to numeric
ukb_clean <- ukb_clean %>%
  mutate(across(any_of(biological_samples_cols), ~ suppressWarnings(as.numeric(as.character(.)))))

sapply(biological_samples_cols, function(v) class(ukb_clean[[v]]))

# 4) Apply range rules
ukb_clean <- ukb_clean %>%
  mutate(
    wbc.0.0 = ifelse(wbc.0.0 < 0 | wbc.0.0 > 389.7, NA, wbc.0.0),
    platelets.0.0 = ifelse(platelets.0.0 < 0.3 | platelets.0.0 > 1821, NA, platelets.0.0),
    albumin.0.0 = ifelse(albumin.0.0 < 17.38 | albumin.0.0 > 59.8, NA, albumin.0.0),
    alp.0.0 = ifelse(alp.0.0 < 8 | alp.0.0 > 1416.7, NA, alp.0.0),
    alt.0.0 = ifelse(alt.0.0 < 3.01 | alt.0.0 > 495.19, NA, alt.0.0),
    apoa1.0.0 = ifelse(apoa1.0.0 < 0.419 | apoa1.0.0 > 2.5, NA, apoa1.0.0),
    apob.0.0 = ifelse(apob.0.0 < 0.4 | apob.0.0 > 2, NA, apob.0.0),
    ast.0.0 = ifelse(ast.0.0 < 3.3 | ast.0.0 > 947.2, NA, ast.0.0),
    bilirubin.0.0 = ifelse(bilirubin.0.0 < 1 | bilirubin.0.0 > 70.06, NA, bilirubin.0.0),
    total_cholesterol.0.0 = ifelse(total_cholesterol.0.0 < 0.601 | total_cholesterol.0.0 > 15.46, NA, total_cholesterol.0.0),
    creatinine.0.0 = ifelse(creatinine.0.0 < 10.7 | creatinine.0.0 > 1499.3, NA, creatinine.0.0),
    crp.0.0 = ifelse(crp.0.0 < 0.08 | crp.0.0 > 79.96, NA, crp.0.0),
    cystatin_c.0.0 = ifelse(cystatin_c.0.0 < 0.295 | cystatin_c.0.0 > 7.487, NA, cystatin_c.0.0),
    ggt.0.0 = ifelse(ggt.0.0 < 5 | ggt.0.0 > 1184.9, NA, ggt.0.0),
    hba1c.0.0 = ifelse(hba1c.0.0 < 15 | hba1c.0.0 > 515.2, NA, hba1c.0.0),
    hdl_cholesterol.0.0 = ifelse(hdl_cholesterol.0.0 < 0.219 | hdl_cholesterol.0.0 > 4.401, NA, hdl_cholesterol.0.0),
    ldl_direct.0.0 = ifelse(ldl_direct.0.0 < 0.266 | ldl_direct.0.0 > 9.797, NA, ldl_direct.0.0),
    lpa.0.0 = ifelse(lpa.0.0 < 3.8 | lpa.0.0 > 189, NA, lpa.0.0),
    triglycerides.0.0 = ifelse(triglycerides.0.0 < 0.231 | triglycerides.0.0 > 11.278, NA, triglycerides.0.0),
    neutrophill_count.0.0 = ifelse(neutrophill_count.0.0 < 0 | neutrophill_count.0.0 > 76.42, NA, neutrophill_count.0.0),
    lymphocyte_count.0.0 = ifelse(lymphocyte_count.0.0 < 0 | lymphocyte_count.0.0 > 196.41, NA, lymphocyte_count.0.0),
    monocyte_count.0.0 = ifelse(monocyte_count.0.0 < 0 | monocyte_count.0.0 > 113.39, NA, monocyte_count.0.0)
  )

# 5) NA summary for these columns
cols_to_check <- biological_samples_cols

na_summary <- data.frame(
  variable = cols_to_check,
  n_na = sapply(cols_to_check, function(v) sum(is.na(ukb_clean[[v]]))),
  pct_na = sapply(cols_to_check, function(v) mean(is.na(ukb_clean[[v]])) * 100)
)

na_summary

# Save RDS

saveRDS(ukb_clean, here("00_data", "01_processed", "bio_samples_cleaned.rds"))
