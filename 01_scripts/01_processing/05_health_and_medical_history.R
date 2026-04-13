rm(list = ls())

# Import packages

library(here)
library(dplyr)
library(ggplot2)
library(purrr)


# Load data

ukb_clean <- readRDS(here("00_data", "01_processed", "ukb_initial_cleaned.rds"))

# Extract relevant columns

pattern <- paste(c(
  "n_noncancer_ill",
  "long_ill",
  "diabetes_doctor_diagnosed",
  "serious_dx",
  "age_high_bp_dx",
  "age_diabetes_dx",
  # "early_insulin",
  "pregnancy",
  "angina_dx_age",
  "mi_dx_age",
  "dvt_dx_age",
  "pe_dx_age",
  "stroke_dx_age",
  "heart_problems",
  "meds_chol_bp_diabetes_hormone",
  "meds_chol_bp_diabetes",
  # "n_previous_pregnancies",
  "date_mi",
  "stroke_date",
  "cm_dx_date",
  "af_dx_date",
  "hf_dx_date",
  "pvd_dx_date"
), collapse = "|")

matching_cols <- grep(pattern, colnames(ukb_clean), value = TRUE)

ukb_clean <- ukb_clean %>% select(eid, any_of(matching_cols))

# Data cleaning

drop_prefixes <- c(
  "n_noncancer_ill",
  "date_mi",
  "stroke_date",
  "cm_dx_date",
  "af_dx_date",
  "hf_dx_date",
  "pvd_dx_date"
)

drop_pattern <- paste0("^(", paste(drop_prefixes, collapse = "|"), ")")
cols_to_drop <- grep(drop_pattern, colnames(ukb_clean), value = TRUE)

ukb_clean <- ukb_clean %>% select(-any_of(cols_to_drop))

# Turn na_strings values into NAs

na_strings <- c("Prefer not to answer", "Do not know", "Unsure")

base_fields <- c(
  "long_ill",
  "diabetes_doctor_diagnosed",
  "serious_dx",
  # "early_insulin",
  "pregnancy",
  "heart_problems",
  "meds_chol_bp_diabetes_hormone",
  "meds_chol_bp_diabetes"
)

cols_to_recode <- grep(paste0("^(", paste(base_fields, collapse = "|"), ")"), names(ukb_clean), value = TRUE)

ukb_clean <- ukb_clean %>%
  mutate(
    across(
      any_of(cols_to_recode),
      ~ {
        x <- .
        if (is.factor(x)) x <- as.character(x)
        if (is.character(x)) x[x %in% na_strings] <- NA_character_
        x
      }
    )
  ) %>%
  mutate(across(any_of(cols_to_recode), ~ if (is.character(.x)) as.factor(.x) else .x))

# Collapse variables with multiple columns

heart_cols <- c("heart_problems.0.0", "heart_problems.0.1", "heart_problems.0.2", "heart_problems.0.3")
heart_yes_vals <- c("Heart attack", "Angina", "Stroke", "High blood pressure")

meds_horm_cols <- c(
  "meds_chol_bp_diabetes_hormone.0.0",
  "meds_chol_bp_diabetes_hormone.0.1",
  "meds_chol_bp_diabetes_hormone.0.2",
  "meds_chol_bp_diabetes_hormone.0.3"
)
meds_horm_yes_vals <- c(
  "Cholesterol lowering medication",
  "Blood pressure medication",
  "Insulin",
  "Hormone replacement therapy",
  "Oral contraceptive pill or minipill"
)

meds_cols <- c("meds_chol_bp_diabetes.0.0", "meds_chol_bp_diabetes.0.1", "meds_chol_bp_diabetes.0.2")
meds_yes_vals <- c("Cholesterol lowering medication", "Blood pressure medication", "Insulin")

preg_cols <- c(
  "n_previous_pregnancies.0.0",
  "n_previous_pregnancies.0.1",
  "n_previous_pregnancies.0.2",
  "n_previous_pregnancies.0.3",
  "n_previous_pregnancies.0.4",
  "n_previous_pregnancies.0.5"
)

ukb_clean <- ukb_clean %>%
  mutate(across(any_of(preg_cols), ~ suppressWarnings(as.numeric(as.character(.))))) %>%
  mutate(
    heart_problems_collapsed = case_when(
      if_any(any_of(heart_cols), ~ .x %in% heart_yes_vals) ~ "Yes",
      if_all(any_of(heart_cols), ~ .x == "None of the above") ~ "No",
      if_all(any_of(heart_cols), ~ .x == "Prefer not to say") ~ "Prefer not to say",
      TRUE ~ NA_character_
    ),
    
    meds_chol_bp_diabetes_hormone_collapsed = case_when(
      if_any(any_of(meds_horm_cols), ~ .x %in% meds_horm_yes_vals) ~ "Yes",
      if_all(any_of(meds_horm_cols), ~ .x == "None of the above") ~ "No",
      if_all(any_of(meds_horm_cols), ~ .x == "Do not know") ~ "Do not know",
      if_all(any_of(meds_horm_cols), ~ .x == "Prefer not to say") ~ "Prefer not to say",
      TRUE ~ NA_character_
    ),
    
    meds_chol_bp_diabetes_collapsed = case_when(
      if_any(any_of(meds_cols), ~ .x %in% meds_yes_vals) ~ "Yes",
      if_all(any_of(meds_cols), ~ .x == "None of the above") ~ "No",
      if_all(any_of(meds_cols), ~ .x == "Do not know") ~ "Do not know",
      if_all(any_of(meds_cols), ~ .x == "Prefer not to say") ~ "Prefer not to say",
      TRUE ~ NA_character_
    ),
    
    # n_previous_pregnancies_collapsed = pmax(!!!rlang::syms(preg_cols), na.rm = TRUE),
    # n_previous_pregnancies_collapsed = ifelse(
    #   is.infinite(n_previous_pregnancies_collapsed),
    #   NA_real_,
    #   n_previous_pregnancies_collapsed
    # )
  ) %>%
  select(
    -any_of(heart_cols),
    -any_of(meds_horm_cols),
    -any_of(meds_cols),
    -any_of(preg_cols)
  )

# Convert implausible ages to NA

ukb_clean <- ukb_clean %>%
  mutate(
    age_high_bp_dx.0.0 = ifelse(age_high_bp_dx.0.0 < 18 | age_high_bp_dx.0.0 > 87, NA, age_high_bp_dx.0.0),
    age_diabetes_dx.0.0 = ifelse(age_diabetes_dx.0.0 < 1 | age_diabetes_dx.0.0 > 87, NA, age_diabetes_dx.0.0),
    angina_dx_age.0.0 = ifelse(angina_dx_age.0.0 < 18 | angina_dx_age.0.0 > 83, NA, angina_dx_age.0.0),
    mi_dx_age.0.0 = ifelse(mi_dx_age.0.0 < 18 | mi_dx_age.0.0 > 84, NA, mi_dx_age.0.0),
    dvt_dx_age.0.0 = ifelse(dvt_dx_age.0.0 < 18 | dvt_dx_age.0.0 > 84, NA, dvt_dx_age.0.0),
    pe_dx_age.0.0 = ifelse(pe_dx_age.0.0 < 1 | pe_dx_age.0.0 > 84, NA, pe_dx_age.0.0),
    stroke_dx_age.0.0 = ifelse(stroke_dx_age.0.0 < 18 | stroke_dx_age.0.0 > 83, NA, stroke_dx_age.0.0)
  )

# NA counts (ages)

cols_to_check <- c(
  "age_high_bp_dx.0.0",
  "age_diabetes_dx.0.0",
  "angina_dx_age.0.0",
  "mi_dx_age.0.0",
  "dvt_dx_age.0.0",
  "pe_dx_age.0.0",
  "stroke_dx_age.0.0"
)

na_summary <- data.frame(
  variable = cols_to_check,
  n_na = sapply(cols_to_check, function(v) sum(is.na(ukb_clean[[v]]))),
  pct_na = sapply(cols_to_check, function(v) mean(is.na(ukb_clean[[v]])) * 100)
)

na_summary

# # Value counts for all columns (including NAs)
# 
# for (v in names(ukb_clean)) {
#   cat("\n=====================================\n")
#   cat("Counts for:", v, "\n")
#   cat("=====================================\n")
#   
#   x <- ukb_clean[[v]]
#   
#   if (is.numeric(x) || is.integer(x)) {
#     print(summary(x))
#     cat("NA count:", sum(is.na(x)), "\n")
#     cat("NA pct:", mean(is.na(x)) * 100, "\n")
#   } else {
#     print(table(x, useNA = "always"))
#     cat("NA pct:", mean(is.na(x)) * 100, "\n")
#   }
# }

# Check value counts for the collapsed columns

collapsed_cols <- c(
  "heart_problems_collapsed",
  "meds_chol_bp_diabetes_hormone_collapsed",
  "meds_chol_bp_diabetes_collapsed"
  # ,
  # "n_previous_pregnancies_collapsed"
)

for (v in collapsed_cols) {
  cat("\n=====================================\n")
  cat("Counts for:", v, "\n")
  cat("=====================================\n")
  print(table(ukb_clean[[v]], useNA = "always"))
}

# Save RDS

saveRDS(ukb_clean, here("00_data", "01_processed", "medical_history_cleaned.rds"))
