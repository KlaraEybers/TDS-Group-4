# ==============================================================================
# Variable Display Name Dictionary
# Source this file in any script: source(here("01_scripts", "00_variable_dictionary.R"))
# ==============================================================================

# This file is to standardise formatting of variable names across plots, tables, etc.
 
display_names <- c(
  # Confounders
  "age" = "Age",
  "sex" = "Sex",
  "ethnicity" = "Ethnicity",

  # Exposures: Sociodemographic
  "household_income" = "Household Income",
  "highest_qualification" = "Education",
  "employment" = "Employment",
  "housing" = "Housing",
  "household_size" = "Household Size",
  "has_vehicle" = "Owns a vehicle",

  # Exposures: Lifestyle
  "diet_score" = "Diet Score",
  "alcohol_consumption" = "Alcohol Consumption",
  "salt_added" = "Salt Intake",
  "sleep_hrs" = "Sleep Hours",
  "sleep_cat" = "Sleep Category",
  "insomnia" = "Insomnia",
  "smoking_status_refined" = "Smoking Status",
  "pack_yrs" = "Pack Years",
  "passive_total_hrs" = "Passive Smoking Hours",
  "passive_cat" = "Passive Smoking",
  "met_min_week_all" = "Physical Activity",

  # Exposures: Psychosocial
  "mh_sought_help" = "Mental Health Help",
  "satisfaction_mean" = "Satisfaction",
  "neuroticism_score" = "Neuroticism Score",

  # Exposures: Early life
  "breastfed" = "Breastfed",
  "comparative_body_size_10" = "Body Size Age 10",
  "maternal_smoking_birth" = "Maternal Smoking",
  "birth_weight" = "Birth Weight",

  # Exposures: Environment
  "no2_2010" = "NO2",
  "pm2.5_2010" = "PM2.5",
  "pm10_2010" = "PM10",
  "urban_rural_3lvl" = "Urban/Rural",

  # Anthropometric
  "sbp" = "Systolic BP",
  "dbp" = "Diastolic BP",
  "grip_strength" = "Grip Strength",
  "whr" = "Waist-Hip Ratio",
  "heel_bone_density" = "Bone Density",
  "body_fat_mass" = "Body Fat Mass",

  # Biomarkers
  "wbc" = "WBC",
  "platelets" = "Platelets",
  "albumin" = "Albumin",
  "alp" = "ALP",
  "alt" = "ALT",
  "apoa1" = "ApoA1",
  "apob" = "ApoB",
  "ast" = "AST",
  "bilirubin" = "Bilirubin",
  "total_cholesterol" = "Total Cholesterol",
  "creatinine" = "Creatinine",
  "crp" = "CRP",
  "cystatin_c" = "Cystatin C",
  "ggt" = "GGT",
  "hba1c" = "HbA1c",
  "hdl_cholesterol" = "HDL Cholesterol",
  "ldl_direct" = "LDL Direct",
  "lpa" = "Lp(a)",
  "triglycerides" = "Triglycerides",
  "neutrophill_count" = "Neutrophils",
  "lymphocyte_count" = "Lymphocytes",
  "monocyte_count" = "Monocytes",

  # Biomarkers: Medical history
  # "long_ill" = "Long-standing Illness",
  # "diabetes_doctor_diagnosed" = "Diabetes Diagnosed",
  # "serious_dx" = "Serious Diagnosis",
  # "age_high_bp_dx" = "Age at High BP Diagnosis",
  # "age_diabetes_dx" = "Age at Diabetes",
  # "angina_dx_age" = "Age at Angina",
  # "mi_dx_age" = "Age at MI",
  # "dvt_dx_age" = "Age at DVT",
  # "pe_dx_age" = "Age at PE",
  # "stroke_dx_age" = "Age at Stroke",
  # "heart_problems_collapsed" = "Heart Problems",
  # "meds_chol_bp_diabetes_hormone_collapsed" = "Chol/BP/Diabetes/Hormone Meds",
  # "meds_chol_bp_diabetes_collapsed" = "Chol/BP/Diabetes Meds",
  "early_insulin" = "Early Insulin",
  "snoring" = "Snoring",
  "recruitment_date" = "Recruitment Date",
  "death_date" = "Death Date",
  "pregnancy" = "Pregnancy",
  "visits_cat" = "Visits frequency",
  
  # Outcome
  "has_cvd" = "CVD"
  
)

# Helper function: look up display name, return original if not found
get_display_name <- function(var) {
  sapply(var, function(v) {
    
    # 1) Exact match
    if (v %in% names(display_names)) {
      return(display_names[[v]])
    }
    
    # 2) Handle dummy variables (prefix match)
    matched_base <- names(display_names)[startsWith(v, names(display_names))]
    
    if (length(matched_base) > 0) {
      # choose longest match
      base_var <- matched_base[which.max(nchar(matched_base))]
      
      level_part <- sub(paste0("^", base_var), "", v)
      
      # clean level name
      level_part <- gsub("^[_\\.]", "", level_part)
      level_part <- gsub("_", " ", level_part)
      level_part <- gsub("\\.", " ", level_part)
      level_part <- gsub("([a-z])([A-Z])", "\\1 \\2", level_part)
      level_part <- trimws(level_part)
      
      if (level_part == "") {
        return(display_names[[base_var]])
      } else {
        return(paste0(display_names[[base_var]], ": ", level_part))
      }
    }
    
    # 3) fallback
    v <- gsub("_", " ", v)
    v <- gsub("\\.", " ", v)
    trimws(v)
    
  }, USE.NAMES = FALSE)
}

# Old Helper function: look up display name, return original if not found
# get_display_name <- function(var) {
#   ifelse(var %in% names(display_names), display_names[var], var)
# }


# ==============================================================================
# Variable → Domain mapping (aligned with dictionary groups)
# ==============================================================================

var_domains <- c(
  # Sociodemographic
  "age" = "Socio-demographic",
  "sex" = "Socio-demographic",
  "ethnicity" = "Socio-demographic",
  "household_income" = "Socio-demographic",
  "highest_qualification" = "Socio-demographic",
  "employment" = "Socio-demographic",
  "housing" = "Socio-demographic",
  "household_size" = "Socio-demographic",
  "has_vehicle" = "Socio-demographic",

  # Alcohol and diet 
  "diet_score" = "Diet",
  "alcohol_consumption" = "Diet",
  "salt_added" = "Diet",
  
  # Lifestyle
  "sleep_hrs" = "Lifestyle",
  "sleep_cat" = "Lifestyle",
  "insomnia" = "Lifestyle",
  "smoking_status_refined" = "Lifestyle",
  "pack_yrs" = "Lifestyle",
  "passive_total_hrs" = "Lifestyle",
  "passive_cat" = "Lifestyle",
  "met_min_week_all" = "Lifestyle",
  "snoring" = "Lifestyle",

  # Mental Health
  "mh_sought_help" = "Mental health",
  "satisfaction_mean" = "Mental health",
  "neuroticism_score" = "Mental health",
  "visits_cat" = "Mental health",

  # Early life
  "breastfed" = "Early life factors",
  "comparative_body_size_10" = "Early life factors",
  "maternal_smoking_birth" = "Early life factors",
  "birth_weight" = "Early life factors",

  # Environment
  "no2_2010" = "Environment",
  "pm2.5_2010" = "Environment",
  "pm10_2010" = "Environment",
  "urban_rural_3lvl" = "Environment",

  # Anthropometric
  "sbp" = "Anthropometrics",
  "dbp" = "Anthropometrics",
  "grip_strength" = "Anthropometrics",
  "whr" = "Anthropometrics",
  "heel_bone_density" = "Anthropometrics",
  "body_fat_mass" = "Anthropometrics",

  # Blood
  "wbc" = "Biological samples",
  "platelets" = "Biological samples",
  "albumin" = "Biological samples",
  "alp" = "Biological samples",
  "alt" = "Biological samples",
  "apoa1" = "Biological samples",
  "apob" = "Biological samples",
  "ast" = "Biological samples",
  "bilirubin" = "Biological samples",
  "total_cholesterol" = "Biological samples",
  "creatinine" = "Biological samples",
  "crp" = "Biological samples",
  "cystatin_c" = "Biological samples",
  "ggt" = "Biological samples",
  "hba1c" = "Biological samples",
  "hdl_cholesterol" = "Biological samples",
  "ldl_direct" = "Biological samples",
  "lpa" = "Biological samples",
  "triglycerides" = "Biological samples",
  "neutrophill_count" = "Biological samples",
  "lymphocyte_count" = "Biological samples",
  "monocyte_count" = "Biological samples",

  # Medical history
  # "long_ill" = "Medical history",
  # "diabetes_doctor_diagnosed" = "Medical history",
  # "serious_dx" = "Medical history",
  # "age_high_bp_dx" = "Medical history",
  # "age_diabetes_dx" = "Medical history",
  # "angina_dx_age" = "Medical history",
  # "mi_dx_age" = "Medical history",
  # "dvt_dx_age" = "Medical history",
  # "pe_dx_age" = "Medical history",
  # "stroke_dx_age" = "Medical history",
  # "heart_problems_collapsed" = "Medical history",
  # "meds_chol_bp_diabetes_hormone_collapsed" = "Medical history",
  # "meds_chol_bp_diabetes_collapsed" = "Medical history",
  "early_insulin" = "Medical history",
  "snoring" = "Lifestyle",
  # "recruitment_date" = "Medical history",
  # "death_date" = "Medical history",
  "pregnancy" = "Medical history",
  
  # Outcome
  "has_cvd" = "Outcome"
)

get_domain_from_variable <- function(var) {
  unname(ifelse(var %in% names(var_domains), var_domains[var], "Other"))
}

