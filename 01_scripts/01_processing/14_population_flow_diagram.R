# ==============================================================================
# POPULATION FLOW DIAGRAM CALCULATOR
# ==============================================================================
rm(list = ls())

# 1. Import Packages -----------------------------------------------------------
library(here)
library(dplyr)
library(tibble)
library(tidyr)

# Read in config file
source(here("01_scripts", "00_variable_config.R"))

cat("=======================================================\n")
cat("Extracting Population Flow Diagram Values...\n")
cat("=======================================================\n")

# 2. DYNAMICALLY CALCULATE BASELINE EXCLUSIONS (CVD & Pregnancy)
cat("Scanning baseline data to calculate exact exclusions...\n")

# Load sociodemographic to get the initial cohort size
socio <- readRDS(here("00_data", "01_processed", "sociodemographic_cleaned.rds"))
n_initial <- nrow(socio)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  cvd_events_path <- args[1]
} else {
  cvd_events_path <- here("00_data", "00_extracted", "cvd_events.rds")
}

cvd_events <- readRDS(cvd_events_path)

# Load and process CVD events
cvd <- cvd_events %>%
  mutate(date = as.Date(date)) %>%
  group_by(eid) %>% arrange(date) %>% slice(1) %>% ungroup() %>%
  rename(first_cvd_date = date) %>% select(eid, first_cvd_date)

proc_dir <- here("00_data", "01_processed")
domain_files <- list.files(proc_dir, pattern = "_cleaned\\.rds$", full.names = TRUE)
preg_data <- NULL
for (f in domain_files) {
  temp_cols <- colnames(readRDS(f))
  if ("pregnancy.0.0" %in% temp_cols) {
    preg_data <- readRDS(f) %>% select(eid, `pregnancy.0.0`)
    break
  }
}

# Merge to simulate exclusion logic
df_calc <- socio %>% select(eid, recruitment_date) %>%
  left_join(cvd, by = "eid") %>%
  left_join(preg_data, by = "eid")

# A. Calculate Prevalent CVD drop
cvd_drop_mask <- !is.na(df_calc$first_cvd_date) & df_calc$first_cvd_date < df_calc$recruitment_date
n_cvd_dropped <- sum(cvd_drop_mask, na.rm = TRUE)

# B. Calculate CVD-free cohort
df_no_cvd <- df_calc[!cvd_drop_mask, ]
n_cvd_free <- nrow(df_no_cvd)

# C. Calculate Pregnant women drop
preg_drop_mask <- !is.na(df_no_cvd$`pregnancy.0.0`) & df_no_cvd$`pregnancy.0.0` == "Yes"
n_preg_dropped <- sum(preg_drop_mask, na.rm = TRUE)

n_merged <- n_cvd_free - n_preg_dropped

# 3. DYNAMICALLY CALCULATE MISSINGNESS AND CONFOUNDERS
cat("Calculating missingness exclusions directly from final_cleaned.rds...\n")
dat_cleaned <- readRDS(here("00_data", "01_processed", "ukb_final_cleaned.rds"))
n_cleaned <- nrow(dat_cleaned)

n_missing_dropped <- n_merged - n_cleaned

# Simulate dropping missing confounders (Age, Sex, Ethnicity)
dat_full <- dat_cleaned %>%
  drop_na(all_of(intersect(confounder_vars, names(dat_cleaned))))

n_full <- nrow(dat_full)
n_confounders_dropped <- n_cleaned - n_full

events_full <- sum(dat_full$has_cvd %in% c("Event", 1), na.rm = TRUE)
no_events_full <- sum(dat_full$has_cvd %in% c("No Event", 0), na.rm = TRUE)

n_selection <- round(n_full * 0.5)
events_sel <- round(events_full * 0.5)
no_events_sel <- n_selection - events_sel

n_prediction <- n_full - n_selection
events_pred <- events_full - events_sel
no_events_pred <- n_prediction - events_pred

n_train <- round(n_prediction * 0.75)
n_test <- n_prediction - n_train

# 5. BUILD THE FINAL TABLE
pop_flow <- tibble(
  Diagram_Box = c(
    "1. UK Biobank participants",
    "2. Excluded: prevalent CVD at baseline",
    "3. CVD-free at baseline",
    "4. Excluded: pregnant women at baseline",
    "5. Excluded: high missingness variables",
    paste0("6. Excluded: missing confounders (", 
           paste(confounder_vars, collapse = ", "), ")"),
    "7. Full analysis cohort",
    "   -> CVD events (Full Cohort Corrected)",
    "   -> No event (Full Cohort Corrected)",
    "8. Univariate & descriptive analysis cohort",
    "9. Variable selection cohort (50%)",
    "   -> CVD events (Selection)",
    "   -> No event (Selection)",
    "10. CVD prediction cohort (50%)",
    "   -> CVD events (Prediction)",
    "   -> No event (Prediction)",
    "11. Training Set (75% split)",
    "12. Test Set (25% split)"
  ),
  Actual_Number_N = c(
    n_initial,
    n_cvd_dropped,
    n_cvd_free,
    n_preg_dropped,      
    n_missing_dropped,
    n_confounders_dropped, 
    n_full,
    events_full,
    no_events_full,
    n_full,
    n_selection,
    events_sel,
    no_events_sel,
    n_prediction,
    events_pred,
    no_events_pred,
    n_train,
    n_test
  )
)

# 6. EXPORT & PRINT
out_dir <- here("02_results", "01_tables", "00_flow_diagrams")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write.csv(pop_flow, here("02_results", "01_tables", "00_flow_diagrams", "population_flow_numbers.csv"), row.names = FALSE)

print(pop_flow, n = 25)
cat("\n=======================================================\n")
cat("SUCCESS! Updated flow diagram numbers exported to:\n", out_dir, "\n")
cat("=======================================================\n")