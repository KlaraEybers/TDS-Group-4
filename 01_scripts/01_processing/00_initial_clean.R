rm(list = ls())

# Import packages ---------------------------------------------------------

library(dplyr)
library(here)
library(stringr)
library(readxl)
library(tibble)
library(here)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

# Import data -------------------------------------------------------------

# Import UK Biobank data either from job script (real) or synthetic (path)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  ukb_path <- args[1]
} else {
  ukb_path <- here("extraction_and_recoding", "outputs", "ukb_recoded.rds")
}

ukb <- readRDS(ukb_path)
annot <- readRDS(here("extraction_and_recoding", "outputs", "annot.rds"))


# Make eid column ---------------------------------------------------------

if (!("eid" %in% names(ukb))) {
  ukb <- ukb %>% rownames_to_column("eid")
}


# Separating illness onset age columns --------------------------------------

# Identify all ill_onset_age columns

ill_onset_cols <- grep("^ill_onset_age", names(ukb), value = TRUE)

# Create separate dataset containing: eid + all ill_onset_age columns

ill_onset_data <- ukb %>%
  select(eid, all_of(ill_onset_cols))

# Remove ill_onset_age columns from main ukb dataset

ukb <- ukb %>%
  select(-all_of(ill_onset_cols))

# ukb now no longer contains ill_onset_age columns
# ill_onset_data contains them separately



# Remove non-baseline incidences ------------------------------------------

# ---------------------------------------------------------
# Keep:
#   - Columns that do NOT follow the ".i.j" UKB instance pattern
#   - Columns that follow the pattern AND have instance = 0
#
# Drop:
#   - Columns like variable.1.n, variable.2.n, etc.
# ---------------------------------------------------------

cn <- names(ukb)

# Match columns ending in ".i.j" (e.g., grip_strength_l.0.0, .1.0, .2.3, etc.)
is_instance_pattern <- str_detect(cn, "\\.[0-9]+\\.[0-9]+$")

# Extract the instance number (the first number after last dot group)
instance_number <- str_extract(cn, "(?<=\\.)[0-9]+(?=\\.[0-9]+$)")

# Keep columns that:
#  - are NOT instance-style columns
#  - OR are instance-style AND instance == "0"
keep_cols <- !is_instance_pattern | instance_number == "0"

# Subset dataset
ukb <- ukb[, keep_cols, drop = FALSE]

# ukb now retains:
#   variable.0.n columns
# and drops:
#   variable.1.n, variable.2.n, etc.


# Export initially cleaned RDS --------------------------------------------

saveRDS(ukb, here("00_data", "01_processed", "ukb_initial_cleaned.rds"))
