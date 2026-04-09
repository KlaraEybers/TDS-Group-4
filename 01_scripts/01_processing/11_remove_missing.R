# ------------------------------------------------------------------------------
# Script: 11_remove_missing.R
# Purpose: Remove variables and participants with high levels of missing data.
#
# Inputs: 00_data/01_processed/ukb_final_merged.rds
#
# Outputs: 00_data/01_processed/ukb_final_cleaned.rds
#
# Upstream: 10_final_clean.R
#
# Downstream: 12_imputation.R
#
# Key steps:
#   - Calculate proportion of missing data for each variable
#   - Remove variables exceeding predefined missingness threshold
#   - Retain key outcome variables regardless of missingness
#   - Calculate proportion of missing data per participant
#   - Remove participants exceeding missingness threshold
#   - Report final dataset dimensions
#
# Notes:
#   - Variable threshold set at 45% missingness
#   - Participant threshold set at 30% missingness
#   - death_date and first_cvd_date are retained regardless of missingness
# ------------------------------------------------------------------------------

rm(list = ls())

# Import packages ---------------------------------------------------------

library(dplyr)
library(here)

# Import data -------------------------------------------------------------

df <- readRDS(here("00_data", "01_processed", "ukb_final_merged.rds"))


# Set missingness thresholds ----------------------------------------------

variable_threshold <- 0.45
participant_threshold <- 0.3


# Remove missingness above thresholds -------------------------------------

# Columns to protect from removal
protected <- c("death_date", "first_cvd_date")

# Remove high-missingness variables (except protected ones)
var_miss <- sapply(df, function(x) mean(is.na(x)))
high_miss_vars <- sort(var_miss[var_miss > variable_threshold & !names(var_miss) %in% protected], decreasing = TRUE)
cat("Variables with high missingness:\n")
print(high_miss_vars)
cat("\nCount:", length(high_miss_vars), "out of", ncol(df), "\n")


drop_vars <- names(var_miss[var_miss > variable_threshold & !names(var_miss) %in% protected])
df <- df[, !names(df) %in% drop_vars]

# Remove high-missingness participants
row_miss <- rowMeans(is.na(df))
drop_rows <- row_miss > participant_threshold
cat("\nDropping participants:", sum(drop_rows), "out of", nrow(df), "\n")

df <- df[!drop_rows, ]


# Final dimensions --------------------------------------------------------

cat("\nFinal dimensions:", nrow(df), "x", ncol(df), "\n")


# Export data -------------------------------------------------------------

saveRDS(df, here("00_data", "01_processed", "ukb_final_cleaned.rds"))


