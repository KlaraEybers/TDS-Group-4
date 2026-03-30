rm(list = ls())
# install.packages("here")

# Import packages ---------------------------------------------------------

library(here)
library(dplyr)

# Import data -------------------------------------------------------------

getwd()
here::here()

ukb_clean <- readRDS(here("00_data", "01_processed", "ukb_initial_cleaned.rds"))

# Extract relevant columns ------------------------------------------------

# Pattern to match early life variables (coding names)
pattern <- paste(c(
  "breastfed",
  "comparative_body_size_10",
  "maternal_smoking_birth",
  "birth_weight"
), collapse = "|")

matching_cols <- grep(pattern, colnames(ukb_clean), value = TRUE)

selected_data <- ukb_clean %>%
  select(eid, all_of(matching_cols))

# # Distribution of variables (JUST FOR ME TO SEE) -----------------------------
# 
# table(ukb_clean$breastfed.0.0, useNA = "ifany")
# table(ukb_clean$comparative_body_size_10.0.0, useNA = "ifany")
# table(ukb_clean$maternal_smoking_birth.0.0, useNA = "ifany")
# table(ukb_clean$birth_weight.0.0, useNA = "ifany")

# Early life cleaning -----------------------------------------------------

ukb_clean %>%
  select(
    breastfed.0.0,
    comparative_body_size_10.0.0,
    maternal_smoking_birth.0.0,
    birth_weight.0.0
  ) %>%
  str()

early_life_clean <- selected_data %>%
  mutate(
    # ----------------------------------------------------------
    # Breastfed as a baby (DNK/PNA -> NA)
    # ----------------------------------------------------------
    breastfed = case_when(
      breastfed.0.0 %in% c("Yes", "No") ~ breastfed.0.0,
      breastfed.0.0 %in% c("Do not know", "Prefer not to answer") ~ NA_character_,
      TRUE ~ NA_character_
    ),
    breastfed = factor(breastfed, levels = c("No", "Yes")),
    
    # ----------------------------------------------------------
    # Comparative body size at age 10 (DNK/PNA -> NA)
    # ----------------------------------------------------------
    comparative_body_size_10 = case_when(
      comparative_body_size_10.0.0 %in% c("Do not know", "Prefer not to answer") ~ NA_character_,
      TRUE ~ comparative_body_size_10.0.0
    ),
    # Ensure proper order of variables is maintained isntead of alphabetical
    comparative_body_size_10 = factor(comparative_body_size_10,
                                      levels = c("Thinner", "About average", "Plumper")),
    
    # ----------------------------------------------------------
    # Maternal smoking around birth (DNK/PNA -> NA)
    # ----------------------------------------------------------
    maternal_smoking_birth = case_when(
      maternal_smoking_birth.0.0 %in% c("Yes", "No") ~ maternal_smoking_birth.0.0,
      maternal_smoking_birth.0.0 %in% c("Do not know", "Prefer not to answer") ~ NA_character_,
      TRUE ~ NA_character_
    ),
    maternal_smoking_birth = factor(maternal_smoking_birth, levels = c("No", "Yes")),
    
    # ----------------------------------------------------------
    # Birth weight: Cap implausible values (DNK/PNA -> NA)
    # ----------------------------------------------------------
    birth_weight = case_when(
      birth_weight.0.0 %in% c("Do not know", "Prefer not to answer") ~ NA_real_,
      TRUE ~ suppressWarnings(as.numeric(as.character(birth_weight.0.0)))
    ),
    
    birth_weight = ifelse(
      !is.na(birth_weight) & (birth_weight < 0.5 | birth_weight > 7.0),
      NA,
      birth_weight
    )
  ) %>%
  select(eid, breastfed, comparative_body_size_10, maternal_smoking_birth, birth_weight)

# Save RDS ----------------------------------------------------------------

saveRDS(early_life_clean, here("00_data", "01_processed", "early_life_cleaned.rds"))

# Sanity Check ----------------------------------------------------------------

# table(early_life_clean$breastfed, useNA="ifany")
# table(early_life_clean$maternal_smoking_birth, useNA="ifany")
summary(early_life_clean$birth_weight)
mean(is.na(early_life_clean$birth_weight))
hist(
  early_life_clean$birth_weight,
  breaks = 40,
  main = "Distribution of birth weight",
  xlab = "Birth weight",
  col = "lightblue",
  border = "white"
)