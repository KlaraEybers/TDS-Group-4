rm(list = ls())

# Import packages ---------------------------------------------------------

library(here)
library(dplyr)


# Import data -------------------------------------------------------------

ukb_clean <- readRDS(here("00_data", "01_processed", "ukb_initial_cleaned.rds"))

# Extract relevant columns ------------------------------------------------

# Pattern to match diet and alcohol variables
pattern <- paste(c(
  "alcohol_intake_freq", "former_alc", "alcohol_status",
  "cooked_veg", "fresh_fruit", "dried_fruit",  "salad_raw_veg",
  "coffee", "tea", "oily_fish", "non_oily_fish",
  "processed_meat", "beef", "lamb_mutton", "pork",
  "bread_slices_week", "bread_type", "cereaL_bowls_week", "cereal_type",
  "salt_added", "water_glasses_week"
), collapse = "|")


# Extract all columns that match any of these patterns
matching_cols <- grep(pattern, colnames(ukb_clean), value = TRUE)

# Select relevant columnns
selected_data <- ukb_clean %>% select(eid, any_of(matching_cols))


# Diet Score --------------------------------------------------------------

# Reference: https://www.mdpi.com/2072-6643/15/6/1344#app1-nutrients-15-01344
# Diet score (0-7): 7 components (fruit, veg, fish, processed meat, red meat, 
# whole grains, refined grains)
# Scored 1 if meeting threshold, 0 otherwise.
# For 7 composite variables, if one component is NA, full composite set to NA
# Requires ≥5/7 components; proportionally scaled to 0-7 if components missing
# If less than 5, set to NA and imputed later

colnames(selected_data)

diet_score_variables <- c("cooked_veg.0.0","salad_raw_veg.0.0", "fresh_fruit.0.0", 
                          "dried_fruit.0.0", "oily_fish.0.0", "non_oily_fish.0.0",
                          "processed_meat.0.0", "beef.0.0", "lamb_mutton.0.0",         
                          "pork.0.0", "bread_slices_week.0.0", "bread_type.0.0",
                          "cereaL_bowls_week.0.0", "cereal_type.0.0")


# ============================================================================
# DIET SCORE CLEANING (7 components, excluding salt)
# ============================================================================

diet_clean <- selected_data %>%
  mutate(
    # -------------------------------------------------------------------------
    # 0. UPFRONT NA CLEANING - Convert all "Prefer not to answer"/"Do not know" to NA
    # -------------------------------------------------------------------------
    across(
      any_of(diet_score_variables),
      ~case_when(
        . %in% c("Do not know", "Prefer not to answer") ~ NA_character_,
        TRUE ~ .
      )
    ),
    
    # -------------------------------------------------------------------------
    # 1. CONVERT TO NUMERIC SERVINGS
    # -------------------------------------------------------------------------
    
    # Vegetables (servings/day)
    cooked_veg_daily = case_when(
      cooked_veg.0.0 == "Less than one" ~ 0.5,
      TRUE ~ suppressWarnings(as.numeric(cooked_veg.0.0))
    ),
    
    salad_raw_veg_daily = case_when(
      salad_raw_veg.0.0 == "Less than one" ~ 0.5,
      TRUE ~ suppressWarnings(as.numeric(salad_raw_veg.0.0))
    ),
    
    # Fresh fruit (servings/day)
    fresh_fruit_daily = case_when(
      fresh_fruit.0.0 == "Less than one" ~ 0.5,
      TRUE ~ suppressWarnings(as.numeric(fresh_fruit.0.0))
    ),
    
    # Dried fruit (servings/day) - 2 pieces = 1 serving
    dried_fruit_daily = case_when(
      dried_fruit.0.0 == "Less than one" ~ 0.5,
      TRUE ~ suppressWarnings(as.numeric(dried_fruit.0.0))
    ) / 2,  # Divide by 2 because 2 pieces = 1 serving
    
    # Bread slices (per week - convert to daily)
    bread_slices_daily = case_when(
      bread_slices_week.0.0 == "Less than one" ~ 0.5 / 7,
      TRUE ~ suppressWarnings(as.numeric(bread_slices_week.0.0)) / 7
    ),
    
    # Cereal bowls (per week - convert to daily)
    cereal_bowls_daily = case_when(
      cereaL_bowls_week.0.0 == "Less than one" ~ 0.5 / 7,
      TRUE ~ suppressWarnings(as.numeric(cereaL_bowls_week.0.0)) / 7
    ),
    
    # Beef (convert frequency to servings/week)
    beef_weekly = case_when(
      beef.0.0 == "Never" ~ 0,
      beef.0.0 == "Less than once a week" ~ 0.5,
      beef.0.0 == "Once a week" ~ 1,
      beef.0.0 == "2-4 times a week" ~ 3,
      beef.0.0 == "5-6 times a week" ~ 5.5,
      beef.0.0 == "Once or more daily" ~ 7,
      TRUE ~ NA_real_
    ),
    
    # Lamb/mutton (convert frequency to servings/week)
    lamb_mutton_weekly = case_when(
      lamb_mutton.0.0 == "Never" ~ 0,
      lamb_mutton.0.0 == "Less than once a week" ~ 0.5,
      lamb_mutton.0.0 == "Once a week" ~ 1,
      lamb_mutton.0.0 == "2-4 times a week" ~ 3,
      lamb_mutton.0.0 == "5-6 times a week" ~ 5.5,
      lamb_mutton.0.0 == "Once or more daily" ~ 7,
      TRUE ~ NA_real_
    ),
    
    # Pork (convert frequency to servings/week)
    pork_weekly = case_when(
      pork.0.0 == "Never" ~ 0,
      pork.0.0 == "Less than once a week" ~ 0.5,
      pork.0.0 == "Once a week" ~ 1,
      pork.0.0 == "2-4 times a week" ~ 3,
      pork.0.0 == "5-6 times a week" ~ 5.5,
      pork.0.0 == "Once or more daily" ~ 7,
      TRUE ~ NA_real_
    ),
    
    # Oily fish (convert frequency to servings/week)
    oily_fish_weekly = case_when(
      oily_fish.0.0 == "Never" ~ 0,
      oily_fish.0.0 == "Less than once a week" ~ 0.5,
      oily_fish.0.0 == "Once a week" ~ 1,
      oily_fish.0.0 == "2-4 times a week" ~ 3,
      oily_fish.0.0 == "5-6 times a week" ~ 5.5,
      oily_fish.0.0 == "Once or more daily" ~ 7,
      TRUE ~ NA_real_
    ),
    
    # Non-oily fish (convert frequency to servings/week)
    non_oily_fish_weekly = case_when(
      non_oily_fish.0.0 == "Never" ~ 0,
      non_oily_fish.0.0 == "Less than once a week" ~ 0.5,
      non_oily_fish.0.0 == "Once a week" ~ 1,
      non_oily_fish.0.0 == "2-4 times a week" ~ 3,
      non_oily_fish.0.0 == "5-6 times a week" ~ 5.5,
      non_oily_fish.0.0 == "Once or more daily" ~ 7,
      TRUE ~ NA_real_
    ),
    
    # Processed meat (convert frequency to servings/week)
    processed_meat_weekly = case_when(
      processed_meat.0.0 == "Never" ~ 0,
      processed_meat.0.0 == "Less than once a week" ~ 0.5,
      processed_meat.0.0 == "Once a week" ~ 1,
      processed_meat.0.0 == "2-4 times a week" ~ 3,
      processed_meat.0.0 == "5-6 times a week" ~ 5.5,
      processed_meat.0.0 == "Once or more daily" ~ 7,
      TRUE ~ NA_real_
    ),
    
    # -------------------------------------------------------------------------
    # 2. CREATE COMPOSITE FOOD GROUPS
    # -------------------------------------------------------------------------
    
    # Total vegetables (servings/day) - need BOTH components, otherwise NA
    total_veg_daily = case_when(
      !is.na(cooked_veg_daily) & !is.na(salad_raw_veg_daily) ~ 
        (cooked_veg_daily / 2) + (salad_raw_veg_daily / 2),  # 2 tablespoons = 1 serving
      TRUE ~ NA_real_
    ),
    
    # Total fruit (servings/day) - need BOTH fresh and dried
    total_fruit_daily = case_when(
      !is.na(fresh_fruit_daily) & !is.na(dried_fruit_daily) ~ 
        fresh_fruit_daily + (dried_fruit_daily / 2), # 2 dried fruit = 1 serving
      TRUE ~ NA_real_
    ),
    
    # Red meat (servings/week) - need ALL THREE meats
    red_meat_weekly = case_when(
      !is.na(beef_weekly) & !is.na(lamb_mutton_weekly) & !is.na(pork_weekly) ~
        beef_weekly + lamb_mutton_weekly + pork_weekly,
      TRUE ~ NA_real_
    ),
    
    # Total fish (servings/week) - need BOTH oily and non-oily
    total_fish_weekly = case_when(
      !is.na(oily_fish_weekly) & !is.na(non_oily_fish_weekly) ~
        oily_fish_weekly + non_oily_fish_weekly,
      TRUE ~ NA_real_
    ),
    
    # Whole grains (servings/day) - need BOTH bread and cereal data
    whole_grain_bread_daily = case_when(
      bread_type.0.0 == "Wholemeal or wholegrain" ~ bread_slices_daily,
      bread_type.0.0 %in% c("Do not know", "Prefer not to answer") ~ NA_real_,
      TRUE ~ 0
    ),
    
    whole_grain_cereal_daily = case_when(
      cereal_type.0.0 %in% c("Bran cereal (e.g. All Bran, Branflakes)",
                             "Oat cereal (e.g. Ready Brek, porridge)",
                             "Muesli") ~ cereal_bowls_daily,
      cereal_type.0.0 %in% c("Do not know", "Prefer not to answer") ~ NA_real_,
      TRUE ~ 0
    ),
    
    whole_grains_daily = case_when(
      !is.na(whole_grain_bread_daily) & !is.na(whole_grain_cereal_daily) ~
        whole_grain_bread_daily + whole_grain_cereal_daily,
      TRUE ~ NA_real_
    ),
    
    # Refined grains (servings/day) - need BOTH bread and cereal data
    refined_grain_bread_daily = case_when(
      bread_type.0.0 %in% c("White", "Brown", "Other type of bread") ~ bread_slices_daily,
      bread_type.0.0 %in% c("Do not know", "Prefer not to answer") ~ NA_real_,
      TRUE ~ 0
    ),
    
    refined_grain_cereal_daily = case_when(
      cereal_type.0.0 %in% c("Biscuit cereal (e.g. Weetabix)",
                             "Other (e.g. Cornflakes, Frosties)") ~ cereal_bowls_daily,
      cereal_type.0.0 %in% c("Do not know", "Prefer not to answer") ~ NA_real_,
      TRUE ~ 0
    ),
    
    refined_grains_daily = case_when(
      !is.na(refined_grain_bread_daily) & !is.na(refined_grain_cereal_daily) ~
        refined_grain_bread_daily + refined_grain_cereal_daily,
      TRUE ~ NA_real_
    ),
    
    # -------------------------------------------------------------------------
    # 3. CREATE DIET SCORE (0-7 POINTS)
    # -------------------------------------------------------------------------
    
    # Individual score components
    score_fruit = case_when(
      is.na(total_fruit_daily) ~ NA_real_,
      total_fruit_daily >= 4 ~ 1,
      TRUE ~ 0
    ),
    score_veg = case_when(
      is.na(total_veg_daily) ~ NA_real_,
      total_veg_daily >= 4 ~ 1,
      TRUE ~ 0
    ),
    score_fish = case_when(
      is.na(total_fish_weekly) ~ NA_real_,
      total_fish_weekly >= 2 ~ 1,
      TRUE ~ 0
    ),
    score_processed_meat = case_when(
      is.na(processed_meat_weekly) ~ NA_real_,
      processed_meat_weekly <= 1 ~ 1,
      TRUE ~ 0
    ),
    score_red_meat = case_when(
      is.na(red_meat_weekly) ~ NA_real_,
      red_meat_weekly <= 1.5 ~ 1,
      TRUE ~ 0
    ),
    score_whole_grains = case_when(
      is.na(whole_grains_daily) ~ NA_real_,
      whole_grains_daily >= 3 ~ 1,
      TRUE ~ 0
    ),
    score_refined_grains = case_when(
      is.na(refined_grains_daily) ~ NA_real_,
      refined_grains_daily <= 1.5 ~ 1,
      TRUE ~ 0
    ),
    
    # Count non-missing components
    n_components = (!is.na(score_fruit)) + (!is.na(score_veg)) + 
      (!is.na(score_fish)) + (!is.na(score_processed_meat)) +
      (!is.na(score_red_meat)) + (!is.na(score_whole_grains)) + 
      (!is.na(score_refined_grains)),
    
    # Raw score (sum of available components)
    diet_score_raw = rowSums(cbind(score_fruit, score_veg, score_fish, 
                                   score_processed_meat, score_red_meat, 
                                   score_whole_grains, score_refined_grains), 
                             na.rm = TRUE),
    
    # Final diet score: scale to 0-7 if at least 5 components available
    diet_score = case_when(
      n_components >= 5 ~ (diet_score_raw / n_components) * 7,  # Scale to 0-7
      TRUE ~ NA_real_
    )
    
  )


# Remove diet processing variables ----------------------------------------

diet_clean <- diet_clean %>%
  select(
    eid,
    diet_score,
    salt_added.0.0,
    alcohol_intake_freq.0.0,
    former_alc.0.0,
    alcohol_status.0.0,
    salt_added.0.0
  )


# Alcohol -----------------------------------------------------------------

# Alcohol: 5 categories based on frequency
# Never, Former, Social, Moderate, Daily
# Reference: https://www.thelancet.com/journals/eclinm/article/PIIS2589-5370(20)30402-8/fulltext
# Hierarchy: (1) Never drinker status, (2) Former drinker status, 
# (3) Current frequency (special occasions/1-3x month = social; 1-2x/week = moderate; 3+x/week = daily)
# NA only if all three source variables are "Prefer not to answer"

# Create alcohol consumption variable (revised hierarchy)
diet_clean <- diet_clean %>%
  mutate(
    alcohol_consumption = case_when(
      # PRIORITY 1: Current drinking frequency
      alcohol_intake_freq.0.0 %in% c("Special occasions only", "One to three times a month") ~ 
        "Social drinker",
      alcohol_intake_freq.0.0 == "Once or twice a week" ~ 
        "Moderate drinker",
      alcohol_intake_freq.0.0 %in% c("Three or four times a week", "Daily or almost daily") ~ 
        "Daily drinker",
      
      # PRIORITY 2: If frequency = "Never", check former_alc to distinguish never vs former
      alcohol_intake_freq.0.0 == "Never" & former_alc.0.0 == "Yes" ~ "Former drinker",
      alcohol_intake_freq.0.0 == "Never" & former_alc.0.0 == "No" ~ "Never drinker",
      alcohol_intake_freq.0.0 == "Never" & is.na(former_alc.0.0) ~ "Never drinker",  # Default to never if former_alc missing
      
      # Handle prefer not to answer
      alcohol_intake_freq.0.0 == "Prefer not to answer" & 
        former_alc.0.0 == "Prefer not to answer" ~ NA_character_,
      
      TRUE ~ NA_character_
    ),
    
    alcohol_consumption = factor(alcohol_consumption, 
                                 levels = c("Never drinker", "Former drinker", 
                                            "Social drinker", "Moderate drinker", 
                                            "Daily drinker"))
  )

# Drop original alcohol columns
diet_clean <- diet_clean %>%
  select(-alcohol_intake_freq.0.0, -former_alc.0.0, -alcohol_status.0.0)



# Convert salt intake to factor -------------------------------------------

diet_clean <- diet_clean %>%
  mutate(
    # Convert "Prefer not to answer" to NA, then factor
    salt_added = case_when(
      salt_added.0.0 == "Prefer not to answer" ~ NA_character_,
      TRUE ~ salt_added.0.0
    ),
    salt_added = factor(salt_added, 
                        levels = c("Never/rarely", "Sometimes", "Usually", "Always"),
                        ordered = TRUE)
  ) %>%
  select(-salt_added.0.0)

# Save RDS ----------------------------------------------------------------

saveRDS(diet_clean, here("00_data", "01_processed", "diet_cleaned.rds"))
