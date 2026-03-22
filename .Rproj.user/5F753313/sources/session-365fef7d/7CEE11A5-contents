rm(list = ls())

# Import packages ---------------------------------------------------------

library(here)
library(dplyr)

# Import data -------------------------------------------------------------

ukb_clean <- readRDS(here("00_data", "01_processed", "ukb_initial_cleaned.rds"))

# Extract relevant columns ------------------------------------------------

pattern <- paste(c(
  "sleep_hrs",
  "insomnia",
  "snoring",
  "tobacco_exp_home",
  "tobacco_exp_out",
  "smoking_status",
  "pack_yrs",
  "met_min_week_all",
  "age_stop_smok",
  "age_at_recruitment"
), collapse = "|")

matching_cols <- grep(pattern, colnames(ukb_clean), value = TRUE)

selected_data <- ukb_clean %>%
  select(eid, all_of(matching_cols))

# Variable Distributions (JUST FOR ME TO LOOK AT) ------------------------------------------------------

summary(selected_data$sleep_hrs.0.0)
head(selected_data$sleep_hrs.0.0, 20)
min(selected_data$sleep_hrs.0.0)

sleep_raw <- selected_data$sleep_hrs.0.0

sleep_num <- suppressWarnings(as.numeric(sleep_raw))

# Remove UKB negative codes
sleep_num[sleep_num %in% c(-1, -3)] <- NA

summary(sleep_num)
hist(sleep_num)

table(selected_data$sleep_hrs.0.0 %in% c(-1, -3))
table(selected_data$sleep_hrs.0.0, useNA = "ifany")[c("-1","-3")]
class(selected_data$sleep_hrs.0.0)

table(selected_data$insomnia.0.0, useNA = "ifany")
table(selected_data$snoring.0.0, useNA = "ifany")
table(selected_data$smoking_status.0.0, useNA = "ifany")
summary(selected_data$pack_yrs.0.0)

py <- as.numeric(selected_data$pack_yrs.0.0)
head(py, 20)
range(py, na.rm = TRUE)

hist(py,
     breaks = 50,
     main = "Raw Pack Years",
     xlab = "Pack years")

table(selected_data$tobacco_exp_home.0.0, useNA = "ifany")

home_raw <- selected_data$tobacco_exp_home.0.0
home_num <- suppressWarnings(as.numeric(as.character(home_raw)))
home_tab <- table(home_num)
home_tab <- home_tab[order(as.numeric(names(home_tab)))]

barplot(home_tab,
        col = "steelblue",
        main = "Tobacco exposure at home (sorted)",
        xlab = "Hours",
        ylab = "Count")
table(selected_data$tobacco_exp_out.0.0, useNA = "ifany")

out_raw <- selected_data$tobacco_exp_out.0.0
out_num <- suppressWarnings(as.numeric(as.character(out_raw)))
out_tab <- table(out_num)
out_tab <- out_tab[order(as.numeric(names(out_tab)))]

barplot(out_tab,
        col = "steelblue",
        main = "Tobacco exposure outside (sorted)",
        xlab = "Hours",
        ylab = "Count")

table(selected_data$met_min_week_all.0.0 %in% c(-1, -3))
table(selected_data$met_min_week_all.0.0)[c("-1","-3")]

met <- selected_data$pack_yrs.0.0
head(met, 100)

hist(selected_data$met_min_week_all.0.0,
     main = "Raw MET-min/week",
     xlab = "MET-min/week")

mean(met, na.rm = TRUE)
sd(met, na.rm = TRUE)

# Lifestyle cleaning ------------------------------------------------------

lifestyle_clean <- selected_data %>%
  mutate(
    # Sleep duration (numeric + Short/Normal/Long)
    sleep_hrs = suppressWarnings(as.numeric(sleep_hrs.0.0)),
    sleep_hrs = ifelse(sleep_hrs %in% c(-1, -3), NA, sleep_hrs),
    sleep_hrs = ifelse(!is.na(sleep_hrs) & (sleep_hrs < 1 | sleep_hrs > 23), NA, sleep_hrs),
    
    sleep_cat = case_when(
      sleep_hrs <= 6 ~ "Short",
      sleep_hrs >= 7 & sleep_hrs <= 8 ~ "Normal",
      sleep_hrs >= 9 ~ "Long",
      TRUE ~ NA_character_
    ),
    sleep_cat = factor(sleep_cat, levels = c("Normal", "Short", "Long")),
    
    # Insomnia (DNK/PNA -> NA)
    insomnia = case_when(
      insomnia.0.0 %in% c("Do not know", "Prefer not to answer") ~ NA_character_,
      TRUE ~ as.character(insomnia.0.0)
    ),
    insomnia = factor(insomnia),
    
    # Snoring (DNK/PNA -> NA)
    snoring = case_when(
      snoring.0.0 %in% c("Do not know", "Prefer not to answer") ~ NA_character_,
      TRUE ~ as.character(snoring.0.0)
    ),
    snoring = factor(snoring),
    
    # Smoking (create refined version)
    # 1) Clean age stopped smoking into numeric 
    age_stop_smok_num = case_when(
      as.character(age_stop_smok.0.0) %in% c("-3", "-1", "Prefer not to answer", "Do not know", "") ~ NA_real_,
      TRUE ~ suppressWarnings(as.numeric(as.character(age_stop_smok.0.0)))
    ),
    age_stop_smok_num = ifelse(!is.na(age_stop_smok_num) & (age_stop_smok_num < 0 | age_stop_smok_num > 120), NA, age_stop_smok_num),
    
    # 2) Years since quit
    years_since_quit = if_else(
      smoking_status.0.0 == "Previous" &
        !is.na(age_stop_smok_num) &
        !is.na(age_at_recruitment.0.0) &
        age_stop_smok_num <= age_at_recruitment.0.0,
      age_at_recruitment.0.0 - age_stop_smok_num,
      NA_real_
    ),
    
    # 3) Refined smoking status (Former_Unknown -> NA)
    smoking_status_refined = case_when(
      smoking_status.0.0 %in% c("Prefer not to answer", "Do not know") ~ NA_character_,
      smoking_status.0.0 == "Never" ~ "Never",
      smoking_status.0.0 == "Current" ~ "Current",
      
      smoking_status.0.0 == "Previous" & !is.na(years_since_quit) & years_since_quit < 10 ~ "Former_Recent",
      smoking_status.0.0 == "Previous" & !is.na(years_since_quit) & years_since_quit >= 10 ~ "Former_Long",
      
      # Former smokers with missing quit time -> now NA
      smoking_status.0.0 == "Previous" & is.na(years_since_quit) ~ NA_character_,
      
      TRUE ~ NA_character_
    ),
    
    smoking_status_refined = factor(
      smoking_status_refined,
      levels = c("Never", "Former_Recent", "Former_Long", "Current")
    ),
    
    # ----------------------------------------------------------
    # Pack years
    # ----------------------------------------------------------
    pack_yrs = suppressWarnings(as.numeric(pack_yrs.0.0)),
    pack_yrs = ifelse(!is.na(pack_yrs) & pack_yrs < 0, NA, pack_yrs),
    pack_yrs = ifelse(!is.na(pack_yrs) & pack_yrs > 200, NA, pack_yrs),
    
    # ----------------------------------------------------------
    # Passive smoking hours (home/outside): -1/-3 -> NA
    # ----------------------------------------------------------
    tobacco_home_hrs = suppressWarnings(as.numeric(tobacco_exp_home.0.0)),
    tobacco_home_hrs = ifelse(tobacco_home_hrs %in% c(-1, -3), NA, tobacco_home_hrs),
    tobacco_home_hrs = ifelse(!is.na(tobacco_home_hrs) & tobacco_home_hrs < 0, NA, tobacco_home_hrs),
    
    tobacco_out_hrs = suppressWarnings(as.numeric(tobacco_exp_out.0.0)),
    tobacco_out_hrs = ifelse(tobacco_out_hrs %in% c(-1, -3), NA, tobacco_out_hrs),
    tobacco_out_hrs = ifelse(!is.na(tobacco_out_hrs) & tobacco_out_hrs < 0, NA, tobacco_out_hrs),
    
    # Total passive smoke hours:
    # - If both missing -> NA
    # - Otherwise sum available (treat missing component as 0)
    passive_total_hrs = ifelse(
      is.na(tobacco_home_hrs) & is.na(tobacco_out_hrs),
      NA,
      ifelse(is.na(tobacco_home_hrs), 0, tobacco_home_hrs) +
        ifelse(is.na(tobacco_out_hrs), 0, tobacco_out_hrs)
    ),
    
    # Flag if total was computed from partial information (one missing, one present)
    passive_partial_flag = ifelse(
      xor(is.na(tobacco_home_hrs), is.na(tobacco_out_hrs)),
      1, 0
    ),
    
    # Hard cap check: >168 hours/week is impossible -> NA
    passive_total_hrs = ifelse(
      !is.na(passive_total_hrs) & passive_total_hrs > 168,
      NA,
      passive_total_hrs
    ),
    
    # ----------------------------------------------------------
    # MET-min/week (keep as-is for now)
    # ----------------------------------------------------------
    met_min_week_all = suppressWarnings(as.numeric(met_min_week_all.0.0)),
    met_min_week_all = ifelse(!is.na(met_min_week_all) & met_min_week_all < 0, NA, met_min_week_all),
    met_min_week_all = ifelse(!is.na(met_min_week_all) & met_min_week_all > 20000, NA, met_min_week_all)
  ) %>%
  mutate(
    passive_cat = case_when(
      is.na(passive_total_hrs) ~ NA_character_,
      passive_total_hrs == 0 ~ "None",
      passive_total_hrs > 0 & passive_total_hrs <= 20 ~ "Low",
      passive_total_hrs > 20 ~ "High"
    ),
    passive_cat = factor(passive_cat, levels = c("None", "Low", "High"))
  ) %>%
  select(
    eid,
    sleep_hrs,
    sleep_cat,
    insomnia,
    snoring,
    smoking_status_refined,
    pack_yrs,
    passive_total_hrs,
    passive_cat,
    met_min_week_all
  )

saveRDS(lifestyle_clean, here("00_data", "01_processed", "lifestyle_cleaned.rds"))



# Sanity Check! ----------------------------------------------------------

str(lifestyle_clean)
colSums(is.na(lifestyle_clean))

table(lifestyle_clean$sleep_cat, useNA = "ifany")
table(lifestyle_clean$insomnia, useNA = "ifany")
table(lifestyle_clean$snoring, useNA = "ifany")

table(selected_data$smoking_status.0.0, useNA = "ifany")
sum(table(ukb_clean$age_stop_smok))
table(ukb_clean$age_stop_smok)
table(lifestyle_clean$smoking_status_refined, useNA = "ifany")

table(ukb_clean$age_stop_smok.0.0, useNA = "ifany")
table(
  ukb_clean$smoking_status,
  is.na(ukb_clean$age_stop_smok.0.0)
)
sum(is.na(ukb_clean$age_stop_smok.0.0))

table(ukb_clean$age_at_recruitment.0.0, useNA = "ifany")
table(lifestyle_clean$smoking_status_refined, useNA = "ifany")
table(lifestyle_clean$passive_cat, useNA = "ifany")
sum(lifestyle_clean$passive_total_hrs > 168, na.rm = TRUE)
hist(lifestyle_clean$passive_total_hrs,
     breaks = 100,
     main = "Passive Smoking Total Hours",
     xlab = "Total Hours (per week)",
     col = "lightblue",
     border = "white")

summary(lifestyle_clean$passive_total_hrs)
quantile(lifestyle_clean$passive_total_hrs, probs = c(0.5,0.75,0.9,0.95,0.99), na.rm=TRUE)
max(lifestyle_clean$passive_total_hrs, na.rm=TRUE)

summary(lifestyle_clean$met_min_week_all)
summary(lifestyle_clean$pack_yrs)
  