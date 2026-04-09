# Note for automation: some variables may not have category mapping. Need to automate this. 

# ============================================================
# Missingness Analysis
# ============================================================
# Investigates missing data across 77 variables in the UK Biobank dataset.
# Analysis proceeds in three stages:
#   1. Variable-level missingness (bar chart + histogram)
#   2. Category-level missingness (boxplot + summary bar chart)
#   3. Structural missingness screen (are high-NA variables explained by gatekeepers?)
# ============================================================

# ============================================================
# 1. Setup
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(ggplot2)
})
set.seed(1)

# ── Paths ─────────────────────────────────────────────────────────────────────
data_path <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group4/00_data/01_processed/ukb_final_merged.rds"
out_dir   <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group4/02_results/02_reports/missingness"
category_map_path <- file.path(out_dir, "variable_category.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Shared palette ────────────────────────────────────────────────────────────
col_low     <- "#a8d5e2"   # light blue  (low missingness)
col_high    <- "#c0392b"   # red         (high missingness / threshold line)
col_fill    <- "#3a6ea5"   # steel blue  (histogram / neutral fill)
col_point   <- "#2c3e50"   # dark navy   (jitter points)
col_box     <- "#ecf0f1"   # light grey  (boxplot fill)

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_report <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      # Titles
      plot.title         = element_text(face = "bold", size = base_size + 2,
                                        colour = "#2c3e50"),
      plot.subtitle      = element_text(size = base_size, colour = "#555555",
                                        margin = margin(b = 8)),
      # Axes
      axis.title         = element_text(face = "bold", colour = "#2c3e50"),
      axis.text          = element_text(colour = "#444444"),
      # Grid: major x only (horizontal reference lines)
      panel.grid.major.x = element_line(colour = "#e0e0e0", linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      # Legend
      legend.title       = element_text(face = "bold"),
      # Margins
      plot.margin        = margin(10, 15, 10, 12)
    )
}

# ── Save helper ───────────────────────────────────────────────────────────────
save_plot <- function(p, filename, w = 7, h = 5, dpi = 300) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, filename)
  message("Saving plot to: ", out_file)
  ggplot2::ggsave(
    filename = out_file, plot = p,
    width = w, height = h, dpi = dpi,
    units = "in", device = "png", bg = "white"
  )
  if (file.exists(out_file)) message("Saved: ", out_file) else stop("Not saved: ", out_file)
}



# ============================================================
# 2. Load and clean data
# ============================================================

ukb <- readRDS(data_path)

# Ensure eid exists as a column (some versions store it as rownames)
if (!("eid" %in% names(ukb))) {
  ukb <- ukb %>% tibble::rownames_to_column("eid")
}

# Recode non-response strings as NA so they are treated as missing
nonresponse_levels <- c("Prefer not to answer", "Do not know", "Unsure")

ukb <- ukb %>%
  mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
  mutate(across(where(is.character),
                ~ if_else(.x %in% nonresponse_levels, NA_character_, .x)))

# Remove eid and outcome variables before missingness calculations.
# death_date and first_cvd_date are excluded because near-complete missingness
# is expected by design (only events that occurred are recorded).
ukb_no_id <- ukb %>%
  select(-eid) %>%
  select(-any_of(c("death_date", "first_cvd_date")))

n_total <- nrow(ukb_no_id)   # number of participants
n_vars  <- ncol(ukb_no_id)   # number of variables after exclusions


# ============================================================
# 3. Variable-level missingness
# ============================================================

# Compute count and percentage missing for each variable, highest first
missing_col <- tibble(
  variable    = names(ukb_no_id),
  n_missing   = map_int(ukb_no_id, ~ sum(is.na(.x))),
  pct_missing = 100 * n_missing / n_total
) %>%
  arrange(desc(pct_missing))

write_csv(missing_col, file.path(out_dir, "missingness_per_variable.csv"))

# ── Figure 1: Horizontal bar chart, variables ordered by % missing ────────────
p_var_bar <- ggplot(
  missing_col,
  aes(x = pct_missing, y = reorder(variable, pct_missing))
) +
  geom_col(aes(fill = pct_missing), width = 0.7) +
  geom_vline(xintercept = 45, linetype = "dashed",
             colour = col_high, linewidth = 0.5) +
  scale_fill_gradient(low = col_low, high = col_high, guide = "none") +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title    = "Percentage missingness per variable",
    subtitle = "Red dashed line = 45% exclusion threshold",
    x = "% Missing", y = NULL
  ) +
  theme_report(base_size = 11) +
  theme(axis.text.y = element_text(size = 8))

save_plot(p_var_bar, "fig_missingness_per_variable_bar.png",
          w = 8, h = max(4, nrow(missing_col) * 0.25))

# ============================================================
# 4. Participant-level missingness
# ============================================================

# Count missing variables per participant to identify data-sparse individuals
missing_row <- tibble(
  eid         = ukb$eid,
  n_missing   = rowSums(is.na(ukb_no_id)),
  pct_missing = 100 * n_missing / n_vars
)

write_csv(missing_row, file.path(out_dir, "missingness_per_participant.csv"))

# ── Figure 2: Histogram ───────────────────────────────────────────────────────
# Now uses theme_report() + matching fill colour + % x-axis labels
p_row <- ggplot(missing_row, aes(x = pct_missing)) +
  geom_histogram(bins = 50, fill = col_fill, colour = "white", linewidth = 0.2) +
  scale_x_continuous(
    breaks = seq(0, 100, 5),
    labels = function(x) paste0(x, "%")
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Missingness per participant",
    subtitle = "Distribution of % variables missing across all participants",
    x = "% of variables missing",
    y = "Number of participants"
  ) +
  theme_report(base_size = 12) +
  # For histogram, vertical grid lines (major.x) are less useful than horizontal ones
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#e0e0e0", linewidth = 0.4)
  )

save_plot(p_row, "fig_missingness_per_participant.png", w = 7, h = 5)

# ============================================================
# 5. Category-level missingness
# ============================================================
# Variables are grouped into thematic categories defined in variable_category.csv.

# --- Auto-create variable_category.csv from hardcoded mapping ---
if (!file.exists(category_map_path)) {
  variable_category_map <- tribble(
    ~variable,                  ~category,
    "recruitment_date",         "study_design",
    "age",                      "demographic",
    "death_date",               "outcomes",
    "sex",                      "demographic",
    "ethnicity",                "demographic",
    "household_income",         "socioeconomic",
    "highest_qualification",             "socioeconomic",
    "employment",               "socioeconomic",
    "housing",                  "socioeconomic",
    "household_size",           "socioeconomic",
    "has_vehicle",              "socioeconomic",
    "sbp",                      "physical_measures",
    "dbp",                      "physical_measures",
    "grip_strength",            "physical_measures",
    "whr",                      "physical_measures",
    "heel_bone_density",        "physical_measures",
    "body_fat_mass",            "physical_measures",
    "diet_score",               "lifestyle",
    "alcohol_consumption",      "lifestyle",
    "salt_added",               "lifestyle",
    "sleep_hrs",                "lifestyle",
    "sleep_cat",                "lifestyle",
    "insomnia",                 "lifestyle",
    "snoring",                  "lifestyle",
    "smoking_status_refined",   "lifestyle",
    "pack_yrs",                 "lifestyle",
    "passive_cat",              "lifestyle",
    "passive_total_hrs", "lifestyle",
    "met_min_week_all",             "lifestyle",
    "long_ill",                 "clinical_history",
    "diabetes_doctor_diagnosed",             "clinical_history",
    "serious_dx",               "clinical_history",
    "age_diabetes_dx",             "clinical_history",
    "age_high_bp_dx",                  "clinical_history",
    "early_insulin",            "clinical_history",
    "pregnancy",                "clinical_history",
    "angina_dx_age",            "clinical_history",
    "mi_dx_age",                "clinical_history",
    "dvt_dx_age",               "clinical_history",
    "pe_dx_age",                "clinical_history",
    "stroke_dx_age",            "clinical_history",
    "heart_problems_collapsed",           "clinical_history",
    "meds_chol_bp_diabetes_collapsed",              "clinical_history",
    "meds_chol_bp_diabetes_hormone_collapsed",             "clinical_history",
    "n_previous_pregnancies_collapsed",           "clinical_history",
    "mh_sought_help",           "clinical_history",
    "work_job_satisfaction_num",    "socioeconomic",
    "neuroticism_score",              "mental_health",
    "breastfed",                "early_life",
    "comparative_body_size_10",          "early_life",
    "maternal_smoking_birth",         "early_life",
    "birth_weight",             "early_life",
    "no2_2010",                 "environment",
    "pm2.5_2010",               "environment",
    "pm10_2010",                "environment",
    "urban_rural_3lvl",              "environment",
    "wbc",                      "biomarkers",
    "platelets",                "biomarkers",
    "albumin",                  "biomarkers",
    "alp",                      "biomarkers",
    "alt",                      "biomarkers",
    "apoa1",                    "biomarkers",
    "apob",                     "biomarkers",
    "ast",                      "biomarkers",
    "bilirubin",                "biomarkers",
    "total_cholesterol",        "biomarkers",
    "creatinine",               "biomarkers",
    "crp",               "biomarkers",
    "ggt",                      "biomarkers",
    "hba1c",                    "biomarkers",
    "hdl_cholesterol",          "biomarkers",
    "ldl_direct",               "biomarkers",
    "lpa",                      "biomarkers",
    "triglycerides",            "biomarkers",
    "first_cvd_date",           "outcomes",
    "has_cvd",                  "outcomes",
    "cystatin_c",               "biomarkers",
    "lymphocyte_count", "biomarkers",
    "monocyte_count", "biomarkers",
    "neutrophill_count", "biomarkers",
    
    
  )
  
  write_csv(variable_category_map, category_map_path)
  message("variable_category.csv created automatically from hardcoded mapping.")
}

variable_category <- readr::read_delim(
  category_map_path, delim = ",", trim_ws = TRUE, show_col_types = FALSE
)

# Handle edge case where file is read as a single column (encoding issues)
if (ncol(variable_category) == 1) {
  variable_category <- variable_category %>%
    tidyr::separate(
      col  = 1, into = c("variable", "category"),
      sep  = ",", remove = TRUE, extra = "merge", fill = "right"
    )
}

stopifnot(all(c("variable", "category") %in% names(variable_category)))

variable_category <- variable_category %>%
  mutate(
    variable = as.character(variable),
    category = stringr::str_to_lower(stringr::str_trim(as.character(category)))
  )

# Validate: every variable in the data must have a category assigned
cat_df <- missing_col %>%
  left_join(variable_category %>% select(variable, category), by = "variable")

if (any(is.na(cat_df$category))) {
  bad_vars <- cat_df %>% filter(is.na(category)) %>% pull(variable)
  stop("These variables have no category mapping: ", paste(bad_vars, collapse = ", "))
}

# ── Figure 3: Boxplot by category ─────────────────────────────────────────────
p_cat_box <- missing_col %>%
  left_join(variable_category, by = "variable") %>%
  ggplot(aes(x = pct_missing,
             y = fct_reorder(category, pct_missing, .fun = median))) +
  geom_boxplot(
    fill    = col_box,
    colour  = col_point,
    outlier.shape = NA,
    linewidth = 0.4
  ) +
  geom_jitter(
    height = 0.18, alpha = 0.6, size = 1.5,
    colour = col_fill
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title    = "Variable missingness by category",
    subtitle = "Boxes show distribution; points are individual variables",
    x = "% Missing", y = NULL
  ) +
  theme_report(base_size = 12) +
  # For this horizontal plot, major.x grids are the vertical reference lines — keep them
  theme(
    panel.grid.major.x = element_line(colour = "#e0e0e0", linewidth = 0.4)
  )

save_plot(p_cat_box, "fig_missingness_by_category_boxplot.png", w = 9, h = 7)

# ── Summary table ─────────────────────────────────────────────────────────────
missing_by_category <- missing_col %>%
  left_join(variable_category, by = "variable") %>%
  group_by(category) %>%
  summarise(
    n_vars         = n(),
    mean_missing   = mean(pct_missing, na.rm = TRUE),
    median_missing = median(pct_missing, na.rm = TRUE),
    n_vars_over_20 = sum(pct_missing > 20, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_missing))

write_csv(missing_by_category, file.path(out_dir, "missingness_by_category.csv"))

write_csv(missing_by_category, file.path(out_dir, "missingness_by_category.csv"))







