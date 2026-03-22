# ============================================================
# Appendix: Missingness Analysis
# ============================================================


# ============================================================
# Setup
# ============================================================

rm(list = ls())

#Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})
set.seed(1)

library(ggplot2)

#Paths
data_path <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group4/00_data/01_processed/ukb_final_merged.rds"
out_dir   <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group4/02_results/02_reports/missingness"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Theme
theme_report <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size, margin = margin(b = 8)),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 12, 10, 12)
    )
}

# Helper to save plots
save_plot <- function(p, filename, w = 7, h = 5, dpi = 300) {
  ggsave(file.path(out_dir, filename), p, width = w, height = h, dpi = dpi)
}


# Load data ---------------------------------------------------------------

ukb <- readRDS(data_path)

if (!("eid" %in% names(ukb))) {
  ukb <- ukb %>% tibble::rownames_to_column("eid")
}

# Data cleaning -----------------------------------------------------------

# Everything non response becomes NA
nonresponse_levels <- c(
  "Prefer not to answer",
  "Do not know",
  "Unsure"
)
ukb <- ukb %>%
  mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
  mutate(across(where(is.character),
                ~ if_else(.x %in% nonresponse_levels, NA_character_, .x)))

ukb_no_id <- ukb %>% select(-eid)
n_total   <- nrow(ukb_no_id)
ukb_no_id <- ukb_no_id %>% select(-any_of(c("death_date", "first_cvd_date")))
n_vars    <- ncol(ukb_no_id)


# ============================================================
# Missingness Analysis
# ============================================================

# Missingness per variable --------------------------------
missing_col <- tibble(
  variable = names(ukb_no_id),
  n_missing = map_int(ukb_no_id, ~ sum(is.na(.x))),
  pct_missing = 100 * n_missing / n_total
) %>%
  arrange(desc(pct_missing))

# ---- Bar chart: % missingness per variable, ordered ----
p_var_bar <- ggplot(
  missing_col,
  aes(
    x = pct_missing,
    y = reorder(variable, pct_missing)
  )
) +
  geom_col(
    aes(fill = pct_missing),
    width = 0.7
  ) +
  geom_vline(xintercept = 45, linetype = "dashed", colour = "red", linewidth = 0.5) +
  scale_fill_gradient(low = "#a8d5e2", high = "#c0392b", guide = "none") +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "Percentage missingness per variable",
    subtitle = "Red dashed line = 45% threshold",
    x = "% Missing",
    y = NULL
  ) +
  theme_report(base_size = 11) +
  theme(
    axis.text.y        = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(10, 15, 10, 10)
  )

save_plot(
  p_var_bar,
  "fig_missingness_per_variable_bar.png",
  w = 8,
  h = max(4, nrow(missing_col) * 0.25)  # auto-scales height with number of variables
)



# -------------------------------------------------
# Read variable_category.csv robustly
# -------------------------------------------------
category_map_path <- file.path(out_dir, "variable_category.csv")

variable_category <- readr::read_delim(
  category_map_path,
  delim = ",",
  trim_ws = TRUE,
  show_col_types = FALSE
)

# If it got read as one column like "recruitment_date,study_design", split it
if (ncol(variable_category) == 1) {
  variable_category <- variable_category %>%
    tidyr::separate(
      col = 1,
      into = c("variable", "category"),
      sep = ",",
      remove = TRUE,
      extra = "merge",
      fill = "right"
    )
}

# Basic checks
stopifnot(all(c("variable", "category") %in% names(variable_category)))

variable_category <- variable_category %>%
  mutate(
    variable = as.character(variable),
    category = stringr::str_to_lower(stringr::str_trim(as.character(category)))
  )

# -------------------------------------------------
# Box plot: variable-level missingness by category
# -------------------------------------------------
cat_df <- missing_col %>%
  left_join(variable_category %>% select(variable, category), by = "variable")

# Catch any variables that failed to match a category
if (any(is.na(cat_df$category))) {
  bad_vars <- cat_df %>% filter(is.na(category)) %>% pull(variable)
  stop("These variables have no category mapping: ", paste(bad_vars, collapse = ", "))
}


write_csv(missing_col,
          file.path(out_dir, "missingness_per_variable.csv"))

p_var <- ggplot(missing_col, aes(x = pct_missing)) +
  geom_histogram(bins = 40) +
  theme_classic(base_size = 12) +
  labs(
    title = "Distribution of missingness per variable",
    x = "Percent missing variables",
    y = "Number of variables"
  )

save_plot(p_var, "fig_missingness_per_variable.png")

# Missingness per participant ---------------------------------------------
missing_row <- tibble(
  eid = ukb$eid,
  n_missing = rowSums(is.na(ukb_no_id)),
  pct_missing = 100 * n_missing / n_vars
)

write_csv(missing_row,
          file.path(out_dir, "missingness_per_participant.csv"))

p_row <- ggplot(missing_row, aes(x = pct_missing)) +
  geom_histogram(bins = 50) +
  theme_classic(base_size = 12) +
  labs(
    title = "Missingness per participant",
    x = "Percent missing variables",
    y = "Number of participants"
  )

# Identify problematic variables 
high_missing <- missing_col %>%
  filter(pct_missing > 20)

write_csv(high_missing,
          file.path(out_dir, "variables_above_20pct_missing.csv"))

# Missingness per category --------------------------------------------------

p_cat_box <- missing_col %>%
  left_join(variable_category, by = "variable") %>%
  ggplot(
    aes(
      x = pct_missing,
      y = fct_reorder(category, pct_missing, .fun = median)
    )
  ) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0.15, alpha = 0.6, size = 1) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "Variable missingness by category",
    subtitle = "Boxplots show distribution; points are individual variables",
    x = "Percent missing variables",
    y = NULL
  ) +
  theme_report(base_size = 12)

save_plot(
  p_cat_box,
  "fig_missingness_by_category_boxplot.png",
  w = 9,
  h = 7
)

# ============================================================
# Missingness by category
# ============================================================

missing_by_category <- missing_col %>%
  left_join(variable_category, by = "variable") %>%
  group_by(category) %>%
  summarise(
    n_vars = n(),
    mean_missing = mean(pct_missing, na.rm = TRUE),
    median_missing = median(pct_missing, na.rm = TRUE),
    n_vars_over_20 = sum(pct_missing > 20, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_missing))

write_csv(missing_by_category,
          file.path(out_dir, "missingness_by_category.csv"))

p_cat <- ggplot(missing_by_category,
                aes(x = reorder(category, mean_missing),
                    y = mean_missing)) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(
    title = "Mean missingness by category",
    x = NULL,
    y = "Mean percent missing"
  )

save_plot(p_cat, "fig_missingness_by_category.png")









# missingness -------------------------------------------------------------

# -------------------------------------------------
# Structural missingness screen: create a results table
# -------------------------------------------------

# 1) Candidates: only check variables above a threshold
threshold <- 20

structural_candidates <- missing_col %>%
  dplyr::filter(pct_missing > threshold) %>%
  dplyr::pull(variable)

message("Variables above ", threshold, "% missing: ", length(structural_candidates))

# If none, stop early
if (length(structural_candidates) == 0) {
  stop("No variables exceed the missingness threshold, so there is nothing to screen.")
}

# 2) Choose gatekeepers that exist in your dataset
gatekeepers_requested <- c("sex", "age", "recruitment_date", "ethnicity", "employment")

gatekeepers <- intersect(gatekeepers_requested, names(ukb))
missing_gates <- setdiff(gatekeepers_requested, names(ukb))

if (length(missing_gates) > 0) {
  message("These requested gatekeepers are not in ukb and will be skipped: ",
          paste(missing_gates, collapse = ", "))
}
if (length(gatekeepers) == 0) {
  stop("None of the requested gatekeepers exist in ukb. Update gatekeepers_requested.")
}

message("Using gatekeepers: ", paste(gatekeepers, collapse = ", "))

# 3) Function to test if missingness is almost perfectly explained by gatekeeper
test_structural <- function(var, gate, data, cutoff = 0.95) {
  # Return NA if var or gate missing for any reason
  if (!(var %in% names(data)) || !(gate %in% names(data))) return(NA)
  
  df <- data %>%
    dplyr::transmute(
      gate = .data[[gate]],
      is_missing = is.na(.data[[var]])
    ) %>%
    dplyr::filter(!is.na(gate))
  
  # Need at least 2 levels to be meaningful
  if (dplyr::n_distinct(df$gate) < 2) return(NA)
  
  max_props <- df %>%
    dplyr::count(gate, is_missing) %>%
    dplyr::group_by(gate) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::summarise(max_prop = max(prop), .groups = "drop")
  
  any(max_props$max_prop >= cutoff)
}

# 4) Build the results table
structural_flags <- tidyr::expand_grid(
  variable = structural_candidates,
  gatekeeper = gatekeepers
) %>%
  dplyr::mutate(
    structural = purrr::map2_lgl(
      variable, gatekeeper,
      ~ test_structural(.x, .y, data = ukb, cutoff = 0.95)
    )
  ) %>%
  dplyr::arrange(dplyr::desc(structural), gatekeeper, variable)

# 5) Save and show a preview
readr::write_csv(
  structural_flags,
  file.path(out_dir, paste0("structural_missingness_flags_over_", threshold, "pct.csv"))
)

print(dplyr::filter(structural_flags, structural) %>% dplyr::distinct(variable, gatekeeper))

