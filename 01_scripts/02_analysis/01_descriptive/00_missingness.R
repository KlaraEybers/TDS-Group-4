# ============================================================
# Missingness Analysis
# ============================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(here)
  library(dplyr)
  library(purrr)
})
set.seed(123)

# ── Load variable dictionary ──────────────────────────────────────────────────
source(here("01_scripts", "00_variable_dictionary.R"))
# Provides: display_names, get_display_name(), var_domains, get_domain_from_variable()

# ── Paths ─────────────────────────────────────────────────────────────────────
data_path <- here("00_data", "01_processed", "ukb_final_merged.rds")
out_dir   <- here("02_results", "02_missingness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ── Missingness thresholds ────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  variable_threshold <- as.numeric(args[1])
  participant_threshold <- as.numeric(args[2])
} else {
  variable_threshold <- 0.3
  participant_threshold <- 0.3
}


# ── Shared palette ────────────────────────────────────────────────────────────
col_low   <- "#a8d5e2"
col_high  <- "#c0392b"
col_fill  <- "#3a6ea5"
col_point <- "#2c3e50"
col_box   <- "#ecf0f1"

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_report <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold", size = base_size + 2,
                                        colour = "#2c3e50"),
      plot.subtitle      = element_text(size = base_size, colour = "#555555",
                                        margin = margin(b = 8)),
      axis.title         = element_text(face = "bold", colour = "#2c3e50"),
      axis.text          = element_text(colour = "#444444"),
      panel.grid.major.x = element_line(colour = "#e0e0e0", linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.title       = element_text(face = "bold"),
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

if (!("eid" %in% names(ukb))) {
  ukb <- ukb %>% tibble::rownames_to_column("eid")
}

nonresponse_levels <- c("Prefer not to answer", "Do not know", "Unsure")
ukb <- ukb %>%
  mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
  mutate(across(where(is.character),
                ~ if_else(.x %in% nonresponse_levels, NA_character_, .x)))

ukb_no_id <- ukb %>%
  select(-eid) %>%
  select(-any_of(c("death_date", "first_cvd_date"))) %>%
  select(any_of(names(var_domains)))           

n_total <- nrow(ukb_no_id)
n_vars  <- ncol(ukb_no_id)

# ============================================================
# 3. Variable-level missingness
# ============================================================
missing_col <- tibble(
  variable    = names(ukb_no_id),
  n_missing   = map_int(ukb_no_id, ~ sum(is.na(.x))),
  pct_missing = 100 * n_missing / n_total
) %>%
  arrange(desc(pct_missing)) %>%
  mutate(
    variable_label = get_display_name(variable),
    domain         = get_domain_from_variable(variable)
  ) %>%
  filter(domain != "Outcome")   

# Warn about any variables not in var_domains (they get "Other")
# unmapped <- missing_col %>%
#   filter(!(variable %in% names(var_domains))) %>%
#   pull(variable)
# 
# unmapped
# 
# if (length(unmapped) > 0) {
#   warning(
#     length(unmapped), " variable(s) not found in var_domains — assigned to 'Other':\n  ",
#     paste(unmapped, collapse = ", ")
#   )
# }

write.csv(missing_col, file.path(out_dir, "missingness_per_variable.csv"))

# ── Figure 1: Horizontal bar chart ───────────────────────────────────────────
p_var_bar <- ggplot(
  missing_col,
  aes(x = pct_missing, y = reorder(variable_label, pct_missing))
) +
  geom_col(aes(fill = pct_missing), width = 0.7) +
  geom_vline(xintercept = variable_threshold * 100, linetype = "dashed",
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
    subtitle = paste0("Red dashed line = ", variable_threshold * 100, "% exclusion threshold"),
    x = "% Missing", y = NULL
  ) +
  theme_report(base_size = 11) +
  theme(axis.text.y = element_text(size = 8))

save_plot(p_var_bar, "fig_missingness_per_variable_bar.png",
          w = 8, h = max(4, nrow(missing_col) * 0.25))

# ============================================================
# 4. Participant-level missingness
# ============================================================
missing_row <- tibble(
  eid         = ukb$eid,
  n_missing   = rowSums(is.na(ukb_no_id)),
  pct_missing = 100 * n_missing / n_vars
)

write.csv(missing_row, file.path(out_dir, "missingness_per_participant.csv"))

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
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#e0e0e0", linewidth = 0.4)
  )

save_plot(p_row, "fig_missingness_per_participant.png", w = 7, h = 5)

# ============================================================
# 5. Domain-level missingness
# ============================================================
# ── Figure 2: Boxplot by domain ───────────────────────────────────────────────
p_cat_box <- missing_col %>%
  ggplot(aes(x = pct_missing,
             y = reorder(domain, pct_missing, FUN = median))) +
  geom_boxplot(
    fill      = col_box,
    colour    = col_point,
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
    title    = "Variable missingness by domain",
    x = "% Missing", y = "Domain",
  ) +
  theme_report(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(colour = "#e0e0e0", linewidth = 0.4)
  )

save_plot(p_cat_box, "fig_missingness_by_domain_boxplot.png", w = 9, h = 7)

# ── Summary table by domain ───────────────────────────────────────────────────
missing_by_domain <- missing_col %>%
  group_by(domain) %>%
  summarise(
    n_vars         = n(),
    mean_missing   = mean(pct_missing, na.rm = TRUE),
    median_missing = median(pct_missing, na.rm = TRUE),
    n_vars_over_20 = sum(pct_missing > 20, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  arrange(desc(mean_missing))

write.csv(missing_by_domain, file.path(out_dir, "missingness_by_domain.csv"))

message("Missingness analysis complete. Outputs saved to: ", out_dir)

