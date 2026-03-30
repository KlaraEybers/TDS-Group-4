rm(list = ls())

# =========================================================
# Incremental AUC analysis
#
# Evaluates how test ROC AUC changes as variables are added
# to the model. Starts with staged confounders, then adds
# LASSO-selected variables sequentially, using repeated
# subsampling to estimate median and percentile AUC.
# =========================================================

library(here)
library(dplyr)
library(pROC)
library(ggplot2)

source(here("01_scripts", "00_variable_config.R"))
source(here("01_scripts", "00_variable_dictionary.R"))

# =========================================================
# 1) Paths
# =========================================================

out_dir      <- lasso_outcome_dir
plot_out_dir <- incremental_fig_dir

# =========================================================
# 2) Load train/test splitted data
# =========================================================

train_df <- readRDS(file.path(data_input_dir, train_file))
test_df  <- readRDS(file.path(data_input_dir, test_file))

# =========================================================
# 3) Load selected variables from stability selection LASSO 1
# =========================================================

selected_tbl <- read.csv(
  file.path(out_dir, "selected_variables_only.csv"),
  stringsAsFactors = FALSE
)

selected_vars <- selected_tbl %>%
  filter(group != "confounder") %>%
  arrange(desc(max_selection_proportion)) %>%
  pull(source_variable)


# Keep only variables present in both datasets
selected_vars <- selected_vars[selected_vars %in% names(train_df) & selected_vars %in% names(test_df)]
confounders   <- confounders[confounders %in% names(train_df) & confounders %in% names(test_df)]

# =========================================================
# 4) Prepare factor variables consistently
# =========================================================

train_df <- train_df %>%
  mutate(across(where(is.character), as.factor))

test_df <- test_df %>%
  mutate(across(where(is.character), as.factor))

for (v in names(train_df)) {
  if (is.factor(train_df[[v]]) && v %in% names(test_df)) {
    test_df[[v]] <- factor(test_df[[v]], levels = levels(train_df[[v]]))
  }
}

# =========================================================
# 5) Settings
# =========================================================

n_resamples <- 100
subsample_frac <- 0.8
set.seed(123)

# =========================================================
# 6) Pre-create fixed train subsamples (without replacement)
# =========================================================

subsample_size <- floor(subsample_frac * nrow(train_df))

train_subsample_indices <- replicate(
  n_resamples,
  sample.int(nrow(train_df), size = subsample_size, replace = FALSE),
  simplify = FALSE
)

# =========================================================
# 7) Helper function: fit one subsample model and get test AUC
# =========================================================

fit_subsample_auc <- function(vars_to_use, train_data, test_data, idx) {
  
  train_sub <- train_data[idx, , drop = FALSE]
  
  formula_str <- paste("has_cvd ~", paste(vars_to_use, collapse = " + "))
  model_formula <- as.formula(formula_str)
  
  fit <- glm(model_formula, data = train_sub, family = binomial())
  
  probs <- predict(fit, newdata = test_data, type = "response")
  
  if (length(unique(test_data$has_cvd)) < 2) {
    return(NA_real_)
  }
  
  as.numeric(auc(test_data$has_cvd, probs))
}

# =========================================================
# 8) Define staged confounder models
# =========================================================

confounder_steps <- list(
  c("sex"),
  c("sex", "ethnicity"),
  c("age", "sex", "ethnicity")
)

confounder_step_labels <- c(
  "Sex",
  "+ Ethnicity",
  "+ Age (Full confounders)"
)

# Keep only available variables in each step
confounder_steps <- lapply(confounder_steps, function(v) v[v %in% names(train_df) & v %in% names(test_df)])

# =========================================================
# 9) Incremental repeated-fit analysis using fixed subsamples
# =========================================================

results_list <- list()

# First: staged confounder models
for (s in seq_along(confounder_steps)) {
  
  vars_now <- confounder_steps[[s]]
  
  for (b in seq_len(n_resamples)) {
    auc_val <- fit_subsample_auc(
      vars_to_use = vars_now,
      train_data = train_df,
      test_data = test_df,
      idx = train_subsample_indices[[b]]
    )
    
    results_list[[length(results_list) + 1]] <- data.frame(
      step = s - 1,
      stage_type = "confounder",
      added_variable = confounder_step_labels[s],
      auc = auc_val,
      resample_id = b,
      stringsAsFactors = FALSE
    )
  }
}

# Then: add selected variables sequentially on top of full baseline
baseline_vars <- confounder_steps[[length(confounder_steps)]]
start_step <- length(confounder_steps)

for (i in seq_along(selected_vars)) {
  
  vars_now <- c(baseline_vars, selected_vars[1:i])
  
  for (b in seq_len(n_resamples)) {
    auc_val <- fit_subsample_auc(
      vars_to_use = vars_now,
      train_data = train_df,
      test_data = test_df,
      idx = train_subsample_indices[[b]]
    )
    
    results_list[[length(results_list) + 1]] <- data.frame(
      step = start_step - 1 + i,
      stage_type = "selected",
      added_variable = selected_vars[i],
      auc = auc_val,
      resample_id = b,
      stringsAsFactors = FALSE
    )
  }
}

all_auc_runs <- bind_rows(results_list)

# =========================================================
# 10) Summarise to median / 5th / 95th percentiles
# =========================================================

incremental_summary <- all_auc_runs %>%
  group_by(step, stage_type, added_variable) %>%
  summarise(
    median_auc = median(auc, na.rm = TRUE),
    p05_auc = quantile(auc, 0.05, na.rm = TRUE),
    p95_auc = quantile(auc, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(step) %>%
  mutate(
    added_variable_label = case_when(
      stage_type == "confounder" ~ added_variable,
      TRUE ~ get_display_name(added_variable)
    ),
    x_label = case_when(
      stage_type == "confounder" ~ added_variable_label,
      TRUE ~ paste0("+ ", added_variable_label)
    ),
    x_label = factor(x_label, levels = x_label)
  )

# Save raw and summary results
write.csv(
  all_auc_runs,
  file.path(plot_out_dir, "incremental_auc_all_runs_staged_confounders.csv"),
  row.names = FALSE
)

write.csv(
  incremental_summary,
  file.path(plot_out_dir, "incremental_auc_summary_staged_confounders.csv"),
  row.names = FALSE
)

# =========================================================
# 11) Plot incremental AUC with percentile bars
# =========================================================

incremental_plot <- ggplot(
  incremental_summary,
  aes(x = x_label, y = median_auc)
) +
  geom_errorbar(
    aes(ymin = p05_auc, ymax = p95_auc),
    width = 0.15,
    linewidth = 0.6,
    colour = "#7f8c8d"
  ) +
  geom_line(
    aes(group = 1),
    linewidth = 0.7,
    colour = "#7f8c8d"
  ) +
  geom_point(
    aes(colour = stage_type),
    size = 2.6
  ) +
  scale_colour_manual(
    values = c(
      "confounder" = "#2c3e50",
      "selected"   = "#2f80c1"
    )
  ) +
  labs(
    title = "Incremental Test AUC for Stable Variables",
    x = NULL,
    y = "Test ROC AUC"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    
    panel.grid.major.x = element_line(colour = "#e0e0e0", linewidth = 0.4),
    panel.grid.major.y = element_line(colour = "#e0e0e0", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    
    plot.title = element_text(hjust = 0.5, face = "bold", colour = "#2c3e50")
  )

print(incremental_plot)

ggsave(
  filename = file.path(plot_out_dir, "incremental_auc_percentile_plot_staged_confounders.png"),
  plot = incremental_plot,
  width = 13,
  height = 7,
  dpi = 300
)

cat("Incremental percentile plot saved to:\n")
cat(file.path(plot_out_dir, "incremental_auc_percentile_plot_staged_confounders.png"), "\n")

