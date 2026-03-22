rm(list = ls())

library(here)
library(dplyr)
library(pROC)
library(ggplot2)

# =========================================================
# 1) Paths
# =========================================================

out_dir <- here(
  "00_data", "02_analysis", "stability_selection_runs",
  "stability_selection_sharp_wo_sbp_FINAL"
)

# =========================================================
# 2) Load train/test data
# =========================================================

train_df <- readRDS(here("00_data", "01_processed", "03_imputation", "df_train_imputed.rds"))
test_df  <- readRDS(here("00_data", "01_processed", "03_imputation", "df_test_imputed.rds"))

# =========================================================
# 3) Load selected variables from stability selection
# =========================================================

selected_tbl <- read.csv(
  file.path(out_dir, "selected_variables_only.csv"),
  stringsAsFactors = FALSE
)

selected_vars <- selected_tbl %>%
  filter(group != "confounder") %>%
  arrange(desc(max_selection_proportion)) %>%
  pull(source_variable)

confounders <- c("age", "sex", "ethnicity")

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

n_boot <- 100
set.seed(123)

# =========================================================
# 6) Helper function: fit one bootstrap model and get test AUC
# =========================================================

fit_boot_auc <- function(vars_to_use, train_data, test_data) {
  
  idx <- sample.int(nrow(train_data), size = nrow(train_data), replace = TRUE)
  train_boot <- train_data[idx, , drop = FALSE]
  
  formula_str <- paste("has_cvd ~", paste(vars_to_use, collapse = " + "))
  model_formula <- as.formula(formula_str)
  
  fit <- glm(model_formula, data = train_boot, family = binomial())
  
  probs <- predict(fit, newdata = test_data, type = "response")
  
  if (length(unique(test_data$has_cvd)) < 2) {
    return(NA_real_)
  }
  
  as.numeric(auc(test_data$has_cvd, probs))
}

# =========================================================
# 7) Define staged confounder models
# =========================================================

confounder_steps <- list(
  c("sex"),
  c("sex", "ethnicity"),
  c("age", "sex", "ethnicity")
)

confounder_step_labels <- c(
  "Sex only",
  "Sex + Ethnicity",
  "Baseline"
)

# Keep only available variables in each step
confounder_steps <- lapply(confounder_steps, function(v) v[v %in% names(train_df) & v %in% names(test_df)])

# =========================================================
# 8) Incremental repeated-fit analysis
# =========================================================

results_list <- list()

# First: staged confounder models
for (s in seq_along(confounder_steps)) {
  
  vars_now <- confounder_steps[[s]]
  
  for (b in seq_len(n_boot)) {
    auc_val <- fit_boot_auc(
      vars_to_use = vars_now,
      train_data = train_df,
      test_data = test_df
    )
    
    results_list[[length(results_list) + 1]] <- data.frame(
      step = s - 1,
      stage_type = "confounder",
      added_variable = confounder_step_labels[s],
      auc = auc_val,
      stringsAsFactors = FALSE
    )
  }
}

# Then: add selected variables sequentially on top of full baseline
baseline_vars <- confounder_steps[[length(confounder_steps)]]
start_step <- length(confounder_steps)

for (i in seq_along(selected_vars)) {
  
  vars_now <- c(baseline_vars, selected_vars[1:i])
  
  for (b in seq_len(n_boot)) {
    auc_val <- fit_boot_auc(
      vars_to_use = vars_now,
      train_data = train_df,
      test_data = test_df
    )
    
    results_list[[length(results_list) + 1]] <- data.frame(
      step = start_step - 1 + i,
      stage_type = "selected",
      added_variable = selected_vars[i],
      auc = auc_val,
      stringsAsFactors = FALSE
    )
  }
}

all_auc_runs <- bind_rows(results_list)

# =========================================================
# 9) Summarise to median / 5th / 95th percentiles
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
    x_label = case_when(
      stage_type == "confounder" ~ added_variable,
      TRUE ~ paste0("+ ", added_variable)
    ),
    x_label = factor(x_label, levels = x_label)
  )

# Save raw and summary results
write.csv(
  all_auc_runs,
  file.path(out_dir, "incremental_auc_all_runs_staged_confounders.csv"),
  row.names = FALSE
)

write.csv(
  incremental_summary,
  file.path(out_dir, "incremental_auc_summary_staged_confounders.csv"),
  row.names = FALSE
)

# =========================================================
# 10) Plot incremental AUC with percentile bars
# =========================================================

incremental_plot <- ggplot(
  incremental_summary,
  aes(x = x_label, y = median_auc)
) +
  geom_errorbar(
    aes(ymin = p05_auc, ymax = p95_auc),
    width = 0.15,
    linewidth = 0.6,
    colour = "grey40"
  ) +
  geom_line(aes(group = 1), linewidth = 0.7, colour = "grey50") +
  geom_point(
    aes(colour = stage_type),
    size = 2.6
  ) +
  scale_colour_manual(values = c("confounder" = "black", "selected" = "red")) +
  labs(
    title = paste0(
      "Incremental Test AUC with Staged Confounders and Sequentially Added Stable Variables\n",
      "Points = median across ", n_boot, " bootstrap refits; bars = 5th-95th percentiles"
    ),
    x = NULL,
    y = "Test ROC AUC"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(incremental_plot)

ggsave(
  filename = file.path(out_dir, "incremental_auc_percentile_plot_staged_confounders.png"),
  plot = incremental_plot,
  width = 13,
  height = 7,
  dpi = 300
)

cat("Incremental percentile plot saved to:\n")
cat(file.path(out_dir, "incremental_auc_percentile_plot_staged_confounders.png"), "\n")