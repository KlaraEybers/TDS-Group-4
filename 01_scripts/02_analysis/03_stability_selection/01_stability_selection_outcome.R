rm(list = ls())

# =========================================================
# First-step stability selection for outcome model
# Model:
#     Y (CVD) ~ exposures + biomarkers + confounders
#
# Purpose:
#     Identify variables with stable direct associations
#     with the outcome (CVD). Selected exposures and
#     biomarkers will be used in second-step mediator models.
#
# Notes:
#     - Confounders are unpenalised
#     - Exposures and biomarkers are penalised
# =========================================================

# Packages -------------------------------------------------

library(here)
library(dplyr)
library(purrr)
library(ggplot2)
library(sharp)

source(here("01_scripts", "00_variable_config.R"))
source(here("01_scripts", "00_variable_dictionary.R"))

# Load data
df <- readRDS(file.path(data_input_dir, selection_file))

# Basic checks before starting 
stopifnot("has_cvd" %in% names(df))
stopifnot(all(df$has_cvd %in% c(0, 1)))

na_total <- sum(is.na(df))
cat("Total NA values in dataset:", na_total, "\n")

# Define output paths -----------------------------------------------------

out_dir       <- lasso_outcome_dir
fig_out_dir   <- lasso_fig_dir
tbl_out_dir   <- lasso_tbl_dir

# Helper functions --------------------------------------------------------

write_csv_safe <- function(x, path, row.names = FALSE) {
  if (is.null(x)) return(invisible(FALSE))
  if (!is.data.frame(x)) x <- as.data.frame(x)
  write.csv(x, path, row.names = row.names)
  invisible(TRUE)
}

write_csv_both <- function(x, file_name, row.names = FALSE) {
  write_csv_safe(x, file.path(out_dir, file_name), row.names = row.names)
  write_csv_safe(x, file.path(tbl_out_dir, file_name), row.names = row.names)
}

write_csv_table_only <- function(x, file_name, row.names = FALSE) {
  write_csv_safe(x, file.path(tbl_out_dir, file_name), row.names = row.names)
}

write_lines_both <- function(x, file_name) {
  writeLines(x, file.path(out_dir, file_name))
  writeLines(x, file.path(tbl_out_dir, file_name))
}

write_component_csv <- function(obj, comp_name, file_name, table_only = TRUE) {
  if (comp_name %in% names(obj)) {
    comp <- obj[[comp_name]]
    if (is.data.frame(comp) || is.matrix(comp) || is.vector(comp)) {
      comp_df <- as.data.frame(comp)
      if (table_only) {
        write_csv_table_only(comp_df, file_name)
      } else {
        write_csv_both(comp_df, file_name)
      }
    }
  }
}

# Variable wrangling ------------------------------------------------------

# Keep only variables that actually exist in data
confounder_vars <- intersect(confounder_vars, names(df))
exposure_vars   <- intersect(exposure_vars, names(df))
bio_vars        <- intersect(bio_vars, names(df))

# Removes any duplicates
candidate_predictors <- unique(c(confounder_vars, exposure_vars, bio_vars))

# Keep only required columns 
df <- df %>%
  select(all_of(c(intersect("eid", names(.)), "has_cvd", candidate_predictors)))

# Save IDs from data set
if ("eid" %in% names(df)) {
  write_csv_safe(
    df %>% select(eid, has_cvd),
    file.path(out_dir, "selection_dataset_ids.csv"),
    row.names = FALSE
  )
}

# Remove eid before modelling
if ("eid" %in% names(df)) {
  df <- df %>% select(-eid)
}

# Remove excluded variables (e.g. anthropometrics in main analysis)
if (length(exclude_vars) > 0) {
  exclude_in_data <- intersect(exclude_vars, names(df))
  if (length(exclude_in_data) > 0) {
    cat("Removing excluded variables:", paste(exclude_in_data, collapse = ", "), "\n")
    df <- df %>% select(-all_of(exclude_in_data))
  }
}

# Remove sbp in all scenarios — dominates selection and masks other associations
if ("sbp" %in% names(df)) {
  cat("Removing sbp (modelling decision: dominates variable selection)\n")
  df <- df %>% select(-sbp)
}

# Re-sync variable groups after exclusions/removals so downstream logic stays correct
confounder_vars <- intersect(confounder_vars, names(df))
exposure_vars   <- intersect(exposure_vars, names(df))
bio_vars        <- intersect(bio_vars, names(df))
candidate_predictors <- unique(c(confounder_vars, exposure_vars, bio_vars))

# =========================================================
# 2) Prepare data types
# =========================================================

# Convert character -> factor
df <- df %>%
  mutate(across(where(is.character), as.factor))

# Convert logical -> factor
df <- df %>%
  mutate(across(where(is.logical), as.factor))

# Convert date -> numeric
date_cols <- names(df)[vapply(df, inherits, logical(1), what = "Date")]
if (length(date_cols) > 0) {
  df <- df %>%
    mutate(across(all_of(date_cols), ~ as.numeric(.)))
}

# Drop unused factor levels to reduce encoding surprises
df <- df %>%
  mutate(across(where(is.factor), droplevels))

# Ensure outcome is numeric 0/1
df$has_cvd <- as.numeric(df$has_cvd)

# Final sanity checks before encoding
stopifnot("has_cvd" %in% names(df))
stopifnot(all(df$has_cvd %in% c(0, 1)))
stopifnot(length(unique(df$has_cvd)) >= 2)

# =========================================================
# 3) One-hot encoding for sharp
# =========================================================

# Exposure and outcome dataframes
x_df <- df %>% select(-has_cvd)
y <- df$has_cvd

if (ncol(x_df) == 0) {
  stop("No predictor columns remain after variable filtering and exclusions.")
}

# model.matrix one-hot encodes factors automatically
x_mat <- model.matrix(~ ., data = x_df)

# Remove intercept if present
if ("(Intercept)" %in% colnames(x_mat)) {
  x_mat <- x_mat[, colnames(x_mat) != "(Intercept)", drop = FALSE]
}

if (ncol(x_mat) == 0) {
  stop("No encoded features were created by model.matrix().")
}

# Remove zero-variance columns (all same value)
keep_cols <- apply(x_mat, 2, function(z) {
  z_num <- as.numeric(z)
  stats::sd(z_num, na.rm = TRUE) > 0
})

if (length(keep_cols) == 0 || !any(keep_cols)) {
  stop("All encoded features were removed as zero-variance.")
}

x_mat <- x_mat[, keep_cols, drop = FALSE]

# Defensive check
if (ncol(x_mat) == 0) {
  stop("No encoded features remain after zero-variance filtering.")
}

# =========================================================
# 4) Penalty factors: Confounders unpenalised
# =========================================================

penalty_factor <- rep(1, ncol(x_mat))
names(penalty_factor) <- colnames(x_mat)

# Any dummy column derived from confounders gets penalty set to 0
for (cf in confounder_vars) {
  idx <- grep(paste0("^", cf), colnames(x_mat))
  if (length(idx) > 0) penalty_factor[idx] <- 0
}

# =========================================================
# 5) Run sharp variable selection
# =========================================================

# n_cat = 3:
#   1) stably included
#   2) stably excluded
#   3) unstable

set.seed(123)

# 100 iterations should suffice
out <- sharp::VariableSelection(
  xdata = x_mat,
  ydata = y,
  family = "binomial",
  penalty.factor = penalty_factor,
  verbose = FALSE,
  n_cat = 3,
  K = 300
  # Below is for confining the grid to particular region
  # ,pi_list = seq(0.65, 0.79, by = 0.01),
  # Lambda = seq(4e-04, 8e-04, length.out = 60)
)

# Save full sharp object
saveRDS(out, file.path(out_dir, "sharp_variable_selection_object.rds"))

# Save object names and structure for debugging
capture.output(
  {
    cat("Object names:\n")
    print(names(out))
    cat("\nStructure:\n")
    str(out, max.level = 2)
  },
  file = file.path(out_dir, "sharp_object_structure.txt")
)

# =========================================================
# 6) Save object outputs for debugging
# =========================================================

# Save full sharp object
saveRDS(out, file.path(out_dir, "sharp_variable_selection_object.rds"))

# Save object names / structure for debugging / inspection
capture.output(
  {
    cat("Object names:\n")
    print(names(out))
    cat("\nStructure:\n")
    str(out, max.level = 2)
  },
  file = file.path(out_dir, "sharp_object_structure.txt")
)

write_component_csv(out, "summary_variable", "summary_variable.csv", table_only = TRUE)
write_component_csv(out, "summary_model", "summary_model.csv", table_only = TRUE)
write_component_csv(out, "summary_selection", "summary_selection.csv", table_only = TRUE)
write_component_csv(out, "summary_cat", "summary_cat.csv", table_only = TRUE)
write_component_csv(out, "Table", "sharp_table.csv", table_only = TRUE)
write_component_csv(out, "PFER", "sharp_pfer.csv", table_only = TRUE)
write_component_csv(out, "selprop", "sharp_selprop_matrix.csv", table_only = TRUE)
write_component_csv(out, "Lambda", "sharp_lambda.csv", table_only = TRUE)
write_component_csv(out, "Q", "sharp_q.csv", table_only = TRUE)
write_component_csv(out, "P", "sharp_p.csv", table_only = TRUE)
write_component_csv(out, "FDP", "sharp_fdp.csv", table_only = TRUE)

# Calibration plot for viewing lambda and pi patterns
png(
  filename = file.path(fig_out_dir, "sharp_calibration_plot.png"),
  width = 2300,
  height = 1300,
  res = 150
)

par(mar = c(6.5, 6.5, 5.0, 7.0))

try(sharp::CalibrationPlot(out), silent = TRUE)

dev.off()

# =========================================================
# 7) Practical calibrated outputs for debugging
# =========================================================

selprop <- SelectionProportions(out)
hat_params <- Argmax(out)

if (is.null(names(selprop))) {
  names(selprop) <- colnames(x_mat)
}

# Defensive fallback in case Argmax output is unnamed or shorter than expected
sel_threshold <- NA_real_
if (length(hat_params) >= 2) {
  sel_threshold <- as.numeric(hat_params[2])
}

if (is.na(sel_threshold) || !is.finite(sel_threshold)) {
  stop("Could not recover a valid calibrated selection threshold from Argmax(out).")
}

cat("Calibrated selection threshold from sharp:", sel_threshold, "\n")

selection_proportions_raw <- data.frame(
  feature = names(selprop),
  selection_proportion = as.numeric(selprop),
  stringsAsFactors = FALSE
)

write_csv_table_only(
  selection_proportions_raw,
  "selection_proportions_raw.csv",
  row.names = FALSE
)

sharp_argmax_params <- data.frame(
  parameter = c("lambda", "pi"),
  value = as.numeric(hat_params),
  stringsAsFactors = FALSE
)

write_csv_table_only(
  sharp_argmax_params,
  "sharp_argmax_params.csv",
  row.names = FALSE
)

# =========================================================
# 8) Feature and source variable summaries
# =========================================================

selection_table <- data.frame(
  feature = names(selprop),
  selection_proportion = as.numeric(selprop),
  selected = as.numeric(selprop) >= sel_threshold,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(selection_proportion), feature)

get_source_variable <- function(feature_name, original_vars) {
  hits <- original_vars[startsWith(feature_name, original_vars)]
  if (length(hits) == 0) return(NA_character_)
  hits[which.max(nchar(hits))]
}

selection_table$source_variable <- sapply(
  selection_table$feature,
  get_source_variable,
  original_vars = c(confounder_vars, exposure_vars, bio_vars)
)

selection_table$group <- dplyr::case_when(
  selection_table$source_variable %in% confounder_vars ~ "confounder",
  selection_table$source_variable %in% bio_vars        ~ "biological",
  selection_table$source_variable %in% exposure_vars   ~ "exposure",
  TRUE ~ "unknown"
)

write_csv_table_only(
  selection_table,
  "selection_table_feature_level.csv",
  row.names = FALSE
)

selected_variables_final <- selection_table %>%
  filter(!is.na(source_variable)) %>%
  group_by(source_variable, group) %>%
  summarise(
    max_selection_proportion = max(selection_proportion, na.rm = TRUE),
    selected = any(selected),
    .groups = "drop"
  ) %>%
  arrange(desc(selected), desc(max_selection_proportion), group, source_variable)

# Outputs selection proportions for all variables
write_csv_both(
  selected_variables_final,
  "selected_variables_final.csv",
  row.names = FALSE
)

selected_variables_only <- selected_variables_final %>%
  filter(selected)

# Outputs selection proportions for selected variables
write_csv_both(
  selected_variables_only,
  "selected_variables_only.csv",
  row.names = FALSE
)

# Adjust group names so more apropriate for mediation analysis
selected_variable_names <- selected_variables_final %>%
  filter(selected) %>%
  transmute(
    variable = source_variable,
    type = case_when(
      group == "biological" ~ "mediator",
      group == "exposure" ~ "exposure",
      group == "confounder" ~ "confounder",
      TRUE ~ group
    )
  ) %>%
  distinct()

# Clean csv with only selected variable names
write_csv_both(
  selected_variable_names,
  "selected_variable_names.csv",
  row.names = FALSE
)

# Same as above just .txt format
write_lines_both(
  paste(selected_variable_names$variable, selected_variable_names$type, sep = ","),
  "selected_variable_names.txt"
)

# Optional useful mapping table for rebuilding plots later
feature_to_source_map <- selection_table %>%
  select(feature, source_variable, group, selection_proportion, selected)

write_csv_table_only(
  feature_to_source_map,
  "feature_to_source_map.csv",
  row.names = FALSE
)

# =========================================================
# 9) Plots
# =========================================================

# Map coding names to display names
selprop_labels <- get_display_name(names(selprop))

# This plot shows selection proportion of all variables
png(
  filename = file.path(fig_out_dir, "selection_proportion_plot_calibrated.png"),
  width = 1800,
  height = 1100,
  res = 150
)

par(mar = c(14, 5, 4, 1))

plot(
  selprop,
  type = "h",
  lwd = 3,
  las = 1,
  xlab = "",
  ylab = "Selection Proportion",
  main = "Outcome Model: Selection Proportions",
  xaxt = "n",
  col = ifelse(selprop >= sel_threshold, "red", "#7f8c8d"),
  cex.lab = 1.5,
  cex.main = 1.4
)

abline(h = sel_threshold, lty = 2, col = "darkred")

for (i in seq_along(selprop)) {
  axis(
    side = 1,
    at = i,
    labels = selprop_labels[i],
    las = 2,
    cex.axis = 0.7,
    col = ifelse(selprop[i] >= sel_threshold, "red", "#7f8c8d"),
    col.axis = ifelse(selprop[i] >= sel_threshold, "red", "#7f8c8d")
  )
}

dev.off()

# Get top 20 variables by selection proportion
top_n <- min(20, length(selprop))
top_idx <- order(selprop, decreasing = TRUE)[seq_len(top_n)]

selprop_top <- selprop[top_idx]
selprop_top_labels <- selprop_labels[top_idx]

# Plot the top 20 selection proportions for more concise plot than above
png(
  filename = file.path(fig_out_dir, "selection_proportion_plot_top20.png"),
  width = 1800,
  height = 1100,
  res = 150
)

par(mar = c(14, 5, 4, 1))

plot(
  selprop_top,
  type = "h",
  lwd = 4,
  las = 1,
  xlab = "",
  ylab = "Selection Proportion",
  main = "Outcome Model: Top 20 Selection Proportions",
  xaxt = "n",
  col = ifelse(selprop_top >= sel_threshold, "red", "#7f8c8d"),
  cex.lab = 1.5,
  cex.main = 1.4
)

abline(h = sel_threshold, lty = 2, col = "darkred")

for (i in seq_along(selprop_top)) {
  axis(
    side = 1,
    at = i,
    labels = selprop_top_labels[i],
    las = 2,
    cex.axis = 0.9,
    col = ifelse(selprop_top[i] >= sel_threshold, "red", "#7f8c8d"),
    col.axis = ifelse(selprop_top[i] >= sel_threshold, "red", "#7f8c8d")
  )
}

dev.off()

# This plot shows selection proportion of all variables ordered
# Order by decreasing selection proportion
ord <- order(selprop, decreasing = TRUE)

selprop_ord <- selprop[ord]
selprop_labels_ord <- selprop_labels[ord]

# Ordered selection proportion plot with display names
png(
  filename = file.path(fig_out_dir, "selection_proportion_plot_ordered.png"),
  width = 1800,
  height = 1100,
  res = 150
)

par(mar = c(14, 5, 4, 1))

plot(
  selprop_ord,
  type = "h",
  lwd = 3,
  las = 1,
  xlab = "",
  ylab = "Selection Proportion",
  main = "Outcome Model: Ordered Selection Proportions",
  xaxt = "n",
  col = ifelse(selprop_ord >= sel_threshold, "red", "#7f8c8d"),
  cex.lab = 1.5,
  cex.main = 1.4
)

abline(h = sel_threshold, lty = 2, col = "darkred")

for (i in seq_along(selprop_ord)) {
  axis(
    side = 1,
    at = i,
    labels = selprop_labels_ord[i],
    las = 2,
    cex.axis = 0.7,
    col = ifelse(selprop_ord[i] >= sel_threshold, "red", "#7f8c8d"),
    col.axis = ifelse(selprop_ord[i] >= sel_threshold, "red", "#7f8c8d")
  )
}

# =========================================================
# 10) Console summary for us to view after running
# =========================================================

cat("\n====================================================\n")
cat("sharp variable selection completed\n")
cat("====================================================\n")
cat("Selection dataset rows:", nrow(df), "\n")
cat("Number of encoded features:", ncol(x_mat), "\n")
cat("n_cat:", 3, "\n")
cat("Calibrated selection threshold:", sel_threshold, "\n")
cat("Output directory:\n", out_dir, "\n")
cat("Figure directory:\n", fig_out_dir, "\n")
cat("Table directory:\n", tbl_out_dir, "\n")
cat("====================================================\n\n")

cat("Confounders used:\n")
print(confounder_vars)

cat("\nExposure variables used:\n")
print(exposure_vars)

cat("\nBiological variables used:\n")
print(bio_vars)

cat("\nTop encoded features by selection proportion:\n")
print(head(selection_table, 20))

cat("\nRecovered selected source variables:\n")
print(selected_variables_only)

cat("Number of original predictors:", length(candidate_predictors), "\n")
cat("Number of predictors in df:", ncol(df) - 1, "\n")  # minus outcome
cat("Number of encoded features:", ncol(x_mat), "\n")

feature_map <- sapply(candidate_predictors, function(var) {
  sum(startsWith(colnames(x_mat), var))
})

feature_summary <- data.frame(
  variable = names(feature_map),
  n_features = as.numeric(feature_map)
)

# Save feature summary for traceability
write_csv_table_only(
  feature_summary,
  "feature_summary_encoded_columns.csv",
  row.names = FALSE
)

# Variables that made it into the model
cat("\nVariables with features:\n")
print(feature_summary %>% filter(n_features > 0))

cat("\nVariables with zero features:\n")
print(feature_summary %>% filter(n_features == 0))

cat("ncol(x_mat):", ncol(x_mat), "\n")
cat("length(selprop):", length(selprop), "\n")
cat("number of names(selprop):", length(names(selprop)), "\n")

summary(df$sex)
summary(df$ethnicity)




# rm(list = ls())
# 
# # =========================================================
# # First-step stability selection for outcome model
# # Model:
# #     Y (CVD) ~ exposures + biomarkers + confounders
# #
# # Purpose:
# #     Identify variables with stable direct associations
# #     with the outcome (CVD). Selected exposures and
# #     biomarkers will be used in second-step mediator models.
# #
# # Notes:
# #     - Confounders are unpenalised
# #     - Exposures and biomarkers are penalised
# # =========================================================
# 
# # Packages -------------------------------------------------
# 
# library(here)
# library(dplyr)
# library(purrr)
# library(ggplot2)
# library(sharp)
# 
# 
# source(here("01_scripts", "00_variable_config.R"))
# source(here("01_scripts", "00_variable_dictionary.R"))
# 
# # Load data
# df <- readRDS(file.path(data_input_dir, selection_file))
# 
# # Basic checks before starting 
# stopifnot("has_cvd" %in% names(df))
# stopifnot(all(df$has_cvd %in% c(0, 1)))
# 
# na_total <- sum(is.na(df))
# cat("Total NA values in dataset:", na_total, "\n")
# 
# 
# 
# # Define output paths -----------------------------------------------------
# 
# out_dir     <- lasso_outcome_dir
# fig_out_dir <- lasso_fig_dir
# 
# 
# 
# # Variable wrangling ------------------------------------------------------
# 
# # Keep only variables that actually exist in data
# confounder_vars <- intersect(confounder_vars, names(df))
# exposure_vars   <- intersect(exposure_vars, names(df))
# bio_vars        <- intersect(bio_vars, names(df))
# 
# # Removes any duplicates
# candidate_predictors <- unique(c(confounder_vars, exposure_vars, bio_vars))
# 
# # Keep only required columns 
# df <- df %>%
#   select(all_of(c(intersect("eid", names(.)), "has_cvd", candidate_predictors)))
# 
# # Save IDs from data set
# if ("eid" %in% names(df)) {
#   write.csv(
#     df %>% select(eid, has_cvd),
#     file.path(out_dir, "selection_dataset_ids.csv"),
#     row.names = FALSE
#   )
# }
# 
# # Remove eid before modelling
# if ("eid" %in% names(df)) {
#   df <- df %>% select(-eid)
# }
# 
# # Remove excluded variables (e.g. anthropometrics in main analysis)
# if (length(exclude_vars) > 0) {
#   exclude_in_data <- intersect(exclude_vars, names(df))
#   if (length(exclude_in_data) > 0) {
#     cat("Removing excluded variables:", paste(exclude_in_data, collapse = ", "), "\n")
#     df <- df %>% select(-all_of(exclude_in_data))
#   }
# }
# 
# # Remove sbp in all scenarios — dominates selection and masks other associations
# if ("sbp" %in% names(df)) {
#   cat("Removing sbp (modelling decision: dominates variable selection)\n")
#   df <- df %>% select(-sbp)
# }
# 
# # =========================================================
# # 2) Prepare data types
# # =========================================================
# 
# # Convert character -> factor
# df <- df %>%
#   mutate(across(where(is.character), as.factor))
# 
# # Convert date -> numeric
# date_cols <- names(df)[vapply(df, inherits, logical(1), what = "Date")]
# if (length(date_cols) > 0) {
#   df <- df %>%
#     mutate(across(all_of(date_cols), ~ as.numeric(.)))
# }
# 
# # Ensure outcome is numeric 0/1
# df$has_cvd <- as.numeric(df$has_cvd)
# 
# # =========================================================
# # 3) One-hot encoding for sharp
# # =========================================================
# 
# # Exposure and outcome dataframes
# x_df <- df %>% select(-has_cvd)
# y <- df$has_cvd
# 
# # model.matrix one-hot encodes factors automatically
# x_mat <- model.matrix(~ ., data = x_df)[, -1, drop = FALSE]
# 
# # Remove zero-variance columns (all same value)
# keep_cols <- apply(x_mat, 2, function(z) stats::sd(z) > 0)
# x_mat <- x_mat[, keep_cols, drop = FALSE]
# 
# # =========================================================
# # 4) Penalty factors: Confounders unpenalised
# # =========================================================
# 
# penalty_factor <- rep(1, ncol(x_mat))
# names(penalty_factor) <- colnames(x_mat)
# 
# # Any dummy column derived from confounders gets penalty set to 0
# for (cf in confounder_vars) {
#   idx <- grep(paste0("^", cf), colnames(x_mat))
#   if (length(idx) > 0) penalty_factor[idx] <- 0
# }
# 
# # =========================================================
# # 5) Run sharp variable selection
# # =========================================================
# 
# # n_cat = 3:
# #   1) stably included
# #   2) stably excluded
# #   3) unstable
# 
# set.seed(123)
# 
# # 100 iterations should suffice
# out <- sharp::VariableSelection(
#   xdata = x_mat,
#   ydata = y,
#   family = "binomial",
#   penalty.factor = penalty_factor,
#   verbose = FALSE,
#   n_cat = 3,
#   K = 100,
#   # Below is for confining the grid to particular region
#   # pi_list = seq(0.65, 0.79, by = 0.01),
#   # Lambda = seq(4e-04, 8e-04, length.out = 60)
# )
# 
# # Save full sharp object
# saveRDS(out, file.path(out_dir, "sharp_variable_selection_object.rds"))
# 
# # Save object names and structure for debugging
# capture.output(
#   {
#     cat("Object names:\n")
#     print(names(out))
#     cat("\nStructure:\n")
#     str(out, max.level = 2)
#   },
#   file = file.path(out_dir, "sharp_object_structure.txt")
# )
# 
# # =========================================================
# # 6) Save object outputs for debugging
# # =========================================================
# 
# # Save full sharp object
# saveRDS(out, file.path(out_dir, "sharp_variable_selection_object.rds"))
# 
# # Save object names / structure for debugging / inspection
# capture.output(
#   {
#     cat("Object names:\n")
#     print(names(out))
#     cat("\nStructure:\n")
#     str(out, max.level = 2)
#   },
#   file = file.path(out_dir, "sharp_object_structure.txt")
# )
# 
# # Helper function for outputting plots
# write_component_csv <- function(obj, comp_name, file_name) {
#   if (comp_name %in% names(obj)) {
#     comp <- obj[[comp_name]]
#     if (is.data.frame(comp) || is.matrix(comp)) {
#       write.csv(as.data.frame(comp), file.path(out_dir, file_name), row.names = FALSE)
#     }
#   }
# }
# 
# write_component_csv(out, "summary_variable", "summary_variable.csv")
# write_component_csv(out, "summary_model", "summary_model.csv")
# write_component_csv(out, "summary_selection", "summary_selection.csv")
# write_component_csv(out, "summary_cat", "summary_cat.csv")
# write_component_csv(out, "Table", "sharp_table.csv")
# write_component_csv(out, "PFER", "sharp_pfer.csv")
# write_component_csv(out, "selprop", "sharp_selprop_matrix.csv")
# write_component_csv(out, "Lambda", "sharp_lambda.csv")
# write_component_csv(out, "Q", "sharp_q.csv")
# write_component_csv(out, "P", "sharp_p.csv")
# write_component_csv(out, "FDP", "sharp_fdp.csv")
# 
# # Calibration plot for viewing lambda and pi patterns
# png(
#   filename = file.path(fig_out_dir, "sharp_calibration_plot.png"),
#   width = 2300,
#   height = 1300,
#   res = 150
# )
# 
# par(mar = c(6.5, 6.5, 5.0, 7.0))
# 
# try(sharp::CalibrationPlot(out), silent = TRUE)
# 
# dev.off()
# 
# # =========================================================
# # 7) Practical calibrated outputs for debugging
# # =========================================================
# 
# selprop <- SelectionProportions(out)
# hat_params <- Argmax(out)
# 
# if (is.null(names(selprop))) {
#   names(selprop) <- colnames(x_mat)
# }
# 
# sel_threshold <- hat_params[2]
# 
# cat("Calibrated selection threshold from sharp:", sel_threshold, "\n")
# 
# write.csv(
#   data.frame(
#     feature = names(selprop),
#     selection_proportion = as.numeric(selprop),
#     stringsAsFactors = FALSE
#   ),
#   file.path(out_dir, "selection_proportions_raw.csv"),
#   row.names = FALSE
# )
# 
# write.csv(
#   data.frame(
#     parameter = c("lambda", "pi"),
#     value = as.numeric(hat_params),
#     stringsAsFactors = FALSE
#   ),
#   file.path(out_dir, "sharp_argmax_params.csv"),
#   row.names = FALSE
# )
# 
# # =========================================================
# # 8) Feature and source variable summaries
# # =========================================================
# 
# selection_table <- data.frame(
#   feature = names(selprop),
#   selection_proportion = as.numeric(selprop),
#   selected = as.numeric(selprop) >= sel_threshold,
#   stringsAsFactors = FALSE
# ) %>%
#   arrange(desc(selection_proportion), feature)
# 
# get_source_variable <- function(feature_name, original_vars) {
#   hits <- original_vars[startsWith(feature_name, original_vars)]
#   if (length(hits) == 0) return(NA_character_)
#   hits[which.max(nchar(hits))]
# }
# 
# selection_table$source_variable <- sapply(
#   selection_table$feature,
#   get_source_variable,
#   original_vars = c(confounder_vars, exposure_vars, bio_vars)
# )
# 
# selection_table$group <- dplyr::case_when(
#   selection_table$source_variable %in% confounder_vars ~ "confounder",
#   selection_table$source_variable %in% bio_vars        ~ "biological",
#   selection_table$source_variable %in% exposure_vars   ~ "exposure",
#   TRUE ~ "unknown"
# )
# 
# write.csv(
#   selection_table,
#   file.path(out_dir, "selection_table_feature_level.csv"),
#   row.names = FALSE
# )
# 
# selected_variables_final <- selection_table %>%
#   filter(!is.na(source_variable)) %>%
#   group_by(source_variable, group) %>%
#   summarise(
#     max_selection_proportion = max(selection_proportion, na.rm = TRUE),
#     selected = any(selected),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(selected), desc(max_selection_proportion), group, source_variable)
# 
# # Outputs selection proportions for all variables
# write.csv(
#   selected_variables_final,
#   file.path(out_dir, "selected_variables_final.csv"),
#   row.names = FALSE
# )
# 
# selected_variables_only <- selected_variables_final %>%
#   filter(selected)
# 
# # Outputs selection proportions for selected variables
# write.csv(
#   selected_variables_only,
#   file.path(out_dir, "selected_variables_only.csv"),
#   row.names = FALSE
# )
# 
# # Adjust group names so more apropriate for mediation analysis
# selected_variable_names <- selected_variables_final %>%
#   filter(selected) %>%
#   transmute(
#     variable = source_variable,
#     type = case_when(
#       group == "biological" ~ "mediator",
#       group == "exposure" ~ "exposure",
#       group == "confounder" ~ "confounder",
#       TRUE ~ group
#     )
#   ) %>%
#   distinct()
# 
# # Clean csv with only selected variable names
# write.csv(
#   selected_variable_names,
#   file.path(out_dir, "selected_variable_names.csv"),
#   row.names = FALSE
# )
# 
# # Same as above just .txt format
# writeLines(
#   paste(selected_variable_names$variable, selected_variable_names$type, sep = ","),
#   file.path(out_dir, "selected_variable_names.txt")
# )
# 
# # =========================================================
# # 9) Plots
# # =========================================================
# 
# # Map coding names to display names
# selprop_labels <- get_display_name(names(selprop))
# 
# # This plot shows selection proportion of all variables
# png(
#   filename = file.path(fig_out_dir, "selection_proportion_plot_calibrated.png"),
#   width = 1800,
#   height = 1100,
#   res = 150
# )
# 
# par(mar = c(14, 5, 2, 1))
# 
# plot(
#   selprop,
#   type = "h",
#   lwd = 3,
#   las = 1,
#   xlab = "",
#   ylab = "Selection Proportion",
#   xaxt = "n",
#   col = ifelse(selprop >= sel_threshold, "red", "#7f8c8d"),
#   cex.lab = 1.5
# )
# 
# abline(h = sel_threshold, lty = 2, col = "darkred")
# 
# for (i in seq_along(selprop)) {
#   axis(
#     side = 1,
#     at = i,
#     labels = selprop_labels[i],
#     las = 2,
#     cex.axis = 0.7,
#     col = ifelse(selprop[i] >= sel_threshold, "red", "#7f8c8d"),
#     col.axis = ifelse(selprop[i] >= sel_threshold, "red", "#7f8c8d")
#   )
# }
# 
# dev.off()
# 
# # Get top 20 variables by selection proportion
# top_n <- 20
# top_idx <- order(selprop, decreasing = TRUE)[1:top_n]
# 
# selprop_top <- selprop[top_idx]
# selprop_top_labels <- selprop_labels[top_idx]
# 
# # Plot the top 20 selection proportions for more concise plot than above
# png(
#   filename = file.path(fig_out_dir, "selection_proportion_plot_top20.png"),
#   width = 1800,
#   height = 1100,
#   res = 150
# )
# 
# par(mar = c(14, 5, 2, 1))
# 
# plot(
#   selprop_top,
#   type = "h",
#   lwd = 4,
#   las = 1,
#   xlab = "",
#   ylab = "Selection Proportion",
#   xaxt = "n",
#   col = ifelse(selprop_top >= sel_threshold, "red", "#7f8c8d"),
#   cex.lab = 1.5
# )
# 
# abline(h = sel_threshold, lty = 2, col = "darkred")
# 
# for (i in seq_along(selprop_top)) {
#   axis(
#     side = 1,
#     at = i,
#     labels = selprop_top_labels[i],
#     las = 2,
#     cex.axis = 0.9,
#     col = ifelse(selprop_top[i] >= sel_threshold, "red", "#7f8c8d"),
#     col.axis = ifelse(selprop_top[i] >= sel_threshold, "red", "#7f8c8d")
#   )
# }
# 
# dev.off()
# 
# # This plot shows selection proportion of all variables ordered
# # Order by decreasing selection proportion
# ord <- order(selprop, decreasing = TRUE)
# 
# selprop_ord <- selprop[ord]
# selprop_labels_ord <- selprop_labels[ord]
# 
# # Ordered selection proportion plot with display names
# png(
#   filename = file.path(fig_out_dir, "selection_proportion_plot_ordered.png"),
#   width = 1800,
#   height = 1100,
#   res = 150
# )
# 
# par(mar = c(14, 5, 2, 1))
# 
# plot(
#   selprop_ord,
#   type = "h",
#   lwd = 3,
#   las = 1,
#   xlab = "",
#   ylab = "Selection Proportion",
#   xaxt = "n",
#   col = ifelse(selprop_ord >= sel_threshold, "red", "#7f8c8d"),
#   cex.lab = 1.5
# )
# 
# abline(h = sel_threshold, lty = 2, col = "darkred")
# 
# for (i in seq_along(selprop_ord)) {
#   axis(
#     side = 1,
#     at = i,
#     labels = selprop_labels_ord[i],
#     las = 2,
#     cex.axis = 0.7,
#     col = ifelse(selprop_ord[i] >= sel_threshold, "red", "#7f8c8d"),
#     col.axis = ifelse(selprop_ord[i] >= sel_threshold, "red", "#7f8c8d")
#   )
# }
# 
# dev.off()
# 
# # =========================================================
# # 10) Console summary for us to view after running
# # =========================================================
# 
# cat("\n====================================================\n")
# cat("sharp variable selection completed\n")
# cat("====================================================\n")
# cat("Selection dataset rows:", nrow(df), "\n")
# cat("Number of encoded features:", ncol(x_mat), "\n")
# cat("n_cat:", 3, "\n")
# cat("Calibrated selection threshold:", sel_threshold, "\n")
# cat("Output directory:\n", out_dir, "\n")
# cat("====================================================\n\n")
# 
# cat("Confounders used:\n")
# print(confounder_vars)
# 
# cat("\nExposure variables used:\n")
# print(exposure_vars)
# 
# cat("\nBiological variables used:\n")
# print(bio_vars)
# 
# cat("\nTop encoded features by selection proportion:\n")
# print(head(selection_table, 20))
# 
# cat("\nRecovered selected source variables:\n")
# print(selected_variables_only)
# 
# cat("Number of original predictors:", length(candidate_predictors), "\n")
# cat("Number of predictors in df:", ncol(df) - 1, "\n")  # minus outcome
# cat("Number of encoded features:", ncol(x_mat), "\n")
# 
# feature_map <- sapply(candidate_predictors, function(var) {
#   sum(startsWith(colnames(x_mat), var))
# })
# 
# feature_summary <- data.frame(
#   variable = names(feature_map),
#   n_features = as.numeric(feature_map)
# )
# 
# # Variables that made it into the model
# cat("\nVariables with features:\n")
# print(feature_summary %>% filter(n_features > 0))
# 
# cat("\nVariables with zero features:\n")
# print(feature_summary %>% filter(n_features == 0))
# 
# cat("ncol(x_mat):", ncol(x_mat), "\n")
# cat("length(selprop):", length(selprop), "\n")
# cat("number of names(selprop):", length(names(selprop)), "\n")
# 
# summary(df$sex)
# summary(df$ethnicity)
