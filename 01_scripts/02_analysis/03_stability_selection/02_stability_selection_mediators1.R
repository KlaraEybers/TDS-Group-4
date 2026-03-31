rm(list = ls())

# =========================================================
# Second-step stability selection for mediator models
# For each selected biomarker M_j (biomarker or anthropometric):
#     M_j ~ exposures + confounders
#
# Purpose:
#     Identify exposures that are stably associated with each
#     selected mediator (biomarker). These associations represent the
#     'a' paths in the mediation framework.
#
# Notes:
#     - Confounders are unpenalised
#     - Exposures are penalised
#     - Used in combination with LASSO 1 to estimate indirect effects
# =========================================================

library(here)
library(dplyr)
library(sharp)

source(here("01_scripts", "00_variable_config.R"))
source(here("01_scripts", "00_variable_dictionary.R"))

# =========================================================
# 0) Paths
# =========================================================

first_step_dir <- lasso_outcome_dir
out_base_dir   <- lasso_mediator_dir
fig_base_dir   <- file.path(lasso_fig_dir, "stability_selection_runs_mediator")
tbl_base_dir   <- file.path(
  here("02_results", "01_tables", "08_lasso", "stability_selection_lasso"),
  "stability_selection_runs_mediator"
)

dir.create(out_base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_base_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 1) Load selection dataset
# =========================================================

df <- readRDS(file.path(data_input_dir, selection_file))

# =========================================================
# 2) Safety checks
# =========================================================

confounder_vars <- intersect(confounder_vars, names(df))
exposure_vars   <- intersect(exposure_vars, names(df))
bio_vars        <- intersect(bio_vars, names(df))

# =========================================================
# 3) Load first-step selected variables and extract mediators
# =========================================================

# Load variabels stably selected from LASSO 1
selected_tbl <- read.csv(
  file.path(first_step_dir, "selected_variables_only.csv"),
  stringsAsFactors = FALSE
)

selected_mediators <- selected_tbl %>%
  filter(group == "biological") %>%
  pull(source_variable) %>%
  unique()

selected_mediators <- intersect(selected_mediators, bio_vars)

# Remove sbp — dominates selection, excluded from all modelling scenarios
selected_mediators <- setdiff(selected_mediators, "sbp")

selected_exposures <- selected_tbl %>%
  filter(group == "exposure") %>%
  pull(source_variable) %>%
  unique()

selected_exposures <- intersect(selected_exposures, exposure_vars)

cat("Selected mediators from first-step stability selection:\n")
print(selected_mediators)

cat("Selected exposures from first-step stability selection:\n")
print(selected_exposures)

if (length(exposure_vars) == 0) {
  stop("No exposure variables found in df")
}

if (length(selected_exposures) == 0) {
  stop("No selected exposure variables found in selected_variables_only.csv")
}

if (length(selected_mediators) == 0) {
  stop("No selected biological mediators found in selected_variables_only.csv")
}

# =========================================================
# 4) Type prep
# =========================================================

df <- df %>%
  mutate(across(where(is.character), as.factor))

df <- df %>%
  mutate(across(where(is.logical), as.factor))

date_cols <- names(df)[vapply(df, inherits, logical(1), what = "Date")]
if (length(date_cols) > 0) {
  df <- df %>%
    mutate(across(all_of(date_cols), ~ as.numeric(.)))
}

df <- df %>%
  mutate(across(where(is.factor), droplevels))

# =========================================================
# 5) Helper functions
# =========================================================

write_csv_safe <- function(x, path, row.names = FALSE) {
  if (is.null(x)) return(invisible(FALSE))
  if (!is.data.frame(x)) x <- as.data.frame(x)
  write.csv(x, path, row.names = row.names)
  invisible(TRUE)
}

write_csv_both <- function(x, file_name, out_dir, tbl_out_dir, row.names = FALSE) {
  write_csv_safe(x, file.path(out_dir, file_name), row.names = row.names)
  write_csv_safe(x, file.path(tbl_out_dir, file_name), row.names = row.names)
}

write_csv_table_only <- function(x, file_name, tbl_out_dir, row.names = FALSE) {
  write_csv_safe(x, file.path(tbl_out_dir, file_name), row.names = row.names)
}

write_component_csv <- function(obj, comp_name, file_name, tbl_out_dir) {
  if (comp_name %in% names(obj)) {
    comp <- obj[[comp_name]]
    if (is.data.frame(comp) || is.matrix(comp) || is.vector(comp)) {
      write_csv_table_only(as.data.frame(comp), file_name, tbl_out_dir)
    }
  }
}

is_binary_numeric <- function(x) {
  ux <- unique(stats::na.omit(x))
  is.numeric(x) && length(ux) <= 2 && all(ux %in% c(0, 1))
}

infer_mediator_family <- function(x, mediator_name) {
  # Factor binary mediator
  if (is.factor(x)) {
    if (nlevels(x) == 2) {
      return("binomial")
    } else {
      stop(
        paste0(
          "Mediator '", mediator_name, "' is a factor with ",
          nlevels(x), " levels. This script only supports continuous or binary mediators."
        )
      )
    }
  }
  
  # Treats mediators differently depending on it's data type
  if (is.numeric(x)) {
    if (is_binary_numeric(x)) {
      return("binomial")
    } else {
      return("gaussian")
    }
  }
  
  stop(
    paste0(
      "Mediator '", mediator_name, "' is neither numeric nor a supported binary factor."
    )
  )
}

get_source_variable <- function(feature_name, original_vars) {
  hits <- original_vars[startsWith(feature_name, original_vars)]
  if (length(hits) == 0) return(NA_character_)
  hits[which.max(nchar(hits))]
}

all_mediator_exposure_pairs <- list()

# =========================================================
# 6) Run separate stability selection for each mediator
# =========================================================

# Loop stability selection for each mediator (we have to run separate LASSO for each mediator)
for (med in selected_mediators) {
  
  cat("Running mediator stability selection for:", med, "\n")
  
  med_family <- infer_mediator_family(df[[med]], med)
  cat("Detected family:", med_family, "\n")
  
  out_dir <- file.path(out_base_dir, paste0("mediator_", med))
  fig_out_dir <- file.path(fig_base_dir, paste0("mediator_", med))
  tbl_out_dir <- file.path(tbl_base_dir, paste0("mediator_", med))
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tbl_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------------------------------------
  # Keep mediator + confounders + selected exposures only
  # -------------------------------------------------------
  
  vars_to_keep <- unique(c(med, confounder_vars, exposure_vars))
  vars_to_keep <- intersect(vars_to_keep, names(df))
  
  df_med <- df %>%
    select(all_of(vars_to_keep))
  
  # Sharp requires complete data
  df_med <- df_med %>%
    filter(if_all(everything(), ~ !is.na(.)))
  
  # Outputs for debugging
  cat("Rows retained after complete-case filtering:", nrow(df_med), "\n")
  
  if (nrow(df_med) == 0) {
    stop(paste("No complete rows available for mediator:", med))
  }
  
  # -------------------------------------------------------
  # Build mediator ande exposure dataframes
  # -------------------------------------------------------
  
  y <- df_med[[med]]
  x_df <- df_med %>% select(-all_of(med))
  
  if (ncol(x_df) == 0) {
    stop(paste("No predictor columns remain for mediator:", med))
  }
  
  # Convert binary factor mediator to 0/1 if needed
  if (med_family == "binomial") {
    if (is.factor(y)) {
      y <- as.numeric(y) - 1
    } else {
      y <- as.numeric(y)
    }
  } else {
    y <- as.numeric(y)
  }
  
  if (length(unique(stats::na.omit(y))) < 2 && med_family == "binomial") {
    stop(paste("Binary mediator has fewer than 2 observed classes after filtering:", med))
  }
  
  # One-hot encode predictors
  x_mat <- model.matrix(~ ., data = x_df)
  
  if ("(Intercept)" %in% colnames(x_mat)) {
    x_mat <- x_mat[, colnames(x_mat) != "(Intercept)", drop = FALSE]
  }
  
  if (ncol(x_mat) == 0) {
    stop(paste("No encoded predictor features created for mediator:", med))
  }
  
  # Remove zero-variance columns (columns with all the same value)
  keep_cols <- apply(x_mat, 2, function(z) {
    z_num <- as.numeric(z)
    stats::sd(z_num, na.rm = TRUE) > 0
  })
  
  if (length(keep_cols) == 0 || !any(keep_cols)) {
    stop(paste("All encoded predictor features removed as zero-variance for mediator:", med))
  }
  
  x_mat <- x_mat[, keep_cols, drop = FALSE]
  
  if (ncol(x_mat) == 0) {
    stop(paste("No encoded predictor features remain after zero-variance filtering for mediator:", med))
  }
  
  cat("Encoded predictor features:", ncol(x_mat), "\n")
  
  # -------------------------------------------------------
  # Penalty factors: confounders unpenalised
  # -------------------------------------------------------
  
  penalty_factor <- rep(1, ncol(x_mat))
  names(penalty_factor) <- colnames(x_mat)
  
  for (cf in confounder_vars) {
    idx <- grep(paste0("^", cf), colnames(x_mat))
    if (length(idx) > 0) penalty_factor[idx] <- 0
  }
  
  cat("Unpenalised features:", sum(penalty_factor == 0), "\n")
  cat("Penalised features:", sum(penalty_factor == 1), "\n")
  
  # -------------------------------------------------------
  # Run sharp stability selection
  # -------------------------------------------------------
  
  set.seed(123)
  
  out <- sharp::VariableSelection(
    xdata = x_mat,
    ydata = y,
    family = med_family,
    penalty.factor = penalty_factor,
    verbose = FALSE,
    n_cat = 3,
    K = 300
  )
  
  # ----------------------------------------------------------------------------------------
  
  # -------------------------------------------------------
  # Save raw object and structure
  # -------------------------------------------------------
  
  saveRDS(out, file.path(out_dir, "sharp_variable_selection_object.rds"))
  
  capture.output(
    {
      cat("Object names:\n")
      print(names(out))
      cat("\nStructure:\n")
      str(out, max.level = 2)
    },
    file = file.path(out_dir, "sharp_object_structure.txt")
  )
  
  # Save common components if present
  write_component_csv(out, "summary_variable", "summary_variable.csv", tbl_out_dir)
  write_component_csv(out, "summary_model", "summary_model.csv", tbl_out_dir)
  write_component_csv(out, "summary_selection", "summary_selection.csv", tbl_out_dir)
  write_component_csv(out, "summary_cat", "summary_cat.csv", tbl_out_dir)
  write_component_csv(out, "Table", "sharp_table.csv", tbl_out_dir)
  write_component_csv(out, "PFER", "sharp_pfer.csv", tbl_out_dir)
  write_component_csv(out, "selprop", "sharp_selprop_matrix.csv", tbl_out_dir)
  write_component_csv(out, "Lambda", "sharp_lambda.csv", tbl_out_dir)
  write_component_csv(out, "Q", "sharp_q.csv", tbl_out_dir)
  write_component_csv(out, "P", "sharp_p.csv", tbl_out_dir)
  write_component_csv(out, "FDP", "sharp_fdp.csv", tbl_out_dir)
  
  # -------------------------------------------------------
  # Calibration and selection summaries
  # -------------------------------------------------------
  
  selprop <- SelectionProportions(out)
  hat_params <- Argmax(out)
  
  if (is.null(names(selprop))) {
    names(selprop) <- colnames(x_mat)
  }
  
  sel_threshold <- NA_real_
  if (length(hat_params) >= 2) {
    sel_threshold <- as.numeric(hat_params[2])
  }
  
  if (is.na(sel_threshold) || !is.finite(sel_threshold)) {
    stop(paste("Could not recover a valid calibrated selection threshold from Argmax(out) for mediator:", med))
  }
  
  cat("Calibrated selection threshold:", sel_threshold, "\n")
  
  selection_proportions_raw <- data.frame(
    feature = names(selprop),
    selection_proportion = as.numeric(selprop),
    stringsAsFactors = FALSE
  )
  
  write_csv_table_only(
    selection_proportions_raw,
    "selection_proportions_raw.csv",
    tbl_out_dir,
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
    tbl_out_dir,
    row.names = FALSE
  )
  
  selection_table <- data.frame(
    feature = names(selprop),
    selection_proportion = as.numeric(selprop),
    selected = as.numeric(selprop) >= sel_threshold,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(selection_proportion), feature)
  
  selection_table$source_variable <- sapply(
    selection_table$feature,
    get_source_variable,
    original_vars = c(confounder_vars, exposure_vars)
  )
  
  selection_table$group <- dplyr::case_when(
    selection_table$source_variable %in% confounder_vars ~ "confounder",
    selection_table$source_variable %in% exposure_vars   ~ "exposure",
    TRUE ~ "unknown"
  )
  
  write_csv_table_only(
    selection_table,
    "selection_table_feature_level.csv",
    tbl_out_dir,
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
  
  # Outputs all variables with selection proportions
  write_csv_both(
    selected_variables_final,
    "selected_variables_final.csv",
    out_dir,
    tbl_out_dir,
    row.names = FALSE
  )
  
  selected_variables_only <- selected_variables_final %>%
    filter(selected)
  
  # Outputs selected variables with selection proportions
  write_csv_both(
    selected_variables_only,
    "selected_variables_only.csv",
    out_dir,
    tbl_out_dir,
    row.names = FALSE
  )
  
  selected_variable_names <- selected_variables_only %>%
    pull(source_variable) %>%
    unique()
  
  # Cleaner file with only selected variable names
  write_csv_both(
    data.frame(variable = selected_variable_names),
    "selected_variable_names.csv",
    out_dir,
    tbl_out_dir,
    row.names = FALSE
  )
  
  # Output each exposure with the mediator it is associated with
  
  selected_variables_only <- selected_variables_final %>%
    filter(selected)
  
  mediator_exposure_pairs <- selected_variables_only %>%
    filter(group == "exposure") %>%
    transmute(
      mediator = med,
      exposure = source_variable
    ) %>%
    distinct()
  
  # Collects all mediator-exposure associations for one mediator
  write_csv_both(
    mediator_exposure_pairs,
    "selected_mediator_exposure_pairs.csv",
    out_dir,
    tbl_out_dir,
    row.names = FALSE
  )
  
  all_mediator_exposure_pairs[[med]] <- mediator_exposure_pairs
  
  feature_to_source_map <- selection_table %>%
    select(feature, source_variable, group, selection_proportion, selected)
  
  write_csv_table_only(
    feature_to_source_map,
    "feature_to_source_map.csv",
    tbl_out_dir,
    row.names = FALSE
  )
  
  feature_map <- sapply(unique(c(confounder_vars, exposure_vars)), function(var) {
    sum(startsWith(colnames(x_mat), var))
  })
  
  feature_summary <- data.frame(
    variable = names(feature_map),
    n_features = as.numeric(feature_map)
  )
  
  write_csv_table_only(
    feature_summary,
    "feature_summary_encoded_columns.csv",
    tbl_out_dir,
    row.names = FALSE
  )
  
  write_csv_table_only(
    data.frame(
      mediator = med,
      mediator_family = med_family,
      rows_retained_complete_case = nrow(df_med),
      predictor_columns_pre_encoding = ncol(x_df),
      encoded_features = ncol(x_mat),
      unpenalised_features = sum(penalty_factor == 0),
      penalised_features = sum(penalty_factor == 1),
      selection_threshold = sel_threshold,
      stringsAsFactors = FALSE
    ),
    "mediator_model_summary.csv",
    tbl_out_dir,
    row.names = FALSE
  )
  
  # Calibration plot for viewing lambda and pi patterns
  png(
    filename = file.path(fig_out_dir, paste0(med, "_sharp_calibration_plot.png")),
    width = 2300,
    height = 1300,
    res = 150
  )
  
  par(mar = c(6.5, 6.5, 5.0, 7.0))
  
  try(sharp::CalibrationPlot(out), silent = TRUE)
  
  dev.off()
  
  # -------------------------------------------------------
  # Plot of selection proportion of exposures for respective mediator
  # -------------------------------------------------------
  
  # Display variable names for plotting instead of coding names
  selprop_labels <- get_display_name(names(selprop))
  
  png(
    filename = file.path(fig_out_dir, paste0(med, "_sharp_plot_calibrated.png")),
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
    main = paste("Mediator Model:", get_display_name(med)),
    xaxt = "n",
    col = ifelse(selprop >= sel_threshold, "red", "#7f8c8d"),
    cex.lab = 1.5,
    cex.main = 1.3
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
  
  # -------------------------------------------------------
  # Ordered plot of all selection proportions
  # -------------------------------------------------------
  
  ord_all <- order(selprop, decreasing = TRUE)
  selprop_ord <- selprop[ord_all]
  selprop_ord_labels <- get_display_name(names(selprop_ord))
  
  png(
    filename = file.path(fig_out_dir, paste0(med, "_sharp_plot_ordered.png")),
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
    main = paste("Ordered Selection Proportions:", get_display_name(med)[1]),
    xaxt = "n",
    col = ifelse(selprop_ord >= sel_threshold, "red", "#7f8c8d"),
    cex.lab = 1.5,
    cex.main = 1.2
  )
  
  abline(h = sel_threshold, lty = 2, col = "darkred")
  
  for (i in seq_along(selprop_ord)) {
    axis(
      side = 1,
      at = i,
      labels = selprop_ord_labels[i],
      las = 2,
      cex.axis = 0.7,
      col = ifelse(selprop_ord[i] >= sel_threshold, "red", "#7f8c8d"),
      col.axis = ifelse(selprop_ord[i] >= sel_threshold, "red", "#7f8c8d")
    )
  }
  
  dev.off()
  
  # -------------------------------------------------------
  # Ordered plot of top 20 selection proportions
  # -------------------------------------------------------
  
  top_n <- min(20, length(selprop))
  
  ord <- order(selprop, decreasing = TRUE)
  selprop_top <- selprop[ord][seq_len(top_n)]
  selprop_top_labels <- get_display_name(names(selprop_top))
  
  png(
    filename = file.path(fig_out_dir, paste0(med, "_sharp_plot_top20_ordered.png")),
    width = 1600,
    height = 1000,
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
    main = paste("Top 20 Selection Proportions:", get_display_name(med)),
    xaxt = "n",
    col = ifelse(selprop_top >= sel_threshold, "red", "#7f8c8d"),
    cex.lab = 1.4,
    cex.main = 1.2
  )
  
  abline(h = sel_threshold, lty = 2, col = "darkred")
  
  for (i in seq_along(selprop_top)) {
    axis(
      side = 1,
      at = i,
      labels = selprop_top_labels[i],
      las = 2,
      cex.axis = 0.85,
      col = ifelse(selprop_top[i] >= sel_threshold, "red", "#7f8c8d"),
      col.axis = ifelse(selprop_top[i] >= sel_threshold, "red", "#7f8c8d")
    )
  }
  
  dev.off()
  
  cat("\n====================================================\n")
  cat("Mediator stability selection completed for:", med, "\n")
  cat("====================================================\n")
  cat("Rows retained:", nrow(df_med), "\n")
  cat("Encoded features:", ncol(x_mat), "\n")
  cat("Mediator family:", med_family, "\n")
  cat("Calibrated selection threshold:", sel_threshold, "\n")
  cat("Data directory:\n", out_dir, "\n")
  cat("Figure directory:\n", fig_out_dir, "\n")
  cat("Table directory:\n", tbl_out_dir, "\n")
  cat("====================================================\n\n")
  
}  

# -------------------------------------------------------
# Combine all mediator-exposure associations into one master file
# -------------------------------------------------------

all_mediator_exposure_pairs_df <- bind_rows(all_mediator_exposure_pairs) %>%
  distinct() %>%
  arrange(mediator, exposure)

# Final mediator-exposure associations combined file
write.csv(
  all_mediator_exposure_pairs_df,
  file.path(out_base_dir, "all_selected_mediator_exposure_pairs.csv"),
  row.names = FALSE
)

write.csv(
  all_mediator_exposure_pairs_df,
  file.path(tbl_base_dir, "all_selected_mediator_exposure_pairs.csv"),
  row.names = FALSE
)


# rm(list = ls())
# 
# # =========================================================
# # Second-step stability selection for mediator models
# # For each selected biomarker M_j (biomarker or anthropometric):
# #     M_j ~ exposures + confounders
# #
# # Purpose:
# #     Identify exposures that are stably associated with each
# #     selected mediator (biomarker). These associations represent the
# #     'a' paths in the mediation framework.
# #
# # Notes:
# #     - Confounders are unpenalised
# #     - Exposures are penalised
# #     - Used in combination with LASSO 1 to estimate indirect effects
# # =========================================================
# 
# library(here)
# library(dplyr)
# library(sharp)
# 
# source(here("01_scripts", "00_variable_config.R"))
# source(here("01_scripts", "00_variable_dictionary.R"))
# 
# # =========================================================
# # 0) Paths
# # =========================================================
# 
# first_step_dir <- lasso_outcome_dir
# out_base_dir   <- lasso_mediator_dir
# 
# # =========================================================
# # 1) Load selection dataset
# # =========================================================
# 
# df <- readRDS(file.path(data_input_dir, selection_file))
# 
# 
# # =========================================================
# # 2) Safety checks
# # =========================================================
# 
# confounder_vars <- intersect(confounder_vars, names(df))
# exposure_vars   <- intersect(exposure_vars, names(df))
# bio_vars        <- intersect(bio_vars, names(df))
# 
# # =========================================================
# # 3) Load first-step selected variables and extract mediators
# # =========================================================
# 
# # Load variabels stably selected from LASSO 1
# selected_tbl <- read.csv(
#   file.path(first_step_dir, "selected_variables_only.csv"),
#   stringsAsFactors = FALSE
# )
# 
# selected_mediators <- selected_tbl %>%
#   filter(group == "biological") %>%
#   pull(source_variable) %>%
#   unique()
# 
# selected_mediators <- intersect(selected_mediators, bio_vars)
# 
# # Remove sbp — dominates selection, excluded from all modelling scenarios
# selected_mediators <- setdiff(selected_mediators, "sbp")
# 
# selected_exposures <- selected_tbl %>%
#   filter(group == "exposure") %>%
#   pull(source_variable) %>%
#   unique()
# 
# selected_exposures <- intersect(selected_exposures, exposure_vars)
# 
# cat("Selected mediators from first-step stability selection:\n")
# print(selected_mediators)
# 
# cat("Selected exposures from first-step stability selection:\n")
# print(selected_exposures)
# 
# if (length(exposure_vars) == 0) {
#   stop("No exposure variables found in df")
# }
# 
# if (length(selected_exposures) == 0) {
#   stop("No selected exposure variables found in selected_variables_only.csv")
# }
# 
# # =========================================================
# # 4) Type prep
# # =========================================================
# 
# df <- df %>%
#   mutate(across(where(is.character), as.factor))
# 
# date_cols <- names(df)[vapply(df, inherits, logical(1), what = "Date")]
# if (length(date_cols) > 0) {
#   df <- df %>%
#     mutate(across(all_of(date_cols), ~ as.numeric(.)))
# }
# 
# # =========================================================
# # 5) Helper functions
# # =========================================================
# 
# is_binary_numeric <- function(x) {
#   ux <- unique(stats::na.omit(x))
#   is.numeric(x) && length(ux) <= 2 && all(ux %in% c(0, 1))
# }
# 
# infer_mediator_family <- function(x, mediator_name) {
#   # Factor binary mediator
#   if (is.factor(x)) {
#     if (nlevels(x) == 2) {
#       return("binomial")
#     } else {
#       stop(
#         paste0(
#           "Mediator '", mediator_name, "' is a factor with ",
#           nlevels(x), " levels. This script only supports continuous or binary mediators."
#         )
#       )
#     }
#   }
#   
#   # Treats mediators differently depending on it's data type
#   if (is.numeric(x)) {
#     if (is_binary_numeric(x)) {
#       return("binomial")
#     } else {
#       return("gaussian")
#     }
#   }
#   
#   stop(
#     paste0(
#       "Mediator '", mediator_name, "' is neither numeric nor a supported binary factor."
#     )
#   )
# }
# 
# get_source_variable <- function(feature_name, original_vars) {
#   hits <- original_vars[startsWith(feature_name, original_vars)]
#   if (length(hits) == 0) return(NA_character_)
#   hits[which.max(nchar(hits))]
# }
# 
# write_component_csv <- function(obj, comp_name, file_name, out_dir) {
#   if (comp_name %in% names(obj)) {
#     comp <- obj[[comp_name]]
#     if (is.data.frame(comp) || is.matrix(comp)) {
#       write.csv(as.data.frame(comp), file.path(out_dir, file_name), row.names = FALSE)
#     }
#   }
# }
# 
# # =========================================================
# # 6) Run separate stability selection for each mediator
# # =========================================================
# 
# # Loop stability selection for each mediator (we have to run separate LASSO for each mediator)
# for (med in selected_mediators) {
#   
#   cat("Running mediator stability selection for:", med, "\n")
#   
#   med_family <- infer_mediator_family(df[[med]], med)
#   cat("Detected family:", med_family, "\n")
#   
#   out_dir <- file.path(out_base_dir, paste0("mediator_", med))
#   dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#   
#   # -------------------------------------------------------
#   # Keep mediator + confounders + selected exposures only
#   # -------------------------------------------------------
#   
#   vars_to_keep <- unique(c(med, confounder_vars, exposure_vars))
#   vars_to_keep <- intersect(vars_to_keep, names(df))
#   
#   df_med <- df %>%
#     select(all_of(vars_to_keep))
#   
#   # Sharp requires complete data
#   df_med <- df_med %>%
#     filter(if_all(everything(), ~ !is.na(.)))
#   
#   # Outputs for debugging
#   cat("Rows retained after complete-case filtering:", nrow(df_med), "\n")
#   
#   if (nrow(df_med) == 0) {
#     stop(paste("No complete rows available for mediator:", med))
#   }
#   
#   # -------------------------------------------------------
#   # Build mediator ande exposure dataframes
#   # -------------------------------------------------------
#   
#   y <- df_med[[med]]
#   x_df <- df_med %>% select(-all_of(med))
#   
#   # Convert binary factor mediator to 0/1 if needed
#   if (med_family == "binomial") {
#     if (is.factor(y)) {
#       y <- as.numeric(y) - 1
#     } else {
#       y <- as.numeric(y)
#     }
#   } else {
#     y <- as.numeric(y)
#   }
#   
#   # One-hot encode predictors
#   x_mat <- model.matrix(~ ., data = x_df)[, -1, drop = FALSE]
#   
#   # Remove zero-variance columns (columns with all the same value)
#   keep_cols <- apply(x_mat, 2, function(z) stats::sd(z) > 0)
#   x_mat <- x_mat[, keep_cols, drop = FALSE]
#   
#   cat("Encoded predictor features:", ncol(x_mat), "\n")
#   
#   # -------------------------------------------------------
#   # Penalty factors: confounders unpenalised
#   # -------------------------------------------------------
#   
#   penalty_factor <- rep(1, ncol(x_mat))
#   names(penalty_factor) <- colnames(x_mat)
#   
#   for (cf in confounder_vars) {
#     idx <- grep(paste0("^", cf), colnames(x_mat))
#     if (length(idx) > 0) penalty_factor[idx] <- 0
#   }
#   
#   cat("Unpenalised features:", sum(penalty_factor == 0), "\n")
#   cat("Penalised features:", sum(penalty_factor == 1), "\n")
#   
#   # -------------------------------------------------------
#   # Run sharp stability selection
#   # -------------------------------------------------------
#   
#   set.seed(123)
#   
#   out <- sharp::VariableSelection(
#     xdata = x_mat,
#     ydata = y,
#     family = med_family,
#     penalty.factor = penalty_factor,
#     verbose = FALSE,
#     n_cat = 3,
#     K = 100
#   )
#   
#   # ----------------------------------------------------------------------------------------
#   
#   # -------------------------------------------------------
#   # Save raw object and structure
#   # -------------------------------------------------------
#   
#   saveRDS(out, file.path(out_dir, "sharp_variable_selection_object.rds"))
#   
#   capture.output(
#     {
#       cat("Object names:\n")
#       print(names(out))
#       cat("\nStructure:\n")
#       str(out, max.level = 2)
#     },
#     file = file.path(out_dir, "sharp_object_structure.txt")
#   )
#   
#   # Save common components if present
#   write_component_csv(out, "summary_variable", "summary_variable.csv", out_dir)
#   write_component_csv(out, "summary_model", "summary_model.csv", out_dir)
#   write_component_csv(out, "summary_selection", "summary_selection.csv", out_dir)
#   write_component_csv(out, "summary_cat", "summary_cat.csv", out_dir)
#   write_component_csv(out, "Table", "sharp_table.csv", out_dir)
#   write_component_csv(out, "PFER", "sharp_pfer.csv", out_dir)
#   write_component_csv(out, "selprop", "sharp_selprop_matrix.csv", out_dir)
#   write_component_csv(out, "Lambda", "sharp_lambda.csv", out_dir)
#   write_component_csv(out, "Q", "sharp_q.csv", out_dir)
#   write_component_csv(out, "P", "sharp_p.csv", out_dir)
#   write_component_csv(out, "FDP", "sharp_fdp.csv", out_dir)
#   
#   # -------------------------------------------------------
#   # Calibration and selection summaries
#   # -------------------------------------------------------
#   
#   selprop <- SelectionProportions(out)
#   hat_params <- Argmax(out)
#   
#   if (is.null(names(selprop))) {
#     names(selprop) <- colnames(x_mat)
#   }
#   
#   sel_threshold <- hat_params[2]
#   
#   cat("Calibrated selection threshold:", sel_threshold, "\n")
#   
#   write.csv(
#     data.frame(
#       feature = names(selprop),
#       selection_proportion = as.numeric(selprop),
#       stringsAsFactors = FALSE
#     ),
#     file.path(out_dir, "selection_proportions_raw.csv"),
#     row.names = FALSE
#   )
#   
#   write.csv(
#     data.frame(
#       parameter = c("lambda", "pi"),
#       value = as.numeric(hat_params),
#       stringsAsFactors = FALSE
#     ),
#     file.path(out_dir, "sharp_argmax_params.csv"),
#     row.names = FALSE
#   )
#   
#   selection_table <- data.frame(
#     feature = names(selprop),
#     selection_proportion = as.numeric(selprop),
#     selected = as.numeric(selprop) >= sel_threshold,
#     stringsAsFactors = FALSE
#   ) %>%
#     arrange(desc(selection_proportion), feature)
#   
#   selection_table$source_variable <- sapply(
#     selection_table$feature,
#     get_source_variable,
#     original_vars = c(confounder_vars, exposure_vars)
#   )
#   
#   selection_table$group <- dplyr::case_when(
#     selection_table$source_variable %in% confounder_vars ~ "confounder",
#     selection_table$source_variable %in% exposure_vars   ~ "exposure",
#     TRUE ~ "unknown"
#   )
#   
#   write.csv(
#     selection_table,
#     file.path(out_dir, "selection_table_feature_level.csv"),
#     row.names = FALSE
#   )
#   
#   selected_variables_final <- selection_table %>%
#     filter(!is.na(source_variable)) %>%
#     group_by(source_variable, group) %>%
#     summarise(
#       max_selection_proportion = max(selection_proportion, na.rm = TRUE),
#       selected = any(selected),
#       .groups = "drop"
#     ) %>%
#     arrange(desc(selected), desc(max_selection_proportion), group, source_variable)
#   
#   # Outputs all variables with selection proportions
#   write.csv(
#     selected_variables_final,
#     file.path(out_dir, "selected_variables_final.csv"),
#     row.names = FALSE
#   )
#   
#   selected_variables_only <- selected_variables_final %>%
#     filter(selected)
#   
#   # Outputs selected variables with selection proportions
#   write.csv(
#     selected_variables_only,
#     file.path(out_dir, "selected_variables_only.csv"),
#     row.names = FALSE
#   )
#   
#   selected_variable_names <- selected_variables_only %>%
#     pull(source_variable) %>%
#     unique()
#   
#   # Cleaner file with only selected variable names
#   write.csv(
#     data.frame(variable = selected_variable_names),
#     file.path(out_dir, "selected_variable_names.csv"),
#     row.names = FALSE
#   )
#   
#   # Output each exposure with the mediator it is associated with
#   
#   selected_variables_only <- selected_variables_final %>%
#     filter(selected)
#   
#   mediator_exposure_pairs <- selected_variables_only %>%
#     filter(group == "exposure") %>%
#     transmute(
#       mediator = med,
#       exposure = source_variable
#     ) %>%
#     distinct()
#   
#   # Collects all mediator-exposure associations for one mediator
#   write.csv(
#     mediator_exposure_pairs,
#     file.path(out_dir, "selected_mediator_exposure_pairs.csv"),
#     row.names = FALSE
#   )
#   
#   all_mediator_exposure_pairs <- list()
#   all_mediator_exposure_pairs[[med]] <- mediator_exposure_pairs
#   all_mediator_exposure_pairs_df <- bind_rows(all_mediator_exposure_pairs)
#   
#   write.csv(
#     all_mediator_exposure_pairs_df,
#     file.path(out_base_dir, "all_selected_mediator_exposure_pairs.csv"),
#     row.names = FALSE
#   )
#   
#   # -------------------------------------------------------
#   # Plot of selection proportion of exposures for respective mediator
#   # -------------------------------------------------------
#   
#   # Display variable names for plotting instead of coding names
#   selprop_labels <- get_display_name(names(selprop))
#   
#   png(
#     filename = file.path(lasso_fig_dir, paste0(med, "_sharp_plot_calibrated.png")),
#     width = 1800,
#     height = 1100,
#     res = 150
#   )
#   
#   par(mar = c(14, 5, 2, 1))
#   
#   plot(
#     selprop,
#     type = "h",
#     lwd = 3,
#     las = 1,
#     xlab = "",
#     ylab = "Selection Proportion",
#     xaxt = "n",
#     col = ifelse(selprop >= sel_threshold, "red", "#7f8c8d"),
#     cex.lab = 1.5
#   )
#   
#   abline(h = sel_threshold, lty = 2, col = "darkred")
#   
#   for (i in seq_along(selprop)) {
#     axis(
#       side = 1,
#       at = i,
#       labels = selprop_labels[i],
#       las = 2,
#       cex.axis = 0.7,
#       col = ifelse(selprop[i] >= sel_threshold, "red", "#7f8c8d"),
#       col.axis = ifelse(selprop[i] >= sel_threshold, "red", "#7f8c8d")
#     )
#   }
#   
#   dev.off()
#   
#   
#   # -------------------------------------------------------
#   # Ordered plot of top 20 selection proportions
#   # -------------------------------------------------------
#   
#   top_n <- min(20, length(selprop))
#   
#   ord <- order(selprop, decreasing = TRUE)
#   selprop_top <- selprop[ord][1:top_n]
#   selprop_top_labels <- get_display_name(names(selprop_top))
#   
#   png(
#     filename = file.path(lasso_fig_dir, paste0(med, "_sharp_plot_top20_ordered.png")),
#     width = 1600,
#     height = 1000,
#     res = 150
#   )
#   
#   par(mar = c(14, 5, 3, 1))
#   
#   plot(
#     selprop_top,
#     type = "h",
#     lwd = 4,
#     las = 1,
#     xlab = "",
#     ylab = "Selection Proportion",
#     xaxt = "n",
#     col = ifelse(selprop_top >= sel_threshold, "red", "#7f8c8d"),
#     cex.lab = 1.4,
#     main = "Top 20 Variables by Selection Proportion"
#   )
#   
#   abline(h = sel_threshold, lty = 2, col = "darkred")
#   
#   for (i in seq_along(selprop_top)) {
#     axis(
#       side = 1,
#       at = i,
#       labels = selprop_top_labels[i],
#       las = 2,
#       cex.axis = 0.85,
#       col = ifelse(selprop_top[i] >= sel_threshold, "red", "#7f8c8d"),
#       col.axis = ifelse(selprop_top[i] >= sel_threshold, "red", "#7f8c8d")
#     )
#   }
#   
#   dev.off()
# 
# }  
# 
# # -------------------------------------------------------
# # Combine all mediator-exposure associations into one master file
# # -------------------------------------------------------
# 
# mediator_dirs <- list.dirs(out_base_dir, recursive = FALSE, full.names = TRUE)
# 
# all_mediator_exposure_pairs <- lapply(mediator_dirs, function(dir_path) {
#   
#   file_path <- file.path(dir_path, "selected_mediator_exposure_pairs.csv")
#   
#   if (file.exists(file_path)) {
#     read.csv(file_path, stringsAsFactors = FALSE)
#   } else {
#     NULL
#   }
# })
# 
# all_mediator_exposure_pairs_df <- bind_rows(all_mediator_exposure_pairs) %>%
#   distinct() %>%
#   arrange(mediator, exposure)
# 
# # Final mediator-exposure associations combined file
# write.csv(
#   all_mediator_exposure_pairs_df,
#   file.path(out_base_dir, "all_selected_mediator_exposure_pairs.csv"),
#   row.names = FALSE
# )
# 
