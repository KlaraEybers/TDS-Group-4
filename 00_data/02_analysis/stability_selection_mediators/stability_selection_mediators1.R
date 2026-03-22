rm(list = ls())

# =========================================================
# Second-step stability selection for mediator models
# For each selected biomarker M_j:
#     M_j ~ exposures + confounders
# =========================================================

library(here)
library(dplyr)
library(sharp)

# =========================================================
# 0) Paths
# =========================================================

first_step_dir <- here(
  "00_data", "02_analysis", "stability_selection_runs",
  "stability_selection_sharp_wo_sbp_FINAL"
)

out_base_dir <- here(
  "00_data", "02_analysis", "stability_selection_mediators_run2"
)
dir.create(out_base_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 1) Load selection dataset
# =========================================================

df <- readRDS(here("00_data", "01_processed", "03_imputation", "df_selection_imputed.rds"))

# =========================================================
# 2) Define variable groups
# Must match your modelling setup
# =========================================================

confounder_vars <- c(
  "age", "sex", "ethnicity"
)

exposure_vars <- c(
  "household_income", "highest_qualification", "employment", "housing", "household_size",
  "diet_score", "alcohol_consumption", "salt_added",
  "sleep_hrs", "insomnia",
  "smoking_status_refined", "pack_yrs", "passive_total_hrs",
  "met_min_week_all",
  "mh_sought_help", "work_job_satisfaction_num", "neuroticism_score",
  "breastfed", "comparative_body_size_10", "maternal_smoking_birth", "birth_weight",
  "no2_2010", "pm2.5_2010", "pm10_2010", "urban_rural_3lvl"
)

bio_vars <- c(
  "sbp", "dbp", "grip_strength", "whr", "heel_bone_density", "body_fat_mass",
  "wbc", "platelets", "albumin", "alp", "alt", "apoa1", "apob", "ast",
  "bilirubin", "total_cholesterol", "creatinine", "crp", "cystatin_c",
  "ggt", "hba1c", "hdl_cholesterol", "ldl_direct", "lpa", "triglycerides",
  "long_ill", "diabetes_doctor_diagnosed", "serious_dx",
  "age_high_bp_dx", "age_diabetes_dx", "angina_dx_age", "mi_dx_age",
  "dvt_dx_age", "pe_dx_age", "stroke_dx_age",
  "heart_problems_collapsed", "meds_chol_bp_diabetes_hormone_collapsed",
  "meds_chol_bp_diabetes_collapsed", "neutrophill_count", "lymphocyte_count",
  "monocyte_count"
)

confounder_vars <- intersect(confounder_vars, names(df))
exposure_vars   <- intersect(exposure_vars, names(df))
bio_vars        <- intersect(bio_vars, names(df))

# =========================================================
# 3) Load first-step selected variables and extract mediators
# =========================================================

selected_tbl <- read.csv(
  file.path(first_step_dir, "selected_variables_only.csv"),
  stringsAsFactors = FALSE
)

selected_mediators <- selected_tbl %>%
  filter(group == "biological") %>%
  pull(source_variable) %>%
  unique()

selected_mediators <- intersect(selected_mediators, bio_vars)

cat("Selected mediators from first-step stability selection:\n")
print(selected_mediators)

if (length(selected_mediators) == 0) {
  stop("No selected biological mediators found in selected_variables_only.csv")
}

# =========================================================
# 4) Type prep
# =========================================================

df <- df %>%
  mutate(across(where(is.character), as.factor))

date_cols <- names(df)[vapply(df, inherits, logical(1), what = "Date")]
if (length(date_cols) > 0) {
  df <- df %>%
    mutate(across(all_of(date_cols), ~ as.numeric(.)))
}

# =========================================================
# 5) Helper functions
# =========================================================

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
  
  # Numeric mediator
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

write_component_csv <- function(obj, comp_name, file_name, out_dir) {
  if (comp_name %in% names(obj)) {
    comp <- obj[[comp_name]]
    if (is.data.frame(comp) || is.matrix(comp)) {
      write.csv(as.data.frame(comp), file.path(out_dir, file_name), row.names = FALSE)
    }
  }
}

# =========================================================
# 6) Run separate stability selection for each mediator
# =========================================================

for (med in selected_mediators) {
  
  cat("\n====================================================\n")
  cat("Running mediator stability selection for:", med, "\n")
  cat("====================================================\n")
  
  med_family <- infer_mediator_family(df[[med]], med)
  cat("Detected family:", med_family, "\n")
  
  out_dir <- file.path(out_base_dir, paste0("mediator_", med))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------------------------------------
  # Keep mediator + confounders + exposures
  # -------------------------------------------------------
  
  vars_to_keep <- unique(c(med, confounder_vars, exposure_vars))
  vars_to_keep <- intersect(vars_to_keep, names(df))
  
  df_med <- df %>%
    select(all_of(vars_to_keep))
  
  # sharp requires complete xdata/ydata
  df_med <- df_med %>%
    filter(if_all(everything(), ~ !is.na(.)))
  
  cat("Rows retained after complete-case filtering:", nrow(df_med), "\n")
  
  if (nrow(df_med) == 0) {
    stop(paste("No complete rows available for mediator:", med))
  }
  
  # -------------------------------------------------------
  # Build y and x
  # -------------------------------------------------------
  
  y <- df_med[[med]]
  x_df <- df_med %>% select(-all_of(med))
  
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
  
  # One-hot encode predictors
  x_mat <- model.matrix(~ ., data = x_df)[, -1, drop = FALSE]
  
  # Remove zero-variance columns
  keep_cols <- apply(x_mat, 2, function(z) stats::sd(z) > 0)
  x_mat <- x_mat[, keep_cols, drop = FALSE]
  
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
    K = 100
  )
  
  # -------------------------------------------------------
  # Save raw object + structure
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
  write_component_csv(out, "summary_variable", "summary_variable.csv", out_dir)
  write_component_csv(out, "summary_model", "summary_model.csv", out_dir)
  write_component_csv(out, "summary_selection", "summary_selection.csv", out_dir)
  write_component_csv(out, "summary_cat", "summary_cat.csv", out_dir)
  write_component_csv(out, "Table", "sharp_table.csv", out_dir)
  write_component_csv(out, "PFER", "sharp_pfer.csv", out_dir)
  write_component_csv(out, "selprop", "sharp_selprop_matrix.csv", out_dir)
  write_component_csv(out, "Lambda", "sharp_lambda.csv", out_dir)
  write_component_csv(out, "Q", "sharp_q.csv", out_dir)
  write_component_csv(out, "P", "sharp_p.csv", out_dir)
  write_component_csv(out, "FDP", "sharp_fdp.csv", out_dir)
  
  # -------------------------------------------------------
  # Calibration and selection summaries
  # -------------------------------------------------------
  
  selprop <- SelectionProportions(out)
  hat_params <- Argmax(out)
  
  if (is.null(names(selprop))) {
    names(selprop) <- colnames(x_mat)
  }
  
  sel_threshold <- hat_params[2]
  
  cat("Calibrated selection threshold:", sel_threshold, "\n")
  
  write.csv(
    data.frame(
      feature = names(selprop),
      selection_proportion = as.numeric(selprop),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "selection_proportions_raw.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      parameter = c("lambda", "pi"),
      value = as.numeric(hat_params),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "sharp_argmax_params.csv"),
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
  
  write.csv(
    selection_table,
    file.path(out_dir, "selection_table_feature_level.csv"),
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
  
  write.csv(
    selected_variables_final,
    file.path(out_dir, "selected_variables_final.csv"),
    row.names = FALSE
  )
  
  selected_variables_only <- selected_variables_final %>%
    filter(selected)
  
  write.csv(
    selected_variables_only,
    file.path(out_dir, "selected_variables_only.csv"),
    row.names = FALSE
  )
  
  selected_variable_names <- selected_variables_only %>%
    pull(source_variable) %>%
    unique()
  
  write.csv(
    data.frame(variable = selected_variable_names),
    file.path(out_dir, "selected_variable_names.csv"),
    row.names = FALSE
  )
  
  # -------------------------------------------------------
  # Plot
  # -------------------------------------------------------
  
  png(
    filename = file.path(out_dir, "sharp_plot_calibrated.png"),
    width = 1800,
    height = 1100,
    res = 150
  )
  
  par(mar = c(12, 5, 2, 1))
  
  plot(
    selprop,
    type = "h",
    lwd = 3,
    las = 1,
    xlab = "",
    ylab = "Selection Proportion",
    xaxt = "n",
    col = ifelse(selprop >= sel_threshold, "red", "grey"),
    cex.lab = 1.5
  )
  
  abline(h = sel_threshold, lty = 2, col = "darkred")
  
  for (i in seq_along(selprop)) {
    axis(
      side = 1,
      at = i,
      labels = names(selprop)[i],
      las = 2,
      cex.axis = 0.7,
      col = ifelse(selprop[i] >= sel_threshold, "red", "grey"),
      col.axis = ifelse(selprop[i] >= sel_threshold, "red", "grey")
    )
  }
  
  dev.off()
  
  # -------------------------------------------------------
  # Console summary
  # -------------------------------------------------------
  
  cat("Mediator:", med, "\n")
  cat("Family used:", med_family, "\n")
  cat("Rows used:", nrow(df_med), "\n")
  cat("Selected exposures/confounders:\n")
  print(selected_variables_only)
}





# rm(list = ls())
# 
# # =========================================================
# # Second-step stability selection for mediator models
# # For each selected biomarker M_j:
# #     M_j ~ exposures + confounders
# # =========================================================
# 
# library(here)
# library(dplyr)
# library(sharp)
# 
# # =========================================================
# # 0) Paths
# # =========================================================
# 
# first_step_dir <- here(
#   "00_data", "02_analysis", "stability_selection_runs",
#   "stability_selection_sharp_wo_sbp_FINAL"
# )
# 
# out_base_dir <- here(
#   "00_data", "02_analysis", "stability_selection_mediators_run2"
# )
# dir.create(out_base_dir, recursive = TRUE, showWarnings = FALSE)
# 
# # =========================================================
# # 1) Load selection dataset
# # =========================================================
# 
# df <- readRDS(here("00_data", "01_processed", "03_imputation", "df_selection_imputed.rds"))
# 
# # =========================================================
# # 2) Define variable groups
# # Must match your modelling setup
# # =========================================================
# 
# confounder_vars <- c(
#   "age", "sex", "ethnicity"
# )
# 
# exposure_vars <- c(
#   "household_income", "highest_qualification", "employment", "housing", "household_size",
#   "diet_score", "alcohol_consumption", "salt_added",
#   "sleep_hrs", "insomnia",
#   "smoking_status_refined", "pack_yrs", "passive_total_hrs",
#   "met_min_week_all",
#   "mh_sought_help", "work_job_satisfaction_num", "neuroticism_score",
#   "breastfed", "comparative_body_size_10", "maternal_smoking_birth", "birth_weight",
#   "no2_2010", "pm2.5_2010", "pm10_2010", "urban_rural_3lvl"
# )
# 
# bio_vars <- c(
#   "sbp", "dbp", "grip_strength", "whr", "heel_bone_density", "body_fat_mass",
#   "wbc", "platelets", "albumin", "alp", "alt", "apoa1", "apob", "ast",
#   "bilirubin", "total_cholesterol", "creatinine", "crp", "cystatin_c",
#   "ggt", "hba1c", "hdl_cholesterol", "ldl_direct", "lpa", "triglycerides",
#   "long_ill", "diabetes_doctor_diagnosed", "serious_dx",
#   "age_high_bp_dx", "age_diabetes_dx", "angina_dx_age", "mi_dx_age",
#   "dvt_dx_age", "pe_dx_age", "stroke_dx_age",
#   "heart_problems_collapsed", "meds_chol_bp_diabetes_hormone_collapsed",
#   "meds_chol_bp_diabetes_collapsed", "neutrophill_count", "lymphocyte_count",
#   "monocyte_count"
# )
# 
# confounder_vars <- intersect(confounder_vars, names(df))
# exposure_vars   <- intersect(exposure_vars, names(df))
# bio_vars        <- intersect(bio_vars, names(df))
# 
# # =========================================================
# # 3) Load first-step selected variables and extract mediators
# # =========================================================
# 
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
# cat("Selected mediators from first-step stability selection:\n")
# print(selected_mediators)
# 
# if (length(selected_mediators) == 0) {
#   stop("No selected biological mediators found in selected_variables_only.csv")
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
#   # Numeric mediator
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
# for (med in selected_mediators) {
#   
#   cat("\n====================================================\n")
#   cat("Running mediator stability selection for:", med, "\n")
#   cat("====================================================\n")
#   
#   med_family <- infer_mediator_family(df[[med]], med)
#   cat("Detected family:", med_family, "\n")
#   
#   out_dir <- file.path(out_base_dir, paste0("mediator_", med))
#   dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#   
#   # -------------------------------------------------------
#   # Keep mediator + confounders + exposures
#   # -------------------------------------------------------
#   
#   vars_to_keep <- unique(c(med, confounder_vars, exposure_vars))
#   vars_to_keep <- intersect(vars_to_keep, names(df))
#   
#   df_med <- df %>%
#     select(all_of(vars_to_keep))
#   
#   # sharp requires complete xdata/ydata
#   df_med <- df_med %>%
#     filter(if_all(everything(), ~ !is.na(.)))
#   
#   cat("Rows retained after complete-case filtering:", nrow(df_med), "\n")
#   
#   if (nrow(df_med) == 0) {
#     stop(paste("No complete rows available for mediator:", med))
#   }
#   
#   # -------------------------------------------------------
#   # Build y and x
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
#   # Remove zero-variance columns
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
#   # -------------------------------------------------------
#   # Save raw object + structure
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
#   write.csv(
#     selected_variables_final,
#     file.path(out_dir, "selected_variables_final.csv"),
#     row.names = FALSE
#   )
#   
#   selected_variables_only <- selected_variables_final %>%
#     filter(selected)
#   
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
#   write.csv(
#     data.frame(variable = selected_variable_names),
#     file.path(out_dir, "selected_variable_names.csv"),
#     row.names = FALSE
#   )
#   
#   # -------------------------------------------------------
#   # Plot
#   # -------------------------------------------------------
#   
#   png(
#     filename = file.path(out_dir, "sharp_plot_calibrated.png"),
#     width = 1800,
#     height = 1100,
#     res = 150
#   )
#   
#   par(mar = c(12, 5, 2, 1))
#   
#   plot(
#     selprop,
#     type = "h",
#     lwd = 3,
#     las = 1,
#     xlab = "",
#     ylab = "Selection Proportion",
#     xaxt = "n",
#     col = ifelse(selprop >= sel_threshold, "red", "grey"),
#     cex.lab = 1.5
#   )
#   
#   abline(h = sel_threshold, lty = 2, col = "darkred")
#   
#   for (i in seq_along(selprop)) {
#     axis(
#       side = 1,
#       at = i,
#       labels = names(selprop)[i],
#       las = 2,
#       cex.axis = 0.7,
#       col = ifelse(selprop[i] >= sel_threshold, "red", "grey"),
#       col.axis = ifelse(selprop[i] >= sel_threshold, "red", "grey")
#     )
#   }
#   
#   dev.off()
#   
#   # -------------------------------------------------------
#   # Console summary
#   # -------------------------------------------------------
#   
#   cat("Mediator:", med, "\n")
#   cat("Family used:", med_family, "\n")
#   cat("Rows used:", nrow(df_med), "\n")
#   cat("Selected exposures/confounders:\n")
#   print(selected_variables_only)
# }