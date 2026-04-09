rm(list = ls())

# Packages -------------------------------------------------

library(here)
library(dplyr)
library(purrr)
library(ggplot2)
library(sharp)

# Load data ------------------------------------------------

df <- readRDS(here("00_data", "01_processed", "03_imputation", "df_selection_imputed.rds"))

# Basic checks --------------------------------------------

stopifnot("has_cvd" %in% names(df))
stopifnot(all(df$has_cvd %in% c(0, 1)))

na_total <- sum(is.na(df))
cat("Total NA values in loaded dataset:", na_total, "\n")

# Create output directory ---------------------------------

out_dir <- here("00_data", "02_analysis", "stability_selection_runs", "stability_selection_sharp_wo_sbp_comorb_wo_sleepsmoke_zoomed_100")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 1) Define variable groups
# =========================================================

admin_vars <- c(
  "eid", "recruitment_date", "death_date", "first_cvd_date", "has_cvd"
)

confounder_vars <- c(
  "age", 
  "sex", "ethnicity"
)

exposure_vars <- c(
  "household_income", "highest_qualification", "employment", "housing", "household_size",
  "diet_score", "alcohol_consumption", "salt_added",
  "sleep_hrs", "insomnia",
  # "sleep_cat", "passive_cat",
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

# Keep only variables that actually exist ------------------

admin_vars      <- intersect(admin_vars, names(df))
confounder_vars <- intersect(confounder_vars, names(df))
exposure_vars   <- intersect(exposure_vars, names(df))
bio_vars        <- intersect(bio_vars, names(df))

candidate_predictors <- unique(c(confounder_vars, exposure_vars, bio_vars))

# Keep only required columns -------------------------------

df <- df %>%
  select(all_of(c(intersect("eid", names(.)), "has_cvd", candidate_predictors)))

# Save IDs from this already-defined selection dataset -----

if ("eid" %in% names(df)) {
  write.csv(
    df %>% select(eid, has_cvd),
    file.path(out_dir, "selection_dataset_ids.csv"),
    row.names = FALSE
  )
}

# Remove eid before modelling ------------------------------

if ("eid" %in% names(df)) {
  df <- df %>% select(-eid)
}

# Also remove sbp to try see if it improves variable selection
if ("sbp" %in% names(df)) {
  df <- df %>% select(-sbp)
}

# # Stratifying by sex to see if it picks different variables
# df <- df %>%
#   filter(sex == "Male")
# 
# summary(df$sex)

# =========================================================
# 2) Prepare data types
# =========================================================

# Convert character -> factor
df <- df %>%
  mutate(across(where(is.character), as.factor))

# Convert Date -> numeric
date_cols <- names(df)[vapply(df, inherits, logical(1), what = "Date")]
if (length(date_cols) > 0) {
  df <- df %>%
    mutate(across(all_of(date_cols), ~ as.numeric(.)))
}

# Ensure outcome is numeric 0/1
df$has_cvd <- as.numeric(df$has_cvd)

# Check for remaining missing values
na_after_prep <- sum(is.na(df))
cat("Total NA values after preprocessing:", na_after_prep, "\n")

if (na_after_prep > 0) {
  stop("Dataset still contains missing values. sharp requires no NAs in xdata/ydata.")
}

# =========================================================
# 3) One-hot encoding for sharp
# =========================================================

x_df <- df %>% select(-has_cvd)
y <- df$has_cvd

# model.matrix one-hot encodes factors automatically
x_mat <- model.matrix(~ ., data = x_df)[, -1, drop = FALSE]

# Remove zero-variance columns
keep_cols <- apply(x_mat, 2, function(z) stats::sd(z) > 0)
x_mat <- x_mat[, keep_cols, drop = FALSE]

# =========================================================
# 4) Penalty factors: confounders unpenalised
# =========================================================

penalty_factor <- rep(1, ncol(x_mat))
names(penalty_factor) <- colnames(x_mat)

# Any dummy column derived from confounders gets penalty 0
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

out <- sharp::VariableSelection(
  xdata = x_mat,
  ydata = y,
  family = "binomial",
  penalty.factor = penalty_factor,
  verbose = FALSE,
  n_cat = 3,
  K = 100,
  pi_list = seq(0.65, 0.79, by = 0.01),
  Lambda = seq(8e-04, 1.2e-03, length.out = 60)
)

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

# =========================================================
# 6) Save object / generic outputs for debugging
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

write_component_csv <- function(obj, comp_name, file_name) {
  if (comp_name %in% names(obj)) {
    comp <- obj[[comp_name]]
    if (is.data.frame(comp) || is.matrix(comp)) {
      write.csv(as.data.frame(comp), file.path(out_dir, file_name), row.names = FALSE)
    }
  }
}

write_component_csv(out, "summary_variable", "summary_variable.csv")
write_component_csv(out, "summary_model", "summary_model.csv")
write_component_csv(out, "summary_selection", "summary_selection.csv")
write_component_csv(out, "summary_cat", "summary_cat.csv")
write_component_csv(out, "Table", "sharp_table.csv")
write_component_csv(out, "PFER", "sharp_pfer.csv")
write_component_csv(out, "selprop", "sharp_selprop_matrix.csv")
write_component_csv(out, "Lambda", "sharp_lambda.csv")
write_component_csv(out, "Q", "sharp_q.csv")
write_component_csv(out, "P", "sharp_p.csv")
write_component_csv(out, "FDP", "sharp_fdp.csv")

png(
  filename = file.path(out_dir, "sharp_calibration_plot.png"),
  width = 2300,
  height = 1300,
  res = 150
)

try(sharp::CalibrationPlot(out), silent = TRUE)

par(mar = c(6, 6, 4, 6))

dev.off()



# =========================================================
# 7) Practical-style calibrated outputs
# =========================================================

selprop <- SelectionProportions(out)
hat_params <- Argmax(out)

if (is.null(names(selprop))) {
  names(selprop) <- colnames(x_mat)
}

sel_threshold <- hat_params[2]

cat("Calibrated selection threshold from sharp:", sel_threshold, "\n")

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

# =========================================================
# 8) Feature-level and source-variable-level summaries
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

# Optional: non-confounder selected variables only
selected_variables_nonconf <- selected_variables_only %>%
  filter(group != "confounder")

write.csv(
  selected_variables_nonconf,
  file.path(out_dir, "selected_variables_nonconfounders_only.csv"),
  row.names = FALSE
)
selected_variable_names <- selected_variables_final %>%
  filter(selected) %>%
  pull(source_variable) %>%
  unique()

write.csv(
  data.frame(variable = selected_variable_names),
  file.path(out_dir, "selected_variable_names.csv"),
  row.names = FALSE
)

writeLines(
  selected_variable_names,
  file.path(out_dir, "selected_variable_names.txt")
)

# =========================================================
# 9) Plots
# =========================================================

# Practical-style calibrated plot
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

# Keep the default sharp plot too
png(
  filename = file.path(out_dir, "sharp_plot_default.png"),
  width = 1400,
  height = 1000,
  res = 150
)

try(plot(out), silent = TRUE)

dev.off()

# =========================================================
# 10) Console summary
# =========================================================

cat("\n====================================================\n")
cat("sharp variable selection completed\n")
cat("====================================================\n")
cat("Selection dataset rows:", nrow(df), "\n")
cat("Number of encoded features:", ncol(x_mat), "\n")
cat("n_cat:", 3, "\n")
cat("Calibrated selection threshold:", sel_threshold, "\n")
cat("Output directory:\n", out_dir, "\n")
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

# Variables that made it into the model
cat("\nVariables WITH features:\n")
print(feature_summary %>% filter(n_features > 0))

# Variables that were DROPPED
cat("\nVariables with ZERO features (DROPPED):\n")
print(feature_summary %>% filter(n_features == 0))

cat("ncol(x_mat):", ncol(x_mat), "\n")
cat("length(selprop):", length(selprop), "\n")
cat("number of names(selprop):", length(names(selprop)), "\n")

summary(df$sex)
summary(df$ethnicity)

# # =========================================================
# # 6) Try to export standard outputs if present
# # =========================================================
# 
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
# 
# # =========================================================
# # 7) Build a generic encoded-feature summary table
# # =========================================================
# 
# feature_summary <- NULL
# 
# possible_feature_components <- c(
#   "summary_variable",
#   "summary_selection",
#   "Table"
# )
# 
# for (nm in possible_feature_components) {
#   if (nm %in% names(out)) {
#     tmp <- as.data.frame(out[[nm]])
#     if (nrow(tmp) > 0) {
#       feature_summary <- tmp
#       break
#     }
#   }
# }
# 
# if (!is.null(feature_summary)) {
#   
#   if (!("feature" %in% names(feature_summary))) {
#     if ("variable" %in% names(feature_summary)) {
#       feature_summary$feature <- feature_summary$variable
#     } else if ("covariate" %in% names(feature_summary)) {
#       feature_summary$feature <- feature_summary$covariate
#     } else if ("Var1" %in% names(feature_summary)) {
#       feature_summary$feature <- feature_summary$Var1
#     }
#   }
#   
#   get_source_variable <- function(feature_name, original_vars) {
#     hits <- original_vars[startsWith(feature_name, original_vars)]
#     if (length(hits) == 0) return(NA_character_)
#     hits[which.max(nchar(hits))]
#   }
#   
#   if ("feature" %in% names(feature_summary)) {
#     feature_summary$source_variable <- sapply(
#       feature_summary$feature,
#       get_source_variable,
#       original_vars = c(confounder_vars, exposure_vars, bio_vars)
#     )
#     
#     feature_summary$group <- dplyr::case_when(
#       feature_summary$source_variable %in% confounder_vars ~ "confounder",
#       feature_summary$source_variable %in% bio_vars        ~ "biological",
#       feature_summary$source_variable %in% exposure_vars   ~ "exposure",
#       TRUE ~ "unknown"
#     )
#     
#     write.csv(
#       feature_summary,
#       file.path(out_dir, "feature_summary_mapped.csv"),
#       row.names = FALSE
#     )
#   }
# }
# 
# # =========================================================
# # 8) Save selected source variables if recoverable
# # =========================================================
# 
# selected_variables_final <- NULL
# 
# if (!is.null(feature_summary) && "source_variable" %in% names(feature_summary)) {
#   
#   category_col <- names(feature_summary)[tolower(names(feature_summary)) %in% c(
#     "category", "cat", "selected", "status", "selection_category"
#   )]
#   
#   if (length(category_col) == 1) {
#     selected_variables_final <- feature_summary %>%
#       filter(.data[[category_col]] %in% c(
#         "stably included", "included", "selected", 1
#       )) %>%
#       filter(!is.na(source_variable)) %>%
#       distinct(source_variable, group) %>%
#       arrange(group, source_variable)
#   }
# }
# 
# if (!is.null(selected_variables_final) && nrow(selected_variables_final) > 0) {
#   write.csv(
#     selected_variables_final,
#     file.path(out_dir, "selected_variables_final.csv"),
#     row.names = FALSE
#   )
# }
# 
# # =========================================================
# # 9) Plot using sharp's plot method if available
# # =========================================================
# 
# png(
#   filename = file.path(out_dir, "sharp_plot.png"),
#   width = 1400,
#   height = 1000,
#   res = 150
# )
# 
# try(plot(out), silent = TRUE)
# 
# dev.off()
# 
# # =========================================================
# # 10) Console summary
# # =========================================================
# 
# cat("\n====================================================\n")
# cat("sharp variable selection completed\n")
# cat("====================================================\n")
# cat("Selection dataset rows:", nrow(df), "\n")
# cat("Number of encoded features:", ncol(x_mat), "\n")
# cat("n_cat:", 3, "\n")
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
# if (!is.null(selected_variables_final) && nrow(selected_variables_final) > 0) {
#   cat("\nRecovered selected source variables:\n")
#   print(selected_variables_final)
# } else {
#   cat("\nNo selected source-variable table could be recovered automatically.\n")
#   cat("Check the saved files:\n")
#   cat("- sharp_object_structure.txt\n")
#   cat("- summary_variable.csv / summary_selection.csv / sharp_table.csv\n")
#   cat("- sharp_plot.png\n")
# }