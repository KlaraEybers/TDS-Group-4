# ==============================================================================
# Evaluate Model Performance on Train and Test Sets 
# ==============================================================================

cat("Setting up paths and loading data...\n")
data_dir <- "/rds/general/user/sh425/projects/hda_25-26/live/TDS/TDS_Group4/00_data/01_processed/03_imputation"
results_dir <- "/rds/general/user/sh425/projects/hda_25-26/live/TDS/TDS_Group4/02_results/01_tables"

# ------------------------------------------------------------------------------
# 0. Custom Function to Calculate AUC Natively (Bypassing pROC)
# ------------------------------------------------------------------------------
calculate_auc_native <- function(actual, predicted) {
  ord <- order(predicted, decreasing = TRUE)
  act_ord <- actual[ord]
  
  tp <- cumsum(act_ord == 1)
  fp <- cumsum(act_ord == 0)
  
  P <- sum(act_ord == 1)
  N <- sum(act_ord == 0)
  
  tpr <- c(0, tp / P)
  fpr <- c(0, fp / N)
  
  # Trapezoidal rule for AUC
  auc_val <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
  return(auc_val)
}

# ------------------------------------------------------------------------------
# 1. Load Train and Test datasets
# ------------------------------------------------------------------------------
train_imputed <- readRDS(file.path(data_dir, "df_train_imputed.rds"))
test_imputed <- readRDS(file.path(data_dir, "df_test_imputed.rds"))

# Ensure outcome variable is binary (0/1)
if(is.character(train_imputed$has_cvd) | is.factor(train_imputed$has_cvd)) {
  train_imputed$has_cvd <- ifelse(train_imputed$has_cvd == "Event", 1, 0)
}
if(is.character(test_imputed$has_cvd) | is.factor(test_imputed$has_cvd)) {
  test_imputed$has_cvd <- ifelse(test_imputed$has_cvd == "Event", 1, 0)
}

# ------------------------------------------------------------------------------
# 2. Load the previously trained (locked) models
# ------------------------------------------------------------------------------
cat("Loading locked models...\n")
model_total <- readRDS(file.path(results_dir, "model_total_final1.rds"))
model_direct <- readRDS(file.path(results_dir, "model_direct_final1.rds"))

# ------------------------------------------------------------------------------
# 3. Prediction Phase 
# ------------------------------------------------------------------------------
cat("Predicting probabilities...\n")

prob_train_total <- predict(model_total, newdata = train_imputed, type = "response")
prob_train_direct <- predict(model_direct, newdata = train_imputed, type = "response")

prob_test_total <- predict(model_total, newdata = test_imputed, type = "response")
prob_test_direct <- predict(model_direct, newdata = test_imputed, type = "response")

# ------------------------------------------------------------------------------
# 4. Evaluation Phase - Print results
# ------------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("          MODEL PERFORMANCE (AUC SCORES)               \n")
cat("=======================================================\n")

# Native AUC calculation
auc_train_tot <- calculate_auc_native(train_imputed$has_cvd, prob_train_total)
auc_train_dir <- calculate_auc_native(train_imputed$has_cvd, prob_train_direct)

auc_test_tot <- calculate_auc_native(test_imputed$has_cvd, prob_test_total)
auc_test_dir <- calculate_auc_native(test_imputed$has_cvd, prob_test_direct)

cat(sprintf("1. TOTAL MODEL (No Mediator):\n"))
cat(sprintf("   - Train AUC : %.4f\n", auc_train_tot))
cat(sprintf("   - Test AUC  : %.4f\n\n", auc_test_tot))

cat(sprintf("2. DIRECT MODEL (With Mediator):\n"))
cat(sprintf("   - Train AUC : %.4f\n", auc_train_dir))
cat(sprintf("   - Test AUC  : %.4f\n", auc_test_dir))
cat("=======================================================\n")

