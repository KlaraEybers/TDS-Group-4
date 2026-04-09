# ==============================================================================
# Prediction Analysis: Logistic Regression and AUC Evaluation
# ==============================================================================
rm(list = ls())

# ------------------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------------------
library(here)
library(dplyr)
library(pROC)
library(broom)
library(ggplot2)


# Define file paths -------------------------------------------------------

data_dir <- here("00_data", "01_processed", "03_imputation")
results_dir <- here("02_results")
source(here("01_scripts", "00_variable_dictionary.R"))


# Import data -------------------------------------------------------------

train <- readRDS(file.path(data_dir, "df_train_imputed.rds"))
test  <- readRDS(file.path(data_dir, "df_test_imputed.rds"))


# Data checks -------------------------------------------------------------

# Ensure binary outcome
for (df_name in c("train", "test")) {
  df <- get(df_name)
  if (is.character(df$has_cvd) | is.factor(df$has_cvd)) {
    df$has_cvd <- ifelse(df$has_cvd == "Event", 1, 0)
    assign(df_name, df)
  }
}



# Define variable categories ----------------------------------------------

# Confounders
confounders <- c("age", "sex", "ethnicity")

# ------------------------------------------------------------------------------
# 3. LASSO Selected Variables (placeholder — replace with CSV reading later)
# ------------------------------------------------------------------------------
selected_exposures <- c("smoking_status_refined", "diet_score", "household_income",
                        "met_min_week_all", "alcohol_consumption")

selected_mediators <- c("crp", "hba1c", "hdl_cholesterol")


# Fit models on training set ----------------------------------------------

# Minimal: confounders only
form_minimal <- as.formula(paste("has_cvd ~", paste(confounders, collapse = " + ")))
model_minimal <- glm(form_minimal, data = train, family = binomial)

# Full: confounders + exposures + mediators
form_full <- as.formula(
  paste("has_cvd ~", paste(c(confounders, selected_exposures, selected_mediators), collapse = " + "))
)
model_full <- glm(form_full, data = train, family = binomial)


# Predict on Test Set and Compute AUC -------------------------------------

models <- list(
  "Minimal (Age + Sex + Ethnicity)" = model_minimal,
  "Full (+ Exposures + Biomarkers)" = model_full
)

roc_curves <- list()
auc_results <- data.frame()

for (name in names(models)) {
  pred <- predict(models[[name]], newdata = test, type = "response")
  roc_obj <- roc(test$has_cvd, pred, quiet = TRUE)
  ci_obj <- ci.auc(roc_obj)
  
  roc_curves[[name]] <- roc_obj
  
  auc_results <- rbind(auc_results, data.frame(
    Model = name,
    AUC = round(auc(roc_obj), 4),
    CI_Lower = round(ci_obj[1], 4),
    CI_Upper = round(ci_obj[3], 4)
  ))
}

print(auc_results)



# Train Set AUC (Overfitting check) ---------------------------------------

train_auc <- data.frame()

for (name in names(models)) {
  pred_train <- predict(models[[name]], newdata = train, type = "response")
  roc_train <- roc(train$has_cvd, pred_train, quiet = TRUE)
  
  train_auc <- rbind(train_auc, data.frame(
    Model = name,
    Train_AUC = round(auc(roc_train), 4)
  ))
}

# Combine with test AUC for comparison
overfit_check <- merge(train_auc, auc_results, by = "Model") %>%
  mutate(Difference = Train_AUC - AUC) %>%
  rename(Test_AUC = AUC)

cat("\nOverfitting check:\n")
print(overfit_check %>% select(Model, Train_AUC, Test_AUC, Difference))


# Test ROC comparison plot -------------------------------------------------

colours <- c("#95A5A6", "#12436D")

png(here(results_dir, "00_figures","03_multivariate", "roc_comparison.png"), width = 8, height = 7, 
    units = "in", res = 300)

plot(roc_curves[[1]], col = colours[1], lwd = 2,
     main = "ROC Curves: Model Comparison",
     legacy.axes = TRUE, print.auc = FALSE)

plot(roc_curves[[2]], col = colours[2], lwd = 2, add = TRUE)

legend("bottomright",
       legend = paste0(names(roc_curves), " (AUC: ", auc_results$AUC, ")"),
       col = colours, lwd = 2, cex = 0.9)

dev.off()



# Odds ratio forest plot --------------------------------------------------

# Extract odds ratios with CIs
or_table <- tidy(model_full, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  rowwise() %>%
  mutate(
    # Identify which variable this term belongs to
    base_var = names(display_names)[which(sapply(names(display_names), 
                                                 function(v) grepl(paste0("^", v), term)))][1],
    category = sub(paste0("^", base_var), "", term) %>% gsub("_", " ", .) %>% trimws(),
    var_label = get_display_name(base_var),
    display_label = ifelse(category == "", var_label, paste0(var_label, ": ", category)),
    significant = p.value < 0.05
  ) %>%
  ungroup()

# Plot
or_table %>%
  mutate(display_label = reorder(display_label, estimate)) %>%
  ggplot(aes(x = estimate, y = display_label)) +
  geom_point(aes(colour = significant), size = 3) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high, colour = significant),
                width = 0.2, orientation = "y") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("TRUE" = "#12436D", "FALSE" = "#95A5A6"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05")) +
  scale_x_log10() +
  labs(
    title = "Odds Ratios and 95% CIs from CVD Prediction Model",
    x = "Odds Ratio (log scale)",
    y = NULL,
    colour = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold", hjust = 0.8),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    plot.margin = margin(10, 60, 10, 10)
  )



# Save results ------------------------------------------------------------

write.csv(auc_results, here(results_dir, "01_tables", "01_multivariate", "auc_comparison.csv"), row.names = FALSE)
saveRDS(model_minimal, here(results_dir, "01_tables", "01_multivariate", "model_minimal.rds"))
saveRDS(model_full, here(results_dir, "01_tables", "01_multivariate", "model_full.rds"))
write.csv(or_table %>% select(display_label, estimate, conf.low, conf.high, p.value),
          here(results_dir, "01_tables", "01_multivariate", "odds_ratios.csv"), row.names = FALSE)
ggsave(here(results_dir, "00_figures", "03_multivariate", "prediction_odds_ratio_forest_plot.png"),
       width = 10, height = 8, dpi = 300)



cat("\nDone. AUC results:\n")
print(auc_results)






