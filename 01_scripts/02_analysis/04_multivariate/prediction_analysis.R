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

source(here("01_scripts", "00_variable_config.R"))
source(here("01_scripts", "00_variable_dictionary.R"))

data_dir   <- data_input_dir
fig_dir    <- figures_multi_dir
tab_dir    <- tables_multi_dir
lasso1_dir <- lasso_outcome_dir


# Import data -------------------------------------------------------------

train <- readRDS(file.path(data_dir, train_file))
test  <- readRDS(file.path(data_dir, test_file))
lasso1 <- read.csv(file.path(lasso1_dir, "selected_variable_names.csv"))


# Data checks -------------------------------------------------------------

# Ensure binary outcome
for (df_name in c("train", "test")) {
  df <- get(df_name)
  if (is.character(df$has_cvd) | is.factor(df$has_cvd)) {
    df$has_cvd <- ifelse(df$has_cvd == "Event", 1, 0)
    assign(df_name, df)
  }
}


# Extract LASSO variables -------------------------------------------------

selected_exposures <- lasso1$variable[lasso1$type == "exposure"]
selected_mediators <- lasso1$variable[lasso1$type == "mediator"]

cat("Selected exposures:", length(selected_exposures), "\n")
print(selected_exposures)
cat("Selected mediators:", length(selected_mediators), "\n")
print(selected_mediators)

# Ensure all variables exist in the data
selected_exposures <- intersect(selected_exposures, names(train))
selected_mediators <- intersect(selected_mediators, names(train))

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
  "Minimal (Confounders)" = model_minimal,
  "Full (Confounders + Exposures + Biomarkers)" = model_full
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


# Compute CI Bands for each model -----------------------------------------

ci_bands <- list()
for (name in names(models)) {
  ci_bands[[name]] <- ci.se(roc_curves[[name]], 
                            specificities = seq(0, 1, by = 0.01),
                            quiet = TRUE)
}



# Export ROC plot data ----------------------------------------------------

roc_plot_data <- do.call(rbind, lapply(names(roc_curves), function(name) {
  roc_obj <- roc_curves[[name]]
  data.frame(
    model = name,
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )
}))

write.csv(roc_plot_data, file.path(tab_dir, "roc_plot_data.csv"), row.names = FALSE)


# Test ROC comparison plot with CI bands ----------------------------------

model_colours <- c("#95A5A6", "#12436D")

png(file.path(fig_dir, "roc_comparison.png"),
    width = 8, height = 7, units = "in", res = 300)

plot(roc_curves[[1]], col = model_colours[1], lwd = 2,
     main = "ROC Curves: Model Comparison",
     legacy.axes = TRUE, print.auc = FALSE)
plot(ci_bands[[1]], type = "shape", col = adjustcolor(model_colours[1], alpha.f = 0.15), 
     border = NA)

plot(roc_curves[[2]], col = model_colours[2], lwd = 2, add = TRUE)
plot(ci_bands[[2]], type = "shape", col = adjustcolor(model_colours[2], alpha.f = 0.15), 
     border = NA)

# Replot curves on top so they're not obscured by shading
plot(roc_curves[[1]], col = model_colours[1], lwd = 2, add = TRUE)
plot(roc_curves[[2]], col = model_colours[2], lwd = 2, add = TRUE)

legend("bottomright",
       legend = paste0(names(roc_curves), " (AUC: ", auc_results$AUC, 
                       " [", auc_results$CI_Lower, "-", auc_results$CI_Upper, "])"),
       col = model_colours, lwd = 2, cex = 0.9)

dev.off()



# Odds ratio forest plot --------------------------------------------------

# Extract odds ratios with CIs
or_table <- tidy(model_full, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  filter(!grepl(paste(paste0("^", confounders), collapse = "|"), term)) %>%
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

# Export or_table data
write.csv(or_table %>% select(term, base_var, display_label, estimate, conf.low, conf.high, p.value, significant),
          file.path(tab_dir, "or_plot_data.csv"), row.names = FALSE)

# Plot
or_table %>%
  mutate(display_label = reorder(display_label, estimate)) %>%
  ggplot(aes(x = estimate, y = display_label)) +
  geom_point(aes(colour = significant), size = 3) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high, colour = significant),
                width = 0.2, orientation = "y") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "#95A5A6"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05")) +
  scale_x_log10() +
  labs(
    title = "Odds Ratios and 95% CIs from CVD Prediction Model",
    x = "Odds Ratio",
    y = NULL,
    colour = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.8),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    plot.margin = margin(10, 60, 10, 10),
  )

ggsave(file.path(fig_dir, "prediction_odds_ratio_forest_plot.png"),
       width = 10, height = 8, dpi = 300)


# Confusion Matrix Visualisation -------------------------------------------

pred_prob <- predict(model_full, newdata = test, type = "response")
threshold <- mean(train$has_cvd)
cat("Prevalence-based threshold:", round(threshold, 4), "\n")
pred_class <- ifelse(pred_prob > threshold, 1, 0)

cm <- table(Predicted = pred_class, Actual = test$has_cvd)

# Convert to long format for ggplot
cm_df <- as.data.frame(cm)
cm_df$Predicted <- factor(cm_df$Predicted, levels = c("1", "0"))
cm_df$Actual <- factor(cm_df$Actual, levels = c("0", "1"))

# Export
write.csv(cm_df, file.path(tab_dir, "confusion_matrix_plot_data.csv"), row.names = FALSE)

ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Freq)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = Freq), size = 8, fontface = "bold", colour = "white") +
  scale_fill_gradient(low = "#95A5A6", high = "black") +
  labs(
    title = paste0("Confusion Matrix: Full Model (Threshold = ", round(threshold, 3), ")"),
    x = "Actual",
    y = "Predicted"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

ggsave(file.path(fig_dir, "confusion_matrix.png"),
       width = 6, height = 5, dpi = 300)


# Calibration Plot ---------------------------------------------------------

# Create deciles of predicted probability
cal_df <- data.frame(
  observed = test$has_cvd,
  predicted = pred_prob
) %>%
  mutate(decile = ntile(predicted, 10)) %>%
  group_by(decile) %>%
  summarise(
    mean_predicted = mean(predicted),
    mean_observed = mean(observed),
    n = n(),
    n_events = sum(observed),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    ci_lower = binom.test(n_events, n)$conf.int[1],
    ci_upper = binom.test(n_events, n)$conf.int[2]
  ) %>%
  ungroup()

# Export
write.csv(cal_df, file.path(tab_dir, "calibration_plot_data.csv"), row.names = FALSE)

ggplot(cal_df, aes(x = mean_predicted, y = mean_observed)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, colour = "black", linewidth = 0.6) +
  geom_point(size = 3.5, colour = "black") +
  geom_line(colour = "black", linewidth = 0.8) +
  labs(
    title = "Calibration Plot of Full Model with 95% CIs",
    x = "Mean Predicted Probability",
    y = "Observed Proportion",
    caption = "Dashed line = perfect calibration"
  ) +
  coord_cartesian(xlim = range(cal_df$mean_predicted), ylim = range(c(cal_df$ci_lower, cal_df$ci_upper))) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    plot.caption = element_text(size = 12, colour = "grey40")
  )

# Export
ggsave(file.path(fig_dir, "calibration_plot.png"),
       width = 7, height = 7, dpi = 300)

# Save results ------------------------------------------------------------

write.csv(auc_results, file.path(tab_dir, "auc_comparison.csv"), row.names = FALSE)
saveRDS(model_minimal, file.path(tab_dir, "model_minimal.rds"))
write.csv(tidy(model_minimal, conf.int = TRUE, exponentiate = TRUE),
          file.path(tab_dir, "model_minimal_summary.csv"), row.names = FALSE)
saveRDS(model_full, file.path(tab_dir, "model_full.rds"))
write.csv(tidy(model_full, conf.int = TRUE, exponentiate = TRUE),
          file.path(tab_dir, "model_full_summary.csv"), row.names = FALSE)
write.csv(overfit_check,
          file.path(tab_dir, "overfit_check.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(cm),
          file.path(tab_dir, "confusion_matrix.csv"))


cat("\nDone. AUC results:\n")
print(auc_results)






