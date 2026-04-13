# ==============================================================================
# SUB-ANALYSIS: SEX-STRATIFIED MODELS (PREDICTION & MEDIATION)
# ==============================================================================
rm(list = ls())

# 1. Import Packages -----------------------------------------------------------
library(here)
library(dplyr)
library(pROC)
library(broom)
library(ggplot2)
library(boot)

cat("=======================================================\n")
cat("Starting Sex-Stratified Sub-Analysis\n")
cat("=======================================================\n")

# 2. Define Paths -------------------------------------------------
source(here("01_scripts", "00_variable_config.R"))
source(here("01_scripts", "00_variable_dictionary.R"))

data_dir   <- data_input_dir
lasso1_dir <- lasso_outcome_dir
lasso2_dir <- lasso_mediator_dir
sub_tbl_dir <- sex_tbl_dir
sub_fig_dir <- sex_fig_dir



# 3. Load Variables from LASSO -------------------------------------------------
lasso1 <- read.csv(file.path(lasso1_dir, "selected_variable_names.csv"))
lasso2 <- read.csv(file.path(lasso2_dir, "all_selected_mediator_exposure_pairs.csv"))

selected_exposures <- lasso1$variable[lasso1$type == "exposure"]
selected_mediators <- lasso1$variable[lasso1$type == "mediator"]
mediation_exposures <- unique(lasso2$exposure)

confounders_no_sex <- setdiff(confounders, "sex")


# 4. Load Imputed Data & Stratify by Sex ---------------------------------------
cat("\nLoading imputed datasets and stratifying by Sex...\n")
train <- readRDS(file.path(data_dir, train_file))
test  <- readRDS(file.path(data_dir, test_file))

for (df_name in c("train", "test")) {
  df <- get(df_name)
  if (is.character(df$has_cvd) | is.factor(df$has_cvd)) {
    df$has_cvd <- ifelse(df$has_cvd == "Event", 1, 0)
    assign(df_name, df)
  }
}

# Ensure all variables exist in the data
selected_exposures  <- intersect(selected_exposures, names(train))
selected_mediators  <- intersect(selected_mediators, names(train))
mediation_exposures <- intersect(mediation_exposures, names(train))
all_exposures       <- intersect(union(selected_exposures, mediation_exposures), names(train))

# Also filter lasso2 pairs to only include variables present in data
lasso2 <- lasso2 %>%
  filter(exposure %in% names(train), mediator %in% names(train))

sapply(train[, sapply(train, is.factor)], is.ordered)

# Stratify by sex
male_levels <- c("Male", "male", "M", "1", 1)
train_m <- train %>% filter(sex %in% male_levels)
train_f <- train %>% filter(!sex %in% male_levels)
test_m  <- test %>% filter(sex %in% male_levels)
test_f  <- test %>% filter(!sex %in% male_levels)

# ==============================================================================
# PART A: PREDICTION MODELS & OR PLOTS
# ==============================================================================
cat("\n--- Running Prediction Models (Male & Female) ---\n")

form_full <- as.formula(paste("has_cvd ~", paste(c(confounders_no_sex, selected_exposures, selected_mediators), collapse = " + ")))

mod_full_m <- glm(form_full, data = train_m, family = binomial)
mod_full_f <- glm(form_full, data = train_f, family = binomial)

# 5. Extract ORs and Create Faceted Forest Plot --------------------------------
extract_or <- function(model, cohort_name) {
  tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    rowwise() %>%
    mutate(
      base_var = names(display_names)[which(sapply(names(display_names), function(v) grepl(paste0("^", v), term)))][1],
      category = sub(paste0("^", base_var), "", term) %>% gsub("_", " ", .) %>% trimws(),
      var_label = get_display_name(base_var),
      display_label = ifelse(category == "", var_label, paste0(var_label, ": ", category)),
      significant = p.value < 0.05,
      Sex = cohort_name
    ) %>% ungroup()
}

or_table_m <- extract_or(mod_full_m, "Male")
or_table_f <- extract_or(mod_full_f, "Female")
or_all <- bind_rows(or_table_m, or_table_f) %>% 
  filter(!grepl(paste(paste0("^", confounders_no_sex), collapse = "|"), term)) %>%
  mutate(Sex = factor(Sex, levels = c("Male", "Female")))


write.csv(or_all, file.path(sub_tbl_dir, "sub_prediction_odds_ratios.csv"), row.names = FALSE)

order_df <- or_all %>% filter(Sex == "Male") %>% arrange(estimate)
or_all$display_label <- factor(or_all$display_label, levels = unique(order_df$display_label))

p_or <- ggplot(or_all, aes(x = estimate, y = display_label)) +
  geom_point(aes(colour = significant), size = 3) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high, colour = significant), width = 0.2, orientation = "y") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ Sex) +
 
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "#95A5A6"), labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05")) +
  scale_x_log10() +
  labs(title = "Sub-Analysis: Odds Ratios and 95% CIs by Sex", x = "Odds Ratio", y = NULL, colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "bold", hjust = 0.5), strip.text = element_text(size = 13, face = "bold"), strip.background = element_rect(fill = "gray90", color = NA), axis.text.y = element_text(size = 11), panel.grid.minor = element_blank())

ggsave(file.path(sub_fig_dir, "sub_prediction_odds_ratio_forest_plot.png"), plot = p_or, width = 11, height = 7, dpi = 300, bg="white")


# 6. ROC Curves -------------------------------------------------
cat("\n--- Calculating ROC curves and AUC Confidence Intervals ---\n")
pred_m <- predict(mod_full_m, newdata = test_m, type = "response")
pred_f <- predict(mod_full_f, newdata = test_f, type = "response")

roc_m <- roc(test_m$has_cvd, pred_m, quiet = TRUE)
roc_f <- roc(test_f$has_cvd, pred_f, quiet = TRUE)

roc_data_m <- data.frame(FPR = 1 - roc_m$specificities, TPR = roc_m$sensitivities, Cohort = "Male")
roc_data_f <- data.frame(FPR = 1 - roc_f$specificities, TPR = roc_f$sensitivities, Cohort = "Female")
roc_all <- bind_rows(roc_data_m, roc_data_f) %>% mutate(Cohort = factor(Cohort, levels = c("Male", "Female")))

auc_ci_m <- ci.auc(roc_m)
auc_ci_f <- ci.auc(roc_f)


auc_summary <- data.frame(
  Sex = c("Male", "Female"),
  AUC = c(auc_ci_m[2], auc_ci_f[2]),
  CI_Lower_95 = c(auc_ci_m[1], auc_ci_f[1]),
  CI_Upper_95 = c(auc_ci_m[3], auc_ci_f[3])
)
write.csv(auc_summary, file.path(sub_tbl_dir, "sub_roc_auc_summary.csv"), row.names = FALSE)
write.csv(roc_all, file.path(sub_tbl_dir, "sub_roc_curve_coordinates.csv"), row.names = FALSE)

auc_labels <- data.frame(
  Cohort = factor(c("Male", "Female"), levels = c("Male", "Female")),
  label = c(sprintf("Test AUC: %.3f\n95%% CI: %.3f - %.3f", auc_ci_m[2], auc_ci_m[1], auc_ci_m[3]),
            sprintf("Test AUC: %.3f\n95%% CI: %.3f - %.3f", auc_ci_f[2], auc_ci_f[1], auc_ci_f[3]))
)

p_roc <- ggplot() +
  geom_line(data = roc_all, aes(x = FPR, y = TPR), colour = "#12436D", linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "#95A5A6") +
  facet_wrap(~ Cohort) +
  geom_text(data = auc_labels, aes(x = 0.65, y = 0.20, label = label), color = "black", inherit.aes = FALSE, size = 4) +
  labs(title = "Sub-Analysis: ROC Curves by Sex", x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0), strip.text = element_text(size = 13, face = "bold"), strip.background = element_rect(fill = "gray90", color = NA))

ggsave(file.path(sub_fig_dir, "sub_roc_curves_by_sex.png"), plot = p_roc, width = 10, height = 6, dpi = 300, bg="white")


# ==============================================================================
# PART B: MEDIATION ANALYSIS (BOOTSTRAPPED & EXPONENTIATED TO OR)
# ==============================================================================
cat("\n--- Running Mediation Models (Robust Bootstrapping) ---\n")


n_cores <- min(parallel::detectCores(), 16) 

run_stratified_boot <- function(data_subset, cohort_label, R_boot = 500) {
  form_dir <- as.formula(paste("has_cvd ~", paste(c(confounders_no_sex, all_exposures, selected_mediators), collapse = " + ")))
  
  
  mod_dir_main <- glm(form_dir, data = data_subset, family = binomial)
  b_cfs_main <- coef(mod_dir_main)[selected_mediators]
  
  expected_paths <- list()
  for (med in selected_mediators) {
    m_exposures <- lasso2$exposure[lasso2$mediator == med]
    mod_a_main <- lm(as.formula(paste(med, "~", paste(c(confounders_no_sex, m_exposures), collapse = " + "))), data = data_subset)
    a_cfs_main <- coef(mod_a_main)
    
    for (exp in m_exposures) {
      match_lvl <- setdiff(grep(paste0("^", exp), names(a_cfs_main), value = TRUE), "(Intercept)")
      for (lvl in match_lvl) {
        expected_paths <- append(expected_paths, list(list(exposure = exp, level = lvl, mediator = med)))
      }
    }
  }
  
  compute_indirect_strat <- function(d, indices) {
    dat <- d[indices, ]
    mod_b <- glm(form_dir, data = dat, family = binomial)
    b_cfs <- coef(mod_b)
    
    ab_res <- numeric(length(expected_paths))
    for (i in seq_along(expected_paths)) {
      path <- expected_paths[[i]]
      med <- path$mediator
      lvl <- path$level
      
      b_val <- b_cfs[med]
      m_exposures <- lasso2$exposure[lasso2$mediator == med]
      mod_a <- lm(as.formula(paste(med, "~", paste(c(confounders_no_sex, m_exposures), collapse = " + "))), data = dat)
      a_cfs <- coef(mod_a)
      
      
      a_val <- ifelse(lvl %in% names(a_cfs), a_cfs[lvl], NA)
      ab_res[i] <- a_val * b_val
    }
    return(ab_res)
  }
  
  set.seed(123)
  boot_res <- boot(data = data_subset, statistic = compute_indirect_strat, R = R_boot, parallel = "multicore", ncpus = n_cores)
  
  strat_cis <- data.frame()
  for (i in seq_along(expected_paths)) {
    path <- expected_paths[[i]]
    ci <- tryCatch(boot.ci(boot_res, index = i, type = "perc"), error = function(e) NULL)
    
    strat_cis <- rbind(strat_cis, 
                       data.frame(Sex = cohort_label, 
                                  exposure = path$exposure, 
                                  level = path$level, 
                                  mediator = path$mediator, 
                                  ab = boot_res$t0[i], 
                                  ci_lower = if(is.null(ci)) NA else ci$percent[4], 
                                  ci_upper = if(is.null(ci)) NA else ci$percent[5]))
  }
  return(strat_cis %>% mutate(significant = !is.na(ci_lower) & !is.na(ci_upper) & !(ci_lower <= 0 & ci_upper >= 0)))
}

med <- "triglycerides"
m_exposures <- lasso2$exposure[lasso2$mediator == med]
mod_test <- lm(as.formula(paste(med, "~", paste(c(confounders_no_sex, m_exposures), collapse = " + "))), data = train_m)
names(coef(mod_test))[grep("household", names(coef(mod_test)))]



med_male <- run_stratified_boot(train_m, "Male", R_boot = 500)
med_female <- run_stratified_boot(train_f, "Female", R_boot = 500)


med_all <- bind_rows(med_male, med_female) %>% 
  mutate(
    Sex = factor(Sex, levels = c("Male", "Female")),
    ab = exp(ab),
    ci_lower = exp(ci_lower),
    ci_upper = exp(ci_upper)
  )


write.csv(med_all, file.path(sub_tbl_dir, "sub_mediation_odds_ratios_by_sex.csv"), row.names = FALSE)

# 7. Faceted Forest Plot for Indirect Effects (OR Scale) -----------------------
plot_med_data <- med_all %>%
  rowwise() %>%
  mutate(
    category = sub(paste0("^", exposure), "", level) %>% gsub("_", " ", .) %>% trimws(),
    exp_short = get_display_name(exposure),
    med_short = get_display_name(mediator),
    pathway_label = ifelse(category == "" | category == gsub("_", " ", exposure), paste0(exp_short, "  →  ", med_short), paste0(exp_short, ": ", category, "  →  ", med_short))
  ) %>% ungroup()

if(nrow(plot_med_data) > 0) {
  order_med_df <- plot_med_data %>% filter(Sex == "Male") %>% arrange(ab)
  if(nrow(order_med_df) > 0) {
    plot_med_data$pathway_label <- factor(plot_med_data$pathway_label, levels = unique(order_med_df$pathway_label))
  }
  
  p_med_forest <- ggplot(plot_med_data, aes(y = pathway_label, x = ab)) +
    geom_point(aes(colour = significant), size = 3) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper, colour = significant), width = 0.2, orientation = "y") +
    
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    facet_wrap(~ Sex) +
    
    scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "#95A5A6"), labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
    
    scale_x_log10() +
    labs(title = "Sub-Analysis: Indirect Effects by Sex", x = "Indirect Effect (Exponentiated)", y = NULL, colour = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "bold", hjust = 0.5), strip.text = element_text(size = 13, face = "bold"), strip.background = element_rect(fill = "gray90", color = NA), axis.text.y = element_text(size = 10))
  
  ggsave(file.path(sub_fig_dir, "sub_mediation_forest_plot_by_sex.png"), plot = p_med_forest, width = 12, height = 7, dpi = 300, bg="white")
}

cat("\n=======================================================\n")
cat("SUCCESS! All Sub-Analysis results generated and saved securely.\n")
cat("=======================================================\n")