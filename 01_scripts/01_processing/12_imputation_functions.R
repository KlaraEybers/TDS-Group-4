############################################################################
# DEFINE PIPELINE FUNCTIONS
############################################################################

# Define data pre-processing function --------------------------------------

# Function to process data for each different imputation set
prep_for_imputation <- function(data, vars_to_impute, predictor_vars, min_level_count = 5) {
  df_out <- data %>%
    select(any_of(c(vars_to_impute, predictor_vars))) %>%
    mutate(
      across(c(any_of(c("age_high_bp_dx", "age_diabetes_dx", "angina_dx_age",
                        "mi_dx_age", "dvt_dx_age", "pe_dx_age", "stroke_dx_age"))),
             as.numeric),
      across(any_of(log_vars), \(x) log(x + 1)),
      across(where(is.character), as.factor),
      across(where(is.ordered), \(x) factor(x, ordered = FALSE)),
      across(where(is.factor), function(x) {
        counts <- table(x)
        rare <- names(counts[counts < min_level_count])
        if (length(rare) > 0 & length(rare) < length(counts)) {
          levels(x)[levels(x) %in% rare] <- "Other"
        }
        x
      })
    ) %>%
    droplevels()
  return(df_out)
}


# Define full imputation function ------------------------------------------

run_imputation <- function(data, maxiter = 1, num_trees = 100, return_models = FALSE) {
  n_cores <- parallel::detectCores()
  
  mr <- miceRanger(
    data = data,
    m = 1,
    maxiter = maxiter,
    num.trees = num_trees,
    returnModels = return_models,
    verbose = TRUE,
    num.threads = n_cores
  )
  return(mr)
}

# Define DAG imputation function ------------------------------------------

run_dag_imputation <- function(data, label, return_models = FALSE) {
  
  cat(paste0("\n--- ", label, ": Stage 1 (exposures) ---\n"))
  exp_impute <- prep_for_imputation(
    data = data,
    vars_to_impute = exposure_vars,
    predictor_vars = confounder_vars
  )
  mr_exp <- run_imputation(exp_impute, maxiter = ifelse(TEST_MODE, 1, 5), return_models = return_models)
  exp_filled <- completeData(mr_exp)[[1]]
  
  cat(paste0("\n--- ", label, ": Stage 2 (biomarkers) ---\n"))
  bio_data <- bind_cols(exp_filled, data[, bio_vars])
  bio_impute <- prep_for_imputation(
    data = bio_data,
    vars_to_impute = bio_vars,
    predictor_vars = c(confounder_vars, exposure_vars)
  )
  mr_bio <- run_imputation(bio_impute, maxiter = ifelse(TEST_MODE, 1, 5), return_models = return_models)
  bio_filled <- completeData(mr_bio)[[1]]
  
  # Combine
  df_filled <- bind_cols(
    exp_filled,
    bio_filled %>% select(any_of(bio_vars))
  ) %>% back_transform(log_vars)
  
  # Original for diagnostics
  df_original <- bind_cols(
    data %>% select(any_of(confounder_vars)),
    data %>% select(any_of(exposure_vars)),
    data %>% select(any_of(bio_vars))
) %>%
    mutate(
      across(c(any_of(c("age_high_bp_dx", "age_diabetes_dx", "angina_dx_age",
                        "mi_dx_age", "dvt_dx_age", "pe_dx_age", "stroke_dx_age"))),
             as.numeric),
      across(where(is.character), as.factor)
    )
  
  result <- list(filled = df_filled, original = df_original)
  if (return_models) {
    result$mr_exp <- mr_exp
    result$mr_bio <- mr_bio
  }
  return(result)
}


# Define function to apply train model to test ----------------------------

apply_dag_imputation <- function(data, train_result, label) {
  
  cat(paste0("\n--- ", label, ": Stage 1 (exposures via train models) ---\n"))
  exp_impute <- prep_for_imputation(
    data = data,
    vars_to_impute = exposure_vars,
    predictor_vars = confounder_vars
  )
  mr_exp_test <- miceRanger::impute(data = exp_impute, miceObj = train_result$mr_exp)
  exp_filled <- mr_exp_test$imputedData[[1]]
  
  cat(paste0("\n--- ", label, ": Stage 2 (bio variables via train models) ---\n"))
  bio_data <- bind_cols(exp_filled, data[, bio_vars])
  bio_impute <- prep_for_imputation(
    data = bio_data,
    vars_to_impute = bio_vars,
    predictor_vars = c(confounder_vars, exposure_vars)
  )
  mr_bio_test <- miceRanger::impute(data = bio_impute, miceObj = train_result$mr_bio)
  bio_filled <- mr_bio_test$imputedData[[1]]
  
  # Combine
  df_filled <- bind_cols(
    exp_filled,
    bio_filled %>% select(any_of(bio_vars))
  ) %>% back_transform(log_vars)
  
  # Original for diagnostics
  df_original <- bind_cols(
    data %>% select(any_of(confounder_vars)),
    data %>% select(any_of(exposure_vars)),
    data %>% select(any_of(bio_vars))
) %>%
    mutate(
      across(c(any_of(c("age_high_bp_dx", "age_diabetes_dx", "angina_dx_age",
                        "mi_dx_age", "dvt_dx_age", "pe_dx_age", "stroke_dx_age"))),
             as.numeric),
      across(where(is.character), as.factor)
    )
  
  return(list(filled = df_filled, original = df_original))
}


# Define back transform function ------------------------------------------
back_transform <- function(data, log_vars) {
  data %>%
    mutate(across(any_of(log_vars), \(x) exp(x) - 1))
}


# Define diagnostics function ---------------------------------------------
save_diagnostics <- function(original, imputed, label) {
  
  # Numeric
  num_vars <- names(original)[sapply(original, is.numeric) & colSums(is.na(original)) > 0]
  num_vars <- intersect(num_vars, names(imputed))
  
  pdf(file.path(imputation_fig_dir, 
                paste0("numeric_diagnostics_", label, ".pdf")), width = 14, height = 10)
  pages <- split(num_vars, ceiling(seq_along(num_vars) / 8))
  for (page in pages) {
    plots <- lapply(page, function(var) {
      orig <- data.frame(value = original[[var]], source = "Original")
      imp <- data.frame(value = imputed[[var]], source = "Imputed")
      combined <- rbind(orig, imp) %>% filter(!is.na(value))
      ggplot(combined, aes(x = value, fill = source)) +
        geom_density(alpha = 0.3) +
        labs(title = get_display_name(var), x = "Value", y = "Density") +
        theme_minimal(base_size = 10) +
        theme(legend.position = "bottom")
    })
    print(wrap_plots(plots, ncol = 4, nrow = 2))
  }
  dev.off()
  
  # Categorical
  cat_vars <- names(original)[sapply(original, is.factor) & colSums(is.na(original)) > 0]
  cat_vars <- intersect(cat_vars, names(imputed))
  
  pdf(file.path(imputation_fig_dir, 
                paste0("categorical_diagnostics_", label, ".pdf")), width = 14, height = 10)
  pages <- split(cat_vars, ceiling(seq_along(cat_vars) / 6))
  for (page in pages) {
    plots <- lapply(page, function(var) {
      all_levels <- union(levels(original[[var]]), levels(imputed[[var]]))
      orig <- table(factor(original[[var]], levels = all_levels), useNA = "no") / sum(!is.na(original[[var]]))
      imp <- table(factor(imputed[[var]], levels = all_levels), useNA = "no") / sum(!is.na(imputed[[var]]))
      combined <- data.frame(
        level = rep(names(orig), 2),
        proportion = c(as.numeric(orig), as.numeric(imp)),
        source = rep(c("Original", "Imputed"), each = length(orig))
      )
      ggplot(combined, aes(x = level, y = proportion, fill = source)) +
        geom_col(position = "dodge", alpha = 0.5) +
        labs(title = get_display_name(var), x = "Category", y = "Proportion") +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
    })
    print(wrap_plots(plots, ncol = 3, nrow = 2))
  }
  dev.off()
}






# Define column add function ----------------------------------------------

add_back_columns <- function(original_df, imputed_df) {
  bind_cols(
    original_df %>% select(any_of(admin_vars)),
    imputed_df
  )
}
