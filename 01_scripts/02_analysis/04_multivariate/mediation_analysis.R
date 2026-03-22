# ==============================================================================
# Mediation Analysis: Total, Direct, and Indirect Effects 
# ==============================================================================
rm(list = ls())


# Import packages ---------------------------------------------------------

library(here)
library(dplyr)
library(boot)
library(parallel)


# Paths -------------------------------------------------------------------

data_dir <- here("00_data", "01_processed", "03_imputation")
results_dir <- here("02_results", "01_tables")
lasso_dir <- here("00_data", "02_analysis", "stability_selection_runs", "stability_selection_sharp_wo_sbp_FINAL")


# Load data ---------------------------------------------------------------

train <- readRDS(file.path(data_dir, "df_train_imputed.rds"))
source(here("01_scripts", "00_variable_dictionary.R"))




# Ensure binary outcome (0/1)
if(is.character(train$has_cvd) | is.factor(train$has_cvd)){
  train$has_cvd <- ifelse(train$has_cvd == "Event", 1, 0)
}


# Load variables from LASSO selection -------------------------------------

# lasso_vars <- readLines(file.path(lasso_dir, "selected_variable_names.txt"))
# lasso_vars <- trimws(lasso_vars)
# lasso_vars <- lasso_vars[lasso_vars != ""]

# Confounders
confounders <- c("age", "sex", "ethnicity")

# ============================================================
# PLACEHOLDER: Replace with LASSO output reading when available
# ============================================================

# Selected exposures from LASSO 1
selected_exposures <- c("smoking_status_refined", "diet_score", "household_income", 
                        "met_min_week_all", "alcohol_consumption")

# Selected mediators from LASSO 1
selected_mediators <- c("crp", "hba1c", "hdl_cholesterol")

# Selected a-paths from LASSO 2 (which exposures influence which mediators)
lasso2 <- data.frame(
  mediator = c("crp", "crp", "hba1c", "hba1c", "hdl_cholesterol", "hdl_cholesterol"),
  exposure = c("smoking_status_refined", "diet_score", "diet_score", 
               "alcohol_consumption", "diet_score", "met_min_week_all")
)

# ============================================================
# END PLACEHOLDER
# ============================================================



# Mediation Analysis Models -----------------------------------------------
# References: https://nicolarighetti.github.io/mediation-moderation-and-conditional-process-analysis-with-r/mediation.html
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6341620/
# Analysis assumes parallel mediators


# Total effect model: CVD ~ confounders + exposures (no mediators) --------

form_total <- as.formula(
  paste("has_cvd ~", paste(c(confounders, selected_exposures), collapse = " + "))
)
model_total <- glm(form_total, data = train, family = binomial)
summary(model_total)


# Direct effect model (b-path): CVD ~ confounders + exposures + mediators -

form_direct <- as.formula(
  paste("has_cvd ~", paste(c(confounders, selected_exposures, selected_mediators), collapse = " + "))
)
model_direct <- glm(form_direct, data = train, family = binomial)


# A-path model: mediator ~ confounders + selected exposures ---------------
# Fit linear models as outcomes are continuous

a_path_models <- list()
a_path_coefs <- list()

# Extract models and coefficients for each exposure -> mediator effect 
for (med in selected_mediators) {
  med_exposures <- lasso2$exposure[lasso2$mediator == med]
  
  form_a <- as.formula(
    paste(med, "~", paste(c(confounders, med_exposures), collapse = " + "))
  )
  
  a_path_models[[med]] <- lm(form_a, data = train)
  
  # Extract exposure coefficients by matching prefix (handles categorical variables)
  all_coefs <- coef(a_path_models[[med]])
  exp_coefs <- list()
  
  for (exp in med_exposures) {
    matching <- grep(paste0("^", exp), names(all_coefs), value = TRUE)
    matching <- setdiff(matching, "(Intercept)")
    exp_coefs[[exp]] <- all_coefs[matching]
  }
  
  a_path_coefs[[med]] <- exp_coefs
  cat("  ", med, "~ ", paste(med_exposures, collapse = " + "), "\n")
}


# Extract b coefficients --------------------------------------------------

# Mediators are continuous so each has a single b coefficient
b_coefs <- coef(model_direct)[selected_mediators]
b_coefs


# Compute indirect effects (ab) per pathway -------------------------------
cat("Computing indirect effects...\n")
pathways <- data.frame()

for (med in selected_mediators) {
  med_exposures <- lasso2$exposure[lasso2$mediator == med]
  b <- b_coefs[med]
  
  for (exp in med_exposures) {
    a_vals <- a_path_coefs[[med]][[exp]]
    
    for (level in names(a_vals)) {
      ab <- a_vals[level] * b
      
      pathways <- rbind(pathways, data.frame(
        exposure = exp,
        level = level,
        mediator = med,
        a_coef = a_vals[level],
        b_coef = b,
        ab = ab,
        row.names = NULL
      ))
    }
  }
}

print(pathways)


# Total indirect per exposure level ---------------------------------------

total_indirect <- pathways %>%
  group_by(level) %>%
  summarise(
    exposure = first(exposure),
    total_ab = sum(ab),
    .groups = "drop"
  )


# Get c' (direct effect) per exposure level from direct model -------------

direct_coefs <- coef(model_direct)

# Match each level to its direct effect coefficient
total_indirect <- total_indirect %>%
  mutate(
    c_prime = direct_coefs[level],
    proportion_mediated = total_ab / (total_ab + c_prime) * 100
  )


# Sense check: ab vs c-c' -------------------------------------------------

# Get total effect coefficients
total_coefs <- coef(model_total)

sense_check <- total_indirect %>%
  mutate(
    c_total = total_coefs[level],
    c_minus_c_prime = c_total - c_prime,
    ab_method = total_ab,
    same_sign = sign(c_minus_c_prime) == sign(ab_method)
  ) %>%
  select(level, exposure, ab_method, c_minus_c_prime, same_sign)

print("Sense check whether difference between total and direct effects are loosely equal to calculated indirect effect per exposure")
print(sense_check)



# Bootstrap CI for ab -----------------------------------------------------

# Bootstrap function for indirect effects
compute_indirect <- function(data, indices) {
  d <- data[indices, ]
  
  # Fit b-path model (direct effect model)
  form_direct <- as.formula(
    paste("has_cvd ~", paste(c(confounders, selected_exposures, selected_mediators), collapse = " + "))
  )
  model_b <- glm(form_direct, data = d, family = binomial)
  b_coefs <- coef(model_b)[selected_mediators]
  
  # Fit a-path models and compute ab for each pathway
  ab_results <- c()
  
  for (med in selected_mediators) {
    med_exposures <- lasso2$exposure[lasso2$mediator == med]
    b <- b_coefs[med]
    
    form_a <- as.formula(
      paste(med, "~", paste(c(confounders, med_exposures), collapse = " + "))
    )
    model_a <- lm(form_a, data = d)
    all_coefs <- coef(model_a)
    
    for (exp in med_exposures) {
      matching <- grep(paste0("^", exp), names(all_coefs), value = TRUE)
      matching <- setdiff(matching, "(Intercept)")
      
      for (level in matching) {
        ab <- all_coefs[level] * b
        ab_results <- c(ab_results, ab)
      }
    }
  }
  
  return(ab_results)
}

# Run bootstrap
cat("Running bootstrap (this may take a few minutes)...\n")
set.seed(123)

n_cores <- parallel::detectCores()

# boot_results <- boot(
#   data = train,
#   statistic = compute_indirect,
#   R = 100,
#   parallel = "multicore",
#   ncpus = n_cores
# )

system.time({
  boot_results <- boot(
    data = train,
    statistic = compute_indirect,
    R = 10
  )
})


 # Extract CIs
boot_cis <- data.frame()

for (i in 1:nrow(pathways)) {
  ci <- boot.ci(boot_results, index = i, type = "perc")
  
  boot_cis <- rbind(boot_cis, data.frame(
    exposure = pathways$exposure[i],
    level = pathways$level[i],
    mediator = pathways$mediator[i],
    ab = pathways$ab[i],
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

# Add significance flag (CI doesn't cross zero)
boot_cis <- boot_cis %>%
  mutate(significant = !(ci_lower <= 0 & ci_upper >= 0))

print(boot_cis)



# Results -----------------------------------------------------------------


# Detailed results --------------------------------------------------------

# Full table of all coefficients and confidence intervals
mediation_detail <- boot_cis %>%
  left_join(
    pathways %>% select(exposure, level, mediator, a_coef, b_coef),
    by = c("exposure", "level", "mediator")
  ) %>%
  select(exposure, level, mediator, a_coef, b_coef, ab, ci_lower, ci_upper, significant) %>%
  arrange(exposure, mediator, level)

print(mediation_detail)

# Summary per exposure ----------------------------------------------------

# Get c' from direct model and total effect from total model
direct_coefs <- coef(model_direct)
total_coefs <- coef(model_total)

mediation_summary <- mediation_detail %>%
  group_by(level, exposure) %>%
  summarise(
    n_pathways = n(),
    total_indirect = sum(ab),
    any_significant = any(significant),
    .groups = "drop"
  ) %>%
  mutate(
    c_prime = direct_coefs[level],
    c_total = total_coefs[level],
    proportion_mediated = total_indirect / (total_indirect + c_prime) * 100
  ) %>%
  select(exposure, level, c_total, c_prime, total_indirect, proportion_mediated, n_pathways, any_significant) %>%
  arrange(exposure, level)

print(mediation_summary)

# Forest plot of indirect effects -----------------------------------------

mediation_detail %>%
  rowwise() %>%
  mutate(
    # Row-wise stripping of variable name from coefficient name
    category = sub(paste0("^", exposure), "", level),
    category = gsub("_", " ", category) %>% trimws(),
    
    # Look up short names
    exp_short = get_display_name(exposure),
    med_short = get_display_name(mediator),
    
    # For continuous variables, category will be empty or equal to cleaned variable name
    is_continuous = (category == "" | category == gsub("_", " ", exposure)),
    
    pathway_label = ifelse(
      is_continuous,
      paste0(exp_short, "  →  ", med_short),
      paste0(exp_short, ": ", category, "  →  ", med_short)
    )
  ) %>%
  ungroup() %>%
  mutate(pathway_label = reorder(pathway_label, ab)) %>%
  ggplot(aes(y = pathway_label, x = ab)) +
  geom_point(aes(colour = significant), size = 3) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper, colour = significant),
                width = 0.2, orientation = "y") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6"),
                      labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_x_continuous(labels = scales::scientific) +
  labs(
    title = "Indirect Effects with 95% Bootstrap CIs",
    x = "Indirect Effect",
    y = NULL,
    colour = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    plot.margin = margin(10, 60, 10, 10)
  )


# DAG visualisation -------------------------------------------------------

dot_code <- c('digraph DAG {', 
              '  rankdir = LR',
              '  node [style=filled, shape=box, fontsize=10]')

# Confounder nodes
for (c in confounders) {
  dot_code <- c(dot_code, paste0('  "', c, '" [fillcolor="#D5E8D4"]'))
}
# Exposure nodes
for (e in selected_exposures) {
  dot_code <- c(dot_code, paste0('  "', e, '" [fillcolor="#DAE8FC"]'))
}
# Mediator nodes
for (m in selected_mediators) {
  dot_code <- c(dot_code, paste0('  "', m, '" [fillcolor="#FFE6CC"]'))
}
dot_code <- c(dot_code, '  "CVD" [fillcolor="#F8CECC"]')

# Edges

# Confounder edges to everything
for (c in confounders) {
  dot_code <- c(dot_code, paste0('  "', c, '" -> "CVD" [style=dashed, color=gray]'))
  for (exp in selected_exposures) {
    dot_code <- c(dot_code, paste0('  "', c, '" -> "', exp, '" [style=dashed, color=gray]'))
  }
  for (med in selected_mediators) {
    dot_code <- c(dot_code, paste0('  "', c, '" -> "', med, '" [style=dashed, color=gray]'))
  }
}

for (exp in selected_exposures) dot_code <- c(dot_code, paste0('  "', exp, '" -> "CVD"'))
for (i in 1:nrow(lasso2)) {
  dot_code <- c(dot_code, paste0('  "', lasso2$exposure[i], '" -> "', lasso2$mediator[i], '"'))
}

for (med in selected_mediators) dot_code <- c(dot_code, paste0('  "', med, '" -> "CVD"'))

dot_code <- c(dot_code, '}')



# Export tables and plots -------------------------------------------------

write.csv(mediation_detail, here("02_results", "01_tables", "01_multivariate", "mediation_detail.csv"), row.names = FALSE)
ggsave(here("02_results", "00_figures", "03_multivariate", "mediation_forest_plot.png"), width = 8, height = 6, dpi = 300)
write.csv(mediation_summary, here("02_results", "01_tables", "01_multivariate", "mediation_summary.csv"), row.names = FALSE)
writeLines(dot_code, here("02_results", "00_figures", "03_multivariate", "dag.dot"))

