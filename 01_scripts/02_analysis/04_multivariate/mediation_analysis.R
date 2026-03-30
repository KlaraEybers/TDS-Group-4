# ==============================================================================
# Mediation Analysis: Total, Direct, and Indirect Effects 
# ==============================================================================
rm(list = ls())


# Import packages ---------------------------------------------------------

library(here)
library(dplyr)
library(boot)
library(parallel)
library(ggplot2)
library(scales)


# Paths -------------------------------------------------------------------

source(here("01_scripts", "00_variable_config.R"))
source(here("01_scripts", "00_variable_dictionary.R"))

data_dir   <- data_input_dir
fig_dir    <- figures_multi_dir
tab_dir    <- tables_multi_dir
lasso1_dir <- lasso_outcome_dir
lasso2_dir <- lasso_mediator_dir


# Load data ---------------------------------------------------------------

train <- readRDS(file.path(data_dir, train_file))
lasso1 <- read.csv(file.path(lasso1_dir, "selected_variable_names.csv"))
lasso2 <- read.csv(file.path(lasso2_dir, "all_selected_mediator_exposure_pairs.csv"))



# Data checks and definitions ---------------------------------------------

# Ensure binary outcome (0/1)
if(is.character(train$has_cvd) | is.factor(train$has_cvd)){
  train$has_cvd <- ifelse(train$has_cvd == "Event", 1, 0)
}


# Prepare selected variables ----------------------------------------------

# Extract Lasso 1 results
selected_exposures <- lasso1$variable[lasso1$type == "exposure"]
selected_mediators <- lasso1$variable[lasso1$type == "mediator"]

mediation_exposures <- unique(lasso2$exposure)

# Ensure all variables exist in the data
selected_exposures  <- intersect(selected_exposures, names(train))
selected_mediators  <- intersect(selected_mediators, names(train))
mediation_exposures <- intersect(mediation_exposures, names(train))
all_exposures       <- intersect(union(selected_exposures, mediation_exposures), names(train))

# Also filter lasso2 pairs to only include variables present in data
lasso2 <- lasso2 %>%
  filter(exposure %in% names(train), mediator %in% names(train))


cat("LASSO 1 exposures:", length(selected_exposures), "\n")
cat("LASSO 2 exposures:", length(mediation_exposures), "\n")
cat("Union (for mediation models):", length(all_exposures), "\n")




# Mediation Analysis Models -----------------------------------------------
# References: https://nicolarighetti.github.io/mediation-moderation-and-conditional-process-analysis-with-r/mediation.html
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6341620/
# Analysis assumes parallel mediators


# Total effect model: CVD ~ confounders + exposures (no mediators) --------

form_total <- as.formula(
  paste("has_cvd ~", paste(c(confounders, all_exposures), collapse = " + "))
)
model_total <- glm(form_total, data = train, family = binomial)
summary(model_total)

total_coefs_df <- as.data.frame(cbind(
  coef(model_total),
  confint(model_total)
))
colnames(total_coefs_df) <- c("c_total", "ci_lower_95", "ci_upper_95")
total_coefs_df$exposure_level <- rownames(total_coefs_df)


# Direct effect model (b-path): CVD ~ confounders + exposures + mediators
form_direct <- as.formula(
  paste("has_cvd ~", paste(c(confounders, all_exposures, selected_mediators), collapse = " + "))
)
model_direct <- glm(form_direct, data = train, family = binomial)
summary(model_direct)


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

boot_results <- boot(
  data = train,
  statistic = compute_indirect,
  R = 500,
  parallel = "multicore",
  ncpus = n_cores
)

# system.time({
#   boot_results <- boot(
#     data = train,
#     statistic = compute_indirect,
#     R = 10
#   )
# })


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

# Add exponentiated results 
boot_cis <- boot_cis %>%
  mutate(
    ab_exp = exp(ab),
    ci_lower_exp = exp(ci_lower),
    ci_upper_exp = exp(ci_upper)
  )

print(boot_cis)


# Results -----------------------------------------------------------------


# Detailed results --------------------------------------------------------

# Full table of all coefficients and confidence intervals
mediation_detail <- boot_cis %>%
  left_join(
    pathways %>% select(exposure, level, mediator, a_coef, b_coef),
    by = c("exposure", "level", "mediator")
  ) %>%
  select(exposure, level, mediator, a_coef, b_coef, ab, ab_exp, ci_lower, ci_upper, ci_lower_exp, ci_upper_exp, significant) %>%
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
    total_indirect_exp = exp(total_indirect),
    any_significant = any(significant),
    .groups = "drop"
  ) %>%
  mutate(
    c_prime = direct_coefs[level],
    c_total = total_coefs[level],
    c_prime_exp = exp(c_prime),
    c_total_exp = exp(c_total)
  ) %>%
  left_join(total_coefs_df[, c("exposure_level", "ci_lower_95", "ci_upper_95")],
            by = c("level" = "exposure_level")) %>%
  mutate(
    ci_lower_95_exp = exp(ci_lower_95),
    ci_upper_95_exp = exp(ci_upper_95),
    proportion_mediated = total_indirect / (total_indirect + c_prime) * 100,
    c_total_significant = !(ci_lower_95 < 0 & ci_upper_95 > 0)
  ) %>%
  select(exposure, level, c_total, c_total_exp, ci_lower_95, ci_upper_95, ci_lower_95_exp, ci_upper_95_exp, c_total_significant,
         c_prime, c_prime_exp, total_indirect, total_indirect_exp, proportion_mediated, n_pathways, any_significant) %>%
  arrange(exposure, level)

print(mediation_summary)

# Forest plot of indirect effects -----------------------------------------

# Indirect effects forest plot data
indirect_plot_data <- boot_cis %>%
  rowwise() %>%
  mutate(
    category = sub(paste0("^", exposure), "", level),
    category = gsub("_", " ", category) %>% trimws(),
    exp_short = get_display_name(exposure),
    med_short = get_display_name(mediator),
    is_continuous = (category == "" | category == gsub("_", " ", exposure)),
    pathway_label = ifelse(
      is_continuous,
      paste0(exp_short, "  →  ", med_short),
      paste0(exp_short, ": ", category, "  →  ", med_short)
    )
  ) %>%
  ungroup()

write.csv(indirect_plot_data, file.path(tab_dir, "indirect_effects_plot_data.csv"), row.names = FALSE)

indirect_plot_data %>%
  mutate(pathway_label = reorder(pathway_label, ab_exp)) %>%
  ggplot(aes(y = pathway_label, x = ab_exp)) +
  geom_point(aes(colour = significant), size = 3) +
  geom_errorbar(aes(xmin = ci_lower_exp, xmax = ci_upper_exp, colour = significant),
                width = 0.2, orientation = "y") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "#95A5A6"),
                      labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.001, drop0trailing = TRUE)) +
  labs(
    title = "Indirect Effects with 95% Bootstrap CIs",
    x = "Indirect Effect (exponentiated)",
    y = "Exposure → Mediator Pathway",
    colour = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    plot.margin = margin(10, 60, 10, 10)
  )

ggsave(file.path(fig_dir, "mediation_forest_plot.png"), width = 8, height = 6, dpi = 300)


# Total effects forest plot -----------------------------------------------

total_plot_data <- total_coefs_df %>%
  filter(
    !exposure_level %in% c("(Intercept)", confounders),
    !grepl(paste(paste0("^", confounders), collapse = "|"), exposure_level)
  ) %>%
  mutate(
    significant = !(ci_lower_95 < 0 & ci_upper_95 > 0),
    c_total_exp = exp(c_total),
    ci_lower_95_exp = exp(ci_lower_95),
    ci_upper_95_exp = exp(ci_upper_95),
    base_var = sapply(exposure_level, function(lvl) {
      m <- Filter(function(e) startsWith(lvl, e), all_exposures)
      if (length(m) == 0) return(lvl) else return(m[[1]])
    }),
    category = mapply(function(lvl, bv) {
      gsub("_", " ", sub(paste0("^", bv), "", lvl)) %>% trimws()
    }, exposure_level, base_var),
    display = sapply(base_var, get_display_name),
    exposure_label = ifelse(category == "", display, paste0(display, ": ", category))
  )

write.csv(total_plot_data, file.path(tab_dir, "total_effects_plot_data.csv"), row.names = FALSE)

total_plot_data %>%
  mutate(exposure_label = reorder(exposure_label, c_total_exp)) %>%
  ggplot(aes(y = exposure_label, x = c_total_exp)) +
  geom_point(aes(colour = significant), size = 3) +
  geom_errorbar(aes(xmin = ci_lower_95_exp, xmax = ci_upper_95_exp, colour = significant),
                width = 0.2, orientation = "y") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "#95A5A6"),
                      labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = "Total Effects with 95% CIs",
    x = "Odds Ratio (Total Effect)",
    y = "Exposure",
    colour = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    plot.margin = margin(10, 60, 10, 10)
  )

ggsave(file.path(fig_dir, "total_effects_forest_plot.png"), width = 8, height = 6, dpi = 300)


# DAG visualisation -------------------------------------------------------
dot_code <- c('digraph DAG {', 
              '  rankdir = LR',
              '  node [style=filled, shape=box, fontsize=10]')

# Confounder nodes - Light Grey
for (c in confounders) {
  dot_code <- c(dot_code, paste0('  "', get_display_name(c), '" [fillcolor="#EBEEEE", fontcolor="#000000"]'))
}

# LASSO 1 exposure nodes - Imperial Blue
for (e in selected_exposures) {
  dot_code <- c(dot_code, paste0('  "', get_display_name(e), '" [fillcolor="#006EAF", fontcolor="#FFFFFF"]'))
}

# LASSO 2 only exposure nodes - Light Blue
lasso2_only <- setdiff(mediation_exposures, selected_exposures)
for (e in lasso2_only) {
  dot_code <- c(dot_code, paste0('  "', get_display_name(e), '" [fillcolor="#D4EFFC", fontcolor="#000000"]'))
}

# Mediator nodes - Purple
for (m in selected_mediators) {
  dot_code <- c(dot_code, paste0('  "', get_display_name(m), '" [fillcolor="#653098", fontcolor="#FFFFFF"]'))
}

# CVD node - Navy
dot_code <- c(dot_code, '  "CVD" [fillcolor="#003E74", fontcolor="#FFFFFF"]')

# Confounder edges - to all exposures (LASSO 1 + LASSO 2) and mediators
for (c in confounders) {
  dot_code <- c(dot_code, paste0('  "', get_display_name(c), '" -> "CVD" [style=dashed, color="#CCCCCC"]'))
  for (exp in all_exposures) {
    dot_code <- c(dot_code, paste0('  "', get_display_name(c), '" -> "', get_display_name(exp), '" [style=dashed, color="#CCCCCC"]'))
  }
  for (med in selected_mediators) {
    dot_code <- c(dot_code, paste0('  "', get_display_name(c), '" -> "', get_display_name(med), '" [style=dashed, color="#CCCCCC"]'))
  }
}

# Exposure -> CVD edges (LASSO 1 only)
for (exp in selected_exposures) {
  dot_code <- c(dot_code, paste0('  "', get_display_name(exp), '" -> "CVD"'))
}

# Exposure -> Mediator edges (from lasso2)
for (i in 1:nrow(lasso2)) {
  dot_code <- c(dot_code, paste0('  "', get_display_name(lasso2$exposure[i]), '" -> "', get_display_name(lasso2$mediator[i]), '"'))
}

# Mediator -> CVD edges
for (med in selected_mediators) {
  dot_code <- c(dot_code, paste0('  "', get_display_name(med), '" -> "CVD"'))
}

dot_code <- c(dot_code, '}')

# Heatmap of proportion mediated ------------------------------------------

mediation_heatmap_data <- mediation_detail %>%
  distinct(exposure, level, mediator, .keep_all = TRUE) %>%
  mutate(
    category = mapply(function(e, l) sub(paste0("^", e), "", l), exposure, level),
    category = gsub("_", " ", category) %>% trimws(),
    exp_short = get_display_name(exposure),
    med_short = get_display_name(mediator),
    y_label = ifelse(category == "" | category == gsub("_", " ", exposure),
                     exp_short,
                     paste0(exp_short, ": ", category))
  )

write.csv(mediation_heatmap_data, file.path(tab_dir, "heatmap_plot_data.csv"), row.names = FALSE)

ggplot(mediation_heatmap_data, aes(x = med_short, y = y_label, fill = ab_exp)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = formatC(ab_exp, format = "f", digits = 4)), size = 3) +
  scale_fill_gradient2(
    low = "lightblue", mid = "white", high = "blue",
    midpoint = 1,
    name = "Indirect effect (exponentiated)",
    labels = scales::comma
  ) +
  facet_grid(exp_short ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  labs(x = "Mediator", y = "Exposure level",
       title = "Indirect effects by exposure-mediator pathway") +
  theme(
    strip.text.y = element_text(angle = 0),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

ggsave(file.path(fig_dir, "mediation_heatmap.png"), width = 10, height = 7, dpi = 300)


# Export tables and plots -------------------------------------------------

write.csv(mediation_detail, file.path(tab_dir, "mediation_detail.csv"), row.names = FALSE)
write.csv(mediation_summary, file.path(tab_dir, "mediation_summary.csv"), row.names = FALSE)
writeLines(dot_code, file.path(fig_dir, "dag.dot"))
