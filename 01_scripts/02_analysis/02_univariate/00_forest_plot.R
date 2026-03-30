rm(list = ls())

# Import packages
library(here)
library(dplyr)
library(ggplot2)
library(purrr)

# Load data

# Use sensitivity imputation which includes all variables (including anthropometrics)
df <- readRDS(here("00_data", "01_processed", "03_imputation",
                   "sensitivity_anthro_confounder", "df_full_imputed.rds"))
# Basic checks for CVD

stopifnot("has_cvd" %in% names(df))
stopifnot(all(df$has_cvd %in% c(0, 1)))

# Domains

domains <- list(
  `Socio-demographics` = c(
    "age","sex","ethnicity","household_income","highest_qualification",
    "employment","housing","household_size"
  ),
  `Anthropometrics` = c(
    "sbp","dbp","grip_strength","whr","heel_bone_density","body_fat_mass"
  ),
  `Diet` = c(
    "diet_score","alcohol_consumption","salt_added"
  ),
  `Lifestyle` = c(
    "sleep_hrs","sleep_cat","insomnia","smoking_status_refined","pack_yrs",
    "passive_total_hrs","passive_cat","met_min_week_all"
  ),
  # `Health and medical history` = c(
  #   "long_ill","diabetes_doctor_diagnosed","serious_dx","age_high_bp_dx","age_diabetes_dx",
  #   "angina_dx_age","mi_dx_age","dvt_dx_age","pe_dx_age","stroke_dx_age",
  #   "heart_problems_collapsed","meds_chol_bp_diabetes_hormone_collapsed",
  #   "meds_chol_bp_diabetes_collapsed","n_previous_pregnancies_collapsed"
  # ),
  `Mental health` = c(
    "mh_sought_help", "satisfaction_mean","neuroticism_score"
  ),
  `Early life factors` = c(
    "breastfed","comparative_body_size_10","maternal_smoking_birth","birth_weight"
  ),
  `Environment` = c(
    "no2_2010","pm2.5_2010","pm10_2010","urban_rural_3lvl"
  ),
  `Biological samples` = c(
    "wbc","platelets","albumin","alp","alt","apoa1","apob","ast","bilirubin",
    "total_cholesterol","creatinine","crp","cystatin_c","ggt","hba1c",
    "hdl_cholesterol","ldl_direct","lpa","triglycerides"
  )
)

# Variable -> domain map

var_domain <- bind_rows(lapply(names(domains), function(d) {
  data.frame(variable = domains[[d]], domain = d, stringsAsFactors = FALSE)
}))

# Keep only variables present in dataset

var_domain <- var_domain %>% filter(variable %in% names(df))

# Remove identifier/outcome columns

non_domain <- c("eid", "first_cvd_date", "has_cvd")
var_domain <- var_domain %>% filter(!(variable %in% non_domain))

# Helper to escape regex characters

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

# Univariate logistic regression function

fit_one <- function(dat, v) {
  x <- dat[[v]]
  
  # Date -> numeric
  if (inherits(x, "Date")) {
    dat[[v]] <- as.numeric(x)
    x <- dat[[v]]
  }
  
  # Character -> factor
  if (is.character(x)) {
    dat[[v]] <- as.factor(x)
    x <- dat[[v]]
  }
  
  # Skip constants / all-NA variables
  ux <- unique(x[!is.na(x)])
  if (length(ux) <= 1) {
    return(data.frame(
      variable = v,
      level = NA_character_,
      effect_log_or = NA_real_,
      lower_log_or = NA_real_,
      upper_log_or = NA_real_,
      p_value = NA_real_,
      n_nonmissing = sum(stats::complete.cases(dat[, c("has_cvd", v)])),
      stringsAsFactors = FALSE
    ))
  }
  
  f <- stats::as.formula(paste("has_cvd ~", v))
  
  m <- tryCatch(
    stats::glm(f, data = dat, family = stats::binomial()),
    error = function(e) NULL
  )
  
  if (is.null(m)) {
    return(data.frame(
      variable = v,
      level = NA_character_,
      effect_log_or = NA_real_,
      lower_log_or = NA_real_,
      upper_log_or = NA_real_,
      p_value = NA_real_,
      n_nonmissing = sum(stats::complete.cases(dat[, c("has_cvd", v)])),
      stringsAsFactors = FALSE
    ))
  }
  
  # Overall p-value
  p <- tryCatch({
    d1 <- stats::drop1(m, test = "Chisq")
    as.numeric(d1[v, "Pr(>Chi)"])
  }, error = function(e) NA_real_)
  
  co <- stats::coef(m)
  vc <- tryCatch(stats::vcov(m), error = function(e) NULL)
  
  effect <- NA_real_
  lower  <- NA_real_
  upper  <- NA_real_
  level  <- NA_character_
  
  if (!is.null(vc)) {
    
    if (is.numeric(x) || is.integer(x)) {
      
      if (v %in% names(co) && v %in% rownames(vc)) {
        effect <- unname(co[v])
        se <- sqrt(vc[v, v])
        lower <- effect - 1.96 * se
        upper <- effect + 1.96 * se
        level <- "1-unit increase"
      }
      
    } else if (is.factor(x) || is.ordered(x)) {
      
      idx <- grep(paste0("^", escape_regex(v)), names(co))
      
      if (length(idx) > 0) {
        j <- idx[which.max(abs(co[idx]))]
        effect <- unname(co[j])
        se <- sqrt(vc[j, j])
        lower <- effect - 1.96 * se
        upper <- effect + 1.96 * se
        level <- names(co)[j]
      }
    }
  }
  
  data.frame(
    variable = v,
    level = level,
    effect_log_or = effect,
    lower_log_or = lower,
    upper_log_or = upper,
    p_value = p,
    n_nonmissing = sum(stats::complete.cases(dat[, c("has_cvd", v)])),
    stringsAsFactors = FALSE
  )
}

# Run all univariate models

vars <- var_domain$variable

res <- purrr::map_dfr(vars, ~ fit_one(df, .x)) %>%
  left_join(var_domain, by = c("variable" = "variable")) %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH"),
    or = exp(effect_log_or),
    lower_or = exp(lower_log_or),
    upper_or = exp(upper_log_or)
  )

# Save results

saveRDS(res, here("00_data", "02_analysis", "univariate_logistic_results.rds"))


# Univariate logistic regression plots

# Load variable dictionary
res <- readRDS(here("00_data", "02_analysis", "univariate_logistic_results.rds"))
source(here("01_scripts", "00_variable_dictionary.R"))

# Domain display names
domain_display_names <- c(
  "sociodemographic" = "Socio-demographics",
  "physical_measures" = "Physical Measures",
  "diet" = "Diet",
  "lifestyle" = "Lifestyle",
  # "health_medical" = "Health and Medical History",
  "mental_health" = "Mental Health",
  "early_life" = "Early Life Factors",
  "environment" = "Environment",
  "biological_samples" = "Biological Samples"
)

get_domain_display_name <- function(x) {
  unname(ifelse(x %in% names(domain_display_names), domain_display_names[x], x))
}

# Prepare plotting data
plot_df <- res %>%
  filter(!is.na(or), !is.na(lower_or), !is.na(upper_or), !is.na(domain)) %>%
  mutate(
    sig_bh = factor(p_adj_bh < 0.05, levels = c(FALSE, TRUE), labels = c("No", "Yes")),
    variable_label = get_display_name(variable),
    domain_label = get_domain_display_name(domain),
    effect_direction = factor(
      if_else(or < 1, "Negative OR", "Positive OR"),
      levels = c("Negative OR", "Positive OR")
    )
  )

# Consistent x-axis limits
max_abs_log <- max(
  abs(log(plot_df$lower_or)),
  abs(log(plot_df$upper_or)),
  na.rm = TRUE
)

x_limits <- exp(c(-max_abs_log, max_abs_log))


# Forest plot by domain
plot_df_domain <- plot_df %>%
  group_by(domain_label) %>%
  arrange(or, .by_group = TRUE) %>%
  mutate(
    variable_plot = paste0(variable_label, "___", row_number())
  ) %>%
  ungroup()

p_forest_domain <- ggplot(
  plot_df_domain,
  aes(y = variable_plot, x = or, xmin = lower_or, xmax = upper_or, color = sig_bh)
) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_errorbar(aes(xmin = lower_or, xmax = upper_or), height = 0, linewidth = 0.5, orientation = "y") +
  geom_point(size = 2) +
  scale_x_log10(limits = x_limits) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_color_manual(values = c("No" = "grey70", "Yes" = "red")) +
  facet_wrap(~ domain_label, scales = "free_y") +
  labs(
    title = "Univariate associations with CVD",
    subtitle = "Forest plots by domain: odds ratios with 95% CIs; red = BH-adjusted p < 0.05",
    x = "Odds ratio (log scale)",
    y = NULL,
    color = "BH significant"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  filename = here("02_results", "00_figures", "02_univariate", "forest_univariate_by_domain.png"),
  plot = p_forest_domain,
  width = 16,
  height = 12,
  dpi = 300
)

# Single-panel forest plot
plot_df_single <- plot_df %>%
  arrange(domain_label, or) %>%
  mutate(
    variable_plot = paste0(domain_label, "___", variable_label)
  )

p_forest_single <- ggplot(
  plot_df_single,
  aes(y = variable_plot, x = or, xmin = lower_or, xmax = upper_or, color = sig_bh)
) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_errorbarh(height = 0, linewidth = 0.5) +
  geom_point(size = 2) +
  scale_x_log10(limits = x_limits) +
  scale_y_discrete(labels = function(x) sub(".*___", "", x)) +
  scale_color_manual(values = c("No" = "grey70", "Yes" = "red")) +
  labs(
    title = "Univariate associations with CVD (all domains)",
    subtitle = "Forest plot: odds ratios with 95% CIs; red = BH-adjusted p < 0.05",
    x = "Odds ratio (log scale)",
    y = NULL,
    color = "BH significant"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7)
  )

ggsave(
  filename = here("02_results", "00_figures", "02_univariate", "forest_univariate_singlepanel.png"),
  plot = p_forest_single,
  width = 14,
  height = 16,
  dpi = 300
)

