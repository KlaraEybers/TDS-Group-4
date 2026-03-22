rm(list = ls())

# Import packages ---------------------------------------------------------

library(here)
library(dplyr)
library(ggplot2)
library(purrr)

# Load data ---------------------------------------------------------------

df <- readRDS(here("00_data", "01_processed", "ukb_train_imputed.rds"))

# Basic checks for CVD ----------------------------------------------------

stopifnot("has_cvd" %in% names(df))
stopifnot(all(df$has_cvd %in% c(0, 1)))

# Domains (exclude non-domain columns later) ------------------------------

domains <- list(
  `Socio-demographics` = c(
    "age","sex","ethnicity","household_income","highest_qualification",
    "employment","housing","household_size","recruitment_date","death_date"
  ),
  `Physical measures` = c(
    "sbp","dbp","grip_strength","whr","heel_bone_density","body_fat_mass"
  ),
  `Diet` = c(
    "diet_score","alcohol_consumption","salt_added"
  ),
  `Lifestyle` = c(
    "sleep_hrs","sleep_cat","insomnia","smoking_status_refined","pack_yrs",
    "passive_total_hrs","passive_cat","met_min_week_all"
  ),
  `Health and medical history` = c(
    "long_ill","diabetes_doctor_diagnosed","serious_dx","age_high_bp_dx","age_diabetes_dx",
    "angina_dx_age","mi_dx_age","dvt_dx_age","pe_dx_age","stroke_dx_age",
    "heart_problems_collapsed","meds_chol_bp_diabetes_hormone_collapsed",
    "meds_chol_bp_diabetes_collapsed","n_previous_pregnancies_collapsed"
  ),
  `Mental health` = c(
    "mh_sought_help","work_job_satisfaction_num","neuroticism_score"
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

# Variable -> domain map --------------------------------------------------

var_domain <- bind_rows(lapply(names(domains), function(d) {
  data.frame(variable = domains[[d]], domain = d, stringsAsFactors = FALSE)
}))

# Keep only variables present in dataset ----------------------------------

var_domain <- var_domain %>% filter(variable %in% names(df))

# Remove identifier/outcome columns ---------------------------------------

non_domain <- c("eid", "first_cvd_date", "has_cvd")
var_domain <- var_domain %>% filter(!(variable %in% non_domain))

# Helper to escape regex characters ---------------------------------------

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

# Univariate logistic regression function ---------------------------------
# For numeric variables:
#   OR and 95% CI correspond to a 1-unit increase
# For factor variables:
#   overall p-value is from drop1(..., test = "Chisq")
#   plotted OR and 95% CI correspond to the dummy coefficient with the
#   largest absolute log(OR), matching the previous volcano-plot summary logic

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

# Run all univariate models -----------------------------------------------

vars <- var_domain$variable

res <- purrr::map_dfr(vars, ~ fit_one(df, .x)) %>%
  left_join(var_domain, by = c("variable" = "variable")) %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH"),
    or = exp(effect_log_or),
    lower_or = exp(lower_log_or),
    upper_or = exp(upper_log_or)
  )

# Save results ------------------------------------------------------------

saveRDS(res, here("00_data", "02_analysis", "univariate_logistic_results.rds"))

# Prepare plotting data ---------------------------------------------------

plot_df <- res %>%
  filter(!is.na(or), !is.na(lower_or), !is.na(upper_or), !is.na(domain)) %>%
  mutate(sig_bh = p_adj_bh < 0.05)

# Keep x-axis centred on OR = 1 and identical across all plots
# Done by using log scale and symmetric limits around 1 on the log scale

max_abs_log <- max(
  abs(log(plot_df$lower_or)),
  abs(log(plot_df$upper_or)),
  na.rm = TRUE
)

x_limits <- exp(c(-max_abs_log, max_abs_log))

# Order variables within each domain for forest plot

plot_df_domain <- plot_df %>%
  group_by(domain) %>%
  arrange(or, .by_group = TRUE) %>%
  mutate(variable_plot = paste0(variable, "___", row_number())) %>%
  ungroup()

# Forest plot by domain ---------------------------------------------------
# x-axis = OR with 95% CI
# vertical line at OR = 1
# same x-scale across all facets
# free_y so each domain shows its own variables cleanly
# red = BH-adjusted p < 0.05

p_forest_domain <- ggplot(
  plot_df_domain,
  aes(y = variable_plot, x = or, xmin = lower_or, xmax = upper_or, color = sig_bh)
) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(height = 0) +
  geom_point(size = 2) +
  scale_x_log10(limits = x_limits) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
  facet_wrap(~ domain, scales = "free_y") +
  labs(
    title = "Univariate associations with has_cvd (logistic regression)",
    subtitle = "Forest plots by domain: odds ratios with 95% CIs; red = BH-adjusted p < 0.05",
    x = "Odds ratio (log scale)",
    y = NULL,
    color = "BH significant"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10)
  )

ggsave(
  filename = here("02_results", "00_figures", "02_univariate", "forest_univariate_by_domain.png"),
  plot = p_forest_domain,
  width = 16, height = 12, dpi = 300
)

# Single-panel forest plot ------------------------------------------------

plot_df_single <- plot_df %>%
  arrange(domain, or) %>%
  mutate(variable_plot = paste0(domain, "___", variable))

p_forest_single <- ggplot(
  plot_df_single,
  aes(y = variable_plot, x = or, xmin = lower_or, xmax = upper_or, color = sig_bh)
) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(height = 0) +
  geom_point(size = 2) +
  scale_x_log10(limits = x_limits) +
  scale_y_discrete(labels = function(x) sub(".*___", "", x)) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
  labs(
    title = "Univariate associations with has_cvd (all domains)",
    subtitle = "Forest plot: odds ratios with 95% CIs; red = BH-adjusted p < 0.05",
    x = "Odds ratio (log scale)",
    y = NULL,
    color = "BH significant"
  ) +
  theme_bw()

ggsave(
  filename = here("02_results", "00_figures", "02_univariate", "forest_univariate_singlepanel.png"),
  plot = p_forest_single,
  width = 14, height = 16, dpi = 300
)
