# ==============================================================================
# Script: Full Table 1
# ==============================================================================
rm(list = ls())

library(here)
library(knitr)
suppressWarnings({
  if (requireNamespace("dplyr", quietly = TRUE)) library(dplyr)
})

## ---- load data ----
cat("Loading datasets...\n")
processed_dir <- here("00_data", "01_processed")

base        <- readRDS(file.path(processed_dir, "ukb_final_cleaned.rds"))
socio       <- readRDS(file.path(processed_dir, "sociodemographic_cleaned.rds"))
phys        <- readRDS(file.path(processed_dir, "physical_cleaned.rds"))
diet        <- readRDS(file.path(processed_dir, "diet_cleaned.rds"))
bio         <- readRDS(file.path(processed_dir, "bio_samples_cleaned.rds"))
health      <- readRDS(file.path(processed_dir, "medical_history_cleaned.rds"))
lifestyle   <- readRDS(file.path(processed_dir, "lifestyle_cleaned.rds"))
environment <- readRDS(file.path(processed_dir, "environment_cleaned.rds"))
early_life  <- readRDS(file.path(processed_dir, "early_life_cleaned.rds"))

## ---- merge ----
cat("Merging datasets...\n")
dat <- Reduce(function(x, y) merge(x, y, by = "eid", all.x = TRUE),
              list(base[, c("eid","has_cvd")], socio, phys, diet, bio, health,
                   lifestyle, environment, early_life))

## ---- event group ----
has_cvd_num <- suppressWarnings(as.integer(as.character(dat$has_cvd)))
has_cvd_num[is.na(has_cvd_num) & is.logical(dat$has_cvd)] <- as.integer(dat$has_cvd[is.na(has_cvd_num) & is.logical(dat$has_cvd)])

dat$event <- ifelse(has_cvd_num == 1, "Event",
                    ifelse(has_cvd_num == 0, "No event", NA))
dat$event <- factor(dat$event, levels = c("No event","Event"))

## -----------------------------------------------------------------------
## Health composites
## -----------------------------------------------------------------------
row_any_positive <- function(df, cols) {
  if (length(cols) == 0) return(rep(NA, nrow(df)))
  m <- df[, cols, drop = FALSE]
  m <- as.data.frame(lapply(m, function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
      x_low <- tolower(x)
      return(ifelse(x_low %in% c("yes","y","1","true"), 1,
                    ifelse(x_low %in% c("no","n","0","false"), 0,
                           suppressWarnings(as.numeric(x)))))
    }
    suppressWarnings(as.numeric(x))
  }))
  apply(m, 1, function(r) {
    if (all(is.na(r))) return(NA_integer_)
    as.integer(any(r == 1, na.rm = TRUE))
  })
}

heart_cols <- grep("^heart_problems(\\.0\\.[123])?$", names(dat), value = TRUE)
dat$heart_problems_any <- row_any_positive(dat, heart_cols)

meds_horm_cols <- grep("^meds_chol_bp_diabetes_hormone(\\.0\\.[123])?$", names(dat), value = TRUE)
dat$meds_chol_bp_diabetes_hormone_any <- row_any_positive(dat, meds_horm_cols)

meds_cols <- grep("^meds_chol_bp_diabetes(\\.0\\.[12])?$", names(dat), value = TRUE)
dat$meds_chol_bp_diabetes_any <- row_any_positive(dat, meds_cols)

## ---- helpers: build meta from available columns (auto-skip missing vars) ----
build_meta <- function(domain, candidates, data_names) {
  keep <- intersect(names(candidates), data_names)
  if (length(keep) == 0) return(NULL)
  df <- data.frame(
    domain = domain,
    var = keep,
    label = vapply(keep, function(v) candidates[[v]]$label, character(1)),
    type  = vapply(keep, function(v) candidates[[v]]$type,  character(1)),
    order_in_domain = vapply(keep, function(v) candidates[[v]]$order, numeric(1)),
    stringsAsFactors = FALSE
  )
  df[order(df$order_in_domain), ]
}

## -----------------------------------------------------------------------
## Base meta
## -----------------------------------------------------------------------
meta_base <- data.frame(
  domain = c(
    rep("Sociodemographic", 9),
    rep("Physical measures", 6),
    rep("Diet & alcohol", 2),
    rep("Outcome / follow-up", 1)
  ),
  var = c(
    "age","sex","ethnicity","highest_qualification","employment",
    "household_income","housing","household_size","has_vehicle",
    "sbp","dbp","grip_strength","whr","body_fat_mass","heel_bone_density",
    "diet_score","alcohol_consumption",
    "has_cvd"
  ),
  label = c(
    "Age, years",
    "Sex",
    "Ethnicity",
    "Highest qualification",
    "Employment status",
    "Household income (before tax)",
    "Housing tenure",
    "Household size category",
    "Has vehicle in household",
    "Systolic BP, mmHg",
    "Diastolic BP, mmHg",
    "Grip strength, kg",
    "Waist-to-hip ratio (WHR)",
    "Body fat mass, kg",
    "Heel bone density T-score",
    "Diet score (0–8)",
    "Alcohol consumption category",
    "Has CVD event during follow-up"
  ),
  type = c(
    "continuous",
    rep("categorical", 8),
    rep("continuous", 6),
    "continuous",
    "categorical",
    "categorical"
  ),
  order_in_domain = c(1:9, 1:6, 1:2, 1),
  stringsAsFactors = FALSE
)

## -----------------------------------------------------------------------
## Biological samples candidates
## -----------------------------------------------------------------------
bio_candidates <- list(
  `hba1c.0.0` = list(label="HbA1c", type="continuous", order=1),
  `crp.0.0` = list(label="C-reactive protein (CRP)", type="continuous", order=2),
  `creatinine.0.0` = list(label="Creatinine", type="continuous", order=3),
  `glucose.0.0` = list(label="Glucose", type="continuous", order=4),
  `total_cholesterol.0.0` = list(label="Total cholesterol", type="continuous", order=5),
  `hdl_cholesterol.0.0` = list(label="HDL cholesterol", type="continuous", order=6),
  `ldl_direct.0.0` = list(label="LDL cholesterol (direct)", type="continuous", order=7),
  `triglycerides.0.0` = list(label="Triglycerides", type="continuous", order=8),
  `alt.0.0` = list(label="ALT", type="continuous", order=9),
  `ast.0.0` = list(label="AST", type="continuous", order=10),
  `albumin.0.0` = list(label="Albumin", type="continuous", order=11),
  `wbc.0.0` = list(label="White blood cell count", type="continuous", order=12),
  `platelets.0.0` = list(label="Platelet count", type="continuous", order=13),
  `haemoglobin.0.0` = list(label="Haemoglobin", type="continuous", order=14),
  
  # fallback
  hba1c = list(label="HbA1c", type="continuous", order=1),
  crp = list(label="C-reactive protein (CRP)", type="continuous", order=2),
  creatinine = list(label="Creatinine", type="continuous", order=3),
  glucose = list(label="Glucose", type="continuous", order=4),
  total_cholesterol = list(label="Total cholesterol", type="continuous", order=5),
  hdl_cholesterol = list(label="HDL cholesterol", type="continuous", order=6),
  ldl_direct = list(label="LDL cholesterol (direct)", type="continuous", order=7),
  triglycerides = list(label="Triglycerides", type="continuous", order=8),
  alt = list(label="ALT", type="continuous", order=9),
  ast = list(label="AST", type="continuous", order=10),
  albumin = list(label="Albumin", type="continuous", order=11),
  wbc = list(label="White blood cell count", type="continuous", order=12),
  platelets = list(label="Platelet count", type="continuous", order=13),
  haemoglobin = list(label="Haemoglobin", type="continuous", order=14)
)

## -----------------------------------------------------------------------
## Health & medical history domain
## -----------------------------------------------------------------------
health_candidates <- list(
  long_ill = list(label="Long-standing illness/disability", type="categorical", order=1),
  diabetes_doctor_diagnosed = list(label="Diabetes (doctor-diagnosed)", type="categorical", order=2),
  serious_dx = list(label="Serious medical diagnosis (self-reported)", type="categorical", order=3),
  age_high_bp_dx = list(label="Age at hypertension diagnosis, years", type="continuous", order=4),
  age_diabetes_dx = list(label="Age at diabetes diagnosis, years", type="continuous", order=5),
  angina_dx_age = list(label="Age at angina diagnosis, years", type="continuous", order=6),
  mi_dx_age = list(label="Age at myocardial infarction diagnosis, years", type="continuous", order=7),
  stroke_dx_age = list(label="Age at stroke diagnosis, years", type="continuous", order=8),
  dvt_dx_age = list(label="Age at DVT diagnosis, years", type="continuous", order=9),
  pe_dx_age  = list(label="Age at PE diagnosis, years", type="continuous", order=10),
  heart_problems_any = list(label="Any heart problems", type="categorical", order=11),
  meds_chol_bp_diabetes_any = list(label="Medication for cholesterol/BP/diabetes", type="categorical", order=12),
  meds_chol_bp_diabetes_hormone_any = list(label="Medication for cholesterol/BP/diabetes/hormone therapy", type="categorical", order=13)
)

## -----------------------------------------------------------------------
## Lifestyle, Environment, Early Life candidates
## -----------------------------------------------------------------------
lifestyle_candidates <- list(
  sleep_hrs = list(label="Sleep duration, hours", type="continuous", order=1),
  sleep_cat = list(label="Sleep duration category", type="categorical", order=2),
  insomnia = list(label="Sleeplessness / insomnia", type="categorical", order=3),
  snoring = list(label="Snoring", type="categorical", order=4),
  smoking_status = list(label="Smoking status", type="categorical", order=5),
  pack_yrs_smoking_prop = list(label="Pack-years of smoking", type="continuous", order=6),
  passive_cat = list(label="Passive smoking exposure", type="categorical", order=7),
  met_min_week_all = list(label="Total physical activity (MET-min/week)", type="continuous", order=8)
)

environment_candidates <- list(
  `no2_2010.0.0` = list(label="NO$_2$ (2010)", type="continuous", order=1),
  `pm2.5_2010.0.0` = list(label="PM$_{2.5}$ (2010)", type="continuous", order=2),
  `pm10_2010.0.0` = list(label="PM$_{10}$ (2010)", type="continuous", order=3),
  urban_rural_3lvl = list(label="Urban/rural classification", type="categorical", order=4),
  no2_2010 = list(label="NO$_2$ (2010)", type="continuous", order=1),
  `pm2.5_2010` = list(label="PM$_{2.5}$ (2010)", type="continuous", order=2),
  pm10_2010 = list(label="PM$_{10}$ (2010)", type="continuous", order=3)
)

early_life_candidates <- list(
  birth_weight = list(label="Birth weight, kg", type="continuous", order=1),
  breastfed = list(label="Breastfed as a baby", type="categorical", order=2),
  maternal_smoking_birth = list(label="Maternal smoking around birth", type="categorical", order=3),
  comparative_body_size_10 = list(label="Comparative body size at age 10", type="categorical", order=4)
)

## ---- build meta additions ----
meta_bio       <- build_meta("Biological samples", bio_candidates, names(dat))
meta_health    <- build_meta("Health & medical history", health_candidates, names(dat))
meta_lifestyle <- build_meta("Lifestyle", lifestyle_candidates, names(dat))
meta_env       <- build_meta("Environment", environment_candidates, names(dat))
meta_earlylife <- build_meta("Early life factors", early_life_candidates, names(dat))

meta <- rbind(
  meta_base[meta_base$var %in% names(dat), ],
  if (!is.null(meta_lifestyle)) meta_lifestyle,
  if (!is.null(meta_bio)) meta_bio,
  if (!is.null(meta_env)) meta_env,
  if (!is.null(meta_earlylife)) meta_earlylife,
  if (!is.null(meta_health)) meta_health
)

domain_order <- c(
  "Sociodemographic",
  "Lifestyle",
  "Physical measures",
  "Diet & alcohol",
  "Biological samples",
  "Environment",
  "Early life factors",
  "Health & medical history",
  "Outcome / follow-up"
)
meta$domain <- factor(meta$domain, levels = domain_order)

if ("dplyr" %in% loadedNamespaces()) {
  meta2 <- meta %>% dplyr::arrange(domain, order_in_domain)
} else {
  meta2 <- meta[order(meta$domain, meta$order_in_domain), ]
}

## -----------------------------------------------------------------------
## SMD + p-value helpers
## -----------------------------------------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

smd_cont <- function(x, g) {
  x <- suppressWarnings(as.numeric(x))
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])
  if (nlevels(g) != 2) return(NA_real_)
  x0 <- x[g == levels(g)[1]]
  x1 <- x[g == levels(g)[2]]
  if (length(x0) < 2 || length(x1) < 2) return(NA_real_)
  m0 <- mean(x0); m1 <- mean(x1)
  s0 <- sd(x0);  s1 <- sd(x1)
  sp <- sqrt(((length(x0)-1)*s0^2 + (length(x1)-1)*s1^2) / (length(x0)+length(x1)-2))
  if (is.na(sp) || sp == 0) return(NA_real_)
  (m1 - m0) / sp
}

p_cont <- function(x, g) {
  x <- suppressWarnings(as.numeric(x))
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])
  if (nlevels(g) != 2) return(NA_real_)
  x0 <- x[g == levels(g)[1]]
  x1 <- x[g == levels(g)[2]]
  if (length(x0) < 2 || length(x1) < 2) return(NA_real_)
  tryCatch(t.test(x1, x0)$p.value, error=function(e) NA_real_)
}

smd_cat_max <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- as.character(x[ok]); g <- droplevels(g[ok])
  if (nlevels(g) != 2) return(NA_real_)
  lvls <- sort(unique(x))
  if (length(lvls) < 2) return(NA_real_)
  g0 <- levels(g)[1]; g1 <- levels(g)[2]
  n0 <- sum(g == g0); n1 <- sum(g == g1)
  smds <- sapply(lvls, function(lv) {
    p0 <- mean(x[g == g0] == lv)
    p1 <- mean(x[g == g1] == lv)
    p_pool <- (p0*n0 + p1*n1)/(n0+n1)
    denom <- sqrt(p_pool * (1 - p_pool))
    if (denom == 0) return(0)
    (p1 - p0) / denom
  })
  max(abs(smds), na.rm=TRUE)
}

p_cat <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- as.character(x[ok]); g <- droplevels(g[ok])
  if (nlevels(g) != 2) return(NA_real_)
  tab <- table(x, g)
  tryCatch(chisq.test(tab)$p.value, error=function(e) NA_real_)
}

## -----------------------------------------------------------------------
## Table display helpers
## -----------------------------------------------------------------------
fmt_mean_sd <- function(x, digits=1) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  sprintf(paste0("%.",digits,"f (%.",digits,"f)"), mean(x), sd(x))
}

fmt_n_pct_one_level <- function(x, level, digits=1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  tab <- table(x)
  pct <- 100 * tab / sum(tab)
  if (!level %in% names(tab)) return("0 (0.0%)")
  sprintf("%d (%.1f%%)", as.integer(tab[level]), as.numeric(pct[level]))
}

make_cont_row <- function(df, var, label) {
  smd <- smd_cont(df[[var]], df$event)
  pv  <- p_cont(df[[var]], df$event)
  c(
    Characteristic = label,
    Overall  = fmt_mean_sd(df[[var]]),
    `No event` = fmt_mean_sd(df[df$event=="No event", var]),
    Event    = fmt_mean_sd(df[df$event=="Event", var]),
    SMD      = ifelse(is.na(smd), "", sprintf("%.3f", smd)),
    `P value`= fmt_p(pv)
  )
}

make_cat_block <- function(df, var, label) {
  v_all <- as.character(df[[var]])
  lvls <- sort(unique(v_all[!is.na(v_all)]))
  
  smd <- smd_cat_max(v_all, df$event)
  pv  <- p_cat(v_all, df$event)
  
  block <- list(c(
    Characteristic = label,
    Overall = "", `No event` = "", Event = "",
    SMD = ifelse(is.na(smd), "", sprintf("%.3f", smd)),
    `P value` = fmt_p(pv)
  ))
  
  for (lvl in lvls) {
    block[[length(block)+1]] <- c(
      Characteristic = paste0("  ", lvl),
      Overall   = fmt_n_pct_one_level(v_all, lvl),
      `No event`= fmt_n_pct_one_level(as.character(df[df$event=="No event", var]), lvl),
      Event     = fmt_n_pct_one_level(as.character(df[df$event=="Event", var]), lvl),
      SMD = "", `P value` = ""
    )
  }
  do.call(rbind, block)
}

## ---- build BIG table ----
cat("Building the Table 1 dataframe...\n")
out <- list()

out[[length(out)+1]] <- c(
  Characteristic = "N",
  Overall  = sum(!is.na(dat$event)),
  `No event` = sum(dat$event=="No event", na.rm=TRUE),
  Event    = sum(dat$event=="Event", na.rm=TRUE),
  SMD = "", `P value` = ""
)

current_domain <- NULL
for (i in seq_len(nrow(meta2))) {
  dom <- as.character(meta2$domain[i])
  var <- meta2$var[i]
  lab <- meta2$label[i]
  typ <- meta2$type[i]
  
  if (!identical(dom, current_domain)) {
    out[[length(out)+1]] <- c(
      Characteristic = paste0("<b>", dom, "</b>"),
      Overall = "", `No event` = "", Event = "",
      SMD = "", `P value` = ""
    )
    current_domain <- dom
  }
  
  if (typ == "categorical") out[[length(out)+1]] <- make_cat_block(dat, var, lab)
  else out[[length(out)+1]] <- rbind(make_cont_row(dat, var, lab))
}

table1 <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)

## ---- output and export ----
cat("Processing Outputs (HTML Viewer + CSV file)...\n")

# 1. Output directory
out_dir <- here("02_results", "01_tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 2. Render and Display HTML in RStudio Viewer
html <- knitr::kable(
  table1,
  format = "html",
  escape = FALSE,
  caption = "Table 1. Baseline characteristics of the analysis cohort by event status."
)

html_file <- file.path(out_dir, "Table1_event_status.html")
writeLines(
  c("<html><head><meta charset='utf-8'><style>body{font-family: sans-serif;}</style></head><body>", html, "</body></html>"),
  html_file
)

if (interactive()) {
  tryCatch(
    rstudioapi::viewer(html_file),
    error = function(e) browseURL(html_file)
  )
}

# 3. Clean up formatting and Save as CSV
table1_csv <- table1
table1_csv$Characteristic <- gsub("<b>|</b>", "", table1_csv$Characteristic) # Clean HTML tags

csv_file <- file.path(out_dir, "Table1_event_status.csv")
write.csv(table1_csv, csv_file, row.names = FALSE)

cat("=======================================================\n")
cat("SUCCESS!\n")
cat("HTML displayed in Viewer & saved to:\n", html_file, "\n")
cat("Clean CSV saved to:\n", csv_file, "\n")
cat("=======================================================\n")