# ==============================================================================
# Script: Imputed Cohort Table 1 
# ==============================================================================
rm(list = ls())

library(here)
library(knitr)
suppressWarnings({
  if (requireNamespace("dplyr", quietly = TRUE)) library(dplyr)
})

## =============================================================================
## 0) Paths 
## =============================================================================
cat("Loading imputed data...\n")
train_path <- here("00_data", "01_processed", "03_imputation", "df_train_imputed.rds")
test_path  <- here("00_data", "01_processed", "03_imputation", "df_test_imputed.rds")

train <- readRDS(train_path)
test  <- readRDS(test_path)

## Combine (Table 1 uses full cohort)
dat <- rbind(train, test)

## =============================================================================
## 1) Event grouping
## =============================================================================
cat("Processing event variable...\n")
if ("event" %in% names(dat)) {
  dat$event <- as.character(dat$event)
  dat$event <- ifelse(dat$event %in% c("1","Event","event","Incident CVD"), "Event",
                      ifelse(dat$event %in% c("0","No event","no_event","No incident CVD"), "No event", dat$event))
  dat$event <- factor(dat$event, levels=c("No event","Event"))
} else if ("has_cvd" %in% names(dat)) {
  has_cvd_num <- suppressWarnings(as.integer(as.character(dat$has_cvd)))
  has_cvd_num[is.na(has_cvd_num) & is.logical(dat$has_cvd)] <- as.integer(dat$has_cvd[is.na(has_cvd_num) & is.logical(dat$has_cvd)])
  dat$event <- ifelse(has_cvd_num == 1, "Event",
                      ifelse(has_cvd_num == 0, "No event", NA))
  dat$event <- factor(dat$event, levels=c("No event","Event"))
} else {
  stop("Need `event` or `has_cvd` in imputed dataset.")
}

## =============================================================================
## 2) Helpers: formatting, p, SMD
## =============================================================================
fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_mean_sd <- function(x, digits=1) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  sprintf(paste0("%.",digits,"f (%.",digits,"f)"), mean(x), sd(x))
}

fmt_n_pct_one_level <- function(x, level) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  tab <- table(x)
  pct <- 100 * tab / sum(tab)
  if (!level %in% names(tab)) return("0 (0.0%)")
  sprintf("%d (%.1f%%)", as.integer(tab[level]), as.numeric(pct[level]))
}

smd_cont <- function(x, g) {
  x <- suppressWarnings(as.numeric(x))
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])
  if (nlevels(g) != 2) return(NA_real_)
  x0 <- x[g == levels(g)[1]]
  x1 <- x[g == levels(g)[2]]
  if (length(x0) < 2 || length(x1) < 2) return(NA_real_)
  s0 <- sd(x0); s1 <- sd(x1)
  sp <- sqrt(((length(x0)-1)*s0^2 + (length(x1)-1)*s1^2) / (length(x0)+length(x1)-2))
  if (is.na(sp) || sp == 0) return(NA_real_)
  (mean(x1) - mean(x0)) / sp
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

## =============================================================================
## 3) Define FULL TABLE domains
## =============================================================================
cat("Configuring variables...\n")
meta_base <- data.frame(
  domain = c(
    rep("Sociodemographic", 4),
    rep("Lifestyle", 3),
    rep("Physical measures", 4),
    rep("Diet & alcohol", 2)
  ),
  var = c(
    "age","sex","ethnicity","imd_q_uk",
    "smoking_status","alcohol_status","met_min_week_all",
    "sbp","dbp","whr","body_fat_mass",
    "diet_score","alcohol_consumption"
  ),
  label = c(
    "Age, years","Sex","Ethnicity","Deprivation quintile",
    "Smoking status","Alcohol status","Physical activity (MET-min/week)",
    "Systolic BP, mmHg","Diastolic BP, mmHg","Waist-to-hip ratio (WHR)","Body fat mass, kg",
    "Diet score (0–8)","Alcohol consumption category"
  ),
  type = c(
    "continuous","categorical","categorical","categorical",
    "categorical","categorical","continuous",
    "continuous","continuous","continuous","continuous",
    "continuous","categorical"
  ),
  order_in_domain = c(1:4, 1:3, 1:4, 1:2),
  stringsAsFactors = FALSE
)

bio_candidates <- list(
  hba1c = list(label="HbA1c", type="continuous", order=1),
  crp = list(label="C-reactive protein (CRP)", type="continuous", order=2),
  creatinine = list(label="Creatinine", type="continuous", order=3),
  total_cholesterol = list(label="Total cholesterol", type="continuous", order=4),
  hdl_cholesterol = list(label="HDL cholesterol", type="continuous", order=5),
  ldl_direct = list(label="LDL cholesterol (direct)", type="continuous", order=6),
  triglycerides = list(label="Triglycerides", type="continuous", order=7),
  alt = list(label="ALT", type="continuous", order=8),
  ast = list(label="AST", type="continuous", order=9),
  albumin = list(label="Albumin", type="continuous", order=10),
  wbc = list(label="White blood cell count", type="continuous", order=11),
  platelets = list(label="Platelet count", type="continuous", order=12)
)

early_life_candidates <- list(
  birth_weight = list(label="Birth weight, kg", type="continuous", order=1),
  breastfed = list(label="Breastfed as a baby", type="categorical", order=2),
  maternal_smoking_birth = list(label="Maternal smoking around birth", type="categorical", order=3),
  comparative_body_size_10 = list(label="Comparative body size at age 10", type="categorical", order=4)
)

health_candidates <- list(
  long_ill = list(label="Long-standing illness/disability", type="categorical", order=1),
  diabetes_doctor_diagnosed = list(label="Diabetes (doctor-diagnosed)", type="categorical", order=2),
  serious_dx = list(label="Serious medical diagnosis (self-reported)", type="categorical", order=3),
  heart_problems_any = list(label="Any heart problems", type="categorical", order=4),
  meds_chol_bp_diabetes_collapsed = list(label="Medication for cholesterol/BP/diabetes", type="categorical", order=5),
  meds_chol_bp_diabetes_hormone_collapsed = list(label="Medication for cholesterol/BP/diabetes/hormone therapy", type="categorical", order=6)
)

meta_bio       <- build_meta("Biological samples", bio_candidates, names(dat))
meta_earlylife <- build_meta("Early life factors", early_life_candidates, names(dat))
meta_health    <- build_meta("Health & medical history", health_candidates, names(dat))

meta_full <- rbind(
  meta_base[meta_base$var %in% names(dat), ],
  if (!is.null(meta_bio)) meta_bio,
  if (!is.null(meta_earlylife)) meta_earlylife,
  if (!is.null(meta_health)) meta_health
)

domain_order <- c(
  "Sociodemographic","Lifestyle","Physical measures","Diet & alcohol",
  "Biological samples","Early life factors","Health & medical history"
)
meta_full$domain <- factor(meta_full$domain, levels=domain_order)
meta_full <- meta_full[order(meta_full$domain, meta_full$order_in_domain), ]

## =============================================================================
## 4) Build table function
## =============================================================================
build_table <- function(df, meta_df, caption, out_html, out_csv) {
  out <- list()
  out[[length(out)+1]] <- c(
    Characteristic="N",
    Overall=sum(!is.na(df$event)),
    `No event`=sum(df$event=="No event", na.rm=TRUE),
    Event=sum(df$event=="Event", na.rm=TRUE),
    SMD="", `P value`=""
  )
  
  current_domain <- NULL
  for (i in seq_len(nrow(meta_df))) {
    dom <- as.character(meta_df$domain[i])
    var <- meta_df$var[i]
    lab <- meta_df$label[i]
    typ <- meta_df$type[i]
    
    if (!identical(dom, current_domain)) {
      out[[length(out)+1]] <- c(Characteristic=paste0("<b>", dom, "</b>"),
                                Overall="", `No event`="", Event="",
                                SMD="", `P value`="")
      current_domain <- dom
    }
    
    if (typ=="categorical") {
      out[[length(out)+1]] <- make_cat_block(df, var, lab)
    } else {
      out[[length(out)+1]] <- rbind(make_cont_row(df, var, lab))
    }
  }
  
  tab <- as.data.frame(do.call(rbind, out), stringsAsFactors=FALSE)
  
  # -------------------------
  # 1. Export HTML
  # -------------------------
  html <- knitr::kable(tab, format="html", escape=FALSE, caption=caption)
  writeLines(c("<html><head><meta charset='utf-8'><style>body{font-family: sans-serif;}</style></head><body>", html, "</body></html>"), out_html)
  
  # -------------------------
  # 2. Export Clean CSV
  # -------------------------
  tab_csv <- tab
  tab_csv$Characteristic <- gsub("<b>|</b>", "", tab_csv$Characteristic) 
  write.csv(tab_csv, out_csv, row.names = FALSE)
  
  invisible(tab)
}

## =============================================================================
## 5) Execute the function!
## =============================================================================
cat("Calculating Table 1. This may take a few seconds...\n")
out_dir <- here("02_results", "01_tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

html_path <- file.path(out_dir, "Table1_Imputed_Full_Cohort.html")
csv_path  <- file.path(out_dir, "Table1_Imputed_Full_Cohort.csv")

final_table <- build_table(
  df = dat, 
  meta_df = meta_full, 
  caption = "Baseline characteristics of the complete imputed cohort by event status.", 
  out_html = html_path, 
  out_csv = csv_path
)

if (interactive()) {
  tryCatch(
    rstudioapi::viewer(html_path),
    error = function(e) browseURL(html_path)
  )
}

cat("\n=======================================================\n")
cat("SUCCESS!\n")
cat("HTML displayed in Viewer & saved to:\n-", html_path, "\n\n")
cat("Clean CSV saved to:\n-", csv_path, "\n")
cat("=======================================================\n")