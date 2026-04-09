rm(list = ls())

# Import packages ---------------------------------------------------------

library(here)
library(dplyr)
library(FactoMineR)
library(corrplot)
library(ggplot2)
library(scales)

# Activate config for anthropometry dataset for exploratory
SCENARIO <- "sensitivity_anthro_confounder"

source(here("01_scripts", "00_variable_dictionary.R"))
source(here("01_scripts", "00_variable_config.R"))

# Import data -------------------------------------------------------------

# Use sensitivity imputation which includes all variables (including anthropometrics)
# for descriptive analysis 
# main imputation excludes anthropometrics

df <- readRDS(file.path(data_input_dir, full_imputed_file))

colnames(df)

# Output path -------------------------------------------------------------

fig_dir <- famd_fig_dir
tbl_dir <- famd_tbl_dir

# Remove ID / outcome / dates ------------------------------------------------

df_analysis <- df %>%
  select(
    -any_of(c(
      "eid",
      "has_cvd",
      "recruitment_date",
      "death_date",
      "first_cvd_date",
      "sleep_cat",
      "passive_cat"
      
      # "long_ill",
      # "diabetes_doctor_diagnosed",
      # "serious_dx",
      # "age_high_bp_dx",
      # "age_diabetes_dx",
      # "angina_dx_age",
      # "mi_dx_age",
      # "dvt_dx_age",
      # "pe_dx_age",
      # "stroke_dx_age",
      # "heart_problems_collapsed",
      # "meds_chol_bp_diabetes_hormone_collapsed",
      # "meds_chol_bp_diabetes_collapsed"
    ))
  ) %>%
  mutate(across(where(is.character), as.factor))

# ----------------------------------------------------------
# Run FAMD -------------------------------------------------
# ----------------------------------------------------------

set.seed(1)
famd_res <- FAMD(df_analysis, ncp = 200, graph = FALSE)

# FAMD Plots ----------------------------------------------

num_contrib <- famd_res$quanti.var$contrib
cat_contrib <- famd_res$quali.var$contrib

num_contrib_df <- as.data.frame(num_contrib)
num_contrib_df$variable <- rownames(num_contrib_df)
num_contrib_df$display_name <- get_display_name(num_contrib_df$variable)

write.csv(
  num_contrib_df,
  file.path(tbl_dir, "famd_numeric_contributions.csv"),
  row.names = FALSE
)

cat_contrib_df <- as.data.frame(cat_contrib)
cat_contrib_df$variable <- rownames(cat_contrib_df)
cat_contrib_df$display_name <- get_display_name(cat_contrib_df$variable)

write.csv(
  cat_contrib_df,
  file.path(tbl_dir, "famd_categorical_contributions.csv"),
  row.names = FALSE
)

famd_eig <- as.data.frame(famd_res$eig)
colnames(famd_eig) <- c("eigenvalue", "percent_variance", "cumulative_percent")

famd_eig$dimension <- seq_len(nrow(famd_eig))

write.csv(
  famd_eig,
  file.path(tbl_dir, "famd_eigenvalues.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# Scree Plots
# ----------------------------------------------------------

# Standard Scree Plot -----------------------------------------------

eig_pct <- famd_res$eig[, 2]

png(file.path(fig_dir, "famd_scree_plot.png"), width = 1800, height = 1400, res = 220)

plot(
  eig_pct,
  type = "b",
  pch = 16,
  lwd = 2,
  col = "#2c3e50",
  xlab = "Component",
  ylab = "% variance explained",
  main = "FAMD Scree Plot"
)

dev.off()

famd_res$eig
famd_res$eig[, 3]
sum(famd_res$eig[,2])
tail(famd_res$eig[,3], 1)

colnames(df_analysis)

# Cumulative Scree Plot ------------------------------------------------

cum_var <- famd_res$eig[, 3]

png(file.path(fig_dir, "famd_cumulative_scree_plot.png"), width = 1800, height = 1400, res = 220)

plot(
  cum_var,
  type = "b",
  pch = 16,
  lwd = 2,
  col = "#2c3e50",
  xlab = "FAMD Dimension",
  ylab = "Cumulative variance explained (%)",
  main = "FAMD Cumulative Scree Plot",
  ylim = c(0, 100),
  yaxt = "n"
)

axis(
  2,
  at = seq(0, 100, by = 20),
  labels = seq(0, 100, by = 20),
  las = 1
)

dev.off()

# Interpreting components of the first 7 most informative PCs (Principal Components) ------------------------------------------------

num_contrib <- famd_res$quanti.var$contrib
cat_contrib <- famd_res$quali.var$contrib

sort(num_contrib[,7], decreasing = TRUE)[1:10]
sort(cat_contrib[,7], decreasing = TRUE)[1:10]

# Create plots for a PC

plot_top_contrib <- function(contrib_matrix, pc_num, title_prefix) {
  
  pc_vals <- contrib_matrix[, pc_num]
  
  top_vars <- sort(pc_vals, decreasing = TRUE)[1:6]
  top_labels <- get_display_name(names(top_vars))
  
  par(mar = c(8, 5, 5, 2))
  
  barplot(
    top_vars,
    names.arg = top_labels,
    col = "steelblue",
    border = NA,
    las = 1,
    cex.names = 0.9,
    cex.axis = 1.1,
    cex.lab = 1.3,
    cex.main = 1.5,
    main = paste("Top Numeric Contributors to PC", pc_num),
    xlab = "Variable",
    ylab = "Contribution (%)"
  )
}

# Automatically loop through first 7 PCs to plot

for (i in 1:7) {
  png(file.path(fig_dir, paste0("famd_numeric_contributors_pc", i, ".png")),
      width = 1800, height = 1400, res = 220)
  plot_top_contrib(num_contrib, i, "Numeric Contributors")
  dev.off()
}

# ----------------------------------------------------------
# Scatter plot of first 2 PCs
# ----------------------------------------------------------

# Scatter Plot (coloured by CVD) ------------------------------------------------

scores <- famd_res$ind$coord

col_outcome <- ifelse(
  df$has_cvd == 1,
  adjustcolor("orange", alpha.f = 0.7),
  adjustcolor("#2c3e50", alpha.f = 0.5)
)

png(file.path(fig_dir, "famd_individuals_by_cvd.png"), width = 1800, height = 1400, res = 220)

par(col.lab = "#2c3e50", col.axis = "#2c3e50", col.main = "#2c3e50")

plot(
  scores[,1], scores[,2],
  col = col_outcome,
  pch = 16,
  cex = 0.5,
  xlab = "Dimension 1",
  ylab = "Dimension 2",
  main = "FAMD: Individuals (coloured by CVD)",
  ylim = c(-6, 8),
  yaxt = "n"
)

axis(2, at = seq(-6, 8, by = 2), las = 1, col.axis = "#2c3e50")

legend(
  "topright",
  legend = c("No CVD", "CVD"),
  col = c("#2c3e50", "#3a6ea5"),
  pch = 16,
  bty = "n",
  text.col = "#2c3e50"
)

dev.off()

# Scatter Plot (coloured by sex) ---------------------------

if ("sex" %in% names(df)) {
  
  sex_fac <- as.factor(df$sex)
  sex_levels <- levels(sex_fac)
  
  base_cols <- c(
    adjustcolor("#3a6ea5", alpha.f = 0.5),
    adjustcolor("#c0392b", alpha.f = 0.5),
    adjustcolor("#16a085", alpha.f = 0.5),
    adjustcolor("#8e44ad", alpha.f = 0.5),
    adjustcolor("#f39c12", alpha.f = 0.5)
  )
  
  cols <- setNames(base_cols[seq_along(sex_levels)], sex_levels)
  point_cols <- cols[as.character(sex_fac)]
  point_cols[is.na(point_cols)] <- adjustcolor("grey50", alpha.f = 0.35)
  
  png(file.path(fig_dir, "famd_individuals_by_sex.png"), width = 1800, height = 1400, res = 220)
  
  par(col.lab = "#2c3e50", col.axis = "#2c3e50", col.main = "#2c3e50")
  
  plot(
    scores[,1], scores[,2],
    col = point_cols,
    pch = 16,
    cex = 0.5,
    xlab = "Dimension 1",
    ylab = "Dimension 2",
    main = "FAMD: Individuals (coloured by Sex)",
    ylim = c(-6, 8),
    yaxt = "n"
  )
  
  axis(2, at = seq(-6, 8, by = 2), las = 1, col.axis = "#2c3e50")
  
  legend(
    "topright",
    legend = sex_levels,
    col = cols[sex_levels],
    pch = 16,
    bty = "n",
    text.col = "#2c3e50"
  )
  
  dev.off()
}

# ----------------------------------------------------------
# Calculating Correlations
# ----------------------------------------------------------

# ----------------------------------------------------------
# ----------------------------------------------------------
# Correlation Heatmap
# ----------------------------------------------------------

num_data <- df_analysis %>% select(where(is.numeric))

if (ncol(num_data) >= 2) {
  
  # Ensures we have valid number of columns
  valid_num_cols <- names(num_data)[vapply(num_data, function(x) sum(!is.na(x)) >= 2, logical(1))]
  num_data <- num_data[, valid_num_cols, drop = FALSE]
  
  non_constant_cols <- names(num_data)[vapply(num_data, function(x) {
    x_non_na <- x[!is.na(x)]
    length(unique(x_non_na)) > 1
  }, logical(1))]
  num_data <- num_data[, non_constant_cols, drop = FALSE]
  
  # Calculates pairwise numeric correlations with Pearson
  if (ncol(num_data) >= 2) {
    
    cor_mat <- cor(num_data, use = "pairwise.complete.obs", method = "pearson")
    cor_mat_raw <- cor_mat
    
    # Output the correlation matrix
    write.csv(
      cor_mat,
      file.path(tbl_dir, "correlation_matrix_numeric.csv"),
      row.names = TRUE
    )
    
    colnames(cor_mat) <- get_display_name(colnames(cor_mat))
    rownames(cor_mat) <- get_display_name(rownames(cor_mat))
    
    # Output the correlation heatmap
    png(file.path(fig_dir, "correlation_heatmap.png"), width = 1800, height = 2400, res = 220)
    
    par(mar = c(1, 1, 8, 1), oma = c(0, 0, 3.5, 0))
    
    corrplot(
      cor_mat,
      method = "color",
      type = "upper",
      main = "",
      tl.cex = 0.62,
      tl.col = "black",
      order = "hclust",
      col.lim = c(-1, 1)
    )
    
    mtext(
      "Numeric Variable Correlations",
      side = 3,
      outer = TRUE,
      line = 0.2,
      cex = 1.1,
      font = 2
    )
    
    dev.off()
    
    cor_mat_raw[lower.tri(cor_mat_raw, diag = TRUE)] <- NA
    cor_pairs <- as.data.frame(as.table(cor_mat_raw))
    cor_pairs <- cor_pairs[!is.na(cor_pairs$Freq), , drop = FALSE]
    
    # Output all numeric correlations into a .csv
    if (nrow(cor_pairs) > 0) {
      
      cor_pairs$abs_cor <- abs(cor_pairs$Freq)
      strongest <- cor_pairs[order(-cor_pairs$abs_cor), , drop = FALSE]
      
      all_cor <- strongest
      all_cor$Var1 <- get_display_name(all_cor$Var1)
      all_cor$Var2 <- get_display_name(all_cor$Var2)
      
      write.csv(
        all_cor,
        file.path(tbl_dir, "all_numeric_correlations.csv"),
        row.names = FALSE
      )
      
      top_n <- min(15, nrow(strongest))
      top15_cor <- head(strongest, top_n)
      
      top15_cor$Var1 <- get_display_name(top15_cor$Var1)
      top15_cor$Var2 <- get_display_name(top15_cor$Var2)
      
      # Outputs just the top 15 numeric correlations
      write.csv(
        top15_cor,
        file.path(tbl_dir, "top15_numeric_correlations.csv"),
        row.names = FALSE
      )
      
      top_vars <- unique(c(
        as.character(head(strongest, top_n)$Var1),
        as.character(head(strongest, top_n)$Var2)
      ))
      
      if (length(top_vars) >= 2) {
        cor_mat_top15 <- cor(
          num_data[, top_vars, drop = FALSE],
          use = "pairwise.complete.obs",
          method = "pearson"
        )
        
        colnames(cor_mat_top15) <- get_display_name(colnames(cor_mat_top15))
        rownames(cor_mat_top15) <- get_display_name(rownames(cor_mat_top15))
        
        # Outputs top 15 numeric correlations as a heatmap
        png(file.path(fig_dir, "correlation_heatmap_top15.png"), width = 1800, height = 2000, res = 220)
        
        par(mar = c(1, 1, 8, 1), oma = c(0, 0, 3.5, 0))
        
        corrplot(
          cor_mat_top15,
          method = "color",
          type = "upper",
          main = "",
          tl.cex = 0.85,
          tl.col = "black",
          order = "hclust",
          col.lim = c(-1, 1)
        )
        
        mtext(
          paste("Top", top_n, "Strongest Pairwise Correlations"),
          side = 3,
          outer = TRUE,
          line = 0.2,
          cex = 1.1,
          font = 2
        )
        
        dev.off()
      }
    }
  }
}

# Helper functions -------------------------------------------------------

# Function for Cramer's V (correlation between categoricals)
cramers_v <- function(x, y) {
  if (all(is.na(x)) || all(is.na(y))) return(NA_real_)
  
  # Just some safety checks below
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  
  x <- as.factor(x)
  y <- as.factor(y)
  
  if (nlevels(x) < 2 || nlevels(y) < 2) return(NA_real_)
  
  tab <- table(x, y)
  if (any(dim(tab) < 2)) return(NA_real_)
  
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  n <- sum(tab)
  
  if (n == 0) return(NA_real_)
  
  phi2 <- unname(chi$statistic) / n
  r <- nrow(tab)
  k <- ncol(tab)
  
  denom <- min(k - 1, r - 1)
  if (denom <= 0) return(NA_real_)
  
  sqrt(phi2 / denom)
}

# Calculates correlation ratio (eta) for categorical vs numeric
correlation_ratio <- function(categories, values) {
  if (all(is.na(categories)) || all(is.na(values))) return(NA_real_)
  
  keep <- !is.na(categories) & !is.na(values)
  categories <- categories[keep]
  values <- values[keep]
  
  if (length(values) == 0) return(NA_real_)
  
  categories <- as.factor(categories)
  if (nlevels(categories) < 2) return(NA_real_)
  if (length(unique(values)) < 2) return(NA_real_)
  
  grand_mean <- mean(values)
  
  group_stats <- split(values, categories)
  n_i <- vapply(group_stats, length, numeric(1))
  mean_i <- vapply(group_stats, mean, numeric(1))
  
  ss_between <- sum(n_i * (mean_i - grand_mean)^2)
  ss_total <- sum((values - grand_mean)^2)
  
  if (ss_total == 0) return(NA_real_)
  
  sqrt(ss_between / ss_total)
}

# Wrapper function to generate a formatted correlation heatmap for a matrix and save it as a PNG
safe_corrplot <- function(mat, file_path, title_text, width = 1800, height = 2400,
                          tl_cex = 0.8, mar_top = 8, col_lim = c(0, 1)) {
  
  if (is.null(mat)) return(invisible(NULL))
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (nrow(mat) < 2 || ncol(mat) < 2) return(invisible(NULL))
  
  png(file_path, width = width, height = height, res = 220)
  
  par(mar = c(1, 1, mar_top, 1), oma = c(0, 0, 3.5, 0))
  
  corrplot(
    mat,
    method = "color",
    type = "upper",
    main = "",
    tl.cex = tl_cex,
    tl.col = "black",
    order = "hclust",
    col.lim = col_lim
  )
  
  mtext(
    title_text,
    side = 3,
    outer = TRUE,
    line = -4,
    cex = 1.1,
    font = 2
  )
  
  dev.off()
}

# Define variable sets ---------------------------------------------------

num_vars_all <- names(df_analysis)[vapply(df_analysis, is.numeric, logical(1))]
cat_vars_all <- names(df_analysis)[vapply(df_analysis, function(x) is.factor(x) || is.character(x), logical(1))]

# Convert character categorical vars to factor for association work
df_assoc <- df_analysis %>%
  mutate(across(any_of(cat_vars_all), ~ as.factor(.x)))

# Remove unusable numeric vars
if (length(num_vars_all) > 0) {
  num_vars_all <- num_vars_all[vapply(df_assoc[num_vars_all], function(x) {
    x_non_na <- x[!is.na(x)]
    length(x_non_na) >= 2 && length(unique(x_non_na)) > 1
  }, logical(1))]
}

# Remove unusable categorical vars
if (length(cat_vars_all) > 0) {
  cat_vars_all <- cat_vars_all[vapply(df_assoc[cat_vars_all], function(x) {
    x <- as.factor(x)
    x_non_na <- x[!is.na(x)]
    length(x_non_na) >= 2 && nlevels(droplevels(x_non_na)) > 1
  }, logical(1))]
}

# ----------------------------------------------------------
# 2) Categorical vs categorical: Cramer's V
# ----------------------------------------------------------

if (length(cat_vars_all) >= 2) {
  
  cat_cat_mat <- matrix(
    NA_real_,
    nrow = length(cat_vars_all),
    ncol = length(cat_vars_all),
    dimnames = list(cat_vars_all, cat_vars_all)
  )
  
  for (i in seq_along(cat_vars_all)) {
    for (j in seq_along(cat_vars_all)) {
      if (i == j) {
        cat_cat_mat[i, j] <- 1
      } else {
        cat_cat_mat[i, j] <- cramers_v(df_assoc[[cat_vars_all[i]]], df_assoc[[cat_vars_all[j]]])
      }
    }
  }
  
  # Create the association matrix calculated with Cramer's V (Categoricals)
  write.csv(
    cat_cat_mat,
    file.path(tbl_dir, "association_matrix_categorical_cramers_v.csv"),
    row.names = TRUE
  )
  
  cat_cat_raw <- cat_cat_mat
  colnames(cat_cat_mat) <- get_display_name(colnames(cat_cat_mat))
  rownames(cat_cat_mat) <- get_display_name(rownames(cat_cat_mat))
  
  safe_corrplot(
    mat = cat_cat_mat,
    file_path = file.path(fig_dir, "association_heatmap_categorical_cramers_v.png"),
    title_text = "Association Structure of Categorical Variables (Cramer's V)",
    width = 1800,
    height = 2400,
    tl_cex = 0.72,
    mar_top = 8,
    col_lim = c(0, 1)
  )
  
  cat_cat_raw[lower.tri(cat_cat_raw, diag = TRUE)] <- NA
  cat_cat_pairs <- as.data.frame(as.table(cat_cat_raw))
  cat_cat_pairs <- cat_cat_pairs[!is.na(cat_cat_pairs$Freq), , drop = FALSE]
  
  if (nrow(cat_cat_pairs) > 0) {
    cat_cat_pairs <- cat_cat_pairs[order(-cat_cat_pairs$Freq), , drop = FALSE]
    top_n_cat <- min(15, nrow(cat_cat_pairs))
    top_cat_cat <- head(cat_cat_pairs, top_n_cat)
    top_cat_cat$Var1 <- get_display_name(top_cat_cat$Var1)
    top_cat_cat$Var2 <- get_display_name(top_cat_cat$Var2)
    
    write.csv(
      top_cat_cat,
      file.path(tbl_dir, "top15_categorical_associations_cramers_v.csv"),
      row.names = FALSE
    )
  }
}

# ----------------------------------------------------------
# 3) Numeric vs categorical: correlation ratio (eta)
# ----------------------------------------------------------

if (length(num_vars_all) >= 1 && length(cat_vars_all) >= 1) {
  
  num_cat_mat <- matrix(
    NA_real_,
    nrow = length(num_vars_all),
    ncol = length(cat_vars_all),
    dimnames = list(num_vars_all, cat_vars_all)
  )
  
  for (i in seq_along(num_vars_all)) {
    for (j in seq_along(cat_vars_all)) {
      num_cat_mat[i, j] <- correlation_ratio(
        categories = df_assoc[[cat_vars_all[j]]],
        values = df_assoc[[num_vars_all[i]]]
      )
    }
  }
  
  # Create the .csv with Eta for all associations (Numeric vs Categorical)
  write.csv(
    num_cat_mat,
    file.path(tbl_dir, "association_matrix_numeric_categorical_eta.csv"),
    row.names = TRUE
  )
  
  num_cat_long <- as.data.frame(as.table(num_cat_mat))
  num_cat_long <- num_cat_long[!is.na(num_cat_long$Freq), , drop = FALSE]
  
  if (nrow(num_cat_long) > 0) {
    num_cat_long <- num_cat_long[order(-num_cat_long$Freq), , drop = FALSE]
    
    top_n_numcat <- min(20, nrow(num_cat_long))
    top_num_cat <- head(num_cat_long, top_n_numcat)
    
    top_num_vars <- unique(as.character(top_num_cat$Var1))
    top_cat_vars <- unique(as.character(top_num_cat$Var2))
    
    num_cat_top <- num_cat_mat[top_num_vars, top_cat_vars, drop = FALSE]
    
    top_num_cat_export <- top_num_cat
    top_num_cat_export$Var1 <- get_display_name(top_num_cat_export$Var1)
    top_num_cat_export$Var2 <- get_display_name(top_num_cat_export$Var2)
    
    # Create .csv for top 20 Eta correlations
    write.csv(
      top_num_cat_export,
      file.path(tbl_dir, "top20_numeric_categorical_associations_eta.csv"),
      row.names = FALSE
    )
    
    colnames(num_cat_top) <- get_display_name(colnames(num_cat_top))
    rownames(num_cat_top) <- get_display_name(rownames(num_cat_top))
    
    # Create the association matrix calculated with Eta (Numeric vs Categorical)
    if (nrow(num_cat_top) >= 1 && ncol(num_cat_top) >= 1) {
      png(file.path(fig_dir, "association_heatmap_numeric_categorical_eta.png"),
          width = 1800, height = 2200, res = 220)
      
      par(mar = c(8, 8, 8, 2), oma = c(0, 0, 3.5, 0))
      
      corrplot(
        num_cat_top,
        is.corr = FALSE,
        method = "color",
        main = "",
        tl.cex = 0.78,
        tl.col = "black",
        col.lim = c(0, 1)
      )
      
      mtext(
        "Associations Between Numeric and Categorical Variables (Eta)",
        side = 3,
        outer = TRUE,
        line = 0.2,
        cex = 1.05,
        font = 2
      )
      
      dev.off()
    }
  }
}

# ----------------------------------------------------------
# 4) Overall mixed association heatmap
# ----------------------------------------------------------

mixed_vars <- c(num_vars_all, cat_vars_all)

if (length(mixed_vars) >= 2) {
  
  mixed_mat <- matrix(
    NA_real_,
    nrow = length(mixed_vars),
    ncol = length(mixed_vars),
    dimnames = list(mixed_vars, mixed_vars)
  )
  
  # Nested loop through all mixed variables
  # To create a full correlation heatmap with all types of variables
  for (i in seq_along(mixed_vars)) {
    for (j in seq_along(mixed_vars)) {
      
      var_i <- mixed_vars[i]
      var_j <- mixed_vars[j]
      
      x <- df_assoc[[var_i]]
      y <- df_assoc[[var_j]]
      
      x_is_num <- is.numeric(x)
      y_is_num <- is.numeric(y)
      x_is_cat <- is.factor(x) || is.character(x)
      y_is_cat <- is.factor(y) || is.character(y)
      
      if (i == j) {
        mixed_mat[i, j] <- 1
        
      } else if (x_is_num && y_is_num) {
        keep <- !is.na(x) & !is.na(y)
        x2 <- x[keep]
        y2 <- y[keep]
        
        if (length(x2) >= 2 && length(y2) >= 2 &&
            length(unique(x2)) > 1 && length(unique(y2)) > 1) {
          mixed_mat[i, j] <- abs(cor(x2, y2, method = "pearson"))
        }
        
      } else if (x_is_cat && y_is_cat) {
        mixed_mat[i, j] <- cramers_v(x, y)
        
      } else if (x_is_num && y_is_cat) {
        mixed_mat[i, j] <- correlation_ratio(y, x)
        
      } else if (x_is_cat && y_is_num) {
        mixed_mat[i, j] <- correlation_ratio(x, y)
      }
    }
  }
  
  # Create the .csv matrix with every variable (all types)
  write.csv(
    mixed_mat,
    file.path(tbl_dir, "association_matrix_mixed_variables.csv"),
    row.names = TRUE
  )
  
  # Full mixed plot can get too crowded, so create a top 20 variables plot
  mixed_long <- as.data.frame(as.table(mixed_mat))
  mixed_long <- mixed_long[mixed_long$Var1 != mixed_long$Var2, , drop = FALSE]
  mixed_long <- mixed_long[!is.na(mixed_long$Freq), , drop = FALSE]
  
  if (nrow(mixed_long) > 0) {
    mixed_long$key <- apply(mixed_long[, c("Var1", "Var2")], 1, function(z) {
      paste(sort(z), collapse = "__")
    })
    mixed_long <- mixed_long[!duplicated(mixed_long$key), , drop = FALSE]
    mixed_long <- mixed_long[order(-mixed_long$Freq), , drop = FALSE]
    
    top_n_mixed <- min(20, nrow(mixed_long))
    top_mixed <- head(mixed_long, top_n_mixed)
    
    mixed_top_vars <- unique(c(
      as.character(top_mixed$Var1),
      as.character(top_mixed$Var2)
    ))
    
    top_mixed_export <- top_mixed
    top_mixed_export$Var1 <- get_display_name(top_mixed_export$Var1)
    top_mixed_export$Var2 <- get_display_name(top_mixed_export$Var2)
    
    # Just the top 20 associations (including mixed types) in a .csv
    write.csv(
      top_mixed_export,
      file.path(tbl_dir, "top20_mixed_variable_associations.csv"),
      row.names = FALSE
    )
    
    if (length(mixed_top_vars) >= 2) {
      mixed_top_mat <- mixed_mat[mixed_top_vars, mixed_top_vars, drop = FALSE]
      
      colnames(mixed_top_mat) <- get_display_name(colnames(mixed_top_mat))
      rownames(mixed_top_mat) <- get_display_name(rownames(mixed_top_mat))
      
      # Create a heatmap for mixed variables too
      png(file.path(fig_dir, "association_heatmap_mixed_variables.png"),
          width = 1800, height = 2400, res = 220)
      
      par(mar = c(1, 1, 8, 1), oma = c(0, 0, 3.5, 0))
      
      corrplot(
        mixed_top_mat,
        is.corr = FALSE,
        method = "color",
        type = "upper",
        main = "",
        tl.cex = 0.72,
        tl.col = "black",
        order = "hclust",
        col.lim = c(0, 1)
      )
      
      mtext(
        "Mixed Variable Association Heatmap",
        side = 3,
        outer = TRUE,
        line = -4,
        cex = 1.1,
        font = 2
      )
      
      dev.off()
    }
  }
}
# ----------------------------------------------------------
# Exploratory Plots 
# ----------------------------------------------------------

# We include the top 6 (expcept sbp) Framingham Risk Score variables for predicting CVD:
# Age, Total Chol, HDL chol, Smoking Status, Gender 
# (sbp moved since it was a secondary mediator)

# Add outcome back for EDA
eda_df <- df_analysis
eda_df$has_cvd <- factor(df$has_cvd, levels = c(0,1), labels = c("No CVD","CVD"))

boxplot_vars <- c("age", "ldl_direct", "hdl_cholesterol",
                  "total_cholesterol", "hba1c", "crp")

boxplot_vars <- intersect(boxplot_vars, names(eda_df))
boxplot_vars <- boxplot_vars[vapply(eda_df[boxplot_vars], is.numeric, logical(1))]

# Create summary stats for numeric variables
boxplot_summary <- bind_rows(lapply(boxplot_vars, function(var) {
  eda_df %>%
    group_by(has_cvd) %>%
    summarise(
      variable = var,
      display_name = get_display_name(var),
      n = sum(!is.na(.data[[var]])),
      mean = mean(.data[[var]], na.rm = TRUE),
      sd = sd(.data[[var]], na.rm = TRUE),
      median = median(.data[[var]], na.rm = TRUE),
      q1 = quantile(.data[[var]], 0.25, na.rm = TRUE),
      q3 = quantile(.data[[var]], 0.75, na.rm = TRUE),
      min = min(.data[[var]], na.rm = TRUE),
      max = max(.data[[var]], na.rm = TRUE),
      .groups = "drop"
    )
}))

# Create numeric summary to make plots later if necessary
if (nrow(boxplot_summary) > 0) {
  write.csv(
    boxplot_summary,
    file.path(tbl_dir, "eda_numeric_grouped_summary_by_cvd.csv"),
    row.names = FALSE
  )
}

# Create summary stats for categorical variables
categorical_vars <- c("sex", "smoking_status_refined", "alcohol_consumption")
categorical_vars <- categorical_vars[categorical_vars %in% names(eda_df)]

categorical_summary <- bind_rows(lapply(categorical_vars, function(var) {
  
  if (!var %in% names(eda_df)) return(NULL)
  
  eda_df %>%
    count(.data[[var]], has_cvd, name = "n") %>%
    group_by(.data[[var]]) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(
      variable = var,
      display_name = get_display_name(var)
    ) %>%
    rename(category = 1)
}))

# Create categorical summary to make plots later if necessary
if (nrow(categorical_summary) > 0) {
  write.csv(
    categorical_summary,
    file.path(tbl_dir, "eda_categorical_counts_proportions_by_cvd.csv"),
    row.names = FALSE
  )
}

# ----------------------------------------------------------
# Boxplots
# ----------------------------------------------------------

plot_box <- function(var){
  if (!var %in% names(eda_df)) return(NULL)
  if (!is.numeric(eda_df[[var]])) return(NULL)
  
  var_label <- get_display_name(var)
  
  ggplot(eda_df, aes(x = has_cvd, y = .data[[var]], fill = has_cvd)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    scale_fill_manual(values = c("No CVD" = "#a8d5e2", "CVD" = "#3a6ea5")) +
    labs(title = paste(var_label, "by CVD Status"), x = "", y = var_label) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 13, face = "bold", colour = "#2c3e50", hjust = 0.5),
      axis.title = element_text(face = "bold", colour = "#2c3e50"),
      axis.text = element_text(size = 9, colour = "#2c3e50"),
      axis.line = element_line(colour = "#7f8c8d", linewidth = 0.6),
      axis.ticks = element_line(colour = "#7f8c8d", linewidth = 0.6),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#c7c7c7", linewidth = 0.7),
      panel.grid.minor = element_blank()
    )
}

# For loop to create box plots by CVD status

for (var in boxplot_vars) {
  p_box <- plot_box(var)
  if (!is.null(p_box)) {
    ggsave(
      file.path(fig_dir, paste0("boxplot_", var, "_by_cvd.png")),
      plot = p_box,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
}

# ----------------------------------------------------------
# Density plots (distribution differences) 
# ----------------------------------------------------------

plot_density <- function(var){
  if (!var %in% names(eda_df)) return(NULL)
  if (!is.numeric(eda_df[[var]])) return(NULL)
  
  var_label <- get_display_name(var)
  
  # Density plots for distribution by CVD
  
  ggplot(eda_df, aes(x = .data[[var]], colour = has_cvd)) +
    geom_density(linewidth = 1) +
    scale_colour_manual(values = c("No CVD" = "#1b1b1b", "CVD" = "#2f6db3")) +
    labs(title = paste(var_label, "distribution by CVD"), x = var_label, y = "Density") +
    theme_minimal(base_size = 12) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 13, face = "bold", colour = "#2c3e50", hjust = 0.5),
      axis.title = element_text(face = "bold", colour = "#2c3e50"),
      axis.text = element_text(size = 9, colour = "#2c3e50"),
      axis.line = element_line(colour = "#7f8c8d", linewidth = 0.6),
      axis.ticks = element_line(colour = "#7f8c8d", linewidth = 0.6),
      panel.grid.major.x = element_line(colour = "#d0d0d0", linewidth = 0.6),
      panel.grid.major.y = element_line(colour = "#c7c7c7", linewidth = 0.7),
      panel.grid.minor = element_blank()
    )
}

# For loop to create Density plots for each variable by CVD

for (var in boxplot_vars) {
  p_density <- plot_density(var)
  if (!is.null(p_density)) {
    ggsave(
      file.path(fig_dir, paste0("density_", var, "_by_cvd.png")),
      plot = p_density,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
}

# ----------------------------------------------------------
# Categorical variables vs outcome: Bar plots
# ----------------------------------------------------------

plot_bar <- function(var){
  if (!var %in% names(eda_df)) return(NULL)
  
  var_label <- get_display_name(var)
  
  plot_df <- eda_df %>%
    count(.data[[var]], has_cvd) %>%
    group_by(.data[[var]]) %>%
    mutate(prop = n / sum(n),
           label = scales::percent(prop, accuracy = 1)) %>%
    ungroup()
  
  # Create bar plots for CVD by proportion
  
  ggplot(plot_df, aes(x = .data[[var]], y = prop, fill = has_cvd)) +
    geom_col(position = "stack") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      color = "white",
      size = 4
    ) +
    scale_fill_manual(values = c("No CVD" = "#a8d5e2", "CVD" = "#3a6ea5")) +
    labs(title = paste("CVD proportion by", var_label),
         x = var_label, y = "Proportion") +
    theme_minimal(base_size = 12) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 13, face = "bold", colour = "#2c3e50", hjust = 0.5),
      axis.title = element_text(face = "bold", colour = "#2c3e50"),
      axis.text = element_text(size = 9, colour = "#2c3e50"),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 9, colour = "#2c3e50"),
      axis.line = element_line(colour = "#7f8c8d", linewidth = 0.6),
      axis.ticks = element_line(colour = "#7f8c8d", linewidth = 0.6),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#c7c7c7", linewidth = 0.7),
      panel.grid.minor = element_blank()
    )
}

# For loop to create bar plot for each categorical variable

for (var in categorical_vars) {
  p_bar <- plot_bar(var)
  if (!is.null(p_bar)) {
    ggsave(
      file.path(fig_dir, paste0("barplot_", var, "_by_cvd.png")),
      plot = p_bar,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
}

