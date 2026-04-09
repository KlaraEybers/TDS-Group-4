rm(list = ls())

# Import packages ---------------------------------------------------------

library(here)
library(dplyr)
library(FactoMineR)
library(corrplot)
library(ggplot2)

# Import data -------------------------------------------------------------

df <- readRDS(here("00_data", "01_processed", "03_imputation", "df_train_imputed.rds"))

# Remove ID / outcome / dates ------------------------------------------------

colnames(df)

df_analysis <- df %>%
  select(
    -eid,
    -has_cvd,
    -recruitment_date,
    -death_date,
    -first_cvd_date,
    -sleep_cat,
    -passive_cat,
  ) %>%
  mutate(across(where(is.character), as.factor))

# Run FAMD ------------------------------------------------

set.seed(1)
famd_res <- FAMD(df_analysis, ncp = 200, graph = FALSE)

# Scree Plot ------------------------------------------------

eig_pct <- famd_res$eig[, 2]
plot(eig_pct, type = "b",
     xlab = "Component", ylab = "% variance explained",
     main = "FAMD Scree Plot")

famd_res$eig
famd_res$eig[, 3]
sum(famd_res$eig[,2])
tail(famd_res$eig[,3], 1)

colnames(df_analysis)

# Cumulative Scree ------------------------------------------------

cum_var <- famd_res$eig[, 3]

plot(cum_var, type = "b",
     xlab = "FAMD Dimension",
     ylab = "Cumulative variance explained (%)",
     main = "FAMD Cumulative Scree Plot",
     ylim = c(0, 100),
     yaxt = "n")

axis(2,
     at = seq(0, 100, by = 20),
     labels = seq(0, 100, by = 20),
     las = 1)

# Interpreting the PCs ------------------------------------------------

num_contrib <- famd_res$quanti.var$contrib
cat_contrib <- famd_res$quali.var$contrib

sort(num_contrib[,7], decreasing = TRUE)[1:10]
sort(cat_contrib[,7], decreasing = TRUE)[1:10]

plot_top_contrib <- function(contrib_matrix, pc_num, title_prefix) {
  
  # Extract contributions for chosen PC
  pc_vals <- contrib_matrix[, pc_num]
  
  # Sort and keep top 6
  top_vars <- sort(pc_vals, decreasing = TRUE)[1:6]
  
  # Barplot
  barplot(top_vars,
          las = 0,                     # rotate labels
          col = "steelblue",
          border = NA,
          ylim = c(0, 30), 
          las = 1,
          main = paste(title_prefix, "- PC", pc_num),
          ylab = "Contribution (%)")
}

# for (i in 1:7) {
#   plot_top_contrib(num_contrib, i, "Numeric Contributors")
# }

plot_top_contrib(num_contrib, 6, "Numeric Contributors")

# Plot (coloured by has_cvd) ------------------------------------------------

scores <- famd_res$ind$coord

col_outcome <- ifelse(df$has_cvd == 1,
                      adjustcolor("red", alpha.f = 0.5),
                      adjustcolor("black", alpha.f = 0.35))

plot(scores[,1], scores[,2],
     col = col_outcome,
     pch = 16,
     cex = 0.5,
     xlab = "Dimension 1 (3.15%)",
     ylab = "Dimension 2 (3.13%)",
     main = "FAMD: Individuals (coloured by CVD)",
     ylim = c(-6, 8),
     yaxt = "n")

axis(2, at = seq(-6, 8, by = 2), las = 1)

legend("topright",
       legend = c("No CVD", "CVD"),
       col = c("black", "red"),
       pch = 16,
       bty = "n")

# Plot (coloured by sex) ------------------------------------------------

sex_fac <- as.factor(df$sex)

cols <- c(
  adjustcolor("#106CEF", alpha.f = 0.30), 
  adjustcolor("#EF9310", alpha.f = 0.45)  
)

plot(scores[,1], scores[,2],
     col = cols[sex_fac],
     pch = 16,
     cex = 0.5,
     xlab = "Dimension 1 (3.15%)",
     ylab = "Dimension 2 (3.13%)",
     main = "FAMD: Individuals (coloured by Sex)",
     ylim = c(-6, 8),
     yaxt = "n")

axis(2, at = seq(-6, 8, by = 2), las = 1)

legend("topright",
       legend = levels(sex_fac),
       col = c("#106CEF", "#EF9310"),
       pch = 16,
       bty = "n")

# Correlation Heatmap ------------------------------------------------

num_data <- df_analysis %>% select(where(is.numeric))

# Correlation matrix
cor_mat <- cor(num_data, use = "pairwise.complete.obs", method = "pearson")

# Set lower triangle and diagonal to NA
cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA

# Convert to long format
cor_pairs <- as.data.frame(as.table(cor_mat))

# Remove NAs
cor_pairs <- cor_pairs[!is.na(cor_pairs$Freq), ]

# Rank by absolute correlation
cor_pairs$abs_cor <- abs(cor_pairs$Freq)
strongest <- cor_pairs[order(-cor_pairs$abs_cor), ]

# Top 15
head(strongest, 15)

# Cramer's V for Categorical vs Categorical correlation ----------------

cramers_v <- function(x, y) {
  tab <- table(x, y)
  
  if (min(dim(tab)) < 2) return(NA)
  
  chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  n <- sum(tab)
  
  k <- min(nrow(tab), ncol(tab))
  
  sqrt(chi2 / (n * (k - 1)))
}

cat_data <- df_analysis %>% select(where(is.factor))

vars <- colnames(cat_data)

cat_assoc <- matrix(NA, nrow = length(vars), ncol = length(vars))
colnames(cat_assoc) <- vars
rownames(cat_assoc) <- vars

for (i in seq_along(vars)) {
  for (j in seq_along(vars)) {
    cat_assoc[i, j] <- cramers_v(cat_data[[vars[i]]], cat_data[[vars[j]]])
  }
}

cat_pairs <- as.data.frame(as.table(cat_assoc))
cat_pairs <- cat_pairs[!is.na(cat_pairs$Freq), ]

cat_pairs$abs_assoc <- abs(cat_pairs$Freq)
cat_pairs <- cat_pairs[order(-cat_pairs$abs_assoc), ]

cat_pairs <- as.data.frame(as.table(cat_assoc))

# Remove NA
cat_pairs <- cat_pairs[!is.na(cat_pairs$Freq), ]

# Remove self-comparisons
cat_pairs <- cat_pairs[cat_pairs$Var1 != cat_pairs$Var2, ]

# Remove duplicate pairs (A-B and B-A)
cat_pairs <- cat_pairs %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(Var1, Var2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

# Rank
cat_pairs$abs_assoc <- abs(cat_pairs$Freq)
cat_pairs <- cat_pairs[order(-cat_pairs$abs_assoc), ]

# Top results
head(cat_pairs, 15)

summary(cat_data)

# Anova for Continuous vs Categorical --------------------

# Exploratory Plots --------------------------------------------

# Add outcome back for EDA -----------------------------------------------
eda_df <- df_analysis
eda_df$has_cvd <- factor(df$has_cvd, levels = c(0,1), labels = c("No CVD","CVD"))

# Boxplots for strongest clinical predictors ---------------------------------

plot_box <- function(var){
  ggplot(eda_df, aes(x = has_cvd, y = .data[[var]], fill = has_cvd)) +
    
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    
    labs(title = paste(var, "by CVD Status"), x = "", y = var) +
    
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
}

plot_box("age")
plot_box("sbp")
plot_box("ldl_direct")
plot_box("hdl_cholesterol")
plot_box("total_cholesterol")
plot_box("hba1c")
plot_box("crp")

# Density plots (distribution differences) ----------------------

plot_density <- function(var){
  ggplot(eda_df, aes(x = .data[[var]], colour = has_cvd)) +
    geom_density() +
    labs(title = paste(var, "distribution by CVD"), x = var) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)
          )
}

plot_density("age")
plot_density("sbp")
plot_density("ldl_direct")
plot_density("hdl_cholesterol")
plot_density("total_cholesterol")
plot_density("hba1c")
plot_density("crp")

# Categorical variables vs outcome --------------------------------

plot_bar <- function(var){
  
  plot_df <- eda_df %>%
    count(.data[[var]], has_cvd) %>%
    group_by(.data[[var]]) %>%
    mutate(prop = n / sum(n),
           label = scales::percent(prop, accuracy = 1)) %>%
    ungroup()
  
  ggplot(plot_df, aes(x = .data[[var]], y = prop, fill = has_cvd)) +
    geom_col(position = "stack") +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              color = "white",
              size = 4) +
    labs(title = paste("CVD proportion by", var),
         x = var, y = "Proportion") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_bar("sex")
plot_bar("smoking_status_refined")
# plot_bar("diabetes_doctor_diagnosed")
plot_bar("alcohol_consumption")

-----------------------------------------------------------

# Count numeric contribution: 1 dimension each
num_vars <- sapply(df_analysis, is.numeric)
n_numeric_dims <- sum(num_vars)

# Count categorical contribution: (levels - 1) each
fac_vars <- sapply(df_analysis, is.factor)
factor_levels <- sapply(df_analysis[, fac_vars, drop = FALSE], nlevels)
n_factor_dims <- sum(factor_levels - 1)

# Total effective dimensionality
cd <- n_numeric_dims + n_factor_dims

cat("Estimated maximum FAMD dimensions/components:", max_famd_dims, "\n")
