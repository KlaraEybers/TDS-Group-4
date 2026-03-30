# Finds leaf folders that do not have .gitkeep, hashed section creates the file

# check_gitkeep.R
# Checks all leaf directories for the presence of a .gitkeep file
library(here)

# ── Find leaf dirs within specific folders only ───────────────────────────
search_dirs <- c(here("02_results"), here("00_data"))

leaf_dirs <- unlist(lapply(search_dirs, function(d) {
  system2("find",
          args = c(d, "-type", "d",
                   "-not", "-path", "*/.git/*",
                   "-links", "2"),
          stdout = TRUE)
}))

# ── Check for .gitkeep in each leaf directory ─────────────────────────────
missing <- leaf_dirs[!sapply(leaf_dirs, function(d) {
  file.exists(file.path(d, ".gitkeep"))
})]

# ── Report ────────────────────────────────────────────────────────────────
if (length(missing) == 0) {
  message("All leaf directories contain a .gitkeep file.")
} else {
  message(sprintf("%d leaf director%s missing .gitkeep:\n",
                  length(missing),
                  ifelse(length(missing) == 1, "y is", "ies are")))
  cat(paste0("  ", gsub(here(), ".", missing, fixed = TRUE), "\n"))
}

# dirs <- c(
#   "./02_results/01_tables/04_neural_network/.ipynb_checkpoints",
#   "./02_results/01_tables/archive",
#   "./02_results/00_figures/04_exploratory/lasso_incremental/incremental_analysis_copy",
#   "./02_results/00_figures/04_exploratory/lasso_incremental/incremental_analysis",
#   "./02_results/00_figures/04_exploratory/stability_selection_lasso",
#   "./02_results/00_figures/02_univariate/archive",
#   "./02_results/00_figures/01_imputation/archive",
#   "./02_results/00_figures/07_sensitivity_complete_cases",
#   "./02_results/00_figures/03_multivariate/archive",
#   "./02_results/03_ml/00_neural_network",
#   "./02_results/03_ml/nn_outputs",
#   "./02_results/03_ml/.ipynb_checkpoints",
#   "./02_results/03_ml/tree_models_python_confounders_plus_selected",
#   "./02_results/03_ml/tree_models_python_all/random_forest"
# )
# 
# results <- sapply(dirs, function(d) {
#   full_path <- here(gsub("^\\./", "", d))  # strip leading ./ then anchor with here()
#   file.create(file.path(full_path, ".gitkeep"))
# })
# 
# cat(ifelse(results, "OK ", "FAIL"), names(results), sep = "\n")
