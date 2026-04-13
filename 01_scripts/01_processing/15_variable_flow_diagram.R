# ==============================================================================
# VARIABLE FLOW DIAGRAM CALCULATOR 
# ==============================================================================
rm(list = ls())

library(here)
library(dplyr)
library(tibble)

cat("=======================================================\n")
cat("Extracting Variable Flow Diagram Values...\n")
cat("=======================================================\n")


proc_dir <- here("00_data", "01_processed")

#1.
domain_files <- list(
  "Socio-demographics" = "sociodemographic_cleaned.rds",
  "Lifestyle"          = "lifestyle_cleaned.rds",
  "Diet & Nutrition"   = "diet_cleaned.rds",
  "Environmental"      = "environment_cleaned.rds",
  "Early Life Factors" = "early_life_cleaned.rds",
  "Physical Measures"  = "physical_cleaned.rds",
  "Mental Health"      = "mental_health_cleaned.rds",
  "Medical History"    = "medical_history_cleaned.rds",
  "Biological Samples" = "bio_samples_cleaned.rds"
)

# 2. HELPER FUNCTION: Count variables (excluding the 'eid' ID column)
count_vars <- function(file_name) {
  file_path <- file.path(proc_dir, file_name)
  if (file.exists(file_path)) {
    cols <- colnames(readRDS(file_path))
    # Subtract 'eid' because it's an ID, not a predictor
    var_count <- length(setdiff(cols, c("eid", "ID", "id"))) 
    return(var_count)
  } else {
    warning(paste("Missing file:", file_name))
    return(NA)
  }
}

# 3. CALCULATE DOMAIN VARIABLES (Right Box)
cat("Counting variables in each domain file...\n")

domain_counts <- numeric(length(domain_files))
names(domain_counts) <- names(domain_files)

for (i in seq_along(domain_files)) {
  domain_counts[i] <- count_vars(domain_files[[i]])
}

domain_table <- tibble(
  Domain = names(domain_counts),
  Variables_p = domain_counts
)
total_extracted <- sum(domain_counts, na.rm = TRUE)

# 4. CALCULATE MAIN PIPELINE VARIABLES (Left Boxes)
cat("Counting variables in main pipeline files...\n")

p_merged  <- count_vars("ukb_final_merged.rds")
p_cleaned <- count_vars("ukb_final_cleaned.rds")

# Identify variables removed for high missingness
names_merged  <- names(readRDS(here("00_data", "01_processed", "ukb_final_merged.rds")))
names_cleaned <- names(readRDS(here("00_data", "01_processed", "ukb_final_cleaned.rds")))
excluded_varnames <- setdiff(names_merged, names_cleaned)

# Create Flow Table 
var_flow_table <- tibble(
  Diagram_Box = c(
    "1. Variables extracted from UK Biobank (Sum of Domains)",
    "2. Pre-processed variables (ukb_final_merged.rds)",
    paste0("3. Excluded: high missingness (", paste(excluded_varnames, collapse = ", "), ")"),
    "4. Final analysis variables (ukb_final_cleaned.rds)"
  ),
  Number_of_Variables_p = c(
    total_extracted,
    p_merged,
    ifelse(!is.na(p_merged) & !is.na(p_cleaned), p_merged - p_cleaned, NA),
    p_cleaned
  )
)

# 5. EXPORT AND PRINT
out_dir <- here("02_results", "01_tables", "00_flow_diagrams")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write.csv(domain_table, here("02_results", "01_tables", "00_flow_diagrams", "variable_domains_p.csv"), row.names = FALSE)
write.csv(var_flow_table, here("02_results", "01_tables", "00_flow_diagrams", "variable_flow_p.csv"), row.names = FALSE)

cat("\n--- DOMAINS BOX (Right side of diagram) --- \n")
print(domain_table)
cat(sprintf("Total Extracted 'p' = %d\n", total_extracted))

cat("\n--- VARIABLE FLOW (Left side of diagram) --- \n")
print(var_flow_table)

cat("\n=======================================================\n")
cat("SUCCESS! Variable flow numbers exported to:\n", out_dir, "\n")
cat("=======================================================\n")