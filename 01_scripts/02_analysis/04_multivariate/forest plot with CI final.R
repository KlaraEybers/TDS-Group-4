# ==============================================================================
# Forest Plot: Total vs Direct Effect with 95% CIs (Fully Fixed Layout)
# ==============================================================================
library(ggplot2)
library(dplyr)
library(stringr)

cat("Setting up paths and loading models...\n")
results_dir <- "/rds/general/user/sh425/projects/hda_25-26/live/TDS/TDS_Group4/02_results/01_tables"
figures_dir <- "/rds/general/user/sh425/projects/hda_25-26/live/TDS/TDS_Group4/02_results/00_figures"

# 1. Load locked models
model_total <- readRDS(file.path(results_dir, "model_total_final1.rds"))
model_direct <- readRDS(file.path(results_dir, "model_direct_final1.rds"))

# 2. Function to extract Betas and 95% CIs
get_estimates <- function(model, model_name) {
  est <- summary(model)$coefficients
  vars <- rownames(est)
  
  # Extract Beta and Standard Error
  beta <- est[, "Estimate"]
  se <- est[, "Std. Error"]
  
  # Calculate 95% CI (Beta +/- 1.96 * SE)
  lower <- beta - 1.96 * se
  upper <- beta + 1.96 * se
  
  df <- data.frame(
    Variable = vars,
    Beta_Value = beta,
    CI_Lower = lower,
    CI_Upper = upper,
    Effect_Type = model_name,
    stringsAsFactors = FALSE
  )
  
  # Remove intercept
  df <- df[df$Variable != "(Intercept)", ]
  return(df)
}

# Extract full data from both models
cat("Extracting Betas and 95% CIs...\n")
df_total <- get_estimates(model_total, "Total Effect (Before Mediator)")
df_direct <- get_estimates(model_direct, "Direct Effect (After Mediator)")

# Keep only common variables (exclude mediator from plotting)
common_vars <- intersect(df_total$Variable, df_direct$Variable)
df_total <- df_total[df_total$Variable %in% common_vars, ]
df_direct <- df_direct[df_direct$Variable %in% common_vars, ]

# Combine into long format
plot_data_long <- rbind(df_total, df_direct)

# ------------------------------------------------------------------------------
# 3. Clean variable names for better readability
# ------------------------------------------------------------------------------
plot_data_long$Variable <- gsub("smoking_status_refined", "Smoking: ", plot_data_long$Variable)
plot_data_long$Variable <- gsub("employment", "Employment: ", plot_data_long$Variable)
plot_data_long$Variable <- gsub("housing", "Housing: ", plot_data_long$Variable)
plot_data_long$Variable <- gsub("salt_added", "Added Salt: ", plot_data_long$Variable)
plot_data_long$Variable <- gsub("_", " ", plot_data_long$Variable)

# Order variables by Total Effect size for a staircase look
order_df <- df_total
order_df$Clean_Var <- plot_data_long$Variable[plot_data_long$Effect_Type == "Total Effect (Before Mediator)"]
order_df <- order_df[order(order_df$Beta_Value), ]

plot_data_long$Variable <- factor(plot_data_long$Variable, levels = order_df$Clean_Var)

# ------------------------------------------------------------------------------
# 4. Generate dodged comparative plot with CIs 
# ------------------------------------------------------------------------------
cat("Generating forest plot with fixed title alignment...\n")

dodge_width <- 0.5 

p <- ggplot(plot_data_long, aes(x = Beta_Value, y = Variable, color = Effect_Type)) +
  
  # Reference line at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 1) +
  
  # Horizontal error bars for 95% CIs
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), 
                 height = 0.2, size = 0.8, alpha = 0.7, 
                 position = position_dodge(width = dodge_width)) +
  
  # Points for Beta estimates
  geom_point(size = 3.5, alpha = 0.9, position = position_dodge(width = dodge_width)) +
  
  # Custom colors
  scale_color_manual(values = c("Total Effect (Before Mediator)" = "#1f78b4", 
                                "Direct Effect (After Mediator)" = "#33a02c")) +
  
  # Titles and labels
  labs(title = "Mediation Analysis: Shift in Effect Sizes with 95% CIs",
       subtitle = "Horizontal bars represent 95% Confidence Intervals",
       x = "Beta Coefficient (Log-Odds)",
       y = "",
       color = "Model") +
  
  # Refined theme
  theme_minimal() +
  theme(
    # [CRITICAL FIX 1]: Align title to the entire plot, not just the panel, preventing right cutoff
    plot.title.position = "plot", 
    
    # Increase margins to prevent cutoff (Top, Right, Bottom, Left)
    plot.margin = margin(t = 20, r = 30, b = 20, l = 20),
    
    # Text styling
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray30", margin = margin(b = 15)),
    
    axis.text.y = element_text(size = 11, face = "plain", color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    
    # Legend styling
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

# 5. Export high-res PNG
png_path <- file.path(figures_dir, "Mediation_Forest_Plot_with_CI.png")

# [CRITICAL FIX 2]: Increased width to 11.5 and height to 7 for plenty of horizontal space
ggsave(png_path, plot = p, width = 11.5, height = 7, dpi = 300, bg = "white")

cat("SUCCESS: Perfectly aligned forest plot saved as 'Mediation_Forest_Plot_with_CI.png'!\n")

