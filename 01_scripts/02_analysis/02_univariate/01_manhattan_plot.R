rm(list = ls())

# Import packages ---------------------------------------------------------
library(here)
library(dplyr)
library(ggplot2)

# Load results ------------------------------------------------------------
res <- readRDS(here("00_data", "02_analysis", "univariate_logistic_results.rds"))

# Parameters --------------------------------------------------------------
SIG_RAW <- 0.05
SIG_BH  <- 0.05

# Domain order ------------------------------------------------------------
domain_levels <- c(
  "Socio-demographics",
  "Physical measures",
  "Diet",
  "Lifestyle",
  "Health and medical history",
  "Mental health",
  "Early life factors",
  "Environment",
  "Biological samples"
)

# Prepare plot data -------------------------------------------------------
# neg_log10 columns are computed here from raw p-values
plot_df <- res %>%
  filter(!is.na(effect_log_or), !is.na(p_value), !is.na(p_adj_bh)) %>%
  mutate(
    domain             = factor(domain, levels = domain_levels),
    direction          = ifelse(effect_log_or >= 0, "Positive OR", "Negative OR"),
    neg_log10_p        = -log10(p_value),
    neg_log10_p_adj_bh = -log10(p_adj_bh)
  ) %>%
  arrange(domain, variable) %>%
  mutate(x_pos = row_number())

# Domain midpoints and boundary lines -------------------------------------
domain_mids <- plot_df %>%
  group_by(domain) %>%
  summarise(mid = mean(x_pos), .groups = "drop")

domain_bounds <- plot_df %>%
  group_by(domain) %>%
  summarise(xmin = min(x_pos) - 0.5, xmax = max(x_pos) + 0.5, .groups = "drop")

# Shared theme ------------------------------------------------------------
manhattan_theme <- theme_bw() +
  theme(
    axis.title         = element_text(size = 11, face = "bold", colour = "#2c3e50"),
    plot.title         = element_text(size = 13, face = "bold", colour = "#2c3e50"),
    plot.subtitle      = element_text(size = 10, colour = "grey40",
                                      margin = margin(b = 8)),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "#e0e0e0", linewidth = 0.4),
    plot.margin        = margin(10, 15, 80, 12)
  )

# Helper: build lollipop Manhattan plot -----------------------------------
make_manhattan_lollipop <- function(data, y_var, y_label, sig_line,
                                    sig_col = "firebrick", title_suffix) {
  
  # Handle Inf values (from p = 0)
  y_raw      <- data[[y_var]]
  finite_max <- max(y_raw[is.finite(y_raw)], na.rm = TRUE)
  data <- data %>%
    mutate(across(all_of(y_var), ~ ifelse(is.infinite(.x), finite_max + 10, .x)))
  
  y_raw <- data[[y_var]]
  y_cap <- quantile(y_raw, 0.99, na.rm = TRUE)
  y_max <- y_cap * 1.15
  
  data <- data %>%
    mutate(
      capped      = .data[[y_var]] > y_cap,
      y_plot      = pmin(.data[[y_var]], y_cap),
      significant = .data[[y_var]] >= -log10(sig_line),
      label       = ifelse(capped, paste0(variable, "*"), variable)
    )
  
  sig_y  <- -log10(sig_line)
  y_lbl  <- y_max * 0.01  # labels start just above baseline, inside panel
  
  label_df <- data %>% select(x_pos, label, significant)
  
  ggplot(data, aes(x = x_pos, y = y_plot)) +
    
    # Alternating domain background bands
    geom_rect(
      data        = domain_bounds,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill        = rep(c("grey95", "white"), length.out = nrow(domain_bounds)),
      alpha       = 0.6
    ) +
    
    # Lollipop stem
    geom_segment(
      aes(x = x_pos, xend = x_pos, y = 0, yend = y_plot,
          colour = direction),
      linewidth = 0.45, alpha = 0.75
    ) +
    
    # Lollipop head
    geom_point(
      aes(colour = direction,
          size   = ifelse(significant, 2.4, 1.5),
          alpha  = ifelse(significant, 1, 0.6))
    ) +
    scale_size_identity() +
    scale_alpha_identity() +
    
    scale_colour_manual(
      name   = "Effect direction",
      values = c("Positive OR" = "#3A6EA5", "Negative OR" = "#C0392B")
    ) +
    
    # Significance threshold line
    geom_hline(yintercept = min(sig_y, y_cap * 0.98),
               linetype = "dashed", colour = sig_col, linewidth = 0.65) +
    annotate("text",
             x      = max(data$x_pos) * 0.02,
             y      = min(sig_y, y_cap * 0.98) + y_max * 0.025,
             label  = paste0("p = ", sig_line),
             colour = sig_col, size = 3, hjust = 0) +
    
    # Non-significant variable labels — plain, grey, small, angled inside panel
    geom_text(
      data        = filter(label_df, !significant),
      aes(x = x_pos, y = y_lbl, label = label),
      angle       = 90, hjust = 0, vjust = 0.5,
      size        = 1.9, colour = "grey60", fontface = "plain",
      inherit.aes = FALSE
    ) +
    
    # Significant variable labels — bold, dark, larger
    geom_text(
      data        = filter(label_df, significant),
      aes(x = x_pos, y = y_lbl, label = label),
      angle       = 90, hjust = 0, vjust = 0.5,
      size        = 2.4, colour = "#2c3e50", fontface = "bold",
      inherit.aes = FALSE
    ) +
    
    # X-axis: domain name labels only
    scale_x_continuous(
      breaks = domain_mids$mid,
      labels = domain_mids$domain,
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = c(0, 0),
      breaks = function(lims) pretty(c(0, lims[2]), n = 6),
      labels = function(b) ifelse(b < 0, "", b)
    ) +
    
    labs(
      title    = paste("Manhattan plot \u2014 associations with CVD", title_suffix),
      subtitle = paste0("Bold labels = significant at p < ", sig_line,
                        "  |  * = value capped at 99th percentile"),
      x = "Domain",
      y = y_label
    ) +
    manhattan_theme +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 9,
                                 colour = "#2c3e50", face = "bold")
    )
}

# Plot 1: raw p-value -----------------------------------------------------
p_raw <- make_manhattan_lollipop(
  data         = plot_df,
  y_var        = "neg_log10_p",
  y_label      = expression(-log[10](p~raw)),
  sig_line     = SIG_RAW,
  title_suffix = "(raw p-value)"
)

ggsave(
  filename = here("00_data", "02_analysis", "manhattan_lollipop_raw_p.png"),
  plot     = p_raw,
  width    = 16, height = 9, dpi = 300, bg = "white"
)

# Plot 2: BH-adjusted p-value ---------------------------------------------
p_bh <- make_manhattan_lollipop(
  data         = plot_df,
  y_var        = "neg_log10_p_adj_bh",
  y_label      = expression(-log[10](p~BH-adjusted)),
  sig_line     = SIG_BH,
  title_suffix = "(Benjamini-Hochberg-adjusted p-value)"
)

ggsave(
  filename = here("00_data", "02_analysis", "manhattan_lollipop_bh_p.png"),
  plot     = p_bh,
  width    = 16, height = 9, dpi = 300, bg = "white"
)

message("Lollipop Manhattan plots saved to 00_data/02_analysis/")