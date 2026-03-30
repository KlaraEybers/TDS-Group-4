rm(list = ls())

# Import packages ---------------------------------------------------------
library(here)
library(dplyr)
library(ggplot2)

# Load naming dictionary --------------------------------------------------
source(here("01_scripts", "00_variable_dictionary.R"))


get_domain_display_name <- function(x) {
  unname(ifelse(x %in% names(domain_display_names), domain_display_names[x], x))
}

# Load results ------------------------------------------------------------
res <- readRDS(here("00_data", "02_analysis", "univariate_logistic_results.rds"))

# Parameters --------------------------------------------------------------
SIG_RAW <- 0.05
SIG_BH  <- 0.05

# Domain order ------------------------------------------------------------
domain_levels <- unique(unname(var_domains))

# Prepare plot data -------------------------------------------------------
plot_df <- res %>%
  filter(!is.na(effect_log_or), !is.na(p_value), !is.na(p_adj_bh)) %>%
  mutate(
    p_adj_bonf = p.adjust(p_value, method = "bonferroni"),
    variable_label = get_display_name(variable),
    domain_label   = get_domain_from_variable(variable),
    direction      = ifelse(effect_log_or >= 0, "Positive OR", "Negative OR"),
    neg_log10_p    = -log10(p_value),
    neg_log10_p_adj_bh   = -log10(p_adj_bh),
    neg_log10_p_adj_bonf = -log10(p_adj_bonf),
    domain_label   = factor(domain_label, levels = domain_levels)  
  ) %>%
  arrange(domain_label) %>%           
  mutate(x_pos = row_number())        

# Domain midpoints and boundary lines 
domain_mids <- plot_df %>%
  group_by(domain_label) %>%
  summarise(mid = mean(x_pos), .groups = "drop") %>%
  filter(!is.na(domain_label))        

domain_bounds <- plot_df %>%
  group_by(domain_label) %>%
  summarise(
    xmin = min(x_pos) - 0.5,
    xmax = max(x_pos) + 0.5,
    .groups = "drop"
  ) %>%
  filter(!is.na(domain_label))        

SIG_BONF <- 0.05 / nrow(plot_df)

# Domain midpoints and boundary lines -------------------------------------
domain_mids <- plot_df %>%
  group_by(domain_label) %>%
  summarise(mid = mean(x_pos), .groups = "drop")

domain_bounds <- plot_df %>%
  group_by(domain_label) %>%
  summarise(
    xmin = min(x_pos) - 0.5,
    xmax = max(x_pos) + 0.5,
    .groups = "drop"
  )

# Shared theme ------------------------------------------------------------
theme_cvd <- function() {
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 13, face = "bold", colour = "#2c3e50"),
      plot.subtitle = element_text(size = 10, colour = "grey35", margin = margin(b = 8)),
      axis.title = element_text(size = 11, face = "bold", colour = "#2c3e50"),
      axis.text = element_text(colour = "#2c3e50"),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 9, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "#e0e0e0", linewidth = 0.4),
      plot.margin = margin(10, 15, 80, 12)
    )
}

# Helper: build lollipop Manhattan plot -----------------------------------
make_manhattan_lollipop <- function(data, y_var, y_label, sig_line,
                                    sig_col = "firebrick", title_suffix) {
  
  y_raw <- data[[y_var]]
  finite_max <- max(y_raw[is.finite(y_raw)], na.rm = TRUE)
  
  data <- data %>%
    mutate(
      across(all_of(y_var), ~ ifelse(is.infinite(.x), finite_max + 10, .x))
    )
  
  y_raw <- data[[y_var]]
  y_cap <- quantile(y_raw, 0.99, na.rm = TRUE)
  y_max <- y_cap * 1.15
  
  data <- data %>%
    mutate(
      capped      = .data[[y_var]] > y_cap,
      y_plot      = pmin(.data[[y_var]], y_cap),
      significant = .data[[y_var]] >= -log10(sig_line),
      label       = ifelse(capped, paste0(variable_label, "*"), variable_label)
    )
  
  sig_y <- -log10(sig_line)
  sig_y_plot <- min(sig_y, y_cap * 0.98)
  y_lbl <- y_max * 0.02
  
  label_df <- data %>%
    select(x_pos, label, significant, y_plot)
  
  ggplot(data, aes(x = x_pos, y = y_plot)) +
    
    # Alternating domain background bands
    geom_rect(
      data = domain_bounds,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = rep(c("grey95", "white"), length.out = nrow(domain_bounds)),
      alpha = 0.6
    ) +
    
    # Lollipop stems
    geom_segment(
      aes(
        x = x_pos, xend = x_pos,
        y = 0, yend = y_plot,
        colour = direction
      ),
      linewidth = 0.45,
      alpha = 0.75
    ) +
    
    # Lollipop heads
    geom_point(
      aes(
        colour = direction,
        size = ifelse(significant, 2.4, 1.5),
        alpha = ifelse(significant, 1, 0.6)
      )
    ) +
    scale_size_identity() +
    scale_alpha_identity() +
    
    scale_colour_manual(
      name = "Effect direction",
      values = c(
        "Positive OR" = "#3A6EA5",
        "Negative OR" = "#C0392B"
      )
    ) +
    
    # Significance threshold
    geom_hline(
      yintercept = sig_y_plot,
      linetype = "dashed",
      colour = sig_col,
      linewidth = 0.65
    ) +
    annotate(
      "text",
      x = max(data$x_pos) * 0.02,
      y = sig_y_plot + y_max * 0.025,
      label = paste0("p = ", sig_line),
      colour = sig_col,
      size = 3,
      hjust = 0
    ) +
    
    # Non-significant labels
    geom_text(
      data = filter(label_df, !significant),
      aes(x = x_pos, y = y_lbl, label = label),
      angle = 90,
      hjust = 0,
      vjust = 0.5,
      size = 1.9,
      colour = "grey60",
      fontface = "plain",
      inherit.aes = FALSE
    ) +
    
    # Significant labels
    geom_text(
      data = filter(label_df, significant),
      aes(x = x_pos, y = y_plot + y_max * 0.01, label = label), 
      angle = 90,
      hjust = 0,
      vjust = 0.5,
      size = 2.4,
      colour = "#2c3e50",
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    
    scale_x_continuous(
      breaks = domain_mids$mid,
      labels = domain_mids$domain_label,
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = c(0, 0),
      breaks = function(lims) pretty(c(0, lims[2]), n = 6),
      labels = function(b) ifelse(b < 0, "", b)
    ) +
    
    labs(
      title = paste("Manhattan plot: associations with CVD", title_suffix),
      subtitle = paste0(
        "Bold labels = significant at p < ", sig_line,
        " | * = value capped at 99th percentile"
      ),
      x = "Domain",
      y = y_label
    ) +
    theme_cvd()
}

# Plot 1: raw p-value -----------------------------------------------------
p_raw <- make_manhattan_lollipop(
  data = plot_df,
  y_var = "neg_log10_p",
  y_label = expression(-log[10](p["raw"])),
  sig_line = SIG_RAW,
  title_suffix = "(raw p-value)"
)

ggsave(
  filename = here("02_results", "00_figures", "02_univariate", "manhattan_raw_p.png"),
  plot = p_raw,
  width = 16,
  height = 9,
  dpi = 300,
  bg = "white"
)

# Plot 2: BH-adjusted p-value ---------------------------------------------
p_bh <- make_manhattan_lollipop(
  data = plot_df,
  y_var = "neg_log10_p_adj_bh",
  y_label = expression(-log[10](p["BH-adjusted"])),
  sig_line = SIG_BH,
  title_suffix = "(BH-adjusted p-value)"
)

ggsave(
  filename = here("02_results", "00_figures", "02_univariate", "manhattan_bh_p.png"),
  plot = p_bh,
  width = 16,
  height = 9,
  dpi = 300,
  bg = "white"
)

# Plot 3: Bonferroni-adjusted p-value -------------------------------------
p_bonf <- make_manhattan_lollipop(
  data = plot_df,
  y_var = "neg_log10_p_adj_bonf",
  y_label = expression(-log[10](p["Bonferroni"])),
  sig_line = SIG_BONF,
  title_suffix = "(Bonferroni-adjusted p-value)"
)

ggsave(
  filename = here("02_results", "00_figures", "02_univariate", "manhattan_bonf_p.png"),
  plot = p_bonf,
  width = 16,
  height = 9,
  dpi = 300,
  bg = "white"
)

# Save results table ------------------------------------------------------
write.csv(
  plot_df %>%
    select(
      variable, variable_label, domain_label,
      effect_log_or, p_value, p_adj_bh, p_adj_bonf,
      direction, neg_log10_p, neg_log10_p_adj_bh, neg_log10_p_adj_bonf
    ) %>%
    arrange(domain_label, p_value),
  file = here("02_results", "01_tables", "02_univariate", "univariate_logistic_results.csv"),
  row.names = FALSE
)

message("Results table saved to 02_results/01_tables/02_univariate")

message("Lollipop Manhattan plots saved to 02_results/00_figures/02_univariate")