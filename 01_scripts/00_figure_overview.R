library(png)
library(grid)

load_png <- function(path) {
  img <- png::readPNG(path)
  rasterGrob(img, interpolate = TRUE)
}

# Load figures
figs <- list(
  load_png(here("02_results/00_figures/02_univariate/forest_univariate_by_domain.png")),
  load_png(here("02_results/00_figures/03_multivariate/prediction_odds_ratio_forest_plot.png")),
  load_png(here("02_results/00_figures/03_multivariate/roc_comparison.png")),
  load_png(here("02_results/00_figures/03_multivariate/calibration_plot.png")),
  load_png(here("02_results/00_figures/03_multivariate/mediation_heatmap.png")),
  load_png(here("02_results/00_figures/04_exploratory/FAMD and Correlation/correlation_heatmap.png"))
)

labels <- c("A", "B", "C", "D", "E", "F")
ncol   <- 2
nrow   <- ceiling(length(figs) / ncol)

png(here("02_results/00_figures/overview_panel.png"),
    width = 20, height = nrow * 8, units = "in", res = 150)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow, ncol)))

for (i in seq_along(figs)) {
  row <- ceiling(i / ncol)
  col <- ifelse(i %% ncol == 0, ncol, i %% ncol)
  
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
  grid.draw(figs[[i]])
  grid.text(labels[i], x = 0.02, y = 0.98,
            just = c("left", "top"),
            gp = gpar(fontsize = 16, fontface = "bold", col = "#2c3e50"))
  popViewport()
}

dev.off()