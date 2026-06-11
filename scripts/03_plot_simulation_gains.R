rm(list = ls())

## Plot simulation gains relative to the convex isotonic projection.
## Usage:
##   Rscript scripts/03_plot_simulation_gains.R results/sim_main

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1L) args[1] else Sys.getenv("OUT_DIR", "results/sim_main")

if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

library(data.table)
library(ggplot2)

gain_path <- file.path(OUT_DIR, "sim_gain_summary.csv")
if (!file.exists(gain_path)) stop("Missing file: ", gain_path)
g <- data.table::fread(gain_path)

metric_map <- c(
  grid_mse = "Grid MSE",
  stress_mse = "Stress region MSE",
  jump_err_095 = "Jump error\n0.95 tail",
  jump_err_099 = "Jump error\n0.99 tail",
  ql_950 = "Quantile loss\n0.95 tail",
  ql_990 = "Quantile loss\n0.99 tail"
)
g <- g[metric %in% names(metric_map)]
g[, metric_label := factor(metric_map[metric], levels = metric_map)]
g[, method_label := fifelse(method == "quartic_fixed", "Fixed quartic\n(lambda=0.05)", "Selected quartic")]

gg <- ggplot(g, aes(x = method_label, y = mean_gain_pct)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_gain_pct - 1.96 * se_gain_pct,
                    ymax = mean_gain_pct + 1.96 * se_gain_pct), width = 0.16) +
  geom_text(aes(label = sprintf("%.1f", mean_gain_pct)), vjust = -0.35, size = 3) +
  facet_wrap(~ metric_label, scales = "free_y", nrow = 1) +
  labs(
    title = "Simulation gains from the quartic projection-and-pooling step",
    subtitle = "Positive values favor the quartic estimator; error bars show approximate 95% Monte Carlo intervals",
    x = NULL,
    y = "Mean percentage reduction relative to convex (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8),
    strip.background = element_rect(fill = "grey95"),
    plot.title.position = "plot"
  )

pdf_file <- file.path(OUT_DIR, "sim_gain_plot.pdf")
png_file <- file.path(OUT_DIR, "sim_gain_plot.png")
ggsave(pdf_file, gg, width = 11, height = 4.2)
ggsave(png_file, gg, width = 11, height = 4.2, dpi = 300)
cat("Wrote: ", pdf_file, "\n", sep = "")
cat("Wrote: ", png_file, "\n", sep = "")
