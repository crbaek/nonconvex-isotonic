rm(list = ls())

## Fold-level paired Wilcoxon post-processing for QQQ results.
## Usage:
##   Rscript scripts/06_qqq_paired_wilcoxon_postprocess.R results/qqq

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1L) args[1] else Sys.getenv("OUT_DIR", "results/qqq")

if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
library(data.table)

surface_path <- file.path(OUT_DIR, "surface_eval_by_fold.csv")
if (!file.exists(surface_path)) stop("Missing file: ", surface_path)
surface_eval <- fread(surface_path)

safe_wilcox <- function(d) {
  d <- d[is.finite(d)]
  if (length(d) == 0L) return(NA_real_)
  if (all(abs(d) < .Machine$double.eps^0.5)) return(1.0)
  suppressWarnings(stats::wilcox.test(d, alternative = "greater", exact = FALSE)$p.value)
}

summarize_diff <- function(dt, by_cols, metric_cols, metric_type) {
  out <- list(); rr <- 1L
  keys <- unique(dt[, ..by_cols])
  for (ii in seq_len(nrow(keys))) {
    sub <- dt
    for (cc in by_cols) sub <- sub[get(cc) == keys[[cc]][ii]]
    form <- stats::as.formula(paste(paste(c(by_cols, "fold"), collapse = "+"), "~ method"))
    wide <- data.table::dcast(sub, form, value.var = metric_cols)
    for (mm in metric_cols) {
      conv_col <- paste(mm, "convex", sep = "_")
      quart_col <- paste(mm, "quartic_selected", sep = "_")
      if (!(conv_col %in% names(wide)) || !(quart_col %in% names(wide))) next
      d <- wide[[conv_col]] - wide[[quart_col]]
      pct <- 100 * d / wide[[conv_col]]
      out[[rr]] <- data.table(
        metric_type = metric_type,
        metric = mm,
        keys[ii],
        n_folds = sum(is.finite(d)),
        n_improved = sum(d > 0, na.rm = TRUE),
        n_equal = sum(abs(d) <= .Machine$double.eps^0.5, na.rm = TRUE),
        n_worse = sum(d < 0, na.rm = TRUE),
        mean_convex = mean(wide[[conv_col]], na.rm = TRUE),
        mean_quartic = mean(wide[[quart_col]], na.rm = TRUE),
        mean_difference = mean(d, na.rm = TRUE),
        mean_reduction_pct = mean(pct, na.rm = TRUE),
        median_reduction_pct = stats::median(pct, na.rm = TRUE),
        wilcoxon_p_one_sided = safe_wilcox(d)
      )
      rr <- rr + 1L
    }
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

surf <- surface_eval[method %in% c("convex", "quartic_selected")]
surf_tests <- summarize_diff(
  surf,
  by_cols = c("horizon_minutes", "grid"),
  metric_cols = intersect(c("surf_mse", "stress_surf_mse", "jump_err"), names(surf)),
  metric_type = "surface"
)
fwrite(surf_tests, file.path(OUT_DIR, "paired_wilcoxon_surface.csv"))
cat("Done. Wrote: ", file.path(OUT_DIR, "paired_wilcoxon_surface.csv"), "\n", sep = "")
