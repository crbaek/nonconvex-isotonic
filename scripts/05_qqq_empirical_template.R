rm(list = ls())

## Template for the QQQ fixed-grid upper-tail quantile illustration.
## Requires a CSV with columns: datetime, price.
## Usage:
##   Rscript scripts/05_qqq_empirical_template.R data/intraday_5m_qqq.csv results/qqq

args <- commandArgs(trailingOnly = TRUE)
DATA_FILE <- if (length(args) >= 1L) args[1] else Sys.getenv("QQQ_DATA", "data/intraday_5m_qqq.csv")
OUT_DIR <- if (length(args) >= 2L) args[2] else Sys.getenv("OUT_DIR", "results/qqq")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path("R", "isotonic_quartic_utils.R"))
require_pkg("data.table")
library(data.table)

if (!file.exists(DATA_FILE)) {
  stop("Missing data file: ", DATA_FILE, "\nExpected columns: datetime, price. See data/README.md.")
}

dt <- fread(DATA_FILE)
if (!all(c("datetime", "price") %in% names(dt))) {
  stop("DATA_FILE must contain columns named datetime and price.")
}
dt[, datetime := as.POSIXct(datetime, tz = "America/New_York")]
setorder(dt, datetime)
dt <- dt[is.finite(price) & price > 0]
dt[, logp := log(price)]
dt[, ret5 := logp - shift(logp)]

## Keep within-day calculations only.
dt[, date := as.Date(datetime, tz = "America/New_York")]
for (h in c(30L, 60L, 120L)) {
  steps <- h / 5L
  dt[, paste0("loss_", h) := -(shift(logp, n = steps, type = "lead") - logp), by = date]
}
dt[, rv_1h := sqrt(frollsum(ret5^2, n = 12, align = "right")), by = date]
## Drawdown over the previous hour.
dt[, dd_1h := {
  x <- logp
  out <- rep(NA_real_, .N)
  for (i in seq_along(x)) {
    lo <- max(1L, i - 12L)
    w <- x[lo:i]
    out[i] <- max(outer(w, w, `-`), na.rm = TRUE)
  }
  out
}, by = date]

grid_vec <- as.integer(strsplit(Sys.getenv("QQQ_GRID", "8,10"), ",")[[1]])
horizons <- as.integer(strsplit(Sys.getenv("QQQ_H", "30,60,120"), ",")[[1]])
alphas <- as.numeric(strsplit(Sys.getenv("QQQ_ALPHAS", "0.95,0.975,0.99"), ",")[[1]])
lambda_grid <- as.numeric(strsplit(Sys.getenv("QQQ_LAMBDA_GRID", "0,0.005,0.01,0.02,0.05,0.075,0.10,0.125,0.15"), ",")[[1]])

## Simple rolling folds: 12-month train, 2-month validation, 2-month test.
month_id <- as.integer(format(dt$date, "%Y")) * 12L + as.integer(format(dt$date, "%m"))
months <- sort(unique(month_id))
folds <- data.table()
fold <- 1L
for (start_i in seq(1L, length(months) - 16L + 1L, by = 2L)) {
  train_m <- months[start_i:(start_i + 11L)]
  val_m <- months[(start_i + 12L):(start_i + 13L)]
  test_m <- months[(start_i + 14L):(start_i + 15L)]
  folds <- rbind(folds, data.table(fold = fold, train_start = min(train_m), train_end = max(train_m),
                                   val_start = min(val_m), val_end = max(val_m),
                                   test_start = min(test_m), test_end = max(test_m)))
  fold <- fold + 1L
}
fwrite(folds, file.path(OUT_DIR, "folds.csv"))

assign_grid <- function(x, breaks) pmin(pmax(findInterval(x, breaks, all.inside = TRUE), 1L), length(breaks) - 1L)
cell_quantiles <- function(d, G, h, breaks1, breaks2) {
  d[, c1 := assign_grid(rv_1h, breaks1)]
  d[, c2 := assign_grid(dd_1h, breaks2)]
  loss_col <- paste0("loss_", h)
  out <- CJ(c1 = seq_len(G), c2 = seq_len(G), alpha = alphas)
  vals <- d[, .(q = as.numeric(quantile(get(loss_col), probs = alpha, type = 8, names = FALSE)),
                n = .N), by = .(c1, c2, alpha)]
  out <- merge(out, vals, by = c("c1", "c2", "alpha"), all.x = TRUE)
  out[is.na(q), `:=`(q = median(d[[loss_col]], na.rm = TRUE), n = 0L)]
  out
}

surface_eval <- list(); rr <- 1L
for (h in horizons) for (G in grid_vec) for (ff in folds$fold) {
  f <- folds[fold == ff]
  dtrain <- dt[month_id >= f$train_start & month_id <= f$train_end & is.finite(rv_1h) & is.finite(dd_1h) & is.finite(get(paste0("loss_", h)))]
  dval <- dt[month_id >= f$val_start & month_id <= f$val_end & is.finite(rv_1h) & is.finite(dd_1h) & is.finite(get(paste0("loss_", h)))]
  dtest <- dt[month_id >= f$test_start & month_id <= f$test_end & is.finite(rv_1h) & is.finite(dd_1h) & is.finite(get(paste0("loss_", h)))]
  if (nrow(dtrain) < 1000 || nrow(dval) < 100 || nrow(dtest) < 100) next
  breaks1 <- unique(quantile(dtrain$rv_1h, probs = seq(0, 1, length.out = G + 1), na.rm = TRUE))
  breaks2 <- unique(quantile(dtrain$dd_1h, probs = seq(0, 1, length.out = G + 1), na.rm = TRUE))
  if (length(breaks1) < G + 1L || length(breaks2) < G + 1L) next
  qtrain <- cell_quantiles(copy(dtrain), G, h, breaks1, breaks2)
  qtest <- cell_quantiles(copy(dtest), G, h, breaks1, breaks2)
  Tlev <- length(alphas); K <- G; L <- G
  x_upper <- as.numeric(quantile(dtrain[[paste0("loss_", h)]], probs = 0.999, na.rm = TRUE))
  x_upper <- max(x_upper, max(qtrain$q, na.rm = TRUE), .Machine$double.eps)
  arr_train <- array(qtrain$q / x_upper, dim = c(G, G, Tlev))
  weights <- pmax(qtrain$n, 1)
  edges <- make_product_order_edges(G, G, Tlev)
  g_proj <- solve_isotonic_projection(as.vector(arr_train), weights, edges, solver = Sys.getenv("SIM_SOLVER", "auto"))
  for (method in c("convex", "quartic_selected")) {
    lambda <- if (method == "convex") 0 else 0.10
    xhat <- as.vector(quartic_transform(g_proj, lambda) * x_upper)
    qtest[, xhat := xhat]
    qtest[, e2 := (xhat - q)^2]
    stress <- qtest[c1 >= ceiling(0.67 * G) & c2 >= ceiling(0.67 * G)]
    surface_eval[[rr]] <- data.table(fold = ff, horizon_minutes = h, grid = G,
                                     method = method, lambda = lambda,
                                     surf_mse = weighted.mean(qtest$e2, w = pmax(qtest$n, 1)),
                                     stress_surf_mse = weighted.mean(stress$e2, w = pmax(stress$n, 1)))
    rr <- rr + 1L
  }
}

surface_eval <- rbindlist(surface_eval, fill = TRUE)
fwrite(surface_eval, file.path(OUT_DIR, "surface_eval_by_fold.csv"))
cat("Done. This is a public-data template, not an exact copy of the private Dropbox empirical script.\n")
