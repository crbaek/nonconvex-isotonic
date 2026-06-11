rm(list = ls())

## Optional binary-response logit likelihood example.
## This compact script illustrates the same projection-and-pooling mechanism for
## a binomial logit criterion with a shared quartic component on the logit scale.

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1L) args[1] else Sys.getenv("OUT_DIR", "results/binary_logit")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path("R", "isotonic_quartic_utils.R"))
require_pkg("data.table")
library(data.table)

set.seed(as.integer(Sys.getenv("SIM_SEED", "20260502")))
B <- as.integer(Sys.getenv("BIN_B", "200"))
K <- as.integer(Sys.getenv("BIN_K", "12"))
L <- as.integer(Sys.getenv("BIN_L", "12"))
k0 <- as.integer(Sys.getenv("BIN_K0", "8"))
m0 <- as.integer(Sys.getenv("BIN_M0", "8"))
n_cell <- as.integer(Sys.getenv("BIN_N_CELL", "80"))
lambda_grid <- as.numeric(strsplit(Sys.getenv("BIN_LAMBDA_GRID", "0,0.5,1.0,1.5"), ",")[[1]])
eta_L <- -4; eta_U <- 4
scale_eta <- function(x) (x - eta_L) / (eta_U - eta_L)
inv_logit <- function(x) 1 / (1 + exp(-x))

make_eta_true <- function() {
  eta <- matrix(0, K, L)
  for (k in seq_len(K)) for (m in seq_len(L)) {
    eta[k, m] <- -2.8 + 1.0 * k / K + 0.8 * m / L +
      1.2 * as.numeric(k >= k0) + 0.9 * as.numeric(m >= m0) +
      1.0 * as.numeric(k >= k0 && m >= m0)
  }
  eta
}

logit_quartic_transform_one <- function(u, lambda) {
  ## Maximize u*x - log(1+exp(x)) + lambda*((x-eta_L)/(eta_U-eta_L))^4 over [eta_L,eta_U].
  f <- function(x) u * x - log1p(exp(x)) + lambda * scale_eta(x)^4
  opt <- optimize(f, interval = c(eta_L, eta_U), maximum = TRUE)
  ## Include endpoints to avoid missing boundary maxima.
  cand <- c(eta_L, eta_U, opt$maximum)
  cand[which.max(vapply(cand, f, numeric(1)))]
}
logit_quartic_transform <- function(u, lambda) vapply(u, logit_quartic_transform_one, numeric(1), lambda = lambda)

run_one <- function(b) {
  set.seed(20260502 + b)
  eta0 <- make_eta_true()
  p0 <- inv_logit(eta0)
  y <- matrix(rbinom(K * L, size = n_cell, prob = as.vector(p0)), K, L)
  phat <- pmin(pmax(y / n_cell, 1 / (2 * n_cell)), 1 - 1 / (2 * n_cell))
  g <- qlogis(phat)
  edges <- make_product_order_edges(K, L, 1L)
  g_proj <- solve_isotonic_projection(as.vector(g), rep(n_cell, K * L), edges, solver = Sys.getenv("SIM_SOLVER", "auto"))
  evals <- list()
  for (lam in lambda_grid) {
    eta_hat <- matrix(logit_quartic_transform(g_proj, lam), K, L)
    p_hat <- inv_logit(eta_hat)
    evals[[as.character(lam)]] <- data.table(
      rep = b,
      lambda = lam,
      grid_mse = mean((eta_hat - eta0)^2),
      stress_mse = mean((eta_hat[k0:K, m0:L] - eta0[k0:K, m0:L])^2),
      brier_proxy = mean((p_hat - p0)^2),
      log_loss = -mean(p0 * log(p_hat) + (1 - p0) * log(1 - p_hat))
    )
  }
  rbindlist(evals)
}

res <- rbindlist(lapply(seq_len(B), run_one))
base <- res[lambda == 0]
summary <- res[, .(
  grid_mse = mean(grid_mse),
  stress_mse = mean(stress_mse),
  brier_proxy = mean(brier_proxy),
  log_loss = mean(log_loss)
), by = lambda]
fwrite(res, file.path(OUT_DIR, "binary_logit_by_replication.csv"))
fwrite(summary, file.path(OUT_DIR, "binary_logit_summary.csv"))
cat("Done. Wrote outputs to: ", OUT_DIR, "\n", sep = "")
