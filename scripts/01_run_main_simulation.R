rm(list = ls())

## Main fixed-grid quantile simulation for the manuscript.
## Usage:
##   Rscript scripts/01_run_main_simulation.R results/sim_main

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1L) args[1] else Sys.getenv("OUT_DIR", "results/sim_main")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path("R", "isotonic_quartic_utils.R"))
require_pkg("foreach")
require_pkg("doParallel")
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(data.table) })

cfg <- list(
  B = as.integer(Sys.getenv("SIM_B", "100")),
  n_cores = as.integer(Sys.getenv("SIM_CORES", "2")),
  base_seed = as.integer(Sys.getenv("SIM_SEED", "20260502")),
  K = as.integer(Sys.getenv("SIM_K", "12")),
  L = as.integer(Sys.getenv("SIM_L", "12")),
  alphas = c(0.95, 0.99),
  nu_t = as.numeric(Sys.getenv("SIM_T_DF", "5")),
  n_train = as.integer(Sys.getenv("SIM_N_TRAIN", "60")),
  n_val = as.integer(Sys.getenv("SIM_N_VAL", "120")),
  n_test = as.integer(Sys.getenv("SIM_N_TEST", "200")),
  k0 = as.integer(Sys.getenv("SIM_K0", "8")),
  m0 = as.integer(Sys.getenv("SIM_M0", "8")),
  shock_size = as.numeric(Sys.getenv("SIM_SHOCK_SIZE", "4.0")),
  pi_base = as.numeric(Sys.getenv("SIM_PI_BASE", "0.002")),
  pi_near = as.numeric(Sys.getenv("SIM_PI_NEAR", "0.008")),
  pi_deep = as.numeric(Sys.getenv("SIM_PI_DEEP", "0.020")),
  fixed_lambda = as.numeric(Sys.getenv("SIM_FIXED_LAMBDA", "0.05")),
  lambda_grid = as.numeric(strsplit(Sys.getenv("SIM_LAMBDA_GRID", "0,0.005,0.01,0.02,0.05,0.075,0.10,0.125,0.15"), ",")[[1]]),
  boundary_tol = as.numeric(Sys.getenv("SIM_BOUNDARY_TOL", "1e-8")),
  validation_tol_frac = as.numeric(Sys.getenv("SIM_VALIDATION_TOL_FRAC", "0.005")),
  min_validation_gain_frac = as.numeric(Sys.getenv("SIM_MIN_VALIDATION_GAIN_FRAC", "0.005")),
  x_upper_mult = as.numeric(Sys.getenv("SIM_X_UPPER_MULT", "1.05")),
  solver = Sys.getenv("SIM_SOLVER", "auto"),
  threads_per_worker = as.integer(Sys.getenv("GUROBI_THREADS_PER_WORKER", "1"))
)

settings <- data.table::data.table(name = names(cfg), value = vapply(cfg, function(x) paste(x, collapse = ","), character(1)))
data.table::fwrite(settings, file.path(OUT_DIR, "sim_main_settings.csv"))

n_workers <- max(1L, min(cfg$n_cores, cfg$B))
cl <- parallel::makeCluster(n_workers)
doParallel::registerDoParallel(cl)
on.exit({ parallel::stopCluster(cl) }, add = TRUE)

cat("Running main simulation with B=", cfg$B, ", workers=", n_workers, ", solver=", cfg$solver, "\n", sep = "")
res <- foreach::foreach(
  b = seq_len(cfg$B),
  .combine = data.table::rbindlist,
  .packages = c("Matrix", "data.table")
) %dopar% {
  source(file.path("R", "isotonic_quartic_utils.R"))
  simulate_quartic_replication(b, cfg)
}

data.table::fwrite(res, file.path(OUT_DIR, "sim_main_by_replication.csv"))
summary <- summarize_simulation(res)
gains <- gain_summary(res)
data.table::fwrite(summary, file.path(OUT_DIR, "sim_main_summary.csv"))
data.table::fwrite(gains, file.path(OUT_DIR, "sim_gain_summary.csv"))

lambda_freq <- res[method == "convex", .N, by = selected_lambda][order(selected_lambda)]
data.table::fwrite(lambda_freq, file.path(OUT_DIR, "selected_lambda_frequency.csv"))

cat("Done. Wrote outputs to: ", OUT_DIR, "\n", sep = "")
