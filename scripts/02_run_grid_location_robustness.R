rm(list = ls())

## Robustness simulation across grid sizes and sharp-change locations.
## Usage:
##   Rscript scripts/02_run_grid_location_robustness.R results/sim_robustness

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1L) args[1] else Sys.getenv("OUT_DIR", "results/sim_robustness")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path("R", "isotonic_quartic_utils.R"))
require_pkg("foreach")
require_pkg("doParallel")
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(data.table) })

B <- as.integer(Sys.getenv("ROBUST_B", Sys.getenv("SIM_B", "100")))
N_CORES <- as.integer(Sys.getenv("SIM_CORES", "2"))
base_seed <- as.integer(Sys.getenv("SIM_SEED", "20260502"))
G_grid <- as.integer(strsplit(Sys.getenv("ROBUST_G", "8,12,16"), ",")[[1]])
rho_grid <- as.numeric(strsplit(Sys.getenv("ROBUST_RHO", "0.60,0.67,0.75"), ",")[[1]])

tasks <- data.table::CJ(G = G_grid, rho = rho_grid, rep = seq_len(B))
tasks[, k0 := ceiling(rho * G)]
tasks[, m0 := ceiling(rho * G)]
tasks[, task_id := .I]

data.table::fwrite(tasks, file.path(OUT_DIR, "robustness_tasks.csv"))

n_workers <- max(1L, min(N_CORES, nrow(tasks)))
cl <- parallel::makeCluster(n_workers)
doParallel::registerDoParallel(cl)
on.exit({ parallel::stopCluster(cl) }, add = TRUE)

cat("Running robustness simulation with tasks=", nrow(tasks), ", workers=", n_workers, "\n", sep = "")
res <- foreach::foreach(
  ii = seq_len(nrow(tasks)),
  .combine = data.table::rbindlist,
  .packages = c("Matrix", "data.table")
) %dopar% {
  source(file.path("R", "isotonic_quartic_utils.R"))
  task <- tasks[ii]
  cfg <- list(
    B = B,
    n_cores = 1L,
    base_seed = base_seed + 100000L * task$task_id,
    K = task$G,
    L = task$G,
    alphas = c(0.95, 0.99),
    nu_t = as.numeric(Sys.getenv("SIM_T_DF", "5")),
    n_train = as.integer(Sys.getenv("SIM_N_TRAIN", "60")),
    n_val = as.integer(Sys.getenv("SIM_N_VAL", "120")),
    n_test = as.integer(Sys.getenv("SIM_N_TEST", "200")),
    k0 = task$k0,
    m0 = task$m0,
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
  tmp <- simulate_quartic_replication(task$rep, cfg)
  tmp[, `:=`(G = task$G, rho = task$rho, k0 = task$k0, m0 = task$m0)]
  tmp
}

data.table::fwrite(res, file.path(OUT_DIR, "robustness_by_replication.csv"))

gains <- res[, gain_summary(.SD), by = .(G, rho)]
data.table::fwrite(gains, file.path(OUT_DIR, "robustness_gain_summary.csv"))

selected <- gains[method == "quartic_selected"]
wide <- data.table::dcast(selected, G + rho ~ metric, value.var = "mean_gain_pct")
data.table::fwrite(wide, file.path(OUT_DIR, "robustness_summary.csv"))
cat("Done. Wrote outputs to: ", OUT_DIR, "\n", sep = "")
