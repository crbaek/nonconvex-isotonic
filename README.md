# Projection and pooling for non-convex isotonic regression

This repository contains cleaned R scripts for reproducing the simulation and empirical calculations in the manuscript

**Projection and pooling for non-convex isotonic regression on finite partial orders**.

The code implements the projection-and-pooling estimator for a shared quartic non-convex component.  The order-constrained step is an ordinary weighted isotonic projection on the same finite partial order; the non-convex step is a scalar maximization within each level set of the projection.

## Repository layout

```text
R/
  isotonic_quartic_utils.R          Common routines for partial-order isotonic projection, quartic transforms, simulation DGPs, and evaluation.

scripts/
  01_simulation_main.R              Main fixed-grid quantile simulation used for the manuscript simulation table and figure.
  02_simulation_robustness.R        Grid-size and sharp-change-location robustness simulation.
  03_plot_simulation_figure.R       Plot mean percentage reductions relative to the convex projection.
  04_binary_logit_appendix.R        Optional binary-response logit likelihood example.
  05_qqq_empirical_template.R       Template for the QQQ fixed-grid upper-tail quantile illustration; requires intraday price data.
  06_qqq_paired_wilcoxon_postprocess.R
                                      Fold-level paired Wilcoxon post-processing for QQQ results.

data/
  README.md                         Data instructions.  The large intraday QQQ file is not included.

results/
  Runtime outputs.  Ignored by git except for `.gitkeep`.
```

## R dependencies

Install the required packages:

```r
install.packages(c(
  "Matrix", "data.table", "foreach", "doParallel", "ggplot2", "scales"
))
```

For the partial-order weighted isotonic projection, the scripts use **Gurobi** if the R package `gurobi` is available.  Otherwise they try the open-source `osqp` package.  Install one of these solvers:

```r
# optional open-source QP solver
install.packages("osqp")

# optional commercial solver; install gurobi separately and then its R package
# see Gurobi documentation for R installation
```

The main simulation uses only moderate-size quadratic programs.  Gurobi is faster and was used in the working Dropbox scripts.  The `osqp` fallback is included for GitHub reproducibility.

## Reproduce the main simulation

From the repository root:

```bash
Rscript scripts/01_simulation_main.R results/sim_main
Rscript scripts/03_plot_simulation_figure.R results/sim_main
```

Default settings match the manuscript simulation design: `B=100`, `K=L=12`, quantile levels `0.95` and `0.99`, training/validation/test cell sizes `60/120/200`, and the validation lambda grid

```text
0, 0.005, 0.01, 0.02, 0.05, 0.075, 0.10, 0.125, 0.15
```

You can override the main settings by environment variables:

```bash
SIM_B=500 SIM_CORES=24 SIM_SEED=20260502 Rscript scripts/01_simulation_main.R results/sim_main_B500
```

Important output files:

```text
results/sim_main/sim_main_settings.csv
results/sim_main/sim_main_by_replication.csv
results/sim_main/sim_main_summary.csv
results/sim_main/sim_gain_summary.csv
results/sim_main/sim_gain_plot.pdf
```

## Reproduce the grid/location robustness check

```bash
ROBUST_B=100 SIM_CORES=24 Rscript scripts/02_simulation_robustness.R results/sim_robustness
```

The script runs `G in {8,12,16}` and sharp-change locations `rho in {0.60,0.67,0.75}` where `k0=m0=ceiling(rho*G)`.  Output files:

```text
results/sim_robustness/robustness_by_replication.csv
results/sim_robustness/robustness_summary.csv
```

## Optional binary logit appendix example

```bash
Rscript scripts/04_binary_logit_appendix.R results/binary_logit
```

This is a compact likelihood example showing the same projection-and-pooling mechanism for a fixed-grid binomial logit criterion with a shared quartic component.  It is not required for the main simulation section.

## QQQ empirical template

The QQQ empirical illustration requires intraday prices and is not fully reproducible from the repository unless the user supplies the data file.  Place a CSV file at

```text
data/intraday_5m_qqq.csv
```

with columns

```text
datetime, price
```

Then run

```bash
Rscript scripts/05_qqq_empirical_template.R data/intraday_5m_qqq.csv results/qqq
Rscript scripts/06_qqq_paired_wilcoxon_postprocess.R results/qqq
```

The template implements the same fixed-grid logic as the manuscript: realized-volatility and drawdown bins, grid-cell empirical upper-tail loss quantiles, convex isotonic projection, validation-selected quartic post-projection transform, and fold-level paired diagnostics.

## Notes on reproducibility

Random seeds are fixed by default.  Parallel execution is controlled by `SIM_CORES`.  For strict reproducibility across machines, use one worker (`SIM_CORES=1`) and the same QP solver.

The output values can vary slightly across QP solvers because the projection is computed numerically.  The manuscript comparisons are based on percentage reductions relative to the convex projection, so small numerical differences should not affect the qualitative conclusions.
