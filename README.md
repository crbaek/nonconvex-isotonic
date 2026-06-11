# Replication materials for non-convex isotonic regression simulations

This repository contains R scripts for reproducing the simulation results in

> Baek, C., Cui, Z., Kim, D. W., Lee, C., and Zhu, Y. (2026). *Projection and Pooling for Non-Convex Isotonic Regression on Finite Partial Orders*. Manuscript.

If you use these scripts, please cite the paper above.  A BibTeX entry is provided in [`paper_citation.bib`](paper_citation.bib), and a GitHub citation file is provided in [`CITATION.cff`](CITATION.cff).

The scripts implement the projection-and-pooling estimator for a shared quartic non-convex component.  The order-constrained step is an ordinary weighted isotonic projection on the same finite partial order.  The non-convex step is a scalar maximization within each level set of the projection.

## Recommended repository name

Use

```bash
crbaek/nonconvex-isotonic
```

A longer alternative is `crbaek/nonconvex-isotonic-regression`.  Avoid `nonconvec_isotonic`, which contains a spelling error.

## Repository layout

```text
R/
  isotonic_quartic_utils.R              Common functions for partial-order isotonic projection, quartic transforms, simulation DGPs, and evaluation.

scripts/
  01_run_main_simulation.R              Main fixed-grid quantile simulation used for the manuscript simulation table and figure.
  02_run_grid_location_robustness.R     Robustness simulation over grid sizes and sharp-change locations.
  03_plot_simulation_gains.R            Plot mean percentage reductions relative to the convex projection.

data/
  README.md                             No external data are required for these simulation scripts.

results/
  Runtime outputs.  Ignored by git except for `.gitkeep`.

paper_citation.bib                      BibTeX entry for citing the paper.
CITATION.cff                            GitHub citation metadata for the paper.
MANUSCRIPT_CODE_MAP.md                  Mapping between manuscript components and scripts.
GITHUB_UPLOAD.md                        Minimal commands for creating and pushing the repository.
```

This bundle is limited to simulation replication materials.

## R dependencies

Install the required packages:

```r
install.packages(c(
  "Matrix", "data.table", "foreach", "doParallel", "ggplot2", "scales"
))
```

The partial-order weighted isotonic projection is solved as a quadratic program.  The scripts use **Gurobi** if the R package `gurobi` is available.  Otherwise they try the open-source `osqp` package.

```r
# optional open-source QP solver
install.packages("osqp")

# optional commercial solver; install Gurobi separately and then its R package
# see the Gurobi documentation for R installation
```

Gurobi was used in the original working scripts.  The `osqp` fallback is included to make the repository easier to run on machines without a Gurobi license.

## Reproduce the main simulation

From the repository root:

```bash
Rscript scripts/01_run_main_simulation.R results/sim_main
Rscript scripts/03_plot_simulation_gains.R results/sim_main
```

Default settings match the manuscript simulation design: `B=100`, `K=L=12`, quantile levels `0.95` and `0.99`, training/validation/test cell sizes `60/120/200`, and validation lambda grid

```text
0, 0.005, 0.01, 0.02, 0.05, 0.075, 0.10, 0.125, 0.15
```

You can override major settings by environment variables:

```bash
SIM_B=500 SIM_CORES=24 SIM_SEED=20260502 Rscript scripts/01_run_main_simulation.R results/sim_main_B500
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
ROBUST_B=100 SIM_CORES=24 Rscript scripts/02_run_grid_location_robustness.R results/sim_robustness
```

The script runs `G in {8,12,16}` and sharp-change locations `rho in {0.60,0.67,0.75}`, where `k0=m0=ceiling(rho*G)`.  Important output files:

```text
results/sim_robustness/robustness_by_replication.csv
results/sim_robustness/robustness_summary.csv
```

## Notes on reproducibility

Random seeds are fixed by default.  Parallel execution is controlled by `SIM_CORES`.  For strict reproducibility across machines, use one worker (`SIM_CORES=1`) and the same QP solver.

The output values can vary slightly across QP solvers because the isotonic projection is computed numerically.  The manuscript comparisons are based on percentage reductions relative to the convex projection, so small numerical differences should not affect the qualitative conclusions.
