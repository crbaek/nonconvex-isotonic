# Manuscript-code map

This note maps manuscript components to repository scripts for

> Baek, C., Cui, Z., Kim, D. W., Lee, C., and Zhu, Y. (2026). *Projection and Pooling for Non-Convex Isotonic Regression on Finite Partial Orders*. Manuscript.

## Section 4: Simulation study

Main simulation script:

```text
scripts/01_run_main_simulation.R
```

This script implements the fixed-grid upper-tail quantile simulation with the product partial order over `(volatility bin, drawdown bin, quantile level)`.  It compares

1. the convex isotonic projection (`lambda = 0`),
2. a fixed quartic perturbation (`lambda = 0.05`), and
3. a validation-selected quartic perturbation.

Post-processing and plotting:

```text
scripts/03_plot_simulation_gains.R
```

## Appendix robustness check

```text
scripts/02_run_grid_location_robustness.R
```

This script varies the grid size and the location of the sharp change.  It corresponds to the grid/location robustness table in the appendix.

## Common functions

```text
R/isotonic_quartic_utils.R
```

This file contains the common routines for

- generating the fixed-grid quantile target,
- constructing the product-order constraint matrix,
- solving the weighted isotonic projection,
- applying the scalar quartic transformation,
- selecting the validation lambda, and
- computing MSE, stress-region MSE, jump-error, and quantile-loss diagnostics.

## Scope

This repository is limited to the simulation code used for the manuscript simulation section and the associated robustness check.
