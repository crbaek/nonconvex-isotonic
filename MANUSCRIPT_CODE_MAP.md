# Manuscript-code map

This note maps manuscript components to repository scripts.

## Section 4: Simulation study

Main script:

```text
scripts/01_simulation_main.R
```

This script implements the fixed-grid upper-tail quantile simulation with the product partial order over `(volatility bin, drawdown bin, quantile level)`.  It compares

1. the convex isotonic projection (`lambda = 0`),
2. a fixed quartic perturbation (`lambda = 0.05`), and
3. a validation-selected quartic perturbation.

Post-processing and plot:

```text
scripts/03_plot_simulation_figure.R
```

## Appendix robustness simulation

```text
scripts/02_simulation_robustness.R
```

This script varies the grid size and the location of the sharp change.

## Optional binary-response likelihood example

```text
scripts/04_binary_logit_appendix.R
```

This is included as an additional example of the same projection-and-pooling mechanism for a likelihood criterion.

## QQQ fixed-grid upper-tail quantile illustration

```text
scripts/05_qqq_empirical_template.R
scripts/06_qqq_paired_wilcoxon_postprocess.R
```

The empirical script is a public template.  It does not include the raw intraday QQQ data.  The original Dropbox working script used a private/local data file and is not included directly.

## Dropbox working files inspected

The cleaned scripts were prepared after inspecting the Dropbox working folder `/isotonic/program4`, including the main simulation script, robustness script, binary logit script, QQQ empirical script, and paired Wilcoxon post-processing script.  Large generated output folders and private data are not included in the GitHub bundle.
