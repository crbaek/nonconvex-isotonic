# Non-Convex Isotonic Regression on Finite Partial Orders

Python scripts for the simulation study and QQQ empirical analysis in Sections 4
and 5 of **Projection-and-Pooling for Non-Convex Isotonic Regression on Finite
Partial Orders**.

## Paper and citation

**Changryong Baek, Zhenyu Cui, Dongwoo Kim, Chihoon Lee, and Yunfan Zhu (2026).**
*Projection-and-Pooling for Non-Convex Isotonic Regression on Finite Partial
Orders.* Working paper.

[Paper draft folder](paper/) — reserved for the manuscript PDF, which will be
uploaded separately.

If you use this code in your research, please cite the paper:

```bibtex
@unpublished{baek2026projection,
  author = {Baek, Changryong and Cui, Zhenyu and Kim, Dongwoo
            and Lee, Chihoon and Zhu, Yunfan},
  title  = {Projection-and-Pooling for Non-Convex Isotonic Regression
            on Finite Partial Orders},
  year   = {2026},
  note   = {Working paper}
}
```

## Repository contents

| File | Description |
| --- | --- |
| `section4_simulation.py` | Monte Carlo simulation with three signal designs and Gaussian, t5, and t3 innovations. |
| `section5_qqq_empirical.py` | Rolling QQQ analysis with convex, persistence, and mean-reversion responses. |
| `requirements.txt` | Python dependencies: NumPy, pandas, SciPy, and Numba. |
| `paper/` | Location for the paper draft. |

The two Python scripts and the dependency list are preserved unchanged in this
packaging update. The instructions below describe the settings in those scripts.

## Installation

Create a Python virtual environment and install the dependencies. For example,
in a POSIX shell:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate the environment with `.venv\Scripts\Activate.ps1` in
PowerShell. The run examples below use POSIX environment-variable syntax.

## Paths and data

The supplied scripts retain their original `/mnt/data/` paths:

- Section 4 writes CSV files directly to `/mnt/data/`. This directory must exist
  and be writable before running the simulation.
- Section 5 reads the CSV specified by `DATA` near the top of the script and
  writes to `OUT`, initially `/mnt/data/section5_twosided_output/`.

Use an environment with these paths, or adapt the path strings to your local
machine before running. They have not been changed in this packaging update.

The QQQ analysis requires the licensed Refinitiv/Tick History five-minute price
export. The source data are not included. The script reads `Date-Time` and the
price fields `Last`, `Close Bid`, `Close Ask`, and `Close Mid Price`; its price
selection uses Last, then Close Mid Price, then the bid–ask midpoint.

## Section 4: simulation

To run 500 replications per signal–innovation combination:

```bash
SIM_B=500 OUT_PREFIX=section4_B500 python section4_simulation.py
```

A smaller run can be requested with the existing environment variables:

```bash
SIM_B=2 ONLY_SIGNAL=smooth ONLY_ERROR=gaussian OUT_PREFIX=section4_small python section4_simulation.py
```

| Environment variable | Default | Description |
| --- | --- | --- |
| `SIM_B` | `100` | Replications per signal–innovation combination. |
| `OUT_PREFIX` | `section4_python` | Prefix for output filenames. |
| `ONLY_SIGNAL` | All designs | `smooth`, `weak_transition`, or `pronounced_transition`. |
| `ONLY_ERROR` | All innovations | `gaussian`, `t5`, or `t3`. |

The supplied script uses a 12 × 12 grid, quantile levels 0.95 and 0.99, and
training/validation/test cell sizes of 60/120/200. Its persistence parameters
are `(c_H, w_H) = (0.80, 0.06)`, and its candidate strengths are
`{0, 0.005, 0.0075, 0.010, 0.0125, 0.015}`. These settings are recorded here as
implemented; packaging does not change the experiment configuration.

Outputs, under `/mnt/data/`:

- `<prefix>_raw_metrics.csv`: replication-level metrics.
- `<prefix>_paired_gains.csv`: paired percentage gains relative to the convex fit.
- `<prefix>_gain_summary.csv`: aggregated gains and selection statistics.

## Section 5: QQQ empirical analysis

After supplying the input CSV at the configured `DATA` path, run each horizon
separately:

```bash
H=30 python section5_qqq_empirical.py
H=60 python section5_qqq_empirical.py
H=120 python section5_qqq_empirical.py
```

`H` is the forward horizon in minutes; its default is `30`. The script uses
24 months for training, two months for validation, two months for testing, and
a two-month rolling step. It evaluates grids `G = 8, 10` and quantile levels
`0.95, 0.975, 0.99`, using the response-selection procedure in the supplied code.

For each horizon, the configured output directory receives:

- `raw_h<h>.csv`: fold-level evaluation metrics.
- `selection_h<h>.csv`: selected response families and strengths.
- `paired_h<h>.csv`: paired gains relative to the convex fit.
- `summary_h<h>.csv`: aggregated metrics and selection frequencies.

Runs with the same output filenames overwrite the existing CSV files. Use
distinct simulation prefixes or preserve earlier empirical outputs separately
when retaining multiple runs.

## Data availability

The empirical source data are licensed and are not redistributed with this
repository. Empirical results require access to the corresponding source data.
The `paper/` folder is reserved for the separately uploaded manuscript draft.
