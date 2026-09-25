# DStMM simulation project

This project implements the simulation study for **Deep Skew-t Mixture Models (DStMM)** and compares:

1. **DGMM** — the uploaded Deep Gaussian Mixture Model implementation;
2. **RDMM** — the uploaded Robust Deep Mixture Model implementation;
3. **DStMM** — the proposed generalized-hyperbolic skew-t deep mixture model.

The DStMM implementation follows the manuscript's pathway formulation: a single inverse-gamma scale is shared by all transitions along an observation-path pair, and active layers contain an `H * delta` mean shift. The default active set is `A = {1}` and the default Monte Carlo size is `M = 1`, giving the stochastic-EM version described in the paper.

## Project structure

```text
R/
  dgmm/                   original DGMM source files from the uploaded archive
  robustdeepgmm.R         RDMM source from the uploaded archive
  dstmm.R                 proposed DStMM implementation
  simulation_design.R     manuscript data-generating mechanism and grids
  metrics_and_fitting.R   common fitting, ARI/MR, nu/delta recovery
  simulation_runner.R     resumable simulation runner
  load_all.R              project loader
scripts/
  00_install_packages.R
  01_smoke_test.R
  02_unit_tests.R
  10_experiment1_main.R
  11_experiment1_robustness.R
  20_experiment2.R
  30_summarise_results.R
  40_make_figures.R
  50_export_latex_tables.R
  99_run_all.R
results/
figures/
```

## Required R packages

```r
Rscript scripts/00_install_packages.R
```

Required packages are `mvtnorm`, `GIGrvg`, and `ggplot2`.

## First run: smoke test

Run this before the full Monte Carlo study:

```bash
Rscript scripts/01_smoke_test.R
```

Before the smoke fit, you can also run deterministic formula checks:

```bash
Rscript scripts/02_unit_tests.R
```

The smoke test fits all three models to a small simulated data set and writes:

```text
results/smoke_test_results.csv
```

## Simulation design

### Experiment 1 — when does skewness become necessary?

- `nu = 6`
- `kappa = 0, 0.5, 1, 1.5, 2`
- `n = 1000`
- balanced layer proportions `(0.5, 0.5)`
- full run: `R = 1000` replications

```bash
Rscript scripts/10_experiment1_main.R
```

### Lightweight robustness checks

The revised design keeps the paper compact and adds only:

- sample-size check: `n = 500`, balanced proportions;
- imbalance check: `n = 1000`, layer proportions `(0.7, 0.3)`.

The same `kappa` grid is used.

```bash
Rscript scripts/11_experiment1_robustness.R
```

### Experiment 2 — tail weight × skewness

- `nu = 5, 8, 15`
- `kappa = 0, 1, 2`
- `n = 1000`
- balanced proportions
- full run: `R = 1000`

```bash
Rscript scripts/20_experiment2.R
```

## Fast development mode

Each experiment supports a smoke mode with two replications:

```bash
RUN_MODE=smoke Rscript scripts/10_experiment1_main.R
RUN_MODE=smoke Rscript scripts/11_experiment1_robustness.R
RUN_MODE=smoke Rscript scripts/20_experiment2.R
```

## Splitting a full run across jobs

The full `R = 1000` study is expensive. The scripts support splitting by condition and replication range.

Example:

```bash
CONDITION_INDEX=1 REP_START=1 REP_END=100 Rscript scripts/10_experiment1_main.R
CONDITION_INDEX=1 REP_START=101 REP_END=200 Rscript scripts/10_experiment1_main.R
```

The CSV runner is resumable: completed condition/replicate/model triples are not recomputed.

`MC_SIZE` controls the Monte Carlo size in DStMM:

```bash
MC_SIZE=5 RUN_MODE=smoke Rscript scripts/10_experiment1_main.R
```

The manuscript baseline should use `MC_SIZE=1` unless a larger Monte Carlo size is deliberately studied.

## Model fitting convention

For each replicate, the code fits models in this order:

1. DGMM;
2. RDMM;
3. DStMM.

DStMM uses the fitted RDMM as a warm start whenever the RDMM fit succeeds, exactly matching the practical workflow suggested in the manuscript: first obtain a stable symmetric heavy-tailed fit, then release the skewness parameters. If RDMM fails, DStMM falls back to its own factor-analytic initialization.

All three models use the same known architecture:

```text
p = 20
(r1, r2) = (5, 2)
(K1, K2) = (2, 2)
```

No architecture selection is mixed into the component-shape comparison.

## Raw outputs

The simulation scripts write:

```text
results/experiment1_raw.csv
results/experiment1_robustness_raw.csv
results/experiment2_raw.csv
```

Each model/replication row contains:

- first-layer ARI and misclassification rate;
- complete-path ARI and misclassification rate;
- estimated degrees of freedom where applicable;
- DStMM `delta` RMSE and pathway `alpha` RMSE;
- maximized log-likelihood;
- iteration/convergence information;
- elapsed time;
- failure messages if a fit fails.

## Result tables

After simulations finish:

```bash
Rscript scripts/30_summarise_results.R
```

This produces:

```text
results/model_summary_long.csv
results/clustering_and_fit_summary.csv
results/nu_recovery_summary.csv
results/skewness_recovery_summary.csv
```


To export compact LaTeX result tables after the raw simulations are available:

```bash
Rscript scripts/50_export_latex_tables.R
```

This writes `results/table_experiment1.tex` and `results/table_experiment2.tex`.

## Figures

After result summarization:

```bash
Rscript scripts/40_make_figures.R
```

The script writes both PNG and PDF versions of manuscript-ready figures, including:

- Experiment 1 first-layer ARI versus `kappa`;
- Experiment 1 pathway ARI versus `kappa`;
- DStMM skewness-vector RMSE;
- robustness-check ARI gains;
- Experiment 2 DStMM-minus-RDMM ARI heat map;
- Experiment 2 DStMM-minus-DGMM ARI heat map;
- computational-time comparison.

## Important implementation notes

- The code uses the paper's inverse-gamma convention `H ~ IG(nu/2, nu/2)`.
- For nonzero skewness, `H | y, s` is sampled from the GIG distribution using `GIGrvg::rgig`.
- When the pathway skewness quadratic form is numerically zero, the code switches exactly to the Student-t/inverse-gamma limit rather than evaluating a nearly singular Bessel expression.
- `E(H^-1)` and `E(log H)` are evaluated analytically from GIG moments for the degrees-of-freedom update.
- The same simulated `H` is reused through every layer of a pathway. Using independent layer-specific scales would be a different model.
- `Psi` is kept diagonal, matching the baseline manuscript specification and the uploaded RDMM implementation.

## Reproducibility note

No numerical Monte Carlo results are included in this archive unless the scripts have actually been run. The current project contains the complete code for generating, fitting, summarizing, and plotting the experiments; it does **not** fabricate placeholder performance numbers.

## DGMM visibility in the boxplots

The uploaded Gaussian DGMM implementation calls `is.positive.definite()` and
`make.positive.definite()` from the **corpcor** package without a namespace.
This project now installs/loads `corpcor` explicitly. Older raw CSV files that
contain `Model == "DGMM"` with `Success == FALSE` can be repaired without
rerunning RDMM/DStMM by running `scripts/05_repair_dgmm_results.R`. The current
`main.R` performs that repair automatically before resuming simulations.
The figure script also stops with a diagnostic error if a simulation condition
has no successful DGMM/RDMM/DStMM result, rather than silently dropping a model.

## DGMM two-layer bug fix (2026-08-19)

`R/dgmm/deep.sem.alg.2.R` was patched so the first SEM iteration is allowed when
`ratio` is initialized to `Inf`. The previous finite-ratio guard skipped the
loop entirely and later caused `object 'lik' not found`. The final likelihood
is also now taken safely from `tail(likelihood, 1)`.
