# Convenience launcher. For a quick check run:
#   RUN_MODE=smoke Rscript scripts/99_run_all.R
# Full R=1000 simulations are computationally expensive and are better split by
# CONDITION_INDEX / REP_START / REP_END on a cluster.
root <- normalizePath(if (file.exists("R/load_all.R")) "." else "..")
source(file.path(root, "scripts", "10_experiment1_main.R"), local = new.env(parent = globalenv()))
source(file.path(root, "scripts", "11_experiment1_robustness.R"), local = new.env(parent = globalenv()))
source(file.path(root, "scripts", "20_experiment2.R"), local = new.env(parent = globalenv()))
source(file.path(root, "scripts", "30_summarise_results.R"), local = new.env(parent = globalenv()))
source(file.path(root, "scripts", "40_make_figures.R"), local = new.env(parent = globalenv()))

source(file.path(root, "scripts", "50_export_latex_tables.R"), local = new.env(parent = globalenv()))
