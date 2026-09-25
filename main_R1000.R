# ============================================================
# Main launcher for the DStMM simulation project
# Windows project directory:
# C:/Users/uqjwu15/Desktop/DStMM_simulation_project
#
# Current manuscript output policy:
#   * Figures: BOXPLOTS
#   * BIC: NOT USED / NOT EXPORTED
#
# Recommended first run:
#   RUN_MODE <- "smoke"
# After the smoke run succeeds, change to:
#   RUN_MODE <- "full"
# ============================================================

rm(list = ls())

# -----------------------------
# 1. User configuration
# -----------------------------
PROJECT_ROOT <- "C:/Users/uqjwu15/Desktop/DStMM_simulation_project"

# "smoke" = a very small number of replications for checking the whole pipeline.
# "full"  = the manuscript Monte Carlo study.
RUN_MODE <- "full"

# Replication counts. Change FULL_REPLICATIONS here if you ever want a different
# number; all three experiment scripts inherit this value from main.R.
SMOKE_REPLICATIONS <- 2L
FULL_REPLICATIONS <- 1000L

# Monte Carlo size for DStMM. The manuscript baseline uses M = 1.
MC_SIZE <- 1L

# Manuscript display settings.
FIGURE_STYLE <- "boxplot"
USE_BIC <- FALSE

# Which parts to run.
INSTALL_PACKAGES <- TRUE
RUN_UNIT_TESTS <- TRUE
RUN_MODEL_SMOKE_TEST <- TRUE
RUN_EXPERIMENT1 <- TRUE
RUN_ROBUSTNESS <- TRUE
RUN_EXPERIMENT2 <- TRUE
# Repair missing/failed DGMM rows from older runs without refitting RDMM/DStMM.
RUN_DGMM_REPAIR <- TRUE
RUN_SUMMARY <- TRUE
RUN_FIGURES <- TRUE
RUN_LATEX_TABLES <- TRUE

# Remove obsolete BIC / heatmap outputs from older runs.
CLEAN_OLD_OUTPUTS <- TRUE

# If TRUE, stop immediately when one stage fails.
STOP_ON_ERROR <- TRUE

# -----------------------------
# 2. Project setup
# -----------------------------
if (!dir.exists(PROJECT_ROOT)) {
  stop("Project directory does not exist: ", PROJECT_ROOT)
}

setwd(PROJECT_ROOT)
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)

# The experiment scripts read these settings from environment variables.
Sys.setenv(
  RUN_MODE = RUN_MODE,
  MC_SIZE = as.character(MC_SIZE),
  SMOKE_REPLICATIONS = as.character(SMOKE_REPLICATIONS),
  FULL_REPLICATIONS = as.character(FULL_REPLICATIONS)
)

# Make sure a previous split-job setting does not accidentally restrict this run.
Sys.unsetenv(c("CONDITION_INDEX", "REP_START", "REP_END"))

# Create output directories if needed.
dir.create(file.path(PROJECT_ROOT, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PROJECT_ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PROJECT_ROOT, "logs"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 3. Remove obsolete outputs
# -----------------------------
if (CLEAN_OLD_OUTPUTS) {
  obsolete_results <- c(
    "bic_selection_summary.csv"
  )

  obsolete_figures <- c(
    "exp1_bic_selection.png",
    "exp1_bic_selection.pdf",
    "exp2_dstm_vs_rdmm_heatmap.png",
    "exp2_dstm_vs_rdmm_heatmap.pdf",
    "exp2_dstm_vs_dgmm_heatmap.png",
    "exp2_dstm_vs_dgmm_heatmap.pdf",
    "exp2_ari_heatmap.png",
    "exp2_ari_heatmap.pdf",
    "exp1_robustness_ari_gain_boxplot.png",
    "exp1_robustness_ari_gain_boxplot.pdf"
  )

  old_files <- c(
    file.path(PROJECT_ROOT, "results", obsolete_results),
    file.path(PROJECT_ROOT, "figures", obsolete_figures)
  )

  old_files <- old_files[file.exists(old_files)]
  if (length(old_files) > 0L) {
    unlink(old_files)
    cat("Removed obsolete BIC/heatmap outputs:\n")
    cat(paste0("  - ", old_files), sep = "\n")
    cat("\n\n")
  }
}

cat("============================================================\n")
cat("DStMM simulation project\n")
cat("Project      : ", PROJECT_ROOT, "\n", sep = "")
cat("Mode         : ", RUN_MODE, "\n", sep = "")
cat("Replications : ", if (RUN_MODE == "full") FULL_REPLICATIONS else SMOKE_REPLICATIONS,
    " per condition\n", sep = "")
cat("MC size      : ", MC_SIZE, "\n", sep = "")
cat("Figure style : ", FIGURE_STYLE, "\n", sep = "")
cat("BIC          : disabled\n")
cat("Started      : ", format(Sys.time()), "\n", sep = "")
cat("============================================================\n\n")

# -----------------------------
# 4. Helper to run one script
# -----------------------------
run_stage <- function(label, script) {
  script_path <- file.path(PROJECT_ROOT, "scripts", script)
  if (!file.exists(script_path)) {
    stop("Cannot find script: ", script_path)
  }

  cat("\n------------------------------------------------------------\n")
  cat("START: ", label, "\n", sep = "")
  cat("Script: ", script, "\n", sep = "")
  cat("Time  : ", format(Sys.time()), "\n", sep = "")
  cat("------------------------------------------------------------\n")

  t0 <- proc.time()[[3L]]

  ans <- tryCatch({
    # Use a fresh environment because some component scripts call rm(list = ls()).
    source(script_path, local = new.env(parent = globalenv()), chdir = FALSE)
    TRUE
  }, error = function(e) {
    cat("\nERROR in ", label, ":\n", conditionMessage(e), "\n", sep = "")
    FALSE
  })

  elapsed <- proc.time()[[3L]] - t0

  if (ans) {
    cat("DONE : ", label, "\n", sep = "")
    cat("Elapsed seconds: ", round(elapsed, 2), "\n", sep = "")
  } else if (STOP_ON_ERROR) {
    stop("Pipeline stopped because stage failed: ", label)
  }

  invisible(ans)
}

# -----------------------------
# 5. Dependencies and checks
# -----------------------------
if (INSTALL_PACKAGES) {
  run_stage("Install/check required R packages", "00_install_packages.R")
}

if (RUN_UNIT_TESTS) {
  run_stage("Formula/unit tests", "02_unit_tests.R")
}

if (RUN_MODEL_SMOKE_TEST) {
  run_stage("Three-model smoke test: DGMM / RDMM / DStMM", "01_smoke_test.R")
}

# -----------------------------
# 6. Repair older DGMM failures before continuing simulations
# -----------------------------
# Older raw files may contain DGMM rows with Success=FALSE because the
# original DGMM code requires corpcor. Repair only those DGMM rows first,
# so successful RDMM/DStMM results do not need to be refitted.
if (RUN_DGMM_REPAIR) {
  run_stage("Repair missing/failed DGMM results", "05_repair_dgmm_results.R")
}

# -----------------------------
# 7. Simulation experiments
# -----------------------------
if (RUN_EXPERIMENT1) {
  run_stage("Experiment 1: skewness effect", "10_experiment1_main.R")
}

if (RUN_ROBUSTNESS) {
  run_stage("Experiment 1: lightweight robustness checks", "11_experiment1_robustness.R")
}

if (RUN_EXPERIMENT2) {
  run_stage("Experiment 2: tail weight x skewness", "20_experiment2.R")
}

# -----------------------------
# 8. Results/figures/tables
# -----------------------------
if (RUN_SUMMARY) {
  run_stage("Summarise simulation results (no BIC)", "30_summarise_results.R")
}

if (RUN_FIGURES) {
  run_stage("Generate manuscript BOXPLOTS", "40_make_figures.R")
}

if (RUN_LATEX_TABLES) {
  run_stage("Export LaTeX result tables (no BIC)", "50_export_latex_tables.R")
}

# -----------------------------
# 8. Final report
# -----------------------------
cat("\n============================================================\n")
cat("PIPELINE FINISHED\n")
cat("Finished: ", format(Sys.time()), "\n", sep = "")
cat("\nRaw and summary results:\n  ", file.path(PROJECT_ROOT, "results"), "\n", sep = "")
cat("Boxplot figures:\n  ", file.path(PROJECT_ROOT, "figures"), "\n", sep = "")
cat("BIC outputs: disabled\n")
cat("DGMM dependency: corpcor loaded; failed old DGMM rows repaired when enabled\n")
cat("============================================================\n")

# Show the most important generated result files when they exist.
important_results <- c(
  "results/smoke_test_results.csv",
  "results/experiment1_raw.csv",
  "results/experiment1_robustness_raw.csv",
  "results/experiment2_raw.csv",
  "results/clustering_and_fit_summary.csv",
  "results/nu_recovery_summary.csv",
  "results/skewness_recovery_summary.csv",
  "results/table_experiment1.tex",
  "results/table_experiment2.tex"
)

existing_results <- important_results[
  file.exists(file.path(PROJECT_ROOT, important_results))
]

if (length(existing_results)) {
  cat("\nGenerated key result files:\n")
  cat(paste0("  - ", existing_results, collapse = "\n"), "\n")
}

# Show the current manuscript boxplots when they exist.
important_figures <- c(
  "figures/exp1_first_layer_ari_boxplot.pdf",
  "figures/exp1_first_layer_mr_boxplot.pdf",
  "figures/exp1_pathway_ari_boxplot.pdf",
  "figures/exp1_pathway_mr_boxplot.pdf",
  "figures/exp1_delta_rmse_boxplot.pdf",
  "figures/exp1_nu_recovery_boxplot.pdf",
  "figures/exp1_robustness_first_layer_ari_boxplot.pdf",
  "figures/exp1_robustness_first_layer_mr_boxplot.pdf",
  "figures/exp1_robustness_pathway_ari_boxplot.pdf",
  "figures/exp1_robustness_pathway_mr_boxplot.pdf",
  "figures/exp2_first_layer_ari_boxplot.pdf",
  "figures/exp2_first_layer_mr_boxplot.pdf",
  "figures/exp2_pathway_ari_boxplot.pdf",
  "figures/exp2_pathway_mr_boxplot.pdf",
  "figures/exp2_delta_rmse_boxplot.pdf",
  "figures/exp2_nu_recovery_boxplot.pdf",
  "figures/runtime_comparison_boxplot.pdf"
)

existing_figures <- important_figures[
  file.exists(file.path(PROJECT_ROOT, important_figures))
]

if (length(existing_figures)) {
  cat("\nGenerated manuscript boxplots:\n")
  cat(paste0("  - ", existing_figures, collapse = "\n"), "\n")
}

cat("\nRecommended workflow:\n")
cat("  1. First use RUN_MODE <- \"smoke\".\n")
cat("  2. Check results/ and figures/.\n")
cat("  3. Then change to RUN_MODE <- \"full\" for R = 1000.\n")
cat("\n")
