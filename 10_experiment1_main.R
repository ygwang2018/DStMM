rm(list = ls())
root <- normalizePath(if (file.exists("R/load_all.R")) "." else "..")
source(file.path(root, "R", "load_all.R")); load_dstmm_project(root)
source(file.path(root, "R", "simulation_runner.R"))

RUN_MODE <- Sys.getenv("RUN_MODE", "full")
N_REP <- if (RUN_MODE == "full") {
  read_env_int("FULL_REPLICATIONS", 1000L)
} else {
  read_env_int("SMOKE_REPLICATIONS", 2L)
}
MAX_ITER <- if (RUN_MODE == "full") 300L else 40L
EPS <- if (RUN_MODE == "full") 1e-4 else 5e-3
MIN_ITER <- if (RUN_MODE == "full") 100L else 15L
WINDOW <- if (RUN_MODE == "full") 20L else 5L
M <- read_env_int("MC_SIZE", 1L)
COND <- Sys.getenv("CONDITION_INDEX", "")
COND <- if (nzchar(COND)) as.integer(COND) else NULL
REP_START <- read_env_int("REP_START", 1L)
REP_END <- read_env_int("REP_END", N_REP)

run_simulation_grid(
  grid = experiment1_grid(), experiment_name = "Experiment1",
  output_file = file.path(root, "results", "experiment1_raw.csv"),
  n_replications = N_REP, max_iter = MAX_ITER, eps = EPS,
  min_iter = MIN_ITER, moving_window = WINDOW, M = M,
  base_data_seed = 110000L, base_fit_seed = 310000L,
  condition_index = COND, rep_start = REP_START, rep_end = REP_END
)
