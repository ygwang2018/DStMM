# Source all model and simulation functions from the project root.
load_dstmm_project <- function(root = ".") {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required by the uploaded DGMM implementation. Run scripts/00_install_packages.R")
  }
  if (!requireNamespace("corpcor", quietly = TRUE)) {
    stop("Package 'corpcor' is required by the uploaded DGMM implementation. Run scripts/00_install_packages.R")
  }
  assign("rmvnorm", mvtnorm::rmvnorm, envir = .GlobalEnv)
  assign("dmvnorm", mvtnorm::dmvnorm, envir = .GlobalEnv)
  # The uploaded DGMM code calls these functions without a namespace.
  assign("is.positive.definite", corpcor::is.positive.definite, envir = .GlobalEnv)
  assign("make.positive.definite", corpcor::make.positive.definite, envir = .GlobalEnv)
  # Load the legacy DGMM numerical/SEM helpers first, then the RDMM file that
  # owns the shared high-dimensional initializer.  deepgmm.R is sourced only
  # after those shared helpers exist, so DGMM/RDMM/DStMM use one initializer.
  dgmm_core_files <- c(
    "misc.R", "ginv.R", "adjustedRandIndex.R", "initial_clustering.R",
    "fix_para.R", "factanal_para.R", "ppca_para.R", "valid_args.R",
    "compute_est.R", "compute_lik.R", "deep.sem.alg.1.R",
    "deep.sem.alg.2.R", "deep.sem.alg.3.R"
  )
  for (f in dgmm_core_files) source(file.path(root, "R", "dgmm", f), local = .GlobalEnv)
  source(file.path(root, "R", "robustdeepgmm.R"), local = .GlobalEnv)
  source(file.path(root, "R", "dgmm", "deepgmm.R"), local = .GlobalEnv)
  source(file.path(root, "R", "dgmm", "print.R"), local = .GlobalEnv)
  source(file.path(root, "R", "dstmm.R"), local = .GlobalEnv)
  source(file.path(root, "R", "simulation_design.R"), local = .GlobalEnv)
  source(file.path(root, "R", "metrics_and_fitting.R"), local = .GlobalEnv)
  invisible(TRUE)
}
