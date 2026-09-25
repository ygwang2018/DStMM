# Repair missing/failed DGMM rows in existing raw simulation CSV files.
# This avoids rerunning successful RDMM/DStMM fits from older runs.
rm(list = ls())
root <- normalizePath(if (file.exists("R/load_all.R")) "." else "..")
source(file.path(root, "R", "load_all.R")); load_dstmm_project(root)

RUN_MODE <- Sys.getenv("RUN_MODE", "smoke")
MAX_ITER <- if (RUN_MODE == "full") 300L else 40L
EPS <- if (RUN_MODE == "full") 1e-4 else 5e-3

raw_files <- c(
  "experiment1_raw.csv",
  "experiment1_robustness_raw.csv",
  "experiment2_raw.csv"
)

repair_one <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  d <- read.csv(path, stringsAsFactors = FALSE)
  if (!nrow(d)) return(invisible(NULL))

  key_cols <- c("Experiment", "Design", "N", "Nu", "Kappa", "Pi", "Replication")
  keys <- unique(d[, key_cols, drop = FALSE])
  repaired <- 0L

  for (i in seq_len(nrow(keys))) {
    keep <- rep(TRUE, nrow(d))
    for (g in key_cols) keep <- keep & d[[g]] == keys[[g]][i]
    block <- d[keep, , drop = FALSE]
    dg <- block[block$Model == "DGMM", , drop = FALSE]
    dg_ok <- nrow(dg) > 0L && any(dg$Success %in% TRUE)
    if (dg_ok) next

    # Recreate exactly the same simulated dataset from the stored condition and seed.
    template <- block[1L, , drop = FALSE]
    pis <- pi_from_label(as.character(template$Pi))
    pars <- make_dstmm_parameters(
      nu = as.numeric(template$Nu),
      kappa = as.numeric(template$Kappa),
      pi1 = pis[[1L]], pi2 = pis[[2L]]
    )
    sim <- simulate_dstmm_data(
      n = as.integer(template$N), parameters = pars,
      seed = as.integer(template$Data_seed)
    )

    cat(sprintf("Repairing DGMM: %s | %s | n=%d nu=%g kappa=%g | rep=%d\n",
                template$Experiment, template$Design, template$N,
                template$Nu, template$Kappa, template$Replication))

    t0 <- proc.time()[[3L]]
    fit <- tryCatch(
      fit_dgmm_once(sim$y, seed = as.integer(template$Fit_seed),
                    max_iter = MAX_ITER, eps = EPS),
      error = identity
    )
    elapsed <- proc.time()[[3L]] - t0

    if (inherits(fit, "error")) {
      rr <- failed_result_row("DGMM", conditionMessage(fit), elapsed)
    } else {
      rr <- model_result_row(fit, "DGMM", sim, elapsed, MAX_ITER, EPS)
    }

    # Add metadata columns using the existing row as a template.
    rr$Experiment <- template$Experiment
    rr$Design <- template$Design
    rr$N <- template$N
    rr$Nu <- template$Nu
    rr$Kappa <- template$Kappa
    rr$Pi <- template$Pi
    rr$Replication <- template$Replication
    rr$Data_seed <- template$Data_seed
    rr$Fit_seed <- template$Fit_seed
    rr$MC_size <- if ("MC_size" %in% names(template)) template$MC_size else 1L

    # Remove old DGMM row(s), retain RDMM/DStMM rows, append repaired DGMM.
    remove <- keep & d$Model == "DGMM"
    d <- d[!remove, , drop = FALSE]

    # Match column ordering robustly.
    missing_in_rr <- setdiff(names(d), names(rr))
    for (nm in missing_in_rr) rr[[nm]] <- NA
    extra_in_rr <- setdiff(names(rr), names(d))
    for (nm in extra_in_rr) d[[nm]] <- NA
    rr <- rr[, names(d), drop = FALSE]
    d <- rbind(d, rr)
    repaired <- repaired + 1L
  }

  # Stable ordering for readability.
  if (nrow(d)) {
    model_order <- match(d$Model, c("DGMM", "RDMM", "DStMM"))
    oo <- order(d$Experiment, d$Design, d$N, d$Nu, d$Kappa,
                d$Pi, d$Replication, model_order, na.last = TRUE)
    d <- d[oo, , drop = FALSE]
  }
  write.csv(d, path, row.names = FALSE)
  cat(sprintf("%s: repaired %d DGMM condition-replicate(s).\n", basename(path), repaired))
  invisible(d)
}

for (ff in raw_files) repair_one(file.path(root, "results", ff))
cat("DGMM repair stage completed.\n")
