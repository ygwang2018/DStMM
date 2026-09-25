# Generic resumable simulation runner.

read_env_int <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (!nzchar(x)) return(as.integer(default))
  as.integer(x)
}

read_env_num <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (!nzchar(x)) return(as.numeric(default))
  as.numeric(x)
}

run_simulation_grid <- function(grid, experiment_name, output_file,
                                n_replications = 1000L,
                                max_iter = 300L,
                                eps = 1e-4,
                                min_iter = 100L,
                                moving_window = 20L,
                                M = 1L,
                                base_data_seed = 100000L,
                                base_fit_seed = 300000L,
                                condition_index = NULL,
                                rep_start = 1L,
                                rep_end = n_replications,
                                verbose = TRUE) {
  if (!is.null(condition_index)) {
    condition_index <- as.integer(condition_index)
    if (condition_index < 1L || condition_index > nrow(grid)) stop("Invalid condition_index.")
    grid <- grid[condition_index, , drop = FALSE]
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  old <- if (file.exists(output_file)) read.csv(output_file, stringsAsFactors = FALSE) else NULL
  rep_start <- max(1L, as.integer(rep_start))
  rep_end <- min(as.integer(n_replications), as.integer(rep_end))

  for (gi in seq_len(nrow(grid))) {
    cond <- grid[gi, , drop = FALSE]
    pis <- pi_from_label(cond$Pi)
    pars <- make_dstmm_parameters(
      nu = cond$Nu,
      kappa = cond$Kappa,
      pi1 = pis[[1L]],
      pi2 = pis[[2L]]
    )

    # Stable condition id independent of row ordering.
    condition_key <- paste(cond$Design, cond$N, cond$Nu, cond$Kappa, cond$Pi, sep = "|")
    key_hash <- sum(utf8ToInt(condition_key))

    for (replication in seq.int(rep_start, rep_end)) {
      if (!is.null(old) && nrow(old)) {
        done <- old$Experiment == experiment_name &
          old$Design == cond$Design & old$N == cond$N &
          old$Nu == cond$Nu & old$Kappa == cond$Kappa & old$Pi == cond$Pi &
          old$Replication == replication
        done_rows <- old[done, , drop = FALSE]
        required_models <- c("DGMM", "RDMM", "DStMM")
        completed_ok <- vapply(required_models, function(mm) {
          zz <- done_rows[done_rows$Model == mm, , drop = FALSE]
          nrow(zz) > 0L && any(zz$Success %in% TRUE)
        }, logical(1))
        if (all(completed_ok)) next
      }

      if (verbose) {
        cat(sprintf("\n%s | %s | n=%d nu=%g kappa=%g pi=%s | rep %d/%d\n",
                    experiment_name, cond$Design, cond$N, cond$Nu, cond$Kappa,
                    cond$Pi, replication, n_replications))
      }

      data_seed <- as.integer(base_data_seed + key_hash * 1000L + replication)
      fit_seed <- as.integer(base_fit_seed + key_hash * 1000L + replication)
      sim <- simulate_dstmm_data(n = cond$N, parameters = pars, seed = data_seed)

      fit_out <- fit_three_models(
        sim = sim, seed = fit_seed,
        max_iter = max_iter, eps = eps,
        initial_nu = cond$Nu,
        min_iter = min_iter,
        moving_window = moving_window,
        M = M
      )
      rr <- fit_out$results
      rr$Experiment <- experiment_name
      rr$Design <- cond$Design
      rr$N <- cond$N
      rr$Nu <- cond$Nu
      rr$Kappa <- cond$Kappa
      rr$Pi <- cond$Pi
      rr$Replication <- replication
      rr$Data_seed <- data_seed
      rr$Fit_seed <- fit_seed
      rr$MC_size <- M

      # Replace an incomplete earlier replicate/condition before appending.
      if (!is.null(old) && nrow(old)) {
        replace <- old$Experiment == experiment_name &
          old$Design == cond$Design & old$N == cond$N &
          old$Nu == cond$Nu & old$Kappa == cond$Kappa & old$Pi == cond$Pi &
          old$Replication == replication
        old <- old[!replace, , drop = FALSE]
      }
      old <- if (is.null(old)) rr else rbind(old, rr)
      write.csv(old, output_file, row.names = FALSE)
    }
  }
  invisible(old)
}
