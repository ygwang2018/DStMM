# Metrics, label alignment, common fitting wrapper, and result summaries.

adjusted_rand_index_local <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) stop("x and y must have the same length.")
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  a <- sum(choose2(tab))
  row_pairs <- sum(choose2(rowSums(tab)))
  col_pairs <- sum(choose2(colSums(tab)))
  total <- choose2(sum(tab))
  if (total <= 0) return(0)
  expected <- row_pairs * col_pairs / total
  max_index <- 0.5 * (row_pairs + col_pairs)
  denom <- max_index - expected
  if (abs(denom) < .Machine$double.eps) return(ifelse(a == max_index, 1, 0))
  (a - expected) / denom
}

all_permutations_local <- function(x) {
  x <- as.integer(x)
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    rest <- all_permutations_local(x[-i])
    cbind(x[i], rest)
  }))
}

best_label_mapping <- function(predicted, truth, groups) {
  predicted <- as.integer(predicted)
  truth <- as.integer(truth)
  perms <- all_permutations_local(seq_len(groups))
  errors <- apply(perms, 1L, function(mapping) {
    mean(mapping[predicted] != truth)
  })
  best <- which.min(errors)
  list(mapping = as.integer(perms[best, ]), error = errors[best])
}

misclassification_rate <- function(predicted, truth, groups) {
  best_label_mapping(predicted, truth, groups)$error
}

extract_path_prediction <- function(fit) {
  if (!is.null(fit$ps.y)) {
    return(max.col(as.matrix(fit$ps.y), ties.method = "first"))
  }
  encode_path(fit$s)
}

extract_loglik <- function(fit, model) {
  if (model %in% c("RDMM", "DStMM")) return(as.numeric(fit$loglik)[1L])
  tail(as.numeric(fit$lik), 1L)
}

extract_nu <- function(fit, model) {
  if (!(model %in% c("RDMM", "DStMM"))) return(NA_real_)
  if (!is.null(fit$nu)) return(mean(as.numeric(fit$nu)))
  if (!is.null(fit$nu_path)) return(mean(as.numeric(fit$nu_path)))
  NA_real_
}

extract_delta_rmse <- function(fit, sim) {
  if (is.null(fit$delta)) return(NA_real_)
  pred <- as.integer(fit$s[, 1L])
  truth <- sim$layer1
  align <- best_label_mapping(pred, truth, groups = sim$parameters$k[1L])
  est <- fit$delta[[1L]]
  true <- sim$parameters$delta[[1L]]
  aligned <- matrix(NA_real_, nrow = nrow(est), ncol = ncol(est))
  # mapping[predicted_label] = truth_label
  for (a in seq_len(ncol(est))) aligned[, align$mapping[a]] <- est[, a]
  sqrt(mean((aligned - true)^2))
}

extract_alpha_rmse <- function(fit, sim) {
  if (is.null(fit$path_alpha)) return(NA_real_)
  pred_path <- extract_path_prediction(fit)
  truth_path <- sim$path_id
  g <- prod(sim$parameters$k)
  align <- best_label_mapping(pred_path, truth_path, groups = g)

  # True pathway alphas from the manuscript recursion.
  p <- sim$parameters$p
  true_alpha <- matrix(NA_real_, nrow = p, ncol = g)
  paths <- sim$parameters$path_labels
  for (s in seq_len(g)) {
    a <- paths[s, 1L]
    b <- paths[s, 2L]
    true_alpha[, s] <- sim$parameters$delta[[1L]][, a] +
      sim$parameters$lambda[[1L]][[a]] %*% sim$parameters$delta[[2L]][, b]
  }
  est_alpha <- do.call(cbind, fit$path_alpha)
  aligned <- matrix(NA_real_, nrow = p, ncol = g)
  for (s in seq_len(g)) aligned[, align$mapping[s]] <- est_alpha[, s]
  sqrt(mean((aligned - true_alpha)^2))
}

model_result_row <- function(fit, model, sim, elapsed, max_iter, eps) {
  pred_l1 <- as.integer(fit$s[, 1L])
  pred_path <- extract_path_prediction(fit)
  converged <- if (model == "DStMM" && !is.null(fit$converged)) {
    isTRUE(fit$converged)
  } else if (!is.null(fit$convergence_ratio)) {
    is.finite(fit$convergence_ratio) && fit$convergence_ratio < eps
  } else {
    NA
  }

  data.frame(
    Model = model,
    Layer1_ARI = adjusted_rand_index_local(sim$layer1, pred_l1),
    Layer1_MR = misclassification_rate(pred_l1, sim$layer1,
                                       groups = sim$parameters$k[1L]),
    Pathway_ARI = adjusted_rand_index_local(sim$path_id, pred_path),
    Pathway_MR = misclassification_rate(pred_path, sim$path_id,
                                        groups = prod(sim$parameters$k)),
    Estimated_nu = extract_nu(fit, model),
    Delta_RMSE = if (model == "DStMM") extract_delta_rmse(fit, sim) else NA_real_,
    Alpha_RMSE = if (model == "DStMM") extract_alpha_rmse(fit, sim) else NA_real_,
    Log_likelihood = extract_loglik(fit, model),
    Iterations = if (!is.null(fit$iterations)) as.integer(fit$iterations) else NA_integer_,
    Converged = converged,
    Elapsed_seconds = elapsed,
    Success = TRUE,
    Error = "",
    stringsAsFactors = FALSE
  )
}

failed_result_row <- function(model, error, elapsed) {
  data.frame(
    Model = model,
    Layer1_ARI = NA_real_, Layer1_MR = NA_real_,
    Pathway_ARI = NA_real_, Pathway_MR = NA_real_,
    Estimated_nu = NA_real_, Delta_RMSE = NA_real_, Alpha_RMSE = NA_real_,
    Log_likelihood = NA_real_,
    Iterations = NA_integer_, Converged = FALSE,
    Elapsed_seconds = elapsed, Success = FALSE,
    Error = as.character(error), stringsAsFactors = FALSE
  )
}

fit_dgmm_once <- function(y, seed, max_iter, eps) {
  set.seed(as.integer(seed))
  deepgmm(
    y = y, layers = 2, k = c(2, 2), r = c(5, 2),
    it = as.integer(max_iter), eps = eps,
    init = "kmeans", init_est = "factanal",
    seed = as.integer(seed), scale = FALSE
  )
}

fit_rdmm_once <- function(y, seed, max_iter, eps,
                          initial_nu = 6, min_iter = 100L,
                          moving_window = 20L) {
  set.seed(as.integer(seed))
  robustdeepgmm(
    y = y, layers = 2, k = c(2, 2), r = c(5, 2),
    it = as.integer(max_iter), eps = eps,
    init = "kmeans", init_est = "factanal",
    seed = as.integer(seed), scale = FALSE,
    nu = initial_nu, nu_structure = "common", estimate_nu = TRUE,
    method = "sem", psi_floor = 1e-6,
    nu_bounds = c(2.05, 200),
    min_iter = as.integer(min_iter), moving_window = as.integer(moving_window),
    verbose = FALSE
  )
}

fit_dstmm_once <- function(y, seed, max_iter, eps, rdmm_warm_start = NULL,
                           initial_nu = 6, min_iter = 100L,
                           moving_window = 20L, M = 1L) {
  set.seed(as.integer(seed))
  dstmm(
    y = y, layers = 2, k = c(2, 2), r = c(5, 2),
    it = as.integer(max_iter), eps = eps,
    init = "kmeans", init_est = "factanal",
    seed = as.integer(seed), scale = FALSE,
    nu = initial_nu, nu_structure = "common", estimate_nu = TRUE,
    active_set = 1L, M = as.integer(M),
    warm_start = rdmm_warm_start,
    psi_floor = 1e-6, nu_bounds = c(2.05, 200),
    min_iter = as.integer(min_iter), moving_window = as.integer(moving_window),
    verbose = FALSE
  )
}

fit_three_models <- function(sim, seed, max_iter, eps,
                             initial_nu = 6,
                             min_iter = 100L,
                             moving_window = 20L,
                             M = 1L) {
  rows <- list()
  fits <- list()

  t0 <- proc.time()[[3L]]
  dg <- tryCatch(fit_dgmm_once(sim$y, seed, max_iter, eps), error = identity)
  elapsed <- proc.time()[[3L]] - t0
  if (inherits(dg, "error")) {
    rows[["DGMM"]] <- failed_result_row("DGMM", conditionMessage(dg), elapsed)
  } else {
    fits$DGMM <- dg
    rows[["DGMM"]] <- model_result_row(dg, "DGMM", sim, elapsed, max_iter, eps)
  }

  t0 <- proc.time()[[3L]]
  rd <- tryCatch(
    fit_rdmm_once(sim$y, seed, max_iter, eps, initial_nu, min_iter, moving_window),
    error = identity
  )
  elapsed <- proc.time()[[3L]] - t0
  if (inherits(rd, "error")) {
    rows[["RDMM"]] <- failed_result_row("RDMM", conditionMessage(rd), elapsed)
    rd_warm <- NULL
  } else {
    fits$RDMM <- rd
    rows[["RDMM"]] <- model_result_row(rd, "RDMM", sim, elapsed, max_iter, eps)
    rd_warm <- rd
  }

  t0 <- proc.time()[[3L]]
  ds <- tryCatch(
    fit_dstmm_once(sim$y, seed, max_iter, eps,
                   rdmm_warm_start = rd_warm,
                   initial_nu = initial_nu,
                   min_iter = min_iter,
                   moving_window = moving_window,
                   M = M),
    error = identity
  )
  elapsed <- proc.time()[[3L]] - t0
  if (inherits(ds, "error")) {
    rows[["DStMM"]] <- failed_result_row("DStMM", conditionMessage(ds), elapsed)
  } else {
    fits$DStMM <- ds
    rows[["DStMM"]] <- model_result_row(ds, "DStMM", sim, elapsed, max_iter, eps)
  }

  list(results = do.call(rbind, rows), fits = fits)
}

mean_se <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(Mean = NA_real_, SD = NA_real_, SE = NA_real_, N = 0))
  c(Mean = mean(x), SD = if (length(x) > 1L) sd(x) else NA_real_,
    SE = if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_, N = length(x))
}

summarise_metric <- function(data, group_cols, metric) {
  keys <- unique(data[, group_cols, drop = FALSE])
  out <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    keep <- rep(TRUE, nrow(data))
    for (g in group_cols) keep <- keep & data[[g]] == keys[[g]][i]
    x <- data[[metric]][keep & data$Success]
    ss <- mean_se(x)
    out[[i]] <- cbind(keys[i, , drop = FALSE],
                      Metric = metric,
                      Mean = ss[["Mean"]], SD = ss[["SD"]],
                      SE = ss[["SE"]], N_success = ss[["N"]])
  }
  do.call(rbind, out)
}

build_model_summary <- function(data) {
  metrics <- c("Layer1_ARI", "Layer1_MR", "Pathway_ARI", "Pathway_MR",
               "Estimated_nu", "Delta_RMSE", "Alpha_RMSE",
               "Log_likelihood", "Elapsed_seconds")
  groups <- c("Experiment", "Design", "N", "Nu", "Kappa", "Pi", "Model")
  do.call(rbind, lapply(metrics, function(m) summarise_metric(data, groups, m)))
}
