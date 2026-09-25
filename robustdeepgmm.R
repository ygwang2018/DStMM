# Robust Deep Mixture Model (RDMM)
#
# A pathway-wise Student-t extension of the original DGMM implementation.
# The public interface follows deepgmm() closely. The model retains independent
# layer-specific mixing proportions and introduces one Gamma precision shared
# by every transition in a complete pathway.
#
# Default estimation uses one stochastic draw of the continuous missing data
# per iteration (method = "sem"). Set method = "ecm" for deterministic
# conditional-moment updates.

.rdmm_make_spd <- function(x, floor = 1e-8) {
  x <- (x + t(x)) / 2
  ee <- eigen(x, symmetric = TRUE)
  ee$values <- pmax(ee$values, floor)
  out <- ee$vectors %*% (ee$values * t(ee$vectors))
  (out + t(out)) / 2
}

.rdmm_inverse <- function(x, floor = 1e-8) {
  x <- .rdmm_make_spd(x, floor)
  chol2inv(chol(x))
}

.rdmm_rmvnorm <- function(n, mean, sigma, floor = 1e-8) {
  sigma <- .rdmm_make_spd(sigma, floor)
  q <- length(mean)
  z <- matrix(stats::rnorm(n * q), nrow = n, ncol = q)
  sweep(z %*% chol(sigma), 2L, mean, "+")
}

.rdmm_row_logsumexp <- function(x) {
  m <- apply(x, 1L, max)
  m + log(rowSums(exp(sweep(x, 1L, m, "-"))))
}

.rdmm_logdmvt <- function(y, mu, sigma, nu, floor = 1e-8) {
  y <- as.matrix(y)
  p <- ncol(y)
  sigma <- .rdmm_make_spd(sigma, floor)
  inv_sigma <- chol2inv(chol(sigma))
  centered <- sweep(y, 2L, mu, "-")
  d2 <- rowSums((centered %*% inv_sigma) * centered)
  logdet <- 2 * sum(log(diag(chol(sigma))))
  logdens <- lgamma((nu + p) / 2) - lgamma(nu / 2) -
    0.5 * (p * log(nu * pi) + logdet) -
    0.5 * (nu + p) * log1p(d2 / nu)
  list(logdens = logdens, d2 = d2, inv = inv_sigma)
}

.rdmm_paths <- function(k) {
  grid <- expand.grid(lapply(k, seq_len), KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  out <- as.matrix(grid)
  storage.mode(out) <- "integer"
  colnames(out) <- paste0("layer", seq_along(k))
  out
}

.rdmm_fix_names <- function(init, init_est) {
  ii <- tolower(init)
  if (ii %in% c("kmeans", "k-means", "k")) init <- "kmeans"
  if (ii %in% c("random", "r")) init <- "random"
  if (ii %in% c("hclass", "h")) init <- "hclass"
  if (ii %in% c("mclust", "mclst", "m")) init <- "mclust"

  ie <- tolower(init_est)
  if (ie %in% c("factanal", "factana", "fact", "f")) init_est <- "factanal"
  if (ie %in% c("ppca", "pca", "p")) init_est <- "ppca"
  list(init = init, init_est = init_est)
}

.rdmm_validate <- function(y, layers, k, r, it, eps, init, nu,
                           nu_structure, psi_floor) {
  if (is.data.frame(y)) y <- as.matrix(y)
  if (!is.matrix(y) || !is.numeric(y)) stop("y must be a numeric matrix.")
  if (anyNA(y) || any(!is.finite(y))) stop("y contains missing or non-finite values.")
  if (!(layers %in% 1:3)) stop("layers must be 1, 2, or 3.")
  if (length(k) != layers || any(k < 1) || any(k != as.integer(k))) {
    stop("k must contain one positive integer per layer.")
  }
  if (length(r) != layers || any(r < 1) || any(r != as.integer(r))) {
    stop("r must contain one positive integer per layer.")
  }
  if (any(r >= ncol(y)) || (length(r) > 1L && any(diff(r) >= 0))) {
    stop("The latent dimensions must satisfy p > r[1] > ... > r[layers] >= 1.")
  }
  if (length(it) != 1L || !is.numeric(it) || !is.finite(it) ||
      it < 1L || it != as.integer(it)) {
    stop("it must be one positive integer.")
  }
  if (!is.finite(eps) || eps <= 0) stop("eps must be positive.")
  if (!(init %in% c("kmeans", "random", "hclass", "mclust"))) {
    stop("init must be kmeans, random, hclass, or mclust.")
  }
  if (!(nu_structure %in% c("common", "first_layer", "pathway"))) {
    stop("nu_structure must be common, first_layer, or pathway.")
  }
  if (any(!is.finite(nu)) || any(nu <= 2)) stop("All initial nu values must exceed 2.")
  if (!is.finite(psi_floor) || psi_floor <= 0) stop("psi_floor must be positive.")
  invisible(TRUE)
}

.rdmm_initial_clustering <- function(data, groups, init) {
  n <- nrow(data)
  if (groups == 1L) return(rep.int(1L, n))

  if (init == "kmeans") {
    return(stats::kmeans(data, centers = groups, iter.max = 100,
                         nstart = 30, algorithm = "Hartigan-Wong")$cluster)
  }
  if (init == "hclass") {
    return(stats::cutree(stats::hclust(stats::dist(data), method = "ward.D2"),
                         k = groups))
  }
  if (init == "random") return(sample.int(groups, n, replace = TRUE))
  if (init == "mclust") {
    if (!requireNamespace("mclust", quietly = TRUE)) {
      stop("The mclust package is required when init = 'mclust'.")
    }
    return(mclust::Mclust(data, G = groups, verbose = FALSE)$classification)
  }
  stop("Unknown initialization method.")
}

.rdmm_repair_small_groups <- function(cluster, groups, min_size) {
  n <- length(cluster)
  counts <- tabulate(cluster, nbins = groups)
  for (g in which(counts < min_size)) {
    donor <- which.max(counts)
    donor_idx <- which(cluster == donor)
    take <- min(min_size - counts[g], max(0L, length(donor_idx) - min_size))
    if (take > 0L) {
      move <- sample(donor_idx, take)
      cluster[move] <- g
      counts <- tabulate(cluster, nbins = groups)
    }
  }
  cluster
}

.rdmm_pca_scores <- function(x, q, floor = 1e-8) {
  # High-dimensional-safe PCA fallback shared by DGMM, RDMM and DStMM.
  # The returned scores are approximately whitened so that their sample
  # covariance is I_q, matching the standard-normal latent-factor scale.
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  xc <- sweep(x, 2L, colMeans(x), "-")

  max_rank <- max(0L, min(n - 1L, p))
  q_eff <- min(as.integer(q), max_rank)
  if (q_eff <= 0L) return(matrix(0, nrow = n, ncol = q))

  sv <- svd(xc, nu = q_eff, nv = 0L)
  keep <- seq_len(q_eff)
  z <- sv$u[, keep, drop = FALSE] * sqrt(max(1, n - 1L))

  # Numerical-rank protection: directions with essentially zero singular
  # value carry no usable information and are initialized at zero.
  tol <- max(floor, sqrt(.Machine$double.eps) * max(1, sv$d[1L]))
  bad <- !is.finite(sv$d[keep]) | sv$d[keep] <= tol
  if (any(bad)) z[, bad] <- 0

  if (ncol(z) < q) z <- cbind(z, matrix(0, nrow = n, ncol = q - ncol(z)))
  z[, seq_len(q), drop = FALSE]
}

.rdmm_layer_initialise <- function(data, cluster, groups, q,
                                   init_est = "factanal", psi_floor = 1e-6) {
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  H <- array(0, dim = c(groups, p, q))
  mu <- matrix(0, nrow = p, ncol = groups)
  psi <- psi_inv <- array(0, dim = c(groups, p, p))
  z <- matrix(0, nrow = n, ncol = q)
  weights <- tabulate(cluster, nbins = groups) / n

  for (g in seq_len(groups)) {
    idx <- which(cluster == g)
    xg <- data[idx, , drop = FALSE]
    mu[, g] <- colMeans(xg)
    xc <- sweep(xg, 2L, mu[, g], "-")

    scores <- NULL
    # factanal() is not numerically appropriate when p >= n_g because the
    # within-group sample covariance is singular.  In that case skip the
    # doomed call and go directly to the SVD fallback.  This is important for
    # the real-data cases (e.g. n_g << p at the first layer).
    can_factanal <- init_est == "factanal" &&
      nrow(xg) > max(q + 2L, p + 1L) && p > q
    if (can_factanal) {
      fit <- try(stats::factanal(xg, factors = q, rotation = "none",
                                scores = "regression"), silent = TRUE)
      if (!inherits(fit, "try-error") && !is.null(fit$scores) &&
          ncol(as.matrix(fit$scores)) >= q && all(is.finite(fit$scores))) {
        scores <- as.matrix(fit$scores)[, seq_len(q), drop = FALSE]
      }
    }
    if (is.null(scores)) scores <- .rdmm_pca_scores(xg, q, psi_floor)

    zz <- crossprod(scores)
    H[g, , ] <- crossprod(xc, scores) %*% .rdmm_inverse(zz, psi_floor)
    residual <- xc - scores %*% t(matrix(H[g, , ], nrow = p, ncol = q))
    vv <- if (nrow(residual) > 1L) colMeans(residual^2) else apply(data, 2L, stats::var)
    vv[!is.finite(vv)] <- 1
    vv <- pmax(vv, psi_floor)
    psi[g, , ] <- diag(vv, p)
    psi_inv[g, , ] <- diag(1 / vv, p)
    z[idx, ] <- scores
  }

  list(H = H, mu = mu, psi = psi, psi.inv = psi_inv,
       w = weights, z = z)
}

.rdmm_initialise <- function(y, layers, k, r_full, init, init_est, psi_floor) {
  H <- mu <- psi <- psi_inv <- w <- vector("list", layers)
  current <- y
  n <- nrow(y)

  for (l in seq_len(layers)) {
    cluster <- .rdmm_initial_clustering(current, k[l], init)
    min_size <- min(n, max(2L, r_full[l + 1L] + 2L))
    cluster <- .rdmm_repair_small_groups(cluster, k[l], min_size)
    fit <- .rdmm_layer_initialise(current, cluster, k[l], r_full[l + 1L],
                                  init_est, psi_floor)
    H[[l]] <- fit$H
    mu[[l]] <- fit$mu
    psi[[l]] <- fit$psi
    psi_inv[[l]] <- fit$psi.inv
    w[[l]] <- fit$w
    current <- fit$z
  }

  list(H = H, mu = mu, psi = psi, psi.inv = psi_inv, w = w)
}

.rdmm_joint_path <- function(path, H, mu, psi, r_full, floor = 1e-8) {
  h <- length(path)
  deepest <- r_full[h + 1L]
  mean_stack <- rep(0, deepest)
  cov_stack <- diag(deepest)

  for (l in h:1L) {
    a <- path[l]
    r_prev <- r_full[l]
    r_curr <- r_full[l + 1L]
    Hl <- matrix(H[[l]][a, , ], nrow = r_prev, ncol = r_curr)
    mul <- as.numeric(mu[[l]][, a])
    psil <- matrix(psi[[l]][a, , ], nrow = r_prev, ncol = r_prev)

    idx_curr <- seq_len(r_curr)
    cov_curr_stack <- cov_stack[idx_curr, , drop = FALSE]
    new_mean <- mul + as.numeric(Hl %*% mean_stack[idx_curr])
    cross <- Hl %*% cov_curr_stack
    new_var <- psil + Hl %*% cov_stack[idx_curr, idx_curr, drop = FALSE] %*% t(Hl)
    new_var <- .rdmm_make_spd(new_var, floor)

    cov_stack <- rbind(cbind(new_var, cross), cbind(t(cross), cov_stack))
    cov_stack <- (cov_stack + t(cov_stack)) / 2
    mean_stack <- c(new_mean, mean_stack)
  }

  list(mean = mean_stack, covariance = cov_stack)
}

.rdmm_expand_nu <- function(nu, structure, paths, k) {
  np <- nrow(paths)
  if (structure == "common") return(rep(nu[1L], np))
  if (structure == "first_layer") {
    if (length(nu) == 1L) nu <- rep(nu, k[1L])
    return(as.numeric(nu[paths[, 1L]]))
  }
  if (length(nu) == 1L) nu <- rep(nu, np)
  as.numeric(nu)
}

.rdmm_compact_nu <- function(nu_path, structure, paths, k) {
  if (structure == "common") return(mean(nu_path))
  if (structure == "first_layer") {
    return(vapply(seq_len(k[1L]), function(a) mean(nu_path[paths[, 1L] == a]),
                  numeric(1)))
  }
  nu_path
}

.rdmm_estep <- function(y, H, mu, psi, w, nu_path, paths, r_full,
                        floor = 1e-8) {
  y <- as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  np <- nrow(paths)
  q <- sum(r_full[-1L])

  log_joint <- matrix(NA_real_, nrow = n, ncol = np)
  d2 <- ew <- elogw <- matrix(NA_real_, nrow = n, ncol = np)
  conditional_mean <- vector("list", np)
  conditional_cov <- vector("list", np)
  shape <- rate <- vector("list", np)
  joint_cache <- vector("list", np)

  for (s in seq_len(np)) {
    path <- paths[s, ]
    joint <- .rdmm_joint_path(path, H, mu, psi, r_full, floor)
    joint_cache[[s]] <- joint
    mu_y <- joint$mean[seq_len(p)]
    mu_z <- joint$mean[p + seq_len(q)]
    sigma_yy <- joint$covariance[seq_len(p), seq_len(p), drop = FALSE]
    sigma_yz <- joint$covariance[seq_len(p), p + seq_len(q), drop = FALSE]
    sigma_zz <- joint$covariance[p + seq_len(q), p + seq_len(q), drop = FALSE]
    sigma_yy <- .rdmm_make_spd(sigma_yy, floor)
    inv_yy <- chol2inv(chol(sigma_yy))
    gain <- t(sigma_yz) %*% inv_yy
    centered <- sweep(y, 2L, mu_y, "-")
    mz <- sweep(centered %*% t(gain), 2L, mu_z, "+")
    Cz <- sigma_zz - gain %*% sigma_yz
    Cz <- .rdmm_make_spd(Cz, floor)

    tt <- .rdmm_logdmvt(y, mu_y, sigma_yy, nu_path[s], floor)
    path_prior <- prod(vapply(seq_along(path), function(l) w[[l]][path[l]], numeric(1)))
    path_prior <- max(path_prior, .Machine$double.xmin)
    log_joint[, s] <- log(path_prior) + tt$logdens
    d2[, s] <- tt$d2
    ew[, s] <- (nu_path[s] + p) / (nu_path[s] + tt$d2)
    elogw[, s] <- digamma((nu_path[s] + p) / 2) -
      log((nu_path[s] + tt$d2) / 2)
    shape[[s]] <- rep((nu_path[s] + p) / 2, n)
    rate[[s]] <- (nu_path[s] + tt$d2) / 2
    conditional_mean[[s]] <- mz
    conditional_cov[[s]] <- Cz
  }

  log_py <- .rdmm_row_logsumexp(log_joint)
  tau <- exp(sweep(log_joint, 1L, log_py, "-"))
  tau[!is.finite(tau)] <- 1 / np
  tau <- sweep(tau, 1L, rowSums(tau), "/")

  layers <- ncol(paths)
  layer_prob <- vector("list", layers)
  classification <- matrix(0L, nrow = n, ncol = layers)
  for (l in seq_len(layers)) {
    lp <- matrix(0, nrow = n, ncol = max(paths[, l]))
    for (a in seq_len(ncol(lp))) {
      lp[, a] <- rowSums(tau[, paths[, l] == a, drop = FALSE])
    }
    layer_prob[[l]] <- lp
    classification[, l] <- max.col(lp, ties.method = "first")
  }

  list(loglik = sum(log_py), log_py = log_py, tau = tau,
       layer_prob = layer_prob, s = classification, paths = paths,
       d2 = d2, ew = ew, elogw = elogw, shape = shape, rate = rate,
       conditional_mean = conditional_mean,
       conditional_cov = conditional_cov, joint = joint_cache,
       posterior_w = rowSums(tau * ew))
}

.rdmm_latent_blocks <- function(r_full) {
  dims <- r_full[-1L]
  ends <- cumsum(dims)
  starts <- c(1L, head(ends, -1L) + 1L)
  Map(seq.int, starts, ends)
}

.rdmm_draw_missing <- function(estep, floor = 1e-8) {
  np <- ncol(estep$tau)
  n <- nrow(estep$tau)
  zdraw <- vector("list", np)
  wdraw <- vector("list", np)

  for (s in seq_len(np)) {
    ws <- stats::rgamma(n, shape = estep$shape[[s]], rate = estep$rate[[s]])
    ws <- pmax(ws, .Machine$double.eps)
    C <- .rdmm_make_spd(estep$conditional_cov[[s]], floor)
    q <- ncol(estep$conditional_mean[[s]])
    noise <- matrix(stats::rnorm(n * q), nrow = n, ncol = q) %*% chol(C)
    noise <- noise / sqrt(ws)
    zdraw[[s]] <- estep$conditional_mean[[s]] + noise
    wdraw[[s]] <- ws
  }
  list(z = zdraw, w = wdraw)
}

.rdmm_mstep_sem <- function(y, estep, draws, k, r_full, psi_floor) {
  layers <- length(k)
  paths <- estep$paths
  blocks <- .rdmm_latent_blocks(r_full)
  H <- mu <- psi <- psi_inv <- w <- vector("list", layers)

  for (l in seq_len(layers)) {
    r_prev <- r_full[l]
    r_curr <- r_full[l + 1L]
    H[[l]] <- array(0, c(k[l], r_prev, r_curr))
    mu[[l]] <- matrix(0, r_prev, k[l])
    psi[[l]] <- psi_inv[[l]] <- array(0, c(k[l], r_prev, r_prev))
    w[[l]] <- colMeans(estep$layer_prob[[l]])

    for (a in seq_len(k[l])) {
      D <- matrix(0, 1L + r_curr, 1L + r_curr)
      A <- matrix(0, r_prev, 1L + r_curr)
      G <- matrix(0, r_prev, r_prev)
      N <- 0

      for (s in which(paths[, l] == a)) {
        tau <- estep$tau[, s]
        ww <- draws$w[[s]]
        tw <- tau * ww
        zcurr <- draws$z[[s]][, blocks[[l]], drop = FALSE]
        xprev <- if (l == 1L) y else draws$z[[s]][, blocks[[l - 1L]], drop = FALSE]
        u <- cbind(1, zcurr)
        D <- D + crossprod(u, u * tw)
        A <- A + crossprod(xprev, u * tw)
        G <- G + crossprod(xprev, xprev * tw)
        N <- N + sum(tau)
      }

      D <- .rdmm_make_spd(D, psi_floor)
      B <- A %*% chol2inv(chol(D))
      R <- G - A %*% t(B) - B %*% t(A) + B %*% D %*% t(B)
      R <- .rdmm_make_spd(R / max(N, .Machine$double.eps), psi_floor)
      R <- diag(pmax(diag(R), psi_floor), r_prev)

      mu[[l]][, a] <- B[, 1L]
      H[[l]][a, , ] <- B[, -1L, drop = FALSE]
      psi[[l]][a, , ] <- R
      psi_inv[[l]][a, , ] <- diag(1 / diag(R), r_prev)
    }
  }

  list(H = H, mu = mu, psi = psi, psi.inv = psi_inv, w = w)
}

.rdmm_mstep_ecm <- function(y, estep, k, r_full, psi_floor) {
  layers <- length(k)
  paths <- estep$paths
  blocks <- .rdmm_latent_blocks(r_full)
  H <- mu <- psi <- psi_inv <- w <- vector("list", layers)

  for (l in seq_len(layers)) {
    r_prev <- r_full[l]
    r_curr <- r_full[l + 1L]
    H[[l]] <- array(0, c(k[l], r_prev, r_curr))
    mu[[l]] <- matrix(0, r_prev, k[l])
    psi[[l]] <- psi_inv[[l]] <- array(0, c(k[l], r_prev, r_prev))
    w[[l]] <- colMeans(estep$layer_prob[[l]])

    for (a in seq_len(k[l])) {
      D <- matrix(0, 1L + r_curr, 1L + r_curr)
      A <- matrix(0, r_prev, 1L + r_curr)
      G <- matrix(0, r_prev, r_prev)
      N <- 0

      for (s in which(paths[, l] == a)) {
        tau <- estep$tau[, s]
        ew <- estep$ew[, s]
        alpha <- tau * ew
        mz <- estep$conditional_mean[[s]]
        Cz <- estep$conditional_cov[[s]]
        mcurr <- mz[, blocks[[l]], drop = FALSE]
        Czz <- Cz[blocks[[l]], blocks[[l]], drop = FALSE]

        if (l == 1L) {
          mprev <- y
          Cxx <- matrix(0, r_prev, r_prev)
          Cxz <- matrix(0, r_prev, r_curr)
        } else {
          mprev <- mz[, blocks[[l - 1L]], drop = FALSE]
          Cxx <- Cz[blocks[[l - 1L]], blocks[[l - 1L]], drop = FALSE]
          Cxz <- Cz[blocks[[l - 1L]], blocks[[l]], drop = FALSE]
        }

        nt <- sum(tau)
        D[1L, 1L] <- D[1L, 1L] + sum(alpha)
        d1 <- colSums(mcurr * alpha)
        D[1L, -1L] <- D[1L, -1L] + d1
        D[-1L, 1L] <- D[-1L, 1L] + d1
        D[-1L, -1L] <- D[-1L, -1L] + nt * Czz +
          crossprod(mcurr, mcurr * alpha)

        A[, 1L] <- A[, 1L] + colSums(mprev * alpha)
        A[, -1L] <- A[, -1L] + nt * Cxz +
          crossprod(mprev, mcurr * alpha)
        G <- G + nt * Cxx + crossprod(mprev, mprev * alpha)
        N <- N + nt
      }

      D <- .rdmm_make_spd(D, psi_floor)
      B <- A %*% chol2inv(chol(D))
      R <- G - A %*% t(B) - B %*% t(A) + B %*% D %*% t(B)
      R <- .rdmm_make_spd(R / max(N, .Machine$double.eps), psi_floor)
      R <- diag(pmax(diag(R), psi_floor), r_prev)

      mu[[l]][, a] <- B[, 1L]
      H[[l]][a, , ] <- B[, -1L, drop = FALSE]
      psi[[l]][a, , ] <- R
      psi_inv[[l]][a, , ] <- diag(1 / diag(R), r_prev)
    }
  }

  list(H = H, mu = mu, psi = psi, psi.inv = psi_inv, w = w)
}

.rdmm_update_one_nu <- function(cbar, old, bounds) {
  f <- function(v) log(v / 2) - digamma(v / 2) + 1 + cbar
  lo <- bounds[1L]
  hi <- bounds[2L]
  flo <- f(lo)
  fhi <- f(hi)
  ans <- if (is.finite(flo) && is.finite(fhi) && flo * fhi <= 0) {
    stats::uniroot(f, interval = c(lo, hi), tol = 1e-7)$root
  } else {
    stats::optimize(function(v) f(v)^2, interval = c(lo, hi))$minimum
  }
  if (!is.finite(ans)) old else ans
}

.rdmm_update_nu <- function(estep, nu_path, structure, k, bounds,
                            min_effective = 5) {
  paths <- estep$paths
  np <- nrow(paths)
  updated <- nu_path

  update_group <- function(indices, old) {
    denom <- sum(estep$tau[, indices, drop = FALSE])
    if (!is.finite(denom) || denom < min_effective) return(old)
    cbar <- sum(estep$tau[, indices, drop = FALSE] *
                  (estep$elogw[, indices, drop = FALSE] -
                     estep$ew[, indices, drop = FALSE])) / denom
    .rdmm_update_one_nu(cbar, old, bounds)
  }

  if (structure == "common") {
    v <- update_group(seq_len(np), mean(nu_path))
    updated[] <- v
  } else if (structure == "first_layer") {
    for (a in seq_len(k[1L])) {
      idx <- which(paths[, 1L] == a)
      updated[idx] <- update_group(idx, mean(nu_path[idx]))
    }
  } else {
    for (s in seq_len(np)) updated[s] <- update_group(s, nu_path[s])
  }
  pmin(pmax(updated, bounds[1L]), bounds[2L])
}

.rdmm_parameter_count <- function(k, r_full, nu_structure) {
  layers <- length(k)
  total <- 0
  for (l in seq_len(layers)) {
    r_prev <- r_full[l]
    r_curr <- r_full[l + 1L]
    total <- total + (k[l] - 1L) +
      k[l] * (r_prev * r_curr + r_prev + r_prev) -
      k[l] * r_curr * (r_curr - 1L) / 2
  }
  total + switch(nu_structure,
                 common = 1L,
                 first_layer = k[1L],
                 pathway = prod(k))
}

.rdmm_entropy <- function(tau) {
  z <- pmax(tau, .Machine$double.xmin)
  -sum(z * log(z))
}

.deep_robust_sem <- function(y, numobs, p, r, k, H.list, psi.list,
                             psi.list.inv, mu.list, w.list, it, eps,
                             nu = 10, nu_structure = "common",
                             estimate_nu = TRUE, method = "sem",
                             psi_floor = 1e-6, nu_bounds = c(2.05, 200),
                             min_iter = 10, moving_window = 5,
                             verbose = FALSE) {
  if (length(it) != 1L || !is.numeric(it) || !is.finite(it) ||
      it < 1L || it != as.integer(it)) {
    stop("Internal argument 'it' must be one positive integer.", call. = FALSE)
  }
  it <- as.integer(it)
  min_iter <- as.integer(min_iter[1L])
  moving_window <- as.integer(moving_window[1L])

  paths <- .rdmm_paths(k)
  nu_path <- .rdmm_expand_nu(nu, nu_structure, paths, k)
  likelihood <- numeric(0)
  best <- NULL
  best_loglik <- -Inf
  ratio <- Inf
  iteration <- 0L

  for (iteration in seq_len(it)) {
    estep <- .rdmm_estep(y, H.list, mu.list, psi.list, w.list,
                         nu_path, paths, r, psi_floor)

    if (method == "sem") {
      draws <- .rdmm_draw_missing(estep, psi_floor)
      updated <- .rdmm_mstep_sem(y, estep, draws, k, r, psi_floor)
    } else {
      updated <- .rdmm_mstep_ecm(y, estep, k, r, psi_floor)
    }

    new_nu <- if (estimate_nu) {
      .rdmm_update_nu(estep, nu_path, nu_structure, k, nu_bounds)
    } else nu_path

    new_estep <- .rdmm_estep(y, updated$H, updated$mu, updated$psi,
                             updated$w, new_nu, paths, r, psi_floor)
    likelihood <- c(likelihood, new_estep$loglik)

    H.list <- updated$H
    mu.list <- updated$mu
    psi.list <- updated$psi
    psi.list.inv <- updated$psi.inv
    w.list <- updated$w
    nu_path <- new_nu

    if (new_estep$loglik > best_loglik) {
      best_loglik <- new_estep$loglik
      best <- list(H = H.list, mu = mu.list, psi = psi.list,
                   psi.inv = psi.list.inv, w = w.list,
                   nu_path = nu_path, estep = new_estep)
    }

    if (length(likelihood) >= 2L * moving_window && iteration >= min_iter) {
      old_ma <- mean(likelihood[(length(likelihood) - 2L * moving_window + 1L):
                                  (length(likelihood) - moving_window)])
      new_ma <- mean(tail(likelihood, moving_window))
      ratio <- abs(new_ma - old_ma) / (abs(old_ma) + 1)
    }
    if (verbose) {
      message(sprintf("iteration %d: log-likelihood = %.6f, ratio = %.3g",
                      iteration, tail(likelihood, 1L), ratio))
    }
    if (iteration >= min_iter && is.finite(ratio) && ratio < eps) break
  }

  H.list <- best$H
  mu.list <- best$mu
  psi.list <- best$psi
  psi.list.inv <- best$psi.inv
  w.list <- best$w
  nu_path <- best$nu_path
  final <- best$estep

  hpar <- .rdmm_parameter_count(k, r, nu_structure)
  lik <- final$loglik
  bic <- -2 * lik + hpar * log(numobs)
  aic <- -2 * lik + 2 * hpar
  EN <- .rdmm_entropy(final$tau)
  clc <- -2 * lik + 2 * EN
  icl_bic <- bic + 2 * EN

  list(H = H.list, w = w.list, mu = mu.list, psi = psi.list,
       psi.inv = psi.list.inv, likelihood = likelihood,
       bic = bic, aic = aic, clc = clc, icl_bic = icl_bic,
       s = final$s, h = hpar, ps.y = final$tau,
       ps.y.list = final$layer_prob, paths = paths,
       nu = .rdmm_compact_nu(nu_path, nu_structure, paths, k),
       nu_path = nu_path, posterior_w = final$posterior_w,
       path_w = final$ew, loglik = lik, iterations = iteration,
       convergence_ratio = ratio)
}

# Layer-specific wrappers retain the structure of the original code base.
deep.robust.sem.alg.1 <- function(y, numobs, p, r, k, H.list, psi.list,
                                  psi.list.inv, mu.list, w.list, it, eps,
                                  ...) {
  .deep_robust_sem(y, numobs, p, c(p, r), k, H.list, psi.list,
                   psi.list.inv, mu.list, w.list, it, eps, ...)
}

deep.robust.sem.alg.2 <- function(y, numobs, p, r, k, H.list, psi.list,
                                  psi.list.inv, mu.list, w.list, it, eps,
                                  ...) {
  .deep_robust_sem(y, numobs, p, r, k, H.list, psi.list,
                   psi.list.inv, mu.list, w.list, it, eps, ...)
}

deep.robust.sem.alg.3 <- function(y, numobs, p, r, k, H.list, psi.list,
                                  psi.list.inv, mu.list, w.list, it, eps,
                                  ...) {
  .deep_robust_sem(y, numobs, p, r, k, H.list, psi.list,
                   psi.list.inv, mu.list, w.list, it, eps, ...)
}

robustdeepgmm <- function(y, layers, k, r,
                          it = 250, eps = 0.001,
                          init = "kmeans", init_est = "factanal",
                          seed = NULL, scale = TRUE,
                          nu = 10,
                          nu_structure = c("common", "first_layer", "pathway"),
                          estimate_nu = TRUE,
                          method = c("sem", "ecm"),
                          psi_floor = 1e-6,
                          nu_bounds = c(2.05, 200),
                          min_iter = 10,
                          moving_window = 5,
                          verbose = FALSE) {
  call <- match.call()
  if (is.data.frame(y)) y <- as.matrix(y)
  nu_structure <- match.arg(nu_structure)
  method <- match.arg(method)
  fixed <- .rdmm_fix_names(init, init_est)
  init <- fixed$init
  init_est <- fixed$init_est
  .rdmm_validate(y, layers, k, r, it, eps, init, nu,
                 nu_structure, psi_floor)

  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.finite(seed)) stop("seed must be one finite number.")
    set.seed(as.integer(seed))
  }

  y_original <- y
  center <- rep(0, ncol(y))
  scale_value <- rep(1, ncol(y))
  if (scale) {
    ys <- base::scale(y)
    center <- attr(ys, "scaled:center")
    scale_value <- attr(ys, "scaled:scale")
    scale_value[!is.finite(scale_value) | scale_value == 0] <- 1
    y <- sweep(sweep(y, 2L, center, "-"), 2L, scale_value, "/")
  }

  numobs <- nrow(y)
  p <- ncol(y)
  r_full <- c(p, as.integer(r))
  k <- as.integer(k)
  expected_nu <- switch(nu_structure,
                        common = 1L,
                        first_layer = k[1L],
                        pathway = prod(k))
  if (!(length(nu) %in% c(1L, expected_nu))) {
    stop("nu must have length 1 or the number implied by nu_structure.")
  }
  initial <- .rdmm_initialise(y, layers, k, r_full, init, init_est, psi_floor)

  if (layers == 1L) {
    out <- deep.robust.sem.alg.1(
      y = y, numobs = numobs, p = p, r = r_full[2L], k = k,
      H.list = initial$H, psi.list = initial$psi,
      psi.list.inv = initial$psi.inv, mu.list = initial$mu,
      w.list = initial$w, it = as.integer(it), eps = eps,
      nu = nu, nu_structure = nu_structure,
      estimate_nu = estimate_nu, method = method,
      psi_floor = psi_floor, nu_bounds = nu_bounds,
      min_iter = min_iter, moving_window = moving_window,
      verbose = verbose
    )
  } else if (layers == 2L) {
    out <- deep.robust.sem.alg.2(
      y = y, numobs = numobs, p = p, r = r_full, k = k,
      H.list = initial$H, psi.list = initial$psi,
      psi.list.inv = initial$psi.inv, mu.list = initial$mu,
      w.list = initial$w, it = as.integer(it), eps = eps,
      nu = nu, nu_structure = nu_structure,
      estimate_nu = estimate_nu, method = method,
      psi_floor = psi_floor, nu_bounds = nu_bounds,
      min_iter = min_iter, moving_window = moving_window,
      verbose = verbose
    )
  } else {
    out <- deep.robust.sem.alg.3(
      y = y, numobs = numobs, p = p, r = r_full, k = k,
      H.list = initial$H, psi.list = initial$psi,
      psi.list.inv = initial$psi.inv, mu.list = initial$mu,
      w.list = initial$w, it = as.integer(it), eps = eps,
      nu = nu, nu_structure = nu_structure,
      estimate_nu = estimate_nu, method = method,
      psi_floor = psi_floor, nu_bounds = nu_bounds,
      min_iter = min_iter, moving_window = moving_window,
      verbose = verbose
    )
  }

  out$lik <- out$likelihood
  output <- out[c("H", "w", "mu", "psi", "lik", "bic", "aic", "clc",
                  "icl_bic", "s", "h", "ps.y", "ps.y.list", "paths",
                  "nu", "nu_path", "posterior_w", "path_w", "loglik",
                  "iterations", "convergence_ratio")]
  output <- c(output,
              list(k = k, r = r_full[-1L], numobs = numobs,
                   layers = layers, seed = seed, method = method,
                   nu_structure = nu_structure, estimate_nu = estimate_nu,
                   scaled = scale, center = center, scale = scale_value,
                   original_data = y_original))
  output$call <- call
  class(output) <- "rdmm"
  invisible(output)
}

rdmm <- robustdeepgmm

predict.rdmm <- function(object, newdata = NULL,
                         type = c("path", "layer", "posterior", "weight"),
                         layer = 1L, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    if (type == "path") return(max.col(object$ps.y, ties.method = "first"))
    if (type == "layer") return(object$s[, layer])
    if (type == "posterior") return(object$ps.y)
    return(object$posterior_w)
  }

  y <- as.matrix(newdata)
  if (ncol(y) != length(object$center)) stop("newdata has the wrong number of columns.")
  if (object$scaled) {
    y <- sweep(sweep(y, 2L, object$center, "-"), 2L, object$scale, "/")
  }
  r_full <- c(ncol(y), object$r)
  ee <- .rdmm_estep(y, object$H, object$mu, object$psi, object$w,
                    object$nu_path, object$paths, r_full)
  if (type == "path") return(max.col(ee$tau, ties.method = "first"))
  if (type == "layer") return(ee$s[, layer])
  if (type == "posterior") return(ee$tau)
  ee$posterior_w
}

print.rdmm <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nRobust deep mixture model\n")
  cat("Layers:", x$layers, "\n")
  cat("Components:", paste(x$k, collapse = " x "), "\n")
  cat("Latent dimensions:", paste(x$r, collapse = " > "), "\n")
  cat("Estimation method:", toupper(x$method), "\n")
  cat("Degrees-of-freedom structure:", x$nu_structure, "\n")
  cat("Degrees of freedom:", paste(round(x$nu, 3), collapse = ", "), "\n")
  cat("Log-likelihood:", round(x$loglik, 3), "\n")
  cat("BIC:", round(x$bic, 3), "\n")
  invisible(x)
}

summary.rdmm <- function(object, ...) {
  print(object)
  cat("AIC:", round(object$aic, 3), "\n")
  cat("ICL-BIC:", round(object$icl_bic, 3), "\n")
  cat("Iterations:", object$iterations, "\n")
  cat("Final convergence ratio:", signif(object$convergence_ratio, 4), "\n")
  cat("Posterior precision summary:\n")
  print(summary(object$posterior_w))
  invisible(object)
}
