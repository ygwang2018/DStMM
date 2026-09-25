# Deep Skew-t Mixture Model (DStMM)
#
# Implements the pathway construction in "Deep Skew-t Mixture Models".
# One inverse-gamma scale H is shared by all transitions of an observation-path
# pair. Conditional on H, every transition is Gaussian and an active layer can
# contain a scale-dependent mean shift H * delta.
#
# This file intentionally reuses the stable numerical utilities and
# initialisation code from robustdeepgmm.R. Source robustdeepgmm.R first.

.dstmm_require_rdmm <- function() {
  needed <- c(".rdmm_make_spd", ".rdmm_inverse", ".rdmm_initialise",
              ".rdmm_paths", ".rdmm_fix_names")
  missing <- needed[!vapply(needed, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop("Source R/robustdeepgmm.R before R/dstmm.R. Missing: ",
         paste(missing, collapse = ", "))
  }
}

.dstmm_logsumexp_rows <- function(x) {
  m <- apply(x, 1L, max)
  out <- m + log(rowSums(exp(sweep(x, 1L, m, "-"))))
  out[!is.finite(out)] <- log(.Machine$double.xmin)
  out
}

.dstmm_log_bessel_k <- function(x, order) {
  # Scaled Bessel K is much more stable for large x:
  # besselK(x, nu, expon.scaled=TRUE) = exp(x) K_nu(x).
  x <- pmax(as.numeric(x), 1e-12)
  ord <- abs(as.numeric(order))
  val <- besselK(x, nu = ord, expon.scaled = TRUE)
  val <- pmax(val, .Machine$double.xmin)
  log(val) - x
}

.dstmm_logdmvt <- function(y, mu, sigma, nu, floor = 1e-8) {
  y <- as.matrix(y)
  p <- ncol(y)
  sigma <- .rdmm_make_spd(sigma, floor)
  cc <- chol(sigma)
  inv_sigma <- chol2inv(cc)
  centered <- sweep(y, 2L, mu, "-")
  d2 <- rowSums((centered %*% inv_sigma) * centered)
  logdet <- 2 * sum(log(diag(cc)))
  logdens <- lgamma((nu + p) / 2) - lgamma(nu / 2) -
    0.5 * (p * log(nu * pi) + logdet) -
    0.5 * (nu + p) * log1p(d2 / nu)
  list(logdens = logdens, d2 = d2, inv = inv_sigma, logdet = logdet)
}

.dstmm_logdghst <- function(y, mu, sigma, alpha, nu,
                            floor = 1e-8, alpha_tol = 1e-10) {
  # GHST density induced by
  #   Y | H=h ~ N(mu + h alpha, h Sigma),
  #   H ~ IG(nu/2, nu/2),
  # using the same parameterisation as the manuscript.
  y <- as.matrix(y)
  p <- ncol(y)
  sigma <- .rdmm_make_spd(sigma, floor)
  cc <- chol(sigma)
  inv_sigma <- chol2inv(cc)
  centered <- sweep(y, 2L, mu, "-")
  d2 <- rowSums((centered %*% inv_sigma) * centered)
  alpha <- as.numeric(alpha)
  A <- as.numeric(crossprod(alpha, inv_sigma %*% alpha))
  logdet <- 2 * sum(log(diag(cc)))

  if (!is.finite(A) || A <= alpha_tol) {
    tt <- .dstmm_logdmvt(y, mu, sigma, nu, floor)
    return(list(logdens = tt$logdens, d2 = d2, A = 0,
                inv = inv_sigma, cross = rep(0, nrow(y))))
  }

  cross_term <- as.numeric(centered %*% (inv_sigma %*% alpha))
  chi <- pmax(nu + d2, floor)
  psi <- A
  lambda <- -(nu + p) / 2
  x <- sqrt(psi * chi)

  # Integral of the GIG kernel:
  # 2 (chi/psi)^(lambda/2) K_lambda(sqrt(psi*chi)).
  log_integral <- log(2) + 0.5 * lambda * (log(chi) - log(psi)) +
    .dstmm_log_bessel_k(x, lambda)

  log_prior_const <- (nu / 2) * log(nu / 2) - lgamma(nu / 2)
  log_gaussian_const <- -0.5 * (p * log(2 * pi) + logdet)
  logdens <- log_prior_const + log_gaussian_const + cross_term + log_integral

  list(logdens = logdens, d2 = d2, A = A,
       inv = inv_sigma, cross = cross_term)
}

.dstmm_gig_moments <- function(psi, chi, lambda, alpha_tol = 1e-10) {
  # Returns E(H^-1) and E(log H) for GIG(psi, chi, lambda), where
  # f(h) proportional to h^(lambda-1) exp[-(psi h + chi/h)/2].
  chi <- pmax(as.numeric(chi), 1e-12)
  if (!is.finite(psi) || psi <= alpha_tol) {
    stop("Use the inverse-gamma limiting formulas when psi is zero.")
  }
  psi <- as.numeric(psi)
  x <- sqrt(psi * chi)

  log_k0 <- .dstmm_log_bessel_k(x, lambda)
  log_km1 <- .dstmm_log_bessel_k(x, lambda - 1)
  einvh <- sqrt(psi / chi) * exp(log_km1 - log_k0)

  # Numerical derivative with respect to the Bessel order for E(log H).
  h <- 1e-5
  log_k_plus <- .dstmm_log_bessel_k(x, lambda + h)
  log_k_minus <- .dstmm_log_bessel_k(x, lambda - h)
  dlogk <- (log_k_plus - log_k_minus) / (2 * h)
  elog <- 0.5 * (log(chi) - log(psi)) + dlogk

  list(einvh = einvh, elog = elog)
}

.dstmm_rinv_gamma <- function(n, shape, scale) {
  1 / rgamma(n, shape = shape, rate = scale)
}

.dstmm_rgig <- function(n, psi, chi, lambda, alpha_tol = 1e-10) {
  if (!is.finite(psi) || psi <= alpha_tol) {
    # In the DStMM posterior limit, lambda=-(nu+p)/2 and chi=nu+D.
    # h^(lambda-1) exp(-chi/(2h)) is IG(-lambda, chi/2).
    return(.dstmm_rinv_gamma(n, shape = -lambda, scale = chi / 2))
  }
  if (!requireNamespace("GIGrvg", quietly = TRUE)) {
    stop("Package 'GIGrvg' is required for nonzero skewness. Run scripts/00_install_packages.R")
  }
  as.numeric(GIGrvg::rgig(n, lambda = lambda, chi = chi, psi = psi))
}

.dstmm_expand_nu <- function(nu, structure, paths, k) {
  np <- nrow(paths)
  if (structure == "common") return(rep(nu[1L], np))
  if (structure == "first_layer") {
    if (length(nu) == 1L) nu <- rep(nu, k[1L])
    return(as.numeric(nu[paths[, 1L]]))
  }
  if (length(nu) == 1L) nu <- rep(nu, np)
  as.numeric(nu)
}

.dstmm_compact_nu <- function(nu_path, structure, paths, k) {
  if (structure == "common") return(mean(nu_path))
  if (structure == "first_layer") {
    return(vapply(seq_len(k[1L]), function(a) {
      mean(nu_path[paths[, 1L] == a])
    }, numeric(1)))
  }
  nu_path
}

.dstmm_active_layers <- function(active_set, layers) {
  active_set <- sort(unique(as.integer(active_set)))
  if (!length(active_set) || any(!active_set %in% seq_len(layers))) {
    stop("active_set must contain valid layer indices, e.g. active_set = 1L.")
  }
  active_set
}

.dstmm_zero_delta <- function(k, r_full) {
  lapply(seq_along(k), function(l) {
    matrix(0, nrow = r_full[l], ncol = k[l])
  })
}

.dstmm_collapse_path <- function(path, eta, lambda_list, psi, delta,
                                 r_full, floor = 1e-8) {
  # Lists are indexed so [[l+1]] corresponds to latent level l;
  # [[1]] is the observation level 0.
  h <- length(path)
  mu <- alpha <- Sigma <- vector("list", h + 1L)
  mu[[h + 1L]] <- rep(0, r_full[h + 1L])
  alpha[[h + 1L]] <- rep(0, r_full[h + 1L])
  Sigma[[h + 1L]] <- diag(r_full[h + 1L])

  for (l in h:1L) {
    a <- path[l]
    Lam <- matrix(lambda_list[[l]][a, , ],
                  nrow = r_full[l], ncol = r_full[l + 1L])
    Ps <- matrix(psi[[l]][a, , ], nrow = r_full[l], ncol = r_full[l])
    mu[[l]] <- as.numeric(eta[[l]][, a] + Lam %*% mu[[l + 1L]])
    alpha[[l]] <- as.numeric(delta[[l]][, a] + Lam %*% alpha[[l + 1L]])
    Sigma[[l]] <- .rdmm_make_spd(
      Ps + Lam %*% Sigma[[l + 1L]] %*% t(Lam), floor
    )
  }
  list(mu = mu, alpha = alpha, Sigma = Sigma)
}

.dstmm_estep <- function(y, eta, lambda_list, psi, delta, pi_list,
                         nu_path, paths, r_full, floor = 1e-8,
                         alpha_tol = 1e-10) {
  y <- as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  np <- nrow(paths)
  layers <- ncol(paths)

  log_joint <- matrix(NA_real_, n, np)
  D <- matrix(NA_real_, n, np)
  A <- numeric(np)
  einvh <- elogH <- matrix(NA_real_, n, np)
  collapse <- vector("list", np)

  for (s in seq_len(np)) {
    path <- paths[s, ]
    rec <- .dstmm_collapse_path(path, eta, lambda_list, psi, delta,
                                r_full, floor)
    collapse[[s]] <- rec
    mu0 <- rec$mu[[1L]]
    alpha0 <- rec$alpha[[1L]]
    Sigma0 <- rec$Sigma[[1L]]
    den <- .dstmm_logdghst(y, mu0, Sigma0, alpha0, nu_path[s],
                           floor, alpha_tol)
    D[, s] <- den$d2
    A[s] <- den$A

    prior <- prod(vapply(seq_len(layers), function(l) {
      pi_list[[l]][path[l]]
    }, numeric(1)))
    prior <- max(prior, .Machine$double.xmin)
    log_joint[, s] <- log(prior) + den$logdens

    if (A[s] <= alpha_tol) {
      shape <- (nu_path[s] + p) / 2
      scale_ig <- (nu_path[s] + D[, s]) / 2
      einvh[, s] <- shape / scale_ig
      elogH[, s] <- log(scale_ig) - digamma(shape)
    } else {
      mm <- .dstmm_gig_moments(
        psi = A[s], chi = nu_path[s] + D[, s],
        lambda = -(nu_path[s] + p) / 2,
        alpha_tol = alpha_tol
      )
      einvh[, s] <- mm$einvh
      elogH[, s] <- mm$elog
    }
  }

  log_py <- .dstmm_logsumexp_rows(log_joint)
  tau <- exp(sweep(log_joint, 1L, log_py, "-"))
  tau[!is.finite(tau)] <- 1 / np
  rs <- rowSums(tau)
  bad <- !is.finite(rs) | rs <= 0
  if (any(bad)) tau[bad, ] <- 1 / np
  tau <- sweep(tau, 1L, rowSums(tau), "/")

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
       layer_prob = layer_prob, s = classification,
       paths = paths, D = D, A = A,
       einvh = einvh, elogH = elogH, collapse = collapse)
}

.dstmm_draw_latent_one <- function(y, path_index, estep, eta, lambda_list,
                                   psi, delta, nu_path, r_full,
                                   floor = 1e-8, alpha_tol = 1e-10) {
  y <- as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  path <- estep$paths[path_index, ]
  rec <- estep$collapse[[path_index]]
  A <- estep$A[path_index]
  nu <- nu_path[path_index]
  lambda_gig <- -(nu + p) / 2

  if (A <= alpha_tol) {
    shape <- (nu + p) / 2
    scale_ig <- (nu + estep$D[, path_index]) / 2
    Hdraw <- 1 / rgamma(n, shape = shape, rate = scale_ig)
  } else {
    Hdraw <- vapply(seq_len(n), function(j) {
      .dstmm_rgig(1L, psi = A, chi = nu + estep$D[j, path_index],
                  lambda = lambda_gig, alpha_tol = alpha_tol)
    }, numeric(1))
  }
  Hdraw <- pmax(Hdraw, 1e-10)

  z <- vector("list", length(path))
  zprev <- y
  for (l in seq_along(path)) {
    a <- path[l]
    r_prev <- r_full[l]
    r_curr <- r_full[l + 1L]
    Lam <- matrix(lambda_list[[l]][a, , ], nrow = r_prev, ncol = r_curr)
    Ps <- .rdmm_make_spd(matrix(psi[[l]][a, , ], r_prev, r_prev), floor)
    invPs <- chol2inv(chol(Ps))
    Sig_l <- .rdmm_make_spd(rec$Sigma[[l + 1L]], floor)
    invSig <- chol2inv(chol(Sig_l))

    Omega <- .rdmm_make_spd(
      solve(invSig + t(Lam) %*% invPs %*% Lam), floor
    )

    mu_l <- rec$mu[[l + 1L]]
    alpha_l <- rec$alpha[[l + 1L]]
    eta_l <- eta[[l]][, a]
    delta_l <- delta[[l]][, a]

    base_prior <- as.numeric(invSig %*% mu_l)
    scale_prior <- as.numeric(invSig %*% alpha_l)
    residual <- sweep(zprev, 2L, eta_l, "-") -
      Hdraw * matrix(delta_l, nrow = n, ncol = r_prev, byrow = TRUE)
    M <- t(Lam) %*% invPs
    rhs <- matrix(base_prior, nrow = n, ncol = r_curr, byrow = TRUE) +
      Hdraw * matrix(scale_prior, nrow = n, ncol = r_curr, byrow = TRUE) +
      residual %*% t(M)
    rho <- rhs %*% t(Omega)

    noise <- matrix(rnorm(n * r_curr), nrow = n, ncol = r_curr) %*% chol(Omega)
    noise <- noise * sqrt(Hdraw)
    z[[l]] <- rho + noise
    zprev <- z[[l]]
  }

  list(H = Hdraw, z = z)
}

.dstmm_draw_missing <- function(y, estep, eta, lambda_list, psi, delta,
                                nu_path, r_full, M = 1L,
                                floor = 1e-8, alpha_tol = 1e-10) {
  M <- as.integer(M)
  if (M < 1L) stop("M must be at least 1.")
  np <- nrow(estep$paths)
  out <- vector("list", np)
  for (s in seq_len(np)) {
    out[[s]] <- lapply(seq_len(M), function(m) {
      .dstmm_draw_latent_one(
        y = y, path_index = s, estep = estep,
        eta = eta, lambda_list = lambda_list,
        psi = psi, delta = delta, nu_path = nu_path,
        r_full = r_full, floor = floor, alpha_tol = alpha_tol
      )
    })
  }
  out
}

.dstmm_mstep <- function(y, estep, draws, k, r_full, active_set,
                         psi_floor = 1e-6, ridge = 1e-8) {
  y <- as.matrix(y)
  layers <- length(k)
  paths <- estep$paths
  M <- length(draws[[1L]])

  eta <- lambda_list <- psi <- psi_inv <- pi_list <- delta <- vector("list", layers)

  for (l in seq_len(layers)) {
    r_prev <- r_full[l]
    r_curr <- r_full[l + 1L]
    is_active <- l %in% active_set
    qreg <- 1L + r_curr + if (is_active) 1L else 0L

    eta[[l]] <- matrix(0, nrow = r_prev, ncol = k[l])
    lambda_list[[l]] <- array(0, dim = c(k[l], r_prev, r_curr))
    psi[[l]] <- psi_inv[[l]] <- array(0, dim = c(k[l], r_prev, r_prev))
    delta[[l]] <- matrix(0, nrow = r_prev, ncol = k[l])
    pi_list[[l]] <- colMeans(estep$layer_prob[[l]])
    pi_list[[l]] <- pmax(pi_list[[l]], 1e-12)
    pi_list[[l]] <- pi_list[[l]] / sum(pi_list[[l]])

    for (a in seq_len(k[l])) {
      Q <- matrix(0, qreg, qreg)
      R <- matrix(0, r_prev, qreg)
      N <- 0

      for (s in which(paths[, l] == a)) {
        tau <- estep$tau[, s]
        N <- N + sum(tau)
        for (m in seq_len(M)) {
          dr <- draws[[s]][[m]]
          h <- pmax(dr$H, 1e-10)
          zcurr <- dr$z[[l]]
          xprev <- if (l == 1L) y else dr$z[[l - 1L]]
          u <- if (is_active) cbind(1, zcurr, h) else cbind(1, zcurr)
          wt <- tau / h / M
          Q <- Q + crossprod(u, u * wt)
          R <- R + crossprod(xprev, u * wt)
        }
      }

      Q <- (Q + t(Q)) / 2 + diag(ridge, qreg)
      B <- R %*% solve(Q)
      eta[[l]][, a] <- B[, 1L]
      lambda_list[[l]][a, , ] <- B[, 1L + seq_len(r_curr), drop = FALSE]
      if (is_active) delta[[l]][, a] <- B[, qreg]

      S <- matrix(0, r_prev, r_prev)
      for (s in which(paths[, l] == a)) {
        tau <- estep$tau[, s]
        for (m in seq_len(M)) {
          dr <- draws[[s]][[m]]
          h <- pmax(dr$H, 1e-10)
          zcurr <- dr$z[[l]]
          xprev <- if (l == 1L) y else dr$z[[l - 1L]]
          fitted <- matrix(eta[[l]][, a], nrow = nrow(y), ncol = r_prev, byrow = TRUE) +
            zcurr %*% t(matrix(lambda_list[[l]][a, , ], r_prev, r_curr))
          if (is_active) {
            fitted <- fitted + h * matrix(delta[[l]][, a],
                                           nrow = nrow(y), ncol = r_prev, byrow = TRUE)
          }
          e <- xprev - fitted
          wt <- tau / h / M
          S <- S + crossprod(e, e * wt)
        }
      }
      S <- S / max(N, .Machine$double.eps)
      diagv <- pmax(diag((S + t(S)) / 2), psi_floor)
      psi[[l]][a, , ] <- diag(diagv, r_prev)
      psi_inv[[l]][a, , ] <- diag(1 / diagv, r_prev)
    }
  }

  list(eta = eta, lambda = lambda_list, psi = psi,
       psi.inv = psi_inv, pi = pi_list, delta = delta)
}

.dstmm_update_one_nu <- function(target_mean, old, bounds) {
  f <- function(v) {
    log(v / 2) - digamma(v / 2) + 1 - target_mean
  }
  lo <- bounds[1L]
  hi <- bounds[2L]
  flo <- f(lo)
  fhi <- f(hi)
  ans <- if (is.finite(flo) && is.finite(fhi) && flo * fhi <= 0) {
    uniroot(f, interval = c(lo, hi), tol = 1e-7)$root
  } else {
    optimize(function(v) f(v)^2, interval = c(lo, hi))$minimum
  }
  if (!is.finite(ans)) old else ans
}

.dstmm_update_nu <- function(estep, nu_path, structure, k, bounds,
                             min_effective = 5) {
  paths <- estep$paths
  np <- nrow(paths)
  updated <- nu_path

  update_group <- function(indices, old) {
    ww <- estep$tau[, indices, drop = FALSE]
    denom <- sum(ww)
    if (!is.finite(denom) || denom < min_effective) return(old)
    target <- sum(ww * (estep$elogH[, indices, drop = FALSE] +
                          estep$einvh[, indices, drop = FALSE])) / denom
    .dstmm_update_one_nu(target, old, bounds)
  }

  if (structure == "common") {
    updated[] <- update_group(seq_len(np), mean(nu_path))
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

.dstmm_parameter_count <- function(k, r_full, active_set, nu_structure) {
  layers <- length(k)
  total <- 0
  for (l in seq_len(layers)) {
    r_prev <- r_full[l]
    r_curr <- r_full[l + 1L]
    # Same count as the RDMM implementation: loadings + intercept + diagonal Psi,
    # minus rotational degrees of freedom; plus active skewness vectors.
    total <- total + (k[l] - 1L) +
      k[l] * (r_prev * r_curr + r_prev + r_prev) -
      k[l] * r_curr * (r_curr - 1L) / 2
    if (l %in% active_set) total <- total + k[l] * r_prev
  }
  total + switch(nu_structure,
                 common = 1L,
                 first_layer = k[1L],
                 pathway = prod(k))
}

.dstmm_entropy <- function(tau) {
  z <- pmax(tau, .Machine$double.xmin)
  -sum(z * log(z))
}

.dstmm_init_from_rdmm <- function(warm_start, k, r_full, active_set) {
  if (is.null(warm_start$H) || is.null(warm_start$mu) ||
      is.null(warm_start$psi) || is.null(warm_start$w)) {
    stop("warm_start must be an RDMM fit with H, mu, psi and w.")
  }
  delta <- .dstmm_zero_delta(k, r_full)
  list(lambda = warm_start$H,
       eta = warm_start$mu,
       psi = warm_start$psi,
       pi = warm_start$w,
       delta = delta,
       nu = if (!is.null(warm_start$nu)) warm_start$nu else 10)
}

.dstmm_validate <- function(y, layers, k, r, it, eps, nu, nu_structure,
                            psi_floor, M, nu_bounds) {
  if (is.data.frame(y)) y <- as.matrix(y)
  if (!is.matrix(y) || !is.numeric(y) || anyNA(y) || any(!is.finite(y))) {
    stop("y must be a finite numeric matrix.")
  }
  if (!(layers %in% 1:3)) stop("layers must be 1, 2, or 3.")
  if (length(k) != layers || length(r) != layers) {
    stop("k and r must have one entry per layer.")
  }
  if (any(r >= ncol(y)) || (length(r) > 1L && any(diff(r) >= 0))) {
    stop("Require p > r[1] > ... > r[layers] >= 1.")
  }
  if (it < 1 || it != as.integer(it)) stop("it must be a positive integer.")
  if (!is.finite(eps) || eps <= 0) stop("eps must be positive.")
  if (!(nu_structure %in% c("common", "first_layer", "pathway"))) {
    stop("nu_structure must be common, first_layer, or pathway.")
  }
  if (any(!is.finite(nu)) || any(nu <= 2)) stop("Initial nu must exceed 2.")
  if (any(!is.finite(nu_bounds)) || length(nu_bounds) != 2L ||
      nu_bounds[1L] <= 2 || nu_bounds[2L] <= nu_bounds[1L]) {
    stop("nu_bounds must be c(lower, upper) with 2 < lower < upper.")
  }
  if (!is.finite(psi_floor) || psi_floor <= 0) stop("psi_floor must be positive.")
  if (M < 1 || M != as.integer(M)) stop("M must be a positive integer.")
  invisible(TRUE)
}

#' Fit a Deep Skew-t Mixture Model
#'
#' @param y n x p numeric data matrix.
#' @param layers Number of latent mixture layers.
#' @param k Number of local components per layer.
#' @param r Latent dimensions.
#' @param active_set Layers allowed to have nonzero delta. Manuscript default: 1.
#' @param M Monte Carlo draws per observation-path pair (M=1 gives stochastic EM).
#' @param warm_start Optional rdmm object. Recommended for the manuscript simulations.
#' @return Object of class "dstmm".
dstmm <- function(y, layers, k, r,
                  it = 250, eps = 1e-3,
                  init = "kmeans", init_est = "factanal",
                  seed = NULL, scale = FALSE,
                  nu = 10,
                  nu_structure = c("common", "first_layer", "pathway"),
                  estimate_nu = TRUE,
                  active_set = 1L,
                  M = 1L,
                  warm_start = NULL,
                  psi_floor = 1e-6,
                  nu_bounds = c(2.05, 200),
                  min_iter = 50L,
                  moving_window = 20L,
                  alpha_tol = 1e-10,
                  ridge = 1e-8,
                  verbose = FALSE) {
  .dstmm_require_rdmm()
  call <- match.call()
  if (is.data.frame(y)) y <- as.matrix(y)
  nu_structure <- match.arg(nu_structure)
  fixed <- .rdmm_fix_names(init, init_est)
  init <- fixed$init
  init_est <- fixed$init_est
  active_set <- .dstmm_active_layers(active_set, layers)
  .dstmm_validate(y, layers, k, r, it, eps, nu, nu_structure,
                  psi_floor, M, nu_bounds)

  if (!is.null(seed)) set.seed(as.integer(seed))

  original_data <- y
  center <- rep(0, ncol(y))
  scale_value <- rep(1, ncol(y))
  if (scale) {
    ys <- base::scale(y)
    center <- attr(ys, "scaled:center")
    scale_value <- attr(ys, "scaled:scale")
    scale_value[!is.finite(scale_value) | scale_value == 0] <- 1
    y <- sweep(sweep(y, 2L, center, "-"), 2L, scale_value, "/")
  }

  n <- nrow(y)
  p <- ncol(y)
  k <- as.integer(k)
  r_full <- c(p, as.integer(r))
  paths <- .rdmm_paths(k)

  expected_nu <- switch(nu_structure,
                        common = 1L,
                        first_layer = k[1L],
                        pathway = prod(k))
  if (!(length(nu) %in% c(1L, expected_nu))) {
    stop("nu has incompatible length for nu_structure.")
  }

  if (!is.null(warm_start)) {
    if (scale) stop("For a supplied warm_start use scale = FALSE in dstmm().")
    ini <- .dstmm_init_from_rdmm(warm_start, k, r_full, active_set)
    eta <- ini$eta
    lambda_list <- ini$lambda
    psi <- ini$psi
    pi_list <- ini$pi
    delta <- ini$delta
    if (!is.null(warm_start$nu)) nu <- warm_start$nu
  } else {
    ini0 <- .rdmm_initialise(y, layers, k, r_full, init, init_est, psi_floor)
    eta <- ini0$mu
    lambda_list <- ini0$H
    psi <- ini0$psi
    pi_list <- ini0$w
    delta <- .dstmm_zero_delta(k, r_full)
  }

  nu_path <- .dstmm_expand_nu(nu, nu_structure, paths, k)
  likelihood <- numeric(0)
  best <- NULL
  best_loglik <- -Inf
  ratio <- Inf
  converged <- FALSE
  min_iter <- as.integer(min_iter)
  moving_window <- as.integer(moving_window)

  for (iteration in seq_len(as.integer(it))) {
    estep <- .dstmm_estep(y, eta, lambda_list, psi, delta, pi_list,
                          nu_path, paths, r_full, psi_floor, alpha_tol)
    draws <- .dstmm_draw_missing(y, estep, eta, lambda_list, psi, delta,
                                 nu_path, r_full, M, psi_floor, alpha_tol)
    upd <- .dstmm_mstep(y, estep, draws, k, r_full, active_set,
                        psi_floor, ridge)
    new_nu <- if (estimate_nu) {
      .dstmm_update_nu(estep, nu_path, nu_structure, k, nu_bounds)
    } else nu_path

    new_estep <- .dstmm_estep(y, upd$eta, upd$lambda, upd$psi, upd$delta,
                              upd$pi, new_nu, paths, r_full,
                              psi_floor, alpha_tol)
    likelihood <- c(likelihood, new_estep$loglik)

    eta <- upd$eta
    lambda_list <- upd$lambda
    psi <- upd$psi
    pi_list <- upd$pi
    delta <- upd$delta
    nu_path <- new_nu

    if (is.finite(new_estep$loglik) && new_estep$loglik > best_loglik) {
      best_loglik <- new_estep$loglik
      best <- list(eta = eta, lambda = lambda_list, psi = psi,
                   pi = pi_list, delta = delta, nu_path = nu_path,
                   estep = new_estep, iteration = iteration)
    }

    if (length(likelihood) >= 2L * moving_window && iteration >= min_iter) {
      old_idx <- (length(likelihood) - 2L * moving_window + 1L):
        (length(likelihood) - moving_window)
      old_ma <- mean(likelihood[old_idx])
      new_ma <- mean(tail(likelihood, moving_window))
      ratio <- abs(new_ma - old_ma) / (abs(old_ma) + 1)
    }
    if (verbose) {
      message(sprintf("DStMM iteration %d: loglik=%.6f ratio=%.3g",
                      iteration, tail(likelihood, 1L), ratio))
    }
    if (iteration >= min_iter && is.finite(ratio) && ratio < eps) {
      converged <- TRUE
      break
    }
  }

  if (is.null(best)) stop("DStMM failed to produce a finite likelihood.")
  eta <- best$eta
  lambda_list <- best$lambda
  psi <- best$psi
  pi_list <- best$pi
  delta <- best$delta
  nu_path <- best$nu_path
  final <- best$estep

  hpar <- .dstmm_parameter_count(k, r_full, active_set, nu_structure)
  lik <- final$loglik
  bic <- -2 * lik + hpar * log(n)
  aic <- -2 * lik + 2 * hpar
  entropy <- .dstmm_entropy(final$tau)
  icl_bic <- bic + 2 * entropy

  path_alpha <- lapply(seq_len(nrow(paths)), function(s) {
    final$collapse[[s]]$alpha[[1L]]
  })
  path_mu <- lapply(seq_len(nrow(paths)), function(s) {
    final$collapse[[s]]$mu[[1L]]
  })
  path_sigma <- lapply(seq_len(nrow(paths)), function(s) {
    final$collapse[[s]]$Sigma[[1L]]
  })

  out <- list(
    eta = eta, lambda = lambda_list, psi = psi, delta = delta,
    pi = pi_list, w = pi_list,
    nu = .dstmm_compact_nu(nu_path, nu_structure, paths, k),
    nu_path = nu_path,
    likelihood = likelihood, lik = likelihood, loglik = lik,
    bic = bic, aic = aic, icl_bic = icl_bic,
    h = hpar, entropy = entropy,
    s = final$s, ps.y = final$tau, ps.y.list = final$layer_prob,
    paths = paths, path_alpha = path_alpha,
    path_mu = path_mu, path_sigma = path_sigma,
    iterations = best$iteration, convergence_ratio = ratio,
    converged = converged,
    active_set = active_set, M = as.integer(M),
    k = k, r = r_full[-1L], layers = layers,
    nu_structure = nu_structure, estimate_nu = estimate_nu,
    scaled = scale, center = center, scale = scale_value,
    original_data = original_data, seed = seed, call = call
  )
  class(out) <- "dstmm"
  invisible(out)
}

predict.dstmm <- function(object, newdata = NULL,
                          type = c("path", "layer", "posterior"),
                          layer = 1L, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    if (type == "path") return(max.col(object$ps.y, ties.method = "first"))
    if (type == "layer") return(object$s[, layer])
    return(object$ps.y)
  }
  y <- as.matrix(newdata)
  if (object$scaled) {
    y <- sweep(sweep(y, 2L, object$center, "-"), 2L, object$scale, "/")
  }
  r_full <- c(ncol(y), object$r)
  ee <- .dstmm_estep(y, object$eta, object$lambda, object$psi,
                     object$delta, object$pi, object$nu_path,
                     object$paths, r_full)
  if (type == "path") return(max.col(ee$tau, ties.method = "first"))
  if (type == "layer") return(ee$s[, layer])
  ee$tau
}

print.dstmm <- function(x, ...) {
  cat("Deep Skew-t Mixture Model (DStMM)\n")
  cat("Layers:", x$layers, "\n")
  cat("Components:", paste(x$k, collapse = " x "), "\n")
  cat("Latent dimensions:", paste(x$r, collapse = " > "), "\n")
  cat("Active skewness layers:", paste(x$active_set, collapse = ", "), "\n")
  cat("Monte Carlo size M:", x$M, "\n")
  cat("Degrees of freedom:", paste(round(x$nu, 3), collapse = ", "), "\n")
  cat("Log-likelihood:", round(x$loglik, 3), "\n")
  cat("BIC:", round(x$bic, 3), "\n")
  invisible(x)
}

summary.dstmm <- function(object, ...) {
  print(object)
  cat("AIC:", round(object$aic, 3), "\n")
  cat("ICL-BIC:", round(object$icl_bic, 3), "\n")
  cat("Best iteration:", object$iterations, "\n")
  cat("Final convergence ratio:", signif(object$convergence_ratio, 4), "\n")
  invisible(object)
}
