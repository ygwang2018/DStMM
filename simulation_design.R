# Simulation design for the DStMM manuscript.
# Implements Sections 4.1--4.3 plus two lightweight robustness checks:
# n = 500 versus 1000, and balanced versus (0.7, 0.3) layer proportions.

make_dstmm_parameters <- function(nu = 6, kappa = 0,
                                  pi1 = c(0.5, 0.5),
                                  pi2 = c(0.5, 0.5)) {
  p <- 20L
  r1 <- 5L
  r2 <- 2L
  k1 <- 2L
  k2 <- 2L

  if (nu <= 2) stop("nu must exceed 2 because the centring correction uses E(H).")
  if (length(pi1) != 2L || length(pi2) != 2L ||
      any(pi1 <= 0) || any(pi2 <= 0)) stop("pi1 and pi2 must be positive length-2 vectors.")
  pi1 <- pi1 / sum(pi1)
  pi2 <- pi2 / sum(pi2)

  # Baseline first-layer locations from the manuscript.
  eta1_0 <- matrix(0, nrow = p, ncol = k1)
  eta1_0[1:3, 1L] <- -2.5
  eta1_0[1:3, 2L] <-  2.5

  # Second-layer locations.
  eta2 <- matrix(0, nrow = r1, ncol = k2)
  eta2[1L, 1L] <-  1.25
  eta2[1L, 2L] <- -1.25

  # Layer-1 loading matrices.
  a <- matrix(c(1.0, 0.9, 1.1, 0.8), ncol = 1L)
  A <- kronecker(diag(r1), a)
  D2 <- diag(c(1.10, 0.90, 1.05, 0.95, 1.00))
  lambda1 <- list(A, A %*% D2)

  # Layer-2 loading matrices.
  lambda2 <- list(
    matrix(c(
      1.0,  0.0,
      0.8,  0.2,
      0.0,  1.0,
      0.2,  0.8,
      0.6, -0.6
    ), nrow = r1, byrow = TRUE),
    matrix(c(
       0.9,  0.1,
       0.7, -0.2,
       0.1,  0.9,
      -0.2,  0.7,
       0.5,  0.5
    ), nrow = r1, byrow = TRUE)
  )

  psi1 <- replicate(k1, 0.5 * diag(p), simplify = FALSE)
  psi2 <- replicate(k2, 0.3 * diag(r1), simplify = FALSE)

  # Manuscript skewness direction v = (1,1,1,1,0,...,0)' / 2.
  v <- numeric(p)
  v[1:4] <- 0.5
  delta1 <- matrix(kappa * v, nrow = p, ncol = k1)
  delta2 <- matrix(0, nrow = r1, ncol = k2)

  # Centre correction: E(H) = nu / (nu - 2).
  # Both first-layer components receive the same correction so changing kappa
  # changes shape without intentionally shifting their unconditional centres.
  EH <- nu / (nu - 2)
  eta1 <- eta1_0 - EH * delta1

  path_labels <- as.matrix(expand.grid(
    layer1 = seq_len(k1),
    layer2 = seq_len(k2),
    KEEP.OUT.ATTRS = FALSE
  ))

  list(
    p = p,
    r = c(r1, r2),
    k = c(k1, k2),
    pi = list(pi1, pi2),
    nu = as.numeric(nu),
    kappa = as.numeric(kappa),
    eta0 = list(eta1_0, eta2),
    eta = list(eta1, eta2),
    lambda = list(lambda1, lambda2),
    psi = list(psi1, psi2),
    delta = list(delta1, delta2),
    skew_direction = v,
    path_labels = path_labels,
    path_probabilities = as.vector(outer(pi1, pi2))
  )
}

.sim_make_spd <- function(x, floor = 1e-10) {
  x <- (x + t(x)) / 2
  ee <- eigen(x, symmetric = TRUE)
  ee$values <- pmax(ee$values, floor)
  ee$vectors %*% (ee$values * t(ee$vectors))
}

.sim_rmvnorm <- function(n, mean, covariance) {
  covariance <- .sim_make_spd(covariance)
  z <- matrix(rnorm(n * length(mean)), nrow = n)
  sweep(z %*% chol(covariance), 2L, mean, "+")
}

encode_path <- function(path_matrix, k = apply(path_matrix, 2L, max)) {
  path_matrix <- as.matrix(path_matrix)
  if (ncol(path_matrix) == 1L) return(as.integer(path_matrix[, 1L]))
  # Consistent with expand.grid: first layer varies fastest.
  mult <- cumprod(c(1L, head(as.integer(k), -1L)))
  as.integer(1L + rowSums((path_matrix - 1L) * matrix(mult,
                                                      nrow(path_matrix),
                                                      ncol(path_matrix),
                                                      byrow = TRUE)))
}

simulate_dstmm_data <- function(n = 1000L,
                                parameters = make_dstmm_parameters(),
                                seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  n <- as.integer(n)
  p <- parameters$p
  r1 <- parameters$r[1L]
  r2 <- parameters$r[2L]
  k1 <- parameters$k[1L]
  k2 <- parameters$k[2L]
  nu <- parameters$nu

  s1 <- sample.int(k1, n, replace = TRUE, prob = parameters$pi[[1L]])
  s2 <- sample.int(k2, n, replace = TRUE, prob = parameters$pi[[2L]])
  pathway <- cbind(s1, s2)

  # One H ~ IG(nu/2, nu/2) shared across the whole observation path.
  H <- 1 / rgamma(n, shape = nu / 2, rate = nu / 2)

  z2 <- matrix(rnorm(n * r2), nrow = n, ncol = r2)
  z2 <- z2 * sqrt(H)

  z1 <- matrix(NA_real_, nrow = n, ncol = r1)
  for (b in seq_len(k2)) {
    idx <- which(s2 == b)
    if (!length(idx)) next
    Lam <- parameters$lambda[[2L]][[b]]
    eta <- parameters$eta[[2L]][, b]
    delta <- parameters$delta[[2L]][, b]
    mean_part <- sweep(z2[idx, , drop = FALSE] %*% t(Lam), 2L, eta, "+") +
      H[idx] * matrix(delta, nrow = length(idx), ncol = r1, byrow = TRUE)
    noise <- .sim_rmvnorm(length(idx), rep(0, r1), parameters$psi[[2L]][[b]])
    noise <- noise * sqrt(H[idx])
    z1[idx, ] <- mean_part + noise
  }

  y <- matrix(NA_real_, nrow = n, ncol = p)
  for (a in seq_len(k1)) {
    idx <- which(s1 == a)
    if (!length(idx)) next
    Lam <- parameters$lambda[[1L]][[a]]
    eta <- parameters$eta[[1L]][, a]
    delta <- parameters$delta[[1L]][, a]
    mean_part <- sweep(z1[idx, , drop = FALSE] %*% t(Lam), 2L, eta, "+") +
      H[idx] * matrix(delta, nrow = length(idx), ncol = p, byrow = TRUE)
    noise <- .sim_rmvnorm(length(idx), rep(0, p), parameters$psi[[1L]][[a]])
    noise <- noise * sqrt(H[idx])
    y[idx, ] <- mean_part + noise
  }

  colnames(y) <- paste0("Y", seq_len(p))
  path_id <- encode_path(pathway, parameters$k)

  list(
    y = y, z1 = z1, z2 = z2, H = H,
    layer1 = s1, layer2 = s2,
    pathway = pathway, path_id = path_id,
    parameters = parameters
  )
}

# Main Experiment 1: fixed nu=6, increasing skewness.
experiment1_grid <- function() {
  expand.grid(
    Design = "main",
    N = 1000L,
    Nu = 6,
    Kappa = c(0, 0.5, 1.0, 1.5, 2.0),
    Pi = "balanced",
    stringsAsFactors = FALSE
  )
}

# Two lightweight robustness checks agreed for the revised design.
experiment1_robustness_grid <- function() {
  rbind(
    expand.grid(
      Design = "sample_size",
      N = 500L,
      Nu = 6,
      Kappa = c(0, 0.5, 1.0, 1.5, 2.0),
      Pi = "balanced",
      stringsAsFactors = FALSE
    ),
    expand.grid(
      Design = "unbalanced",
      N = 1000L,
      Nu = 6,
      Kappa = c(0, 0.5, 1.0, 1.5, 2.0),
      Pi = "unbalanced",
      stringsAsFactors = FALSE
    )
  )
}

# Experiment 2: tail weight x skewness factorial design.
experiment2_grid <- function() {
  expand.grid(
    Design = "factorial",
    N = 1000L,
    Nu = c(5, 8, 15),
    Kappa = c(0, 1, 2),
    Pi = "balanced",
    stringsAsFactors = FALSE
  )
}

pi_from_label <- function(label) {
  if (identical(label, "balanced")) {
    list(c(0.5, 0.5), c(0.5, 0.5))
  } else if (identical(label, "unbalanced")) {
    list(c(0.7, 0.3), c(0.7, 0.3))
  } else {
    stop("Unknown Pi label: ", label)
  }
}
