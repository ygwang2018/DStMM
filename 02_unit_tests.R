rm(list = ls())
root <- normalizePath(if (file.exists("R/load_all.R")) "." else "..")
source(file.path(root, "R", "load_all.R")); load_dstmm_project(root)

set.seed(123)
y <- matrix(rnorm(30), nrow = 10, ncol = 3)
mu <- c(0.2, -0.1, 0.3)
Sigma <- crossprod(matrix(rnorm(9), 3, 3)) + diag(3)
nu <- 6

# 1. Exact Student-t limit at alpha = 0.
a <- .dstmm_logdghst(y, mu, Sigma, rep(0, 3), nu)$logdens
b <- .dstmm_logdmvt(y, mu, Sigma, nu)$logdens
stopifnot(max(abs(a - b)) < 1e-10)

# 2. Path skewness vanishes when local deltas vanish.
pars0 <- make_dstmm_parameters(nu = 6, kappa = 0)
# Build arrays in fitted-model format.
lam <- list(
  array(unlist(pars0$lambda[[1L]]), dim = c(2, 20, 5)),
  array(unlist(pars0$lambda[[2L]]), dim = c(2, 5, 2))
)
# Rebuild carefully because list->array memory ordering differs from the desired component slices.
lam[[1L]] <- array(0, c(2, 20, 5)); lam[[2L]] <- array(0, c(2, 5, 2))
for (aidx in 1:2) { lam[[1L]][aidx,,] <- pars0$lambda[[1L]][[aidx]]; lam[[2L]][aidx,,] <- pars0$lambda[[2L]][[aidx]] }
psi <- list(array(0, c(2,20,20)), array(0, c(2,5,5)))
for (aidx in 1:2) { psi[[1L]][aidx,,] <- pars0$psi[[1L]][[aidx]]; psi[[2L]][aidx,,] <- pars0$psi[[2L]][[aidx]] }
rec <- .dstmm_collapse_path(c(1L,1L), pars0$eta, lam, psi, pars0$delta, c(20,5,2))
stopifnot(max(abs(rec$alpha[[1L]])) < 1e-12)

# 3. DStMM adds exactly 2*p active skewness parameters to the matched RDMM here.
paths <- .rdmm_paths(c(2,2))
h_rd <- .rdmm_parameter_count(c(2,2), c(20,5,2), "common")
h_ds <- .dstmm_parameter_count(c(2,2), c(20,5,2), 1L, "common")
stopifnot(h_ds - h_rd == 40)

# 4. Simulated shared H is positive and path encoding has four states.
sim <- simulate_dstmm_data(200, make_dstmm_parameters(nu = 6, kappa = 1), seed = 123)
stopifnot(all(sim$H > 0), all(sim$path_id %in% 1:4))

cat("All deterministic unit tests passed.\n")
