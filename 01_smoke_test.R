rm(list = ls())
root <- normalizePath(if (file.exists("R/load_all.R")) "." else "..")
source(file.path(root, "R", "load_all.R"))
load_dstmm_project(root)

pars <- make_dstmm_parameters(nu = 6, kappa = 1)
sim <- simulate_dstmm_data(n = 150L, parameters = pars, seed = 20260819L)

cat("Smoke test: fitting DGMM, RDMM, DStMM on n=150...\n")
out <- fit_three_models(
  sim = sim, seed = 20260820L,
  max_iter = 25L, eps = 5e-3,
  initial_nu = 6, min_iter = 10L, moving_window = 5L, M = 1L
)
print(out$results)

stopifnot(nrow(out$results) == 3L)
stopifnot(all(out$results$Model == c("DGMM", "RDMM", "DStMM")))
write.csv(out$results, file.path(root, "results", "smoke_test_results.csv"), row.names = FALSE)
if (!all(out$results$Success)) {
  bad <- out$results[!out$results$Success, c("Model", "Error"), drop = FALSE]
  print(bad)
  stop("Smoke test failed: at least one of DGMM/RDMM/DStMM did not fit successfully.")
}
cat("Smoke test completed: DGMM, RDMM and DStMM all succeeded.\n")
