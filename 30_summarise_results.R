rm(list = ls())
root <- normalizePath(if (file.exists("R/load_all.R")) "." else "..")
source(file.path(root, "R", "load_all.R")); load_dstmm_project(root)

files <- c("experiment1_raw.csv", "experiment1_robustness_raw.csv", "experiment2_raw.csv")
paths <- file.path(root, "results", files)
existing <- paths[file.exists(paths)]
if (!length(existing)) stop("No raw result files found. Run simulations first.")
raw <- do.call(rbind, lapply(existing, read.csv, stringsAsFactors = FALSE))

summary_all <- build_model_summary(raw)
write.csv(summary_all, file.path(root, "results", "model_summary_long.csv"), row.names = FALSE)

# Manuscript-friendly compact summaries. BIC is intentionally excluded.
keep_metrics <- c("Layer1_ARI", "Layer1_MR", "Pathway_ARI", "Pathway_MR", "Elapsed_seconds")
write.csv(summary_all[summary_all$Metric %in% keep_metrics, ],
          file.path(root, "results", "clustering_and_fit_summary.csv"), row.names = FALSE)

# nu bias/RMSE from raw replicates.
nu_models <- raw[raw$Model %in% c("RDMM", "DStMM") & raw$Success, ]
if (nrow(nu_models)) {
  keys <- unique(nu_models[, c("Experiment", "Design", "N", "Nu", "Kappa", "Pi", "Model")])
  nr <- lapply(seq_len(nrow(keys)), function(i) {
    keep <- rep(TRUE, nrow(nu_models))
    for (g in names(keys)) keep <- keep & nu_models[[g]] == keys[[g]][i]
    d <- nu_models[keep, ]
    err <- d$Estimated_nu - d$Nu
    cbind(keys[i, ], Nu_Bias = mean(err, na.rm = TRUE),
          Nu_RMSE = sqrt(mean(err^2, na.rm = TRUE)), N_success = sum(is.finite(err)))
  })
  write.csv(do.call(rbind, nr), file.path(root, "results", "nu_recovery_summary.csv"), row.names = FALSE)
}

# DStMM skewness recovery.
ds <- raw[raw$Model == "DStMM" & raw$Success, ]
if (nrow(ds)) {
  groups <- c("Experiment", "Design", "N", "Nu", "Kappa", "Pi", "Model")
  dr <- summarise_metric(ds, groups, "Delta_RMSE")
  ar <- summarise_metric(ds, groups, "Alpha_RMSE")
  write.csv(rbind(dr, ar), file.path(root, "results", "skewness_recovery_summary.csv"), row.names = FALSE)
}

# Remove an old BIC selection summary if it exists, so stale output is not mistaken for current output.
old_bic <- file.path(root, "results", "bic_selection_summary.csv")
if (file.exists(old_bic)) unlink(old_bic)

cat("Summaries written to results/. BIC summaries are not produced.\n")
