rm(list = ls())
root <- normalizePath(if (file.exists("R/load_all.R")) "." else "..")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")
library(ggplot2)

dir.create(file.path(root, "figures"), showWarnings = FALSE, recursive = TRUE)

read_if_exists <- function(name) {
  f <- file.path(root, "results", name)
  if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
}

e1 <- read_if_exists("experiment1_raw.csv")
rob <- read_if_exists("experiment1_robustness_raw.csv")
e2 <- read_if_exists("experiment2_raw.csv")
if (is.null(e1) && is.null(rob) && is.null(e2)) {
  stop("No raw simulation results found. Run the simulation scripts first.")
}

save_both <- function(plot, stem, width = 7.8, height = 5.2) {
  ggsave(file.path(root, "figures", paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 300)
  ggsave(file.path(root, "figures", paste0(stem, ".pdf")), plot,
         width = width, height = height)
}

clean_success <- function(d) {
  if (is.null(d)) return(NULL)
  if ("Success" %in% names(d)) d <- d[d$Success, , drop = FALSE]
  d
}

set_model_factor <- function(d) {
  if (is.null(d)) return(NULL)
  d$Model <- factor(d$Model, levels = c("DGMM", "RDMM", "DStMM"))
  d
}

assert_three_models <- function(d, label, condition_cols) {
  if (is.null(d) || !nrow(d)) return(invisible(TRUE))
  req <- c("DGMM", "RDMM", "DStMM")
  keys <- unique(d[, condition_cols, drop = FALSE])
  bad <- character(0)
  for (i in seq_len(nrow(keys))) {
    keep <- rep(TRUE, nrow(d))
    for (g in condition_cols) keep <- keep & d[[g]] == keys[[g]][i]
    present <- unique(as.character(d$Model[keep]))
    miss <- setdiff(req, present)
    if (length(miss)) {
      desc <- paste(paste0(condition_cols, "=", unlist(keys[i, condition_cols, drop = FALSE])), collapse = ", ")
      bad <- c(bad, paste0(desc, " missing ", paste(miss, collapse = "/")))
    }
  }
  if (length(bad)) {
    stop(label, " is missing successful model results. Run scripts/05_repair_dgmm_results.R first.\n",
         paste(head(bad, 12L), collapse = "\n"))
  }
  invisible(TRUE)
}

three_model_boxplot <- function(data, yvar, ylab, title,
                                facet_var = NULL,
                                width = 7.8, height = 5.2) {
  p <- ggplot(data, aes(x = KappaF, y = .data[[yvar]], fill = Model)) +
    geom_boxplot(position = position_dodge(width = 0.82), width = 0.72,
                 outlier.size = 0.65, na.rm = TRUE) +
    labs(x = expression(kappa), y = ylab, title = title, fill = "Model") +
    theme_bw() +
    theme(legend.position = "top")

  if (!is.null(facet_var)) {
    if (facet_var == "NuF") {
      p <- p + facet_wrap(~ NuF, nrow = 1,
                          labeller = labeller(NuF = function(x) paste0("nu = ", x)))
    } else if (facet_var == "Scenario") {
      p <- p + facet_wrap(~ Scenario, nrow = 1)
    }
  }
  p
}

# ============================================================
# Experiment 1
# ARI and MR: DGMM / RDMM / DStMM shown together in every boxplot
# ============================================================
e1 <- set_model_factor(clean_success(e1))
if (!is.null(e1) && nrow(e1)) {
  assert_three_models(e1, "Experiment 1", c("N", "Nu", "Kappa", "Pi"))
  e1$KappaF <- factor(e1$Kappa, levels = sort(unique(e1$Kappa)))

  p1 <- three_model_boxplot(
    e1, "Layer1_ARI", "First-layer ARI",
    "Experiment 1: first-layer ARI"
  )
  save_both(p1, "exp1_first_layer_ari_boxplot")

  p2 <- three_model_boxplot(
    e1, "Layer1_MR", "First-layer MR",
    "Experiment 1: first-layer misclassification rate"
  )
  save_both(p2, "exp1_first_layer_mr_boxplot")

  p3 <- three_model_boxplot(
    e1, "Pathway_ARI", "Pathway ARI",
    "Experiment 1: pathway ARI"
  )
  save_both(p3, "exp1_pathway_ari_boxplot")

  p4 <- three_model_boxplot(
    e1, "Pathway_MR", "Pathway MR",
    "Experiment 1: pathway misclassification rate"
  )
  save_both(p4, "exp1_pathway_mr_boxplot")

  ds <- e1[e1$Model == "DStMM" & is.finite(e1$Delta_RMSE), , drop = FALSE]
  if (nrow(ds)) {
    p5 <- ggplot(ds, aes(x = KappaF, y = Delta_RMSE)) +
      geom_boxplot(width = 0.65, outlier.size = 0.65, na.rm = TRUE) +
      labs(x = expression(kappa), y = expression(RMSE[delta]),
           title = "Experiment 1: skewness-parameter recovery") +
      theme_bw()
    save_both(p5, "exp1_delta_rmse_boxplot")
  }

  nu_dat <- e1[e1$Model %in% c("RDMM", "DStMM") & is.finite(e1$Estimated_nu), , drop = FALSE]
  if (nrow(nu_dat)) {
    p6 <- ggplot(nu_dat, aes(x = KappaF, y = Estimated_nu, fill = Model)) +
      geom_boxplot(position = position_dodge(width = 0.82), width = 0.72,
                   outlier.size = 0.65, na.rm = TRUE) +
      geom_hline(yintercept = unique(e1$Nu)[1], linetype = 2) +
      labs(x = expression(kappa), y = expression(hat(nu)), fill = "Model",
           title = "Experiment 1: degrees-of-freedom recovery") +
      theme_bw() + theme(legend.position = "top")
    save_both(p6, "exp1_nu_recovery_boxplot")
  }
}

# ============================================================
# Experiment 1 robustness checks
# Again show DGMM / RDMM / DStMM together for ARI and MR.
# ============================================================
rob <- set_model_factor(clean_success(rob))
if (!is.null(rob) && nrow(rob)) {
  assert_three_models(rob, "Robustness experiment", c("Design", "N", "Nu", "Kappa", "Pi"))
  rob$KappaF <- factor(rob$Kappa, levels = sort(unique(rob$Kappa)))
  rob$Scenario <- ifelse(
    rob$Design == "sample_size",
    "n = 500, balanced",
    "n = 1000, unbalanced"
  )
  rob$Scenario <- factor(rob$Scenario,
                         levels = c("n = 500, balanced", "n = 1000, unbalanced"))

  p7 <- three_model_boxplot(
    rob, "Layer1_ARI", "First-layer ARI",
    "Robustness checks: first-layer ARI", facet_var = "Scenario"
  )
  save_both(p7, "exp1_robustness_first_layer_ari_boxplot", width = 10.2, height = 5.1)

  p8 <- three_model_boxplot(
    rob, "Layer1_MR", "First-layer MR",
    "Robustness checks: first-layer misclassification rate", facet_var = "Scenario"
  )
  save_both(p8, "exp1_robustness_first_layer_mr_boxplot", width = 10.2, height = 5.1)

  p9 <- three_model_boxplot(
    rob, "Pathway_ARI", "Pathway ARI",
    "Robustness checks: pathway ARI", facet_var = "Scenario"
  )
  save_both(p9, "exp1_robustness_pathway_ari_boxplot", width = 10.2, height = 5.1)

  p10 <- three_model_boxplot(
    rob, "Pathway_MR", "Pathway MR",
    "Robustness checks: pathway misclassification rate", facet_var = "Scenario"
  )
  save_both(p10, "exp1_robustness_pathway_mr_boxplot", width = 10.2, height = 5.1)
}

# ============================================================
# Experiment 2
# ARI and MR: all three models, faceted by nu.
# ============================================================
e2 <- set_model_factor(clean_success(e2))
if (!is.null(e2) && nrow(e2)) {
  assert_three_models(e2, "Experiment 2", c("N", "Nu", "Kappa", "Pi"))
  e2$KappaF <- factor(e2$Kappa, levels = sort(unique(e2$Kappa)))
  e2$NuF <- factor(e2$Nu, levels = sort(unique(e2$Nu)))

  p11 <- three_model_boxplot(
    e2, "Layer1_ARI", "First-layer ARI",
    "Experiment 2: first-layer ARI", facet_var = "NuF"
  )
  save_both(p11, "exp2_first_layer_ari_boxplot", width = 10.6, height = 5.2)

  p12 <- three_model_boxplot(
    e2, "Layer1_MR", "First-layer MR",
    "Experiment 2: first-layer misclassification rate", facet_var = "NuF"
  )
  save_both(p12, "exp2_first_layer_mr_boxplot", width = 10.6, height = 5.2)

  p13 <- three_model_boxplot(
    e2, "Pathway_ARI", "Pathway ARI",
    "Experiment 2: pathway ARI", facet_var = "NuF"
  )
  save_both(p13, "exp2_pathway_ari_boxplot", width = 10.6, height = 5.2)

  p14 <- three_model_boxplot(
    e2, "Pathway_MR", "Pathway MR",
    "Experiment 2: pathway misclassification rate", facet_var = "NuF"
  )
  save_both(p14, "exp2_pathway_mr_boxplot", width = 10.6, height = 5.2)

  ds2 <- e2[e2$Model == "DStMM" & is.finite(e2$Delta_RMSE), , drop = FALSE]
  if (nrow(ds2)) {
    p15 <- ggplot(ds2, aes(x = KappaF, y = Delta_RMSE)) +
      geom_boxplot(width = 0.65, outlier.size = 0.6, na.rm = TRUE) +
      facet_wrap(~ NuF, nrow = 1,
                 labeller = labeller(NuF = function(x) paste0("nu = ", x))) +
      labs(x = expression(kappa), y = expression(RMSE[delta]),
           title = "Experiment 2: skewness-parameter recovery") +
      theme_bw()
    save_both(p15, "exp2_delta_rmse_boxplot", width = 10.2, height = 5.0)
  }

  nu2 <- e2[e2$Model %in% c("RDMM", "DStMM") & is.finite(e2$Estimated_nu), , drop = FALSE]
  if (nrow(nu2)) {
    p16 <- ggplot(nu2, aes(x = KappaF, y = Estimated_nu, fill = Model)) +
      geom_boxplot(position = position_dodge(width = 0.82), width = 0.72,
                   outlier.size = 0.6, na.rm = TRUE) +
      facet_wrap(~ NuF, nrow = 1, scales = "free_y",
                 labeller = labeller(NuF = function(x) paste0("nu = ", x))) +
      labs(x = expression(kappa), y = expression(hat(nu)), fill = "Model",
           title = "Experiment 2: degrees-of-freedom recovery") +
      theme_bw() + theme(legend.position = "top")
    save_both(p16, "exp2_nu_recovery_boxplot", width = 10.6, height = 5.2)
  }
}

# ============================================================
# Runtime: all three models
# ============================================================
time_parts <- list()
if (!is.null(e1) && nrow(e1)) {
  x <- e1[, c("Model", "Elapsed_seconds"), drop = FALSE]
  x$Experiment <- "Experiment 1"
  time_parts[[length(time_parts) + 1L]] <- x
}
if (!is.null(e2) && nrow(e2)) {
  x <- e2[, c("Model", "Elapsed_seconds"), drop = FALSE]
  x$Experiment <- "Experiment 2"
  time_parts[[length(time_parts) + 1L]] <- x
}
if (length(time_parts)) {
  td <- do.call(rbind, time_parts)
  td <- td[is.finite(td$Elapsed_seconds), , drop = FALSE]
  td$Model <- factor(td$Model, levels = c("DGMM", "RDMM", "DStMM"))
  if (nrow(td)) {
    p17 <- ggplot(td, aes(x = Model, y = Elapsed_seconds, fill = Model)) +
      geom_boxplot(width = 0.65, outlier.size = 0.6, na.rm = TRUE) +
      facet_wrap(~ Experiment, scales = "free_y") +
      labs(x = NULL, y = "Elapsed seconds", title = "Computational cost") +
      theme_bw() + theme(legend.position = "none")
    save_both(p17, "runtime_comparison_boxplot", width = 7.6, height = 4.8)
  }
}

# Remove obsolete two-model gain plot from older runs.
old_gain <- file.path(root, "figures", c(
  "exp1_robustness_ari_gain_boxplot.png",
  "exp1_robustness_ari_gain_boxplot.pdf"
))
unlink(old_gain[file.exists(old_gain)])

cat("Boxplot figures written to figures/ as PNG and PDF.\n")
cat("ARI and MR plots compare DGMM, RDMM, and DStMM side-by-side.\n")
