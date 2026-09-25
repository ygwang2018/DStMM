# ============================================================================
# Case Study 1: UCI Multiple Features (mfeat-fac), digits 0-9
# Ten-class high-dimensional clustering for DGMM / RDMM / DStMM
# Deep-model-only run: K2 is fixed at 1; no external baseline methods are run.
# ============================================================================
# Goal
#   * Observations: all 2000 handwritten-digit patterns from the UCI Multiple
#     Features data set: 200 observations for each digit 0,1,...,9.
#   * Variables: ALL 216 profile-correlation features in mfeat-fac.
#   * No supervised feature selection is performed.
#   * The original digit labels are used ONLY after fitting to compute ARI/MR.
#   * First-layer allocation is the reported clustering result (K1 = 10).
#   * DGMM / RDMM / DStMM use two latent layers and the same high-dimensional
#     initialization framework supplied by the project.
#   * Only DGMM, RDMM and DStMM are fitted in this script.
#   * K2 is fixed at 1. Latent dimensions (r1,r2) can still be compared by BIC.
#   * In DStMM, first-layer skewness vectors delta^(1)_a are estimated
#     separately for a = 1,...,K1; no equality constraint is imposed across
#     digit clusters.
#
# UCI data ordering
#   The mfeat-fac file has 2000 rows. UCI documents that rows 1:200 are digit 0,
#   rows 201:400 are digit 1, ..., rows 1801:2000 are digit 9. Therefore this
#   script uses all ten original digit classes 0,...,9. No class subset is
#   selected using ARI/MR.
#
# Reproducibility
#   * The downloaded mfeat-fac file is cached locally.
#   * The exact all-digit analysis matrix is written to CSV.
#   * Deep-model results are checkpointed after each model/start.
#   * RDMM fit objects are cached so an interrupted DStMM run can resume from
#     the matching RDMM warm start.
# ============================================================================

rm(list = ls())

# ----------------------------------------------------------------------------
# 1. USER CONFIGURATION
# ----------------------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  farg <- grep("^--file=", args, value = TRUE)
  if (length(farg)) {
    return(dirname(normalizePath(sub("^--file=", "", farg[1L]),
                                winslash = "/", mustWork = TRUE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(rstudioapi::isAvailable())) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(dirname(normalizePath(ctx$path, winslash = "/", mustWork = TRUE)))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

PROJECT_ROOT <- get_script_dir()
INSTALL_MISSING_PACKAGES <- TRUE
AUTO_INSTALL_OPTIONAL_PACKAGES <- FALSE
MASTER_SEED <- 20260825L

# Fixed BEFORE fitting. All original UCI digit classes are retained; there is
# no search over class subsets using ARI/MR.
CASE1_DIGITS <- 0:9
K1 <- length(CASE1_DIGITS)

# Run debug first; switch to paper only after DGMM/RDMM/DStMM all succeed.
RUN_MODE <- "paper"  # "debug" or "paper"

# Data download/cache controls.
FORCE_REFRESH_DATA <- FALSE
MFEAT_FAC_URL <- "https://archive.ics.uci.edu/ml/machine-learning-databases/mfeat/mfeat-fac"
MFEAT_TAR_URL <- "https://archive.ics.uci.edu/ml/machine-learning-databases/mfeat/mfeat.tar"

# Deep-model controls shared by DGMM/RDMM/DStMM.
EPS <- 1e-4
MOVING_WINDOW <- 20L
MC_SIZE <- 1L
NU_INIT <- 8
NU_BOUNDS <- c(2.05, 200)
PSI_FLOOR <- 1e-3
DEEP_SCALE <- FALSE  # X is standardized explicitly below
RESUME_FROM_CHECKPOINT <- TRUE
RETRY_FAILED_FITS <- FALSE

# DStMM skewness specification. active_set = 1 means that each first-layer
# local component has its own freely estimated delta^(1)_a. The model does NOT
# force delta^(1)_1 = ... = delta^(1)_K1.
DSTMM_ACTIVE_SET <- 1L

if (identical(RUN_MODE, "debug")) {
  N_DEEP_STARTS <- 1L
  MAX_ITER <- 60L
  MIN_ITER <- 25L
  K2_GRID <- 1L
  LATENT_GRID <- data.frame(r1 = 10L, r2 = 3L)
} else if (identical(RUN_MODE, "paper")) {
  N_DEEP_STARTS <- 3L
  MAX_ITER <- 200L
  MIN_ITER <- 60L
  K2_GRID <- 1L
  LATENT_GRID <- data.frame(
    r1 = c(10L, 15L, 20L, 25L),
    r2 = c( 3L,  5L,  5L,  8L)
  )
} else {
  stop("RUN_MODE must be 'debug' or 'paper'.")
}

# Baseline controls for the genuine ten-class comparison.
RUN_BASELINES <- FALSE
RUN_FACTOR_BASELINES <- FALSE
RUN_SLOW_SKEWT_BASELINES <- FALSE
RUN_SANITY_CHECKS <- FALSE  # do not run K-means/Ward quick checks in deep-only mode
BASELINE_PCA_DIM <- 20L
BASELINE_NSTART <- if (RUN_MODE == "debug") 2L else 3L
BASELINE_MAX_ITER <- if (RUN_MODE == "debug") 80L else 150L
FA_Q_GRID <- 2:4
PGMM_Q_GRID <- 2:4
CFUST_Q_GRID <- 1:2
USTFA_Q_GRID <- 2:4   # unrestricted skew-t factor analyzer (uskewFactors)
USTFA_NSTART <- if (RUN_MODE == "debug") 1L else 2L

# ----------------------------------------------------------------------------
# 2. PROJECT SETUP AND PACKAGES
# ----------------------------------------------------------------------------

if (!dir.exists(PROJECT_ROOT)) stop("PROJECT_ROOT does not exist: ", PROJECT_ROOT)
setwd(PROJECT_ROOT)
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)

RESULT_DIR <- file.path(PROJECT_ROOT, "results",
                        "case1_mfeat_fac_digits_0_to_9_k2_1_component_delta")
DATA_DIR <- file.path(PROJECT_ROOT, "data", "case1_mfeat_fac")
FIT_CACHE_DIR <- file.path(RESULT_DIR, "fit_cache")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIT_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

core_packages <- c("mvtnorm", "corpcor", "GIGrvg", "clue")
optional_packages <- c(
  "mclust", "teigen", "pgmm", "EMMIXmfa", "uskewFactors",
  "EMMIXskew", "EMMIXuskew", "EMMIXcskew"
)

install_cran_if_needed <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && INSTALL_MISSING_PACKAGES) {
    message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
    try(install.packages(missing, dependencies = TRUE), silent = FALSE)
  }
  invisible(NULL)
}
install_cran_if_needed(core_packages)
if (AUTO_INSTALL_OPTIONAL_PACKAGES) install_cran_if_needed(optional_packages)

still_missing_core <- core_packages[
  !vapply(core_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(still_missing_core)) {
  stop("Required packages are missing: ", paste(still_missing_core, collapse = ", "))
}

source(file.path(PROJECT_ROOT, "R", "load_all.R"))
load_dstmm_project(PROJECT_ROOT)
set.seed(MASTER_SEED)

# ----------------------------------------------------------------------------
# 3. DOWNLOAD / LOAD UCI MFEAT-FAC
# ----------------------------------------------------------------------------

MFEAT_FILE <- file.path(DATA_DIR, "mfeat-fac")

obtain_mfeat_fac <- function() {
  if (file.exists(MFEAT_FILE) && !FORCE_REFRESH_DATA) {
    message("Using cached mfeat-fac: ", MFEAT_FILE)
    return(MFEAT_FILE)
  }

  if (file.exists(MFEAT_FILE) && FORCE_REFRESH_DATA) unlink(MFEAT_FILE)

  message("Downloading mfeat-fac from UCI...")
  direct_ok <- tryCatch({
    utils::download.file(MFEAT_FAC_URL, MFEAT_FILE,
                         mode = "wb", quiet = FALSE)
    file.exists(MFEAT_FILE) && file.info(MFEAT_FILE)$size > 100000
  }, error = function(e) FALSE, warning = function(w) FALSE)

  if (!isTRUE(direct_ok)) {
    if (file.exists(MFEAT_FILE)) unlink(MFEAT_FILE)
    message("Direct download failed; trying the UCI mfeat.tar archive...")
    tar_file <- file.path(DATA_DIR, "mfeat.tar")
    tar_ok <- tryCatch({
      utils::download.file(MFEAT_TAR_URL, tar_file,
                           mode = "wb", quiet = FALSE)
      file.exists(tar_file) && file.info(tar_file)$size > 100000
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!isTRUE(tar_ok)) {
      stop(
        "Could not download mfeat-fac from UCI.\n",
        "Manually download the UCI Multiple Features file 'mfeat-fac' and place it at:\n",
        MFEAT_FILE
      )
    }
    untar_ok <- tryCatch({
      utils::untar(tar_file, exdir = DATA_DIR)
      TRUE
    }, error = function(e) FALSE)
    if (!untar_ok || !file.exists(MFEAT_FILE)) {
      # Some archives extract into a nested directory; search once by basename.
      hits <- list.files(DATA_DIR, pattern = "^mfeat-fac$", recursive = TRUE,
                         full.names = TRUE)
      if (length(hits)) file.copy(hits[1L], MFEAT_FILE, overwrite = TRUE)
    }
  }

  if (!file.exists(MFEAT_FILE)) stop("mfeat-fac was not obtained successfully.")
  MFEAT_FILE
}

mfeat_path <- obtain_mfeat_fac()

mfeat_all <- utils::read.table(
  mfeat_path,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = "",
  quote = ""
)
mfeat_all <- as.matrix(mfeat_all)
storage.mode(mfeat_all) <- "double"

# UCI documents 200 patterns per class, ordered 0,1,...,9.
if (nrow(mfeat_all) != 2000L) {
  stop("Expected 2000 rows in mfeat-fac, found ", nrow(mfeat_all), ".")
}
if (ncol(mfeat_all) != 216L) {
  stop("Expected 216 profile-correlation features, found ", ncol(mfeat_all), ".")
}
if (any(!is.finite(mfeat_all))) {
  stop("Non-finite values found in mfeat-fac; UCI documents no missing values.")
}

all_digit <- rep(0:9, each = 200L)
keep <- all_digit %in% CASE1_DIGITS
X_raw_full <- mfeat_all[keep, , drop = FALSE]
digit <- all_digit[keep]

truth_factor <- factor(digit, levels = CASE1_DIGITS,
                       labels = paste0("Digit ", CASE1_DIGITS))
truth <- as.integer(truth_factor)
truth_name <- as.character(truth_factor)

expected_counts <- rep(200L, K1)
if (!identical(as.integer(table(factor(digit, levels = CASE1_DIGITS))), expected_counts)) {
  stop("The expected 200 observations per digit class (0,...,9) were not obtained.")
}

# No feature ranking or supervised screening. Remove only columns that are
# numerically constant within the fixed all-digit study population, because such
# columns cannot be standardized and contain no clustering information.
feature_sd <- apply(X_raw_full, 2L, stats::sd)
good_feature <- is.finite(feature_sd) & feature_sd > 1e-12
if (any(!good_feature)) {
  warning("Dropping ", sum(!good_feature),
          " numerically degenerate mfeat-fac columns; no feature ranking was used.")
}
X_raw <- X_raw_full[, good_feature, drop = FALSE]
retained_feature_index <- which(good_feature)

# Standardize each retained feature across all 2000 observations WITHOUT labels.
X <- scale(X_raw, center = TRUE, scale = TRUE)
X <- as.matrix(X)
storage.mode(X) <- "double"
if (any(!is.finite(X))) stop("Non-finite values remain after standardization.")

feature_map <- data.frame(
  OriginalFeature = seq_len(ncol(X_raw_full)),
  Retained = good_feature,
  AnalysisColumn = NA_integer_,
  stringsAsFactors = FALSE
)
feature_map$AnalysisColumn[good_feature] <- seq_len(sum(good_feature))
utils::write.csv(feature_map,
                 file.path(DATA_DIR, "mfeat_fac_feature_map.csv"), row.names = FALSE)

utils::write.csv(
  data.frame(Digit = digit, X_raw, check.names = FALSE),
  file.path(DATA_DIR, "mfeat_fac_digits_0_to_9_raw_analysis_matrix.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(Digit = digit, X, check.names = FALSE),
  file.path(DATA_DIR, "mfeat_fac_digits_0_to_9_standardized_matrix.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("CASE 1 DATA LOADED: UCI Multiple Features / mfeat-fac\n")
cat("Multiclass task: digits 0-9 (10 original classes)\n")
cat("n =", nrow(X), "; original p =", ncol(X_raw_full),
    "; analysis p =", ncol(X), "; K1 =", K1, "\n")
cat("True class counts (used only for external evaluation):\n")
print(table(truth_name))
cat("============================================================\n\n")

# ----------------------------------------------------------------------------
# 4. UNSUPERVISED DISTRIBUTIONAL DIAGNOSTICS
# ----------------------------------------------------------------------------

sample_skewness <- function(z) {
  z <- z[is.finite(z)]
  if (length(z) < 3L) return(NA_real_)
  s <- stats::sd(z)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  mean((z - mean(z))^3) / s^3
}

sample_excess_kurtosis <- function(z) {
  z <- z[is.finite(z)]
  if (length(z) < 4L) return(NA_real_)
  s <- stats::sd(z)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  mean((z - mean(z))^4) / s^4 - 3
}

feature_diag <- data.frame(
  Feature = retained_feature_index,
  Skewness = apply(X_raw, 2L, sample_skewness),
  ExcessKurtosis = apply(X_raw, 2L, sample_excess_kurtosis),
  OutlierRate3SD = apply(X_raw, 2L, function(z) {
    m <- mean(z); s <- stats::sd(z)
    if (!is.finite(s) || s <= 0) return(NA_real_)
    mean(abs(z - m) > 3 * s)
  }),
  stringsAsFactors = FALSE
)
utils::write.csv(feature_diag,
                 file.path(RESULT_DIR, "distribution_diagnostics_by_feature.csv"),
                 row.names = FALSE)

diag_summary <- data.frame(
  Statistic = c(
    "median_abs_feature_skewness",
    "median_feature_excess_kurtosis",
    "median_feature_3sd_outlier_rate",
    "proportion_abs_skewness_gt_1",
    "proportion_excess_kurtosis_gt_3"
  ),
  Value = c(
    stats::median(abs(feature_diag$Skewness), na.rm = TRUE),
    stats::median(feature_diag$ExcessKurtosis, na.rm = TRUE),
    stats::median(feature_diag$OutlierRate3SD, na.rm = TRUE),
    mean(abs(feature_diag$Skewness) > 1, na.rm = TRUE),
    mean(feature_diag$ExcessKurtosis > 3, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(diag_summary,
                 file.path(RESULT_DIR, "distribution_diagnostics_summary.csv"),
                 row.names = FALSE)
print(diag_summary, row.names = FALSE)

writeLines(c(
  "Dataset: UCI Multiple Features, mfeat-fac representation",
  "UCI dataset ID: 72",
  paste0("Source file: ", MFEAT_FAC_URL),
  "UCI mfeat-fac dimension: 2000 observations x 216 profile-correlation features",
  "UCI ordering: 200 consecutive observations per digit class 0,...,9",
  "All 2000 observations and all 10 digit classes are retained.",
  "Fixed multiclass task: all original digits 0,1,...,9",
  paste0("Analysis n: ", nrow(X)),
  paste0("Original p in mfeat-fac: ", ncol(X_raw_full)),
  paste0("Analysis p after only degenerate-column safeguard: ", ncol(X)),
  paste0("Degenerate columns removed: ", sum(!good_feature)),
  "No supervised feature selection was performed.",
  "No class subset or class pairing was selected using ARI/MR.",
  "All retained variables were standardized without using class labels.",
  "True digit labels were used only after fitting to compute ARI/MR.",
  "The first-layer allocation is the reported clustering result.",
  paste0("First-layer cluster count K1: ", K1),
  paste0("Common deep-model psi floor: ", PSI_FLOOR),
  paste0("Run mode: ", RUN_MODE)
), file.path(RESULT_DIR, "data_provenance.txt"))

# ----------------------------------------------------------------------------
# 5. LABEL-INVARIANT MULTICLASS METRICS AND COMMON HELPERS
# ----------------------------------------------------------------------------

normalize_partition <- function(pred) {
  if (is.factor(pred)) pred <- as.character(pred)
  as.integer(factor(pred))
}

misclassification_rate_hungarian <- function(pred, truth, groups) {
  pred <- normalize_partition(pred)
  truth <- as.integer(truth)
  tab <- table(
    factor(pred, levels = seq_len(groups)),
    factor(truth, levels = seq_len(groups))
  )
  tab <- as.matrix(tab)
  assignment <- clue::solve_LSAP(tab, maximum = TRUE)
  matched <- sum(tab[cbind(seq_len(groups), as.integer(assignment))])
  1 - matched / length(truth)
}

score_partition <- function(pred) {
  pred <- normalize_partition(pred)
  data.frame(
    ARI = adjusted_rand_index_local(truth, pred),
    MR = misclassification_rate_hungarian(pred, truth, groups = K1),
    N_clusters = length(unique(pred)),
    stringsAsFactors = FALSE
  )
}

# Relabel arbitrary cluster IDs to the digit labels with a Hungarian assignment.
# This mapping is used only for reporting confusion matrices/MR after fitting;
# it never feeds back into model estimation or model selection.
hungarian_relabel <- function(pred, truth, groups) {
  pred_n <- normalize_partition(pred)
  truth_i <- as.integer(truth)
  tab <- table(
    factor(pred_n, levels = seq_len(groups)),
    factor(truth_i, levels = seq_len(groups))
  )
  assignment <- clue::solve_LSAP(as.matrix(tab), maximum = TRUE)
  map <- as.integer(assignment)
  mapped <- map[pred_n]
  list(
    mapped = mapped,
    mapping = data.frame(
      Cluster = seq_len(groups),
      TruthIndex = map,
      Digit = CASE1_DIGITS[map],
      stringsAsFactors = FALSE
    ),
    confusion = table(
      PredictedDigit = factor(CASE1_DIGITS[mapped], levels = CASE1_DIGITS),
      TrueDigit = factor(CASE1_DIGITS[truth_i], levels = CASE1_DIGITS)
    )
  )
}

safe_num <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  y <- suppressWarnings(as.numeric(x))
  if (!length(y)) NA_real_ else y[1L]
}

extract_deep_loglik <- function(fit, model) {
  if (model == "DGMM") return(tail(as.numeric(fit$lik), 1L))
  safe_num(fit$loglik)
}

extract_min_psi <- function(fit) {
  if (is.null(fit$psi)) return(NA_real_)
  vals <- unlist(lapply(fit$psi, function(a) {
    if (length(dim(a)) == 3L) {
      unlist(lapply(seq_len(dim(a)[1L]), function(g) {
        diag(a[g, , , drop = FALSE][1L, , ])
      }))
    } else if (is.matrix(a)) {
      diag(a)
    } else {
      as.numeric(a)
    }
  }))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) NA_real_ else min(vals)
}

extract_nu_mean <- function(fit) {
  if (!is.null(fit$nu)) return(mean(as.numeric(fit$nu), na.rm = TRUE))
  if (!is.null(fit$nu_path)) return(mean(as.numeric(fit$nu_path), na.rm = TRUE))
  NA_real_
}

make_architecture_grid <- function() {
  rows <- list()
  ii <- 0L
  for (k2 in K2_GRID) {
    for (jj in seq_len(nrow(LATENT_GRID))) {
      ii <- ii + 1L
      rows[[ii]] <- data.frame(
        K1 = K1,
        K2 = as.integer(k2),
        r1 = as.integer(LATENT_GRID$r1[jj]),
        r2 = as.integer(LATENT_GRID$r2[jj]),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out <- out[out$r1 > out$r2 & out$r1 < ncol(X), , drop = FALSE]
  out$ArchID <- sprintf("K%d-%d_R%d-%d", out$K1, out$K2, out$r1, out$r2)
  rownames(out) <- NULL
  out[, c("ArchID", "K1", "K2", "r1", "r2")]
}

ARCH_GRID <- make_architecture_grid()
utils::write.csv(ARCH_GRID, file.path(RESULT_DIR, "architecture_grid.csv"), row.names = FALSE)
cat("\nArchitecture grid:\n")
print(ARCH_GRID)
cat("Planned deep fits =", nrow(ARCH_GRID) * N_DEEP_STARTS * 3L, "\n")

# Optional quick unsupervised sanity checks. Disabled in this deep-only script.
if (RUN_SANITY_CHECKS) {
  cat("\nQuick unsupervised sanity checks:\n")
  set.seed(MASTER_SEED)
  km_check <- stats::kmeans(X, centers = K1, nstart = 50L, iter.max = 200L)$cluster
  km_score <- score_partition(km_check)
  cat(sprintf("  K-means: ARI=%.4f, MR=%.4f\n", km_score$ARI, km_score$MR))
  wd_check <- stats::cutree(stats::hclust(stats::dist(X), method = "ward.D2"), k = K1)
  wd_score <- score_partition(wd_check)
  cat(sprintf("  Ward.D2: ARI=%.4f, MR=%.4f\n", wd_score$ARI, wd_score$MR))
}

# ----------------------------------------------------------------------------
# 7. DGMM / RDMM / DStMM FITTING
# ----------------------------------------------------------------------------

fit_dgmm_arch <- function(y, k, r, seed) {
  set.seed(seed)
  deepgmm(
    y = y, layers = 2L, k = k, r = r,
    it = MAX_ITER, eps = EPS,
    init = "kmeans", init_est = "factanal",
    seed = seed, scale = DEEP_SCALE,
    psi_floor = PSI_FLOOR
  )
}

fit_rdmm_arch <- function(y, k, r, seed) {
  set.seed(seed)
  robustdeepgmm(
    y = y, layers = 2L, k = k, r = r,
    it = MAX_ITER, eps = EPS,
    init = "kmeans", init_est = "factanal",
    seed = seed, scale = DEEP_SCALE,
    nu = NU_INIT, nu_structure = "common", estimate_nu = TRUE,
    method = "sem", psi_floor = PSI_FLOOR, nu_bounds = NU_BOUNDS,
    min_iter = MIN_ITER, moving_window = MOVING_WINDOW,
    verbose = FALSE
  )
}

fit_dstmm_arch <- function(y, k, r, seed, rdmm_warm = NULL) {
  set.seed(seed)
  dstmm(
    y = y, layers = 2L, k = k, r = r,
    it = MAX_ITER, eps = EPS,
    init = "kmeans", init_est = "factanal",
    seed = seed, scale = DEEP_SCALE,
    nu = NU_INIT, nu_structure = "common", estimate_nu = TRUE,
    # Component-specific first-layer skewness: delta^(1)_a is updated
    # separately for every a = 1,...,K1 inside dstmm().
    active_set = DSTMM_ACTIVE_SET, M = MC_SIZE,
    warm_start = rdmm_warm,
    psi_floor = PSI_FLOOR, nu_bounds = NU_BOUNDS,
    min_iter = MIN_ITER, moving_window = MOVING_WINDOW,
    verbose = FALSE
  )
}

# Extract the first-layer DStMM skewness vectors for diagnostics only.
# The project implementation stores local parameters by layer; this helper is
# intentionally tolerant of common matrix/list orientations and does not alter
# estimation. It lets us verify/report that K1 separate delta vectors were fit.
extract_first_layer_delta_matrix <- function(fit, K_expected = K1) {
  d <- fit$delta
  if (is.null(d)) return(NULL)

  d1 <- if (is.list(d)) d[[1L]] else d
  if (is.null(d1)) return(NULL)

  if (is.vector(d1) && !is.list(d1)) {
    # A single vector would indicate one shared delta, which is not the desired
    # real-data specification when K1 > 1.
    return(matrix(as.numeric(d1), nrow = 1L))
  }

  if (is.matrix(d1)) {
    if (nrow(d1) == K_expected) return(d1)
    if (ncol(d1) == K_expected) return(t(d1))
    return(d1)
  }

  if (is.array(d1) && length(dim(d1)) == 2L) {
    m <- as.matrix(d1)
    if (nrow(m) == K_expected) return(m)
    if (ncol(m) == K_expected) return(t(m))
    return(m)
  }

  if (is.list(d1) && length(d1) == K_expected) {
    lens <- vapply(d1, length, integer(1L))
    if (length(unique(lens)) == 1L) {
      return(do.call(rbind, lapply(d1, as.numeric)))
    }
  }

  NULL
}

save_dstmm_delta_diagnostics <- function(fit, arch, start_id) {
  dm <- extract_first_layer_delta_matrix(fit, K_expected = arch$K1)
  if (is.null(dm)) {
    warning("Could not extract first-layer DStMM delta vectors from the fit object; ",
            "the fit is retained, but no delta diagnostic CSV was written.")
    return(invisible(NULL))
  }

  # For the intended specification there must be one row per first-layer local
  # component. We stop if a single shared vector is detected, because that would
  # contradict the requested component-specific-delta fit.
  if (nrow(dm) != arch$K1) {
    stop("DStMM delta diagnostic found ", nrow(dm),
         " first-layer delta vector(s), but K1 = ", arch$K1,
         ". Expected one separately estimated delta vector per cluster.")
  }

  delta_df <- data.frame(
    Cluster = seq_len(arch$K1),
    DeltaNorm = sqrt(rowSums(dm^2)),
    dm,
    check.names = FALSE
  )
  names(delta_df)[-(1:2)] <- paste0("delta_", seq_len(ncol(dm)))

  fn <- file.path(
    RESULT_DIR,
    sprintf("DStMM_delta_%s_start%d.csv", arch$ArchID, start_id)
  )
  utils::write.csv(delta_df, fn, row.names = FALSE)
  invisible(delta_df)
}


make_deep_result <- function(fit, model, arch, start_id, seed, elapsed) {
  pred <- as.integer(fit$s[, 1L])
  ss <- score_partition(pred)
  min_psi <- extract_min_psi(fit)

  converged <- NA
  if (model == "DStMM" && !is.null(fit$converged)) {
    converged <- isTRUE(fit$converged)
  } else if (!is.null(fit$convergence_ratio)) {
    converged <- is.finite(fit$convergence_ratio) && fit$convergence_ratio < EPS
  }

  data.frame(
    ArchID = as.character(arch$ArchID),
    K1 = arch$K1, K2 = arch$K2, r1 = arch$r1, r2 = arch$r2,
    Start = start_id, Seed = seed, Model = model,
    ARI = ss$ARI, MR = ss$MR, N_clusters = ss$N_clusters,
    BIC = safe_num(fit$bic), ICL_BIC = safe_num(fit$icl_bic),
    LogLik = extract_deep_loglik(fit, model),
    MeanNu = if (model == "DGMM") NA_real_ else extract_nu_mean(fit),
    MinPsi = min_psi,
    AtPsiFloor = is.finite(min_psi) && min_psi <= PSI_FLOOR * 1.001,
    Iterations = if (!is.null(fit$iterations)) as.integer(fit$iterations) else NA_integer_,
    Converged = converged,
    Elapsed_seconds = elapsed,
    Success = TRUE,
    Error = "",
    stringsAsFactors = FALSE
  )
}

make_deep_failure <- function(model, arch, start_id, seed, elapsed, err) {
  data.frame(
    ArchID = as.character(arch$ArchID),
    K1 = arch$K1, K2 = arch$K2, r1 = arch$r1, r2 = arch$r2,
    Start = start_id, Seed = seed, Model = model,
    ARI = NA_real_, MR = NA_real_, N_clusters = NA_integer_,
    BIC = NA_real_, ICL_BIC = NA_real_, LogLik = NA_real_, MeanNu = NA_real_,
    MinPsi = NA_real_, AtPsiFloor = NA,
    Iterations = NA_integer_, Converged = FALSE,
    Elapsed_seconds = elapsed, Success = FALSE,
    Error = conditionMessage(err),
    stringsAsFactors = FALSE
  )
}

CHECKPOINT_FILE <- file.path(RESULT_DIR, "deep_all_starts.csv")

read_checkpoint <- function() {
  if (!RESUME_FROM_CHECKPOINT || !file.exists(CHECKPOINT_FILE)) return(NULL)
  d <- tryCatch(utils::read.csv(CHECKPOINT_FILE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d
}

fit_key_done <- function(chk, model, arch_id, start_id) {
  if (is.null(chk) || !nrow(chk)) return(FALSE)
  hit <- chk$Model == model & chk$ArchID == arch_id & chk$Start == start_id
  if (!any(hit)) return(FALSE)
  if (RETRY_FAILED_FITS) return(any(hit & chk$Success))
  TRUE
}

append_or_replace_checkpoint <- function(chk, row) {
  if (is.null(chk) || !nrow(chk)) chk <- row[0, , drop = FALSE]
  if (nrow(chk)) {
    same <- chk$Model == row$Model[1L] &
      chk$ArchID == row$ArchID[1L] &
      chk$Start == row$Start[1L]
    chk <- chk[!same, , drop = FALSE]
  }
  chk <- rbind(chk, row)
  ord <- order(match(chk$ArchID, ARCH_GRID$ArchID), chk$Start,
               match(chk$Model, c("DGMM", "RDMM", "DStMM")))
  chk <- chk[ord, , drop = FALSE]
  rownames(chk) <- NULL
  utils::write.csv(chk, CHECKPOINT_FILE, row.names = FALSE)
  chk
}

rdmm_cache_file <- function(arch_id, start_id) {
  file.path(FIT_CACHE_DIR, paste0("RDMM_", gsub("[^A-Za-z0-9_-]", "_", arch_id),
                                 "_start", start_id, ".rds"))
}


dstmm_cache_file <- function(arch_id, start_id) {
  file.path(FIT_CACHE_DIR, paste0("DStMM_", gsub("[^A-Za-z0-9_-]", "_", arch_id),
                                 "_start", start_id, ".rds"))
}

run_deep_architecture_grid <- function() {
  chk <- read_checkpoint()
  if (!is.null(chk)) cat("\nResuming from", nrow(chk), "checkpoint row(s).\n")

  for (aa in seq_len(nrow(ARCH_GRID))) {
    arch <- ARCH_GRID[aa, , drop = FALSE]
    k <- c(arch$K1, arch$K2)
    r <- c(arch$r1, arch$r2)

    cat("\n------------------------------------------------------------\n")
    cat("Architecture:", arch$ArchID, "\n")
    cat("k =", paste(k, collapse = " x "), "; r =", paste(r, collapse = " > "), "\n")
    cat("------------------------------------------------------------\n")

    for (ss in seq_len(N_DEEP_STARTS)) {
      seed <- as.integer(MASTER_SEED + aa * 1000L + ss)
      cat("Start", ss, "of", N_DEEP_STARTS, "(seed", seed, ")\n")
      flush.console()

      # DGMM -----------------------------------------------------------------
      if (!fit_key_done(chk, "DGMM", arch$ArchID, ss)) {
        cat("  DGMM ... "); flush.console()
        t0 <- proc.time()[[3L]]
        dg <- tryCatch(fit_dgmm_arch(X, k, r, seed), error = identity)
        elapsed <- proc.time()[[3L]] - t0
        if (inherits(dg, "error")) {
          cat("FAILED:", conditionMessage(dg), "\n")
          row <- make_deep_failure("DGMM", arch, ss, seed, elapsed, dg)
        } else {
          cat(sprintf("done (%.1fs)\n", elapsed))
          row <- make_deep_result(dg, "DGMM", arch, ss, seed, elapsed)
        }
        chk <- append_or_replace_checkpoint(chk, row)
      } else {
        cat("  DGMM: checkpointed\n")
      }

      # RDMM -----------------------------------------------------------------
      rd_file <- rdmm_cache_file(arch$ArchID, ss)
      rd_warm <- NULL
      if (!fit_key_done(chk, "RDMM", arch$ArchID, ss)) {
        cat("  RDMM ... "); flush.console()
        t0 <- proc.time()[[3L]]
        rd <- tryCatch(fit_rdmm_arch(X, k, r, seed), error = identity)
        elapsed <- proc.time()[[3L]] - t0
        if (inherits(rd, "error")) {
          cat("FAILED:", conditionMessage(rd), "\n")
          row <- make_deep_failure("RDMM", arch, ss, seed, elapsed, rd)
        } else {
          cat(sprintf("done (%.1fs)\n", elapsed))
          rd_warm <- rd
          saveRDS(rd, rd_file)
          row <- make_deep_result(rd, "RDMM", arch, ss, seed, elapsed)
        }
        chk <- append_or_replace_checkpoint(chk, row)
      } else {
        cat("  RDMM: checkpointed\n")
        if (file.exists(rd_file)) {
          rd_warm <- tryCatch(readRDS(rd_file), error = function(e) NULL)
        }
      }

      # If DStMM is unfinished but the checkpoint lacks the RDMM object, refit
      # RDMM only to reconstruct the required warm start. Do not add a duplicate
      # result row unless the old row was absent.
      if (!fit_key_done(chk, "DStMM", arch$ArchID, ss) && is.null(rd_warm)) {
        cat("  Reconstructing RDMM warm start ... "); flush.console()
        rd_tmp <- tryCatch(fit_rdmm_arch(X, k, r, seed), error = identity)
        if (!inherits(rd_tmp, "error")) {
          rd_warm <- rd_tmp
          saveRDS(rd_tmp, rd_file)
          cat("done\n")
        } else {
          cat("FAILED\n")
        }
      }

      # DStMM ----------------------------------------------------------------
      if (!fit_key_done(chk, "DStMM", arch$ArchID, ss)) {
        cat("  DStMM ... "); flush.console()
        t0 <- proc.time()[[3L]]
        ds <- tryCatch(fit_dstmm_arch(X, k, r, seed, rd_warm), error = identity)
        elapsed <- proc.time()[[3L]] - t0
        if (inherits(ds, "error")) {
          cat("FAILED:", conditionMessage(ds), "\n")
          row <- make_deep_failure("DStMM", arch, ss, seed, elapsed, ds)
        } else {
          cat(sprintf("done (%.1fs)\n", elapsed))
          saveRDS(ds, dstmm_cache_file(arch$ArchID, ss))
          save_dstmm_delta_diagnostics(ds, arch, ss)
          row <- make_deep_result(ds, "DStMM", arch, ss, seed, elapsed)
        }
        chk <- append_or_replace_checkpoint(chk, row)
      } else {
        cat("  DStMM: checkpointed\n")
      }
    }
  }

  chk
}

select_best_start_by_bic <- function(deep_raw) {
  ok <- deep_raw[deep_raw$Success & is.finite(deep_raw$BIC), , drop = FALSE]
  if (!nrow(ok)) stop("No successful deep-model fits were obtained.")
  key <- interaction(ok$ArchID, ok$Model, drop = TRUE, lex.order = TRUE)
  pieces <- split(ok, key)
  best <- lapply(pieces, function(dd) dd[which.min(dd$BIC), , drop = FALSE])
  ans <- do.call(rbind, best)
  rownames(ans) <- NULL
  ans
}

make_table2 <- function(best) {
  arch_cols <- c("ArchID", "K1", "K2", "r1", "r2")
  arch <- unique(best[, arch_cols, drop = FALSE])
  arch <- arch[order(arch$K2, arch$r1, arch$r2), , drop = FALSE]
  out <- arch
  for (mm in c("DGMM", "RDMM", "DStMM")) {
    dd <- best[best$Model == mm,
               c("ArchID", "ARI", "MR", "BIC", "LogLik", "MeanNu",
                 "MinPsi", "AtPsiFloor", "Start"), drop = FALSE]
    names(dd)[-1L] <- paste0(mm, "_", names(dd)[-1L])
    out <- merge(out, dd, by = "ArchID", all.x = TRUE, sort = FALSE)
  }
  out <- out[match(arch$ArchID, out$ArchID), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ----------------------------------------------------------------------------
# 8. RUN DEEP MODELS + TABLE 2
# ----------------------------------------------------------------------------

deep_raw <- run_deep_architecture_grid()
deep_best <- select_best_start_by_bic(deep_raw)
utils::write.csv(deep_best,
                 file.path(RESULT_DIR, "deep_best_start_by_architecture.csv"),
                 row.names = FALSE)

table2 <- make_table2(deep_best)
utils::write.csv(table2,
                 file.path(RESULT_DIR, "Table2_Deep_Architecture_Comparison.csv"),
                 row.names = FALSE)

best_deep_model <- do.call(rbind, lapply(c("DGMM", "RDMM", "DStMM"), function(mm) {
  dd <- deep_best[deep_best$Model == mm & is.finite(deep_best$BIC), , drop = FALSE]
  if (!nrow(dd)) return(NULL)
  dd[which.min(dd$BIC), , drop = FALSE]
}))
if (!is.null(best_deep_model)) {
  rownames(best_deep_model) <- NULL
  utils::write.csv(best_deep_model,
                   file.path(RESULT_DIR, "deep_best_architecture_by_model.csv"),
                   row.names = FALSE)
}

# ----------------------------------------------------------------------------
# 9. BASELINES FOR TABLE 1
# ----------------------------------------------------------------------------

baseline_row <- function(method, pred, input, family,
                         bic = NA_real_, loglik = NA_real_, detail = "",
                         elapsed = NA_real_, status = "OK") {
  ss <- score_partition(pred)
  data.frame(
    Method = method, Family = family, Input = input,
    ARI = ss$ARI, MR = ss$MR,
    BIC = bic, LogLik = loglik, Elapsed_seconds = elapsed,
    Status = status, Detail = detail,
    stringsAsFactors = FALSE
  )
}

baseline_failure <- function(method, input, family, msg, elapsed = NA_real_) {
  data.frame(
    Method = method, Family = family, Input = input,
    ARI = NA_real_, MR = NA_real_, BIC = NA_real_, LogLik = NA_real_,
    Elapsed_seconds = elapsed, Status = "FAILED/SKIPPED",
    Detail = as.character(msg), stringsAsFactors = FALSE
  )
}

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

run_baselines <- function() {
  rows <- list()
  details <- list()
  ii <- 0L
  add_row <- function(z) { ii <<- ii + 1L; rows[[ii]] <<- z }

  # Common PCA scores for full-covariance baseline families.
  pca_dim <- min(BASELINE_PCA_DIM, nrow(X) - 2L, ncol(X))
  pc <- stats::prcomp(X, center = FALSE, scale. = FALSE, rank. = pca_dim)
  X_pca <- pc$x[, seq_len(min(pca_dim, ncol(pc$x))), drop = FALSE]
  pca_label <- paste0("PCA(", ncol(X_pca), ")")

  cat("\n================ BASELINES ================\n")

  # K-means ------------------------------------------------------------------
  cat("[1] K-means ... "); flush.console()
  t0 <- proc.time()[[3L]]
  km <- tryCatch({
    set.seed(MASTER_SEED)
    stats::kmeans(X, centers = K1, nstart = max(25L, BASELINE_NSTART),
                  iter.max = 200L)$cluster
  }, error = identity)
  el <- proc.time()[[3L]] - t0
  if (inherits(km, "error")) add_row(baseline_failure("K-means", "Standardized mfeat-fac features", "Distance", conditionMessage(km), el))
  else add_row(baseline_row("K-means", km, "Standardized mfeat-fac features", "Distance", elapsed = el))
  cat("done\n")

  # Ward ---------------------------------------------------------------------
  cat("[2] Ward.D2 ... "); flush.console()
  t0 <- proc.time()[[3L]]
  wd <- tryCatch(stats::cutree(stats::hclust(stats::dist(X), method = "ward.D2"), k = K1),
                 error = identity)
  el <- proc.time()[[3L]] - t0
  if (inherits(wd, "error")) add_row(baseline_failure("Ward.D2", "Standardized mfeat-fac features", "Hierarchical", conditionMessage(wd), el))
  else add_row(baseline_row("Ward.D2", wd, "Standardized mfeat-fac features", "Hierarchical", elapsed = el))
  cat("done\n")

  # Mclust -------------------------------------------------------------------
  cat("[3] Mclust on ", pca_label, " ... ", sep = ""); flush.console()
  if (has_pkg("mclust")) {
    t0 <- proc.time()[[3L]]
    mc <- tryCatch(mclust::Mclust(X_pca, G = K1, verbose = FALSE), error = identity)
    el <- proc.time()[[3L]] - t0
    if (inherits(mc, "error")) add_row(baseline_failure("Mclust", pca_label, "Gaussian mixture", conditionMessage(mc), el))
    else add_row(baseline_row("Mclust", mc$classification, pca_label, "Gaussian mixture",
                              bic = safe_num(mc$bic), loglik = safe_num(mc$loglik),
                              detail = paste0("model=", mc$modelName), elapsed = el))
  } else add_row(baseline_failure("Mclust", pca_label, "Gaussian mixture", "Package mclust not installed"))
  cat("done\n")

  # tEIGEN -------------------------------------------------------------------
  cat("[4] tEIGEN on ", pca_label, " ... ", sep = ""); flush.console()
  if (has_pkg("teigen")) {
    t0 <- proc.time()[[3L]]
    tg <- tryCatch({
      set.seed(MASTER_SEED)
      teigen::teigen(
        x = X_pca, Gs = K1, models = "all",
        init = "emem", scale = FALSE, dfstart = 30,
        dfupdate = "approx", eps = c(0.001, 0.01),
        verbose = FALSE, maxit = c(75L, BASELINE_MAX_ITER),
        convstyle = "aitkens", parallel.cores = FALSE,
        ememargs = list(BASELINE_NSTART, 5L, "UUUU", "hard")
      )
    }, error = identity)
    el <- proc.time()[[3L]] - t0
    if (inherits(tg, "error")) add_row(baseline_failure("tEIGEN", pca_label, "t mixture", conditionMessage(tg), el))
    else add_row(baseline_row("tEIGEN", tg$classification, pca_label, "t mixture",
                              bic = safe_num(tg$bic), loglik = safe_num(tg$logl),
                              detail = paste0("model=", tg$modelname), elapsed = el))
  } else add_row(baseline_failure("tEIGEN", pca_label, "t mixture", "Package teigen not installed"))
  cat("done\n")

  # Factor-analytic baselines ------------------------------------------------
  if (RUN_FACTOR_BASELINES) {
    # PGMM
    cat("[5] PGMM (slow) ... "); flush.console()
    if (has_pkg("pgmm")) {
      t0 <- proc.time()[[3L]]
      pg <- tryCatch(pgmm::pgmmEM(
        x = X, rG = K1, rq = PGMM_Q_GRID,
        class = NULL, icl = FALSE, zstart = 1L, cccStart = TRUE,
        loop = BASELINE_NSTART,
        modelSubset = c("CCC", "CCU", "CUC", "CUU", "UCC", "UCU", "UUC", "UUU"),
        seed = MASTER_SEED, tol = 0.01, relax = FALSE
      ), error = identity)
      el <- proc.time()[[3L]] - t0
      if (inherits(pg, "error")) add_row(baseline_failure("PGMM", "Standardized mfeat-fac features", "Gaussian factor mixture", conditionMessage(pg), el))
      else add_row(baseline_row("PGMM", pg$map, "Standardized mfeat-fac features", "Gaussian factor mixture",
                                bic = if (!is.null(pg$bic)) suppressWarnings(max(as.numeric(pg$bic), na.rm = TRUE)) else NA_real_,
                                detail = paste0("best model=", pg$model, ", q=", pg$q), elapsed = el))
    } else add_row(baseline_failure("PGMM", "Standardized mfeat-fac features", "Gaussian factor mixture", "Package pgmm not installed"))
    cat("done\n")

    # EMMIXmfa factor-analyzer families: MFA / MtFA / MCFA / MCtFA
    # These give four directly comparable high-dimensional factor-mixture baselines.
    factor_specs <- list(
      mfa   = list(method = "MFA",   family = "Gaussian factor mixture"),
      mtfa  = list(method = "MtFA",  family = "t factor mixture"),
      mcfa  = list(method = "MCFA",  family = "Gaussian common-factor mixture"),
      mctfa = list(method = "MCtFA", family = "t common-factor mixture")
    )

    for (kind in names(factor_specs)) {
      method_name <- factor_specs[[kind]]$method
      family_name <- factor_specs[[kind]]$family
      cat("[factor] ", method_name, " (slow) ... ", sep = ""); flush.console()

      if (!has_pkg("EMMIXmfa")) {
        add_row(baseline_failure(method_name, "Standardized mfeat-fac features",
                                 family_name, "Package EMMIXmfa not installed"))
      } else {
        fits <- list(); qrows <- list()
        t0 <- proc.time()[[3L]]

        for (q in FA_Q_GRID) {
          set.seed(MASTER_SEED + q + switch(kind,
            mfa = 0L, mtfa = 1000L, mcfa = 2000L, mctfa = 3000L
          ))

          fq <- tryCatch({
            common_args <- list(
              Y = X, g = K1, q = q,
              itmax = BASELINE_MAX_ITER,
              nkmeans = max(3L, BASELINE_NSTART),
              nrandom = max(3L, BASELINE_NSTART),
              tol = 1e-5,
              warn_messages = FALSE
            )

            if (kind == "mfa") {
              do.call(EMMIXmfa::mfa, c(common_args, list(
                sigma_type = "common", D_type = "common"
              )))
            } else if (kind == "mtfa") {
              do.call(EMMIXmfa::mtfa, c(common_args, list(
                df_init = rep(30, K1), df_update = TRUE,
                sigma_type = "common", D_type = "common"
              )))
            } else if (kind == "mcfa") {
              do.call(EMMIXmfa::mcfa, common_args)
            } else {
              do.call(EMMIXmfa::mctfa, c(common_args, list(
                df_init = rep(30, K1), df_update = TRUE
              )))
            }
          }, error = identity)

          if (!inherits(fq, "error")) {
            fits[[as.character(q)]] <- fq
            qrows[[length(qrows) + 1L]] <- data.frame(
              q = q,
              BIC = safe_num(fq$BIC),
              LogLik = safe_num(fq$logL),
              Success = TRUE,
              Error = "",
              stringsAsFactors = FALSE
            )
          } else {
            qrows[[length(qrows) + 1L]] <- data.frame(
              q = q, BIC = NA_real_, LogLik = NA_real_, Success = FALSE,
              Error = conditionMessage(fq), stringsAsFactors = FALSE
            )
          }
        }

        el <- proc.time()[[3L]] - t0
        qdiag <- do.call(rbind, qrows)
        details[[paste0(method_name, "_q_grid")]] <- qdiag
        ok <- qdiag[qdiag$Success & is.finite(qdiag$BIC), , drop = FALSE]

        if (!nrow(ok)) {
          add_row(baseline_failure(method_name, "Standardized mfeat-fac features",
                                   family_name, "No successful q fit", el))
        } else {
          # EMMIXmfa reports BIC as -2 logL + penalty; smaller is preferred.
          best_q <- ok$q[which.min(ok$BIC)]
          ff <- fits[[as.character(best_q)]]
          pred_ff <- if (!is.null(ff$clust)) ff$clust else tryCatch(
            stats::predict(ff, X), error = function(e) NULL
          )

          if (is.null(pred_ff)) {
            add_row(baseline_failure(method_name, "Standardized mfeat-fac features",
                                     family_name, "Fitted model returned no clustering", el))
          } else {
            add_row(baseline_row(
              method_name, pred_ff, "Standardized mfeat-fac features", family_name,
              bic = safe_num(ff$BIC), loglik = safe_num(ff$logL),
              detail = paste0("BIC-selected q=", best_q), elapsed = el
            ))
          }
        }
      }
      cat("done\n")
    }
  } else {
    add_row(baseline_failure("PGMM", "Standardized mfeat-fac features", "Gaussian factor mixture", "Disabled: set RUN_FACTOR_BASELINES=TRUE for final run"))
    add_row(baseline_failure("MFA", "Standardized mfeat-fac features", "Gaussian factor mixture", "Disabled: set RUN_FACTOR_BASELINES=TRUE for final run"))
    add_row(baseline_failure("MtFA", "Standardized mfeat-fac features", "t factor mixture", "Disabled: set RUN_FACTOR_BASELINES=TRUE for final run"))
    add_row(baseline_failure("MCFA", "Standardized mfeat-fac features", "Gaussian common-factor mixture", "Disabled: set RUN_FACTOR_BASELINES=TRUE for final run"))
    add_row(baseline_failure("MCtFA", "Standardized mfeat-fac features", "t common-factor mixture", "Disabled: set RUN_FACTOR_BASELINES=TRUE for final run"))
  }

  # Skew-t baselines ---------------------------------------------------------
  if (RUN_SLOW_SKEWT_BASELINES) {
    # Direct skew-t FACTOR-ANALYZER competitor -------------------------------
    # Murray, Browne & McNicholas: mixture of unrestricted skew-t factor
    # analyzers (uMSTFA). This is fitted on the original standardized 216-D
    # feature space, not on PCA scores.
    cat("[skew-t factor] uMSTFA (uskewFactors) ... "); flush.console()
    if (has_pkg("uskewFactors")) {
      t0 <- proc.time()[[3L]]
      uf_fits <- list(); uf_diag <- list()

      q_grid <- USTFA_Q_GRID[USTFA_Q_GRID >= 1L & USTFA_Q_GRID < ncol(X)]
      for (q in q_grid) {
        for (st in seq_len(USTFA_NSTART)) {
          # init=1 uses k-means; additional starts use random partitions.
          init_code <- if (st == 1L) 1L else 2L
          set.seed(MASTER_SEED + 60000L + 100L * q + st)

          uf <- tryCatch(
            uskewFactors::uskewFA(
              x = X, G = K1, q = q,
              init = init_code,
              max.it = BASELINE_MAX_ITER
            ),
            error = identity
          )

          key <- paste(q, st, sep = "_")
          if (inherits(uf, "error")) {
            uf_diag[[length(uf_diag) + 1L]] <- data.frame(
              q = q, Start = st, Init = init_code,
              BIC = NA_real_, LogLik = NA_real_, Success = FALSE,
              Error = conditionMessage(uf), stringsAsFactors = FALSE
            )
          } else {
            ll <- if (!is.null(uf$likelihood) && length(uf$likelihood)) {
              suppressWarnings(as.numeric(tail(uf$likelihood, 1L)))
            } else NA_real_
            uf_fits[[key]] <- uf
            uf_diag[[length(uf_diag) + 1L]] <- data.frame(
              q = q, Start = st, Init = init_code,
              BIC = safe_num(uf$bic), LogLik = ll,
              Success = !is.null(uf$map), Error = "",
              stringsAsFactors = FALSE
            )
          }
        }
      }

      el <- proc.time()[[3L]] - t0
      ufd <- if (length(uf_diag)) do.call(rbind, uf_diag) else data.frame()
      details[["uMSTFA_q_start_grid"]] <- ufd
      ok <- if (nrow(ufd)) ufd[ufd$Success & is.finite(ufd$BIC), , drop = FALSE] else ufd

      if (!nrow(ok)) {
        add_row(baseline_failure(
          "uMSTFA (uskewFactors)", "Standardized mfeat-fac features",
          "unrestricted skew-t factor mixture",
          "No successful uMSTFA q/start fit", el
        ))
      } else {
        # Murray/McNicholas convention: BIC = 2 log L - penalty; larger is better.
        ibest <- which.max(ok$BIC)
        best_q <- ok$q[ibest]
        best_st <- ok$Start[ibest]
        uf <- uf_fits[[paste(best_q, best_st, sep = "_")]]
        add_row(baseline_row(
          "uMSTFA (uskewFactors)", uf$map,
          "Standardized mfeat-fac features",
          "unrestricted skew-t factor mixture",
          bic = safe_num(uf$bic),
          loglik = if (!is.null(uf$likelihood)) safe_num(tail(uf$likelihood, 1L)) else NA_real_,
          detail = paste0("BIC-selected q=", best_q,
                          "; start=", best_st,
                          "; init=", ifelse(ok$Init[ibest] == 1L, "kmeans", "random")),
          elapsed = el
        ))
      }
    } else {
      add_row(baseline_failure(
        "uMSTFA (uskewFactors)", "Standardized mfeat-fac features",
        "unrestricted skew-t factor mixture",
        "Package uskewFactors not installed"
      ))
    }
    cat("done\n")
    cat("[skew-t] rMST ... "); flush.console()
    if (has_pkg("EMMIXskew")) {
      t0 <- proc.time()[[3L]]
      rs <- tryCatch(EMMIXskew::EmSkew(
        dat = X_pca, g = K1, distr = "mst", ncov = 3L,
        itmax = BASELINE_MAX_ITER, epsilon = 1e-6,
        nkmeans = BASELINE_NSTART, nrandom = BASELINE_NSTART,
        nhclust = FALSE, debug = FALSE, initloop = 10L
      ), error = identity)
      el <- proc.time()[[3L]] - t0
      if (inherits(rs, "error") || is.null(rs$clust)) {
        msg <- if (inherits(rs, "error")) conditionMessage(rs) else "EMMIXskew returned no clustering"
        add_row(baseline_failure("rMST (EMMIXskew)", pca_label, "restricted skew-t mixture", msg, el))
      } else add_row(baseline_row("rMST (EMMIXskew)", rs$clust, pca_label, "restricted skew-t mixture",
                                  bic = safe_num(rs$bic), loglik = safe_num(rs$loglik), elapsed = el))
    } else add_row(baseline_failure("rMST (EMMIXskew)", pca_label, "restricted skew-t mixture", "Package EMMIXskew not installed"))
    cat("done\n")

    cat("[skew-t] uMST ... "); flush.console()
    if (has_pkg("EMMIXuskew")) {
      t0 <- proc.time()[[3L]]
      us <- tryCatch(EMMIXuskew::fmmst(
        g = K1, dat = X_pca, itmax = BASELINE_MAX_ITER, eps = 1e-5,
        nkmeans = BASELINE_NSTART, print = FALSE
      ), error = identity)
      el <- proc.time()[[3L]] - t0
      if (inherits(us, "error") || is.null(us$clusters)) {
        msg <- if (inherits(us, "error")) conditionMessage(us) else "EMMIXuskew returned no clustering"
        add_row(baseline_failure("uMST (EMMIXuskew)", pca_label, "unrestricted skew-t mixture", msg, el))
      } else add_row(baseline_row("uMST (EMMIXuskew)", us$clusters, pca_label, "unrestricted skew-t mixture",
                                  bic = safe_num(us$bic), loglik = safe_num(us$loglik), elapsed = el))
    } else add_row(baseline_failure("uMST (EMMIXuskew)", pca_label, "unrestricted skew-t mixture", "Package EMMIXuskew not installed"))
    cat("done\n")

    cat("[skew-t] CFUST ... "); flush.console()
    if (has_pkg("EMMIXcskew")) {
      t0 <- proc.time()[[3L]]
      cf_fits <- list(); cf_diag <- list()
      for (q in CFUST_Q_GRID[CFUST_Q_GRID <= ncol(X_pca)]) {
        set.seed(MASTER_SEED + 50000L + q)
        cfq <- tryCatch(EMMIXcskew::fmcfust(
          g = K1, dat = X_pca, q = q,
          itmax = BASELINE_MAX_ITER, eps = 1e-5,
          nkmeans = BASELINE_NSTART, print = FALSE
        ), error = identity)
        if (inherits(cfq, "error")) {
          cf_diag[[length(cf_diag) + 1L]] <- data.frame(q = q, BIC = NA_real_, Success = FALSE, Error = conditionMessage(cfq))
        } else {
          cf_fits[[as.character(q)]] <- cfq
          cf_diag[[length(cf_diag) + 1L]] <- data.frame(q = q, BIC = safe_num(cfq$bic), Success = !is.null(cfq$clusters), Error = "")
        }
      }
      el <- proc.time()[[3L]] - t0
      cfd <- do.call(rbind, cf_diag)
      details[["CFUST_q_grid"]] <- cfd
      ok <- cfd[cfd$Success & is.finite(cfd$BIC), , drop = FALSE]
      if (!nrow(ok)) add_row(baseline_failure("CFUST (EMMIXcskew)", pca_label, "CFUST skew-t mixture", "No successful q fit", el))
      else {
        best_q <- ok$q[which.min(ok$BIC)]
        cf <- cf_fits[[as.character(best_q)]]
        add_row(baseline_row("CFUST (EMMIXcskew)", cf$clusters, pca_label, "CFUST skew-t mixture",
                             bic = safe_num(cf$bic), loglik = safe_num(cf$loglik),
                             detail = paste0("BIC-selected q=", best_q), elapsed = el))
      }
    } else add_row(baseline_failure("CFUST (EMMIXcskew)", pca_label, "CFUST skew-t mixture", "Package EMMIXcskew not installed"))
    cat("done\n")
  } else {
    add_row(baseline_failure("uMSTFA (uskewFactors)", "Standardized mfeat-fac features", "unrestricted skew-t factor mixture", "Disabled: set RUN_SLOW_SKEWT_BASELINES=TRUE for final run"))
    add_row(baseline_failure("rMST (EMMIXskew)", pca_label, "restricted skew-t mixture", "Disabled: set RUN_SLOW_SKEWT_BASELINES=TRUE for final run"))
    add_row(baseline_failure("uMST (EMMIXuskew)", pca_label, "unrestricted skew-t mixture", "Disabled: set RUN_SLOW_SKEWT_BASELINES=TRUE for final run"))
    add_row(baseline_failure("CFUST (EMMIXcskew)", pca_label, "CFUST skew-t mixture", "Disabled: set RUN_SLOW_SKEWT_BASELINES=TRUE for final run"))
  }

  list(results = do.call(rbind, rows), details = details)
}

baseline_results <- NULL
if (RUN_BASELINES) {
  bo <- run_baselines()
  baseline_results <- bo$results
  utils::write.csv(baseline_results, file.path(RESULT_DIR, "baseline_diagnostics.csv"), row.names = FALSE)
  if (length(bo$details)) {
    for (nm in names(bo$details)) {
      utils::write.csv(bo$details[[nm]], file.path(RESULT_DIR, paste0(nm, ".csv")), row.names = FALSE)
    }
  }
}

# ----------------------------------------------------------------------------
# 10. TABLE 1: METHOD COMPARISON
# ----------------------------------------------------------------------------

deep_table1 <- NULL
if (!is.null(best_deep_model) && nrow(best_deep_model)) {
  deep_table1 <- do.call(rbind, lapply(c("DGMM", "RDMM", "DStMM"), function(mm) {
    dd <- best_deep_model[best_deep_model$Model == mm, , drop = FALSE]
    if (!nrow(dd)) return(NULL)
    data.frame(
      Method = if (mm == "DStMM") "DStMM (proposed)" else mm,
      Family = switch(mm,
                      DGMM = "deep Gaussian mixture",
                      RDMM = "deep t mixture",
                      DStMM = "deep skew-t mixture"),
      Input = "Standardized mfeat-fac features",
      ARI = dd$ARI, MR = dd$MR,
      BIC = dd$BIC, LogLik = dd$LogLik,
      Elapsed_seconds = dd$Elapsed_seconds,
      Status = "OK",
      Detail = paste0("BIC-selected architecture ", dd$ArchID,
                      "; best start=", dd$Start,
                      "; minPsi=", signif(dd$MinPsi, 4),
                      "; atFloor=", dd$AtPsiFloor),
      stringsAsFactors = FALSE
    )
  }))
}

if (is.null(baseline_results)) {
  table1 <- deep_table1
} else if (is.null(deep_table1)) {
  table1 <- baseline_results
} else {
  table1 <- rbind(baseline_results, deep_table1)
}

if (!is.null(table1) && nrow(table1)) {
  # Paper-facing order: proposed model first, then the closest factor/skew-t
  # competitors, followed by broader mixture and classical clustering methods.
  method_order <- c(
    "DStMM (proposed)",
    "uMSTFA (uskewFactors)",
    "RDMM", "DGMM",
    "MtFA", "MCtFA", "MFA", "MCFA", "PGMM",
    "CFUST (EMMIXcskew)", "rMST (EMMIXskew)", "uMST (EMMIXuskew)",
    "tEIGEN", "Mclust", "K-means", "Ward.D2"
  )
  rank_method <- match(table1$Method, method_order)
  rank_method[is.na(rank_method)] <- length(method_order) + seq_len(sum(is.na(rank_method)))
  table1 <- table1[order(rank_method), , drop = FALSE]
  rownames(table1) <- NULL

  utils::write.csv(table1,
                   file.path(RESULT_DIR, "Table1_Clustering_Methods_ARI_MR.csv"),
                   row.names = FALSE)

  # A compact status file makes it obvious whether the requested comparison
  # set was actually fitted rather than silently skipped because of packages.
  comparison_status <- table1[, c("Method", "Family", "Status", "Detail"), drop = FALSE]
  utils::write.csv(comparison_status,
                   file.path(RESULT_DIR, "comparison_method_status.csv"),
                   row.names = FALSE)

  n_comp_ok <- sum(table1$Method != "DStMM (proposed)" & table1$Status == "OK")
  cat("\nSuccessful comparison methods (excluding DStMM): ", n_comp_ok, "\n", sep = "")
  if (RUN_BASELINES && RUN_MODE == "paper" && n_comp_ok < 8L) {
    warning(
      "Fewer than 8 comparison methods succeeded. Check comparison_method_status.csv ",
      "and install/repair the missing optional packages before using Table 1 in the paper."
    )
  }
}

# ----------------------------------------------------------------------------
# 11. FINAL CONSOLE SUMMARY
# ----------------------------------------------------------------------------

cat("\n\n============================================================\n")
cat("CASE 1 (DIGITS 0-9) FINISHED\n")
cat("n =", nrow(X), "; p =", ncol(X), "; K1 =", K1, "\n")
cat("Results folder:", RESULT_DIR, "\n")
cat("\nBest deep architecture by model (BIC; lower is better):\n")
if (!is.null(best_deep_model)) {
  print(best_deep_model[, c("Model", "ArchID", "ARI", "MR", "BIC", "LogLik",
                            "MeanNu", "MinPsi", "AtPsiFloor", "Start")],
        row.names = FALSE)
}
cat("\nDistribution diagnostic summary:\n")
print(diag_summary, row.names = FALSE)
cat("============================================================\n")

cat("\nMain outputs:\n")
cat("  Table1_Clustering_Methods_ARI_MR.csv\n")
cat("  Table2_Deep_Architecture_Comparison.csv\n")
cat("  deep_all_starts.csv\n")
cat("  distribution_diagnostics_summary.csv\n")
cat("  data_provenance.txt\n")
cat("  comparison_method_status.csv\n")
cat("\nThis script is configured for DGMM / RDMM / DStMM only, with K2 fixed at 1.\n")
cat("DStMM estimates a separate first-layer delta vector for each of the K1 clusters.\n")
cat("Per-start delta estimates are written as DStMM_delta_<ArchID>_start<id>.csv.\n")
