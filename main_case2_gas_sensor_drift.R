# ============================================================================
# Case Study 2: UCI Gas Sensor Array Drift at Different Concentrations
# Six-class high-dimensional clustering for DGMM / RDMM / DStMM
# ============================================================================
# Goal
#   * Data: UCI dataset 270, Gas Sensor Array Drift at Different Concentrations.
#   * Observations: 13,910 measurements from 16 chemical sensors exposed to
#     six gases over ten temporal batches (36 months in total).
#   * Variables: 128 real-valued sensor-response features.
#   * Gas identity is used ONLY after fitting to compute first-layer ARI/MR.
#   * Batch and concentration are NOT used to fit/select the deep models; they
#     are retained only for post-hoc interpretation of the selected DStMM paths.
#   * K1 is fixed at 6 (the six gases); K2 is selected by BIC from a grid that
#     includes K2 > 1 so this case study can test whether a genuine deeper
#     mixture layer is supported by the data.
#   * A fast pilot mode uses a fixed batch-stratified random subsample without
#     looking at gas labels or ARI/MR. Full-data mode retains all 13,910 rows.
#   * Table 1 can include the same broad baseline families as the first case study,
#     including the direct skew-t factor-analyzer competitor uMSTFA.
#
# UCI raw format
#   Each batch*.dat line begins with "gas_id;concentration" and is followed by
#   128 libsvm-style feature:value fields. This script parses all ten batches.
#
# Reproducibility
#   * The UCI archive and extracted batch files are cached locally.
#   * The exact analysis rows and standardized matrix are written to CSV.
#   * Deep-model results are checkpointed after each model/start.
#   * Fitted deep-model objects are cached so interrupted runs can resume and
#     the BIC-selected DStMM can be used for second-layer post-hoc diagnostics.
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
MASTER_SEED <- 20260826L

# UCI gas identities are fixed before fitting. These labels are never used for
# initialization, architecture selection, or parameter estimation.
GAS_LEVELS <- 1:6
GAS_NAMES <- c("Ethanol", "Ethylene", "Ammonia",
               "Acetaldehyde", "Acetone", "Toluene")
K1 <- length(GAS_LEVELS)

# Fast first run: debug + pilot. After the pipeline is verified, switch to
# RUN_MODE="paper" and optionally DATA_MODE="full" for the manuscript run.
RUN_MODE <- "paper"      # "debug" or "paper"
DATA_MODE <- "full"     # "pilot" or "full"
PILOT_N <- 4000L         # used only when DATA_MODE == "pilot"

# UCI download/cache controls. The first URL is the current UCI static archive;
# the second is the historical UCI archive path retained as a fallback.
FORCE_REFRESH_DATA <- FALSE
GAS_ZIP_URLS <- c(
  "https://archive.ics.uci.edu/static/public/270/gas+sensor+array+drift+dataset+at+different+concentrations.zip",
  "https://archive.ics.uci.edu/ml/machine-learning-databases/00270/driftdataset.zip"
)

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

# The debug grid still includes K2>1; otherwise a pilot cannot answer the main
# question of whether a deeper mixture layer is supported.
if (identical(RUN_MODE, "debug")) {
  N_DEEP_STARTS <- 1L
  MAX_ITER <- 100L
  MIN_ITER <- 40L
  K2_GRID <- 1:3
  LATENT_GRID <- data.frame(
    r1 = c(8L, 12L),
    r2 = c(2L,  3L)
  )
} else if (identical(RUN_MODE, "paper")) {
  N_DEEP_STARTS <- 3L
  MAX_ITER <- 200L
  MIN_ITER <- 60L
  K2_GRID <- 1:4
  LATENT_GRID <- data.frame(
    r1 = c(8L, 12L, 16L),
    r2 = c(2L,  3L,  4L)
  )
} else {
  stop("RUN_MODE must be 'debug' or 'paper'.")
}

if (!DATA_MODE %in% c("pilot", "full")) {
  stop("DATA_MODE must be 'pilot' or 'full'.")
}
if (identical(DATA_MODE, "pilot") && PILOT_N < 500L) {
  stop("PILOT_N is too small for a stable six-class high-dimensional case study.")
}

# Baseline controls. Slow factor/skew-t baselines are disabled in debug mode.
RUN_BASELINES <- TRUE
RUN_FACTOR_BASELINES <- (RUN_MODE == "paper")
RUN_SLOW_SKEWT_BASELINES <- (RUN_MODE == "paper")
BASELINE_PCA_DIM <- 20L
BASELINE_NSTART <- if (RUN_MODE == "debug") 2L else 3L
BASELINE_MAX_ITER <- if (RUN_MODE == "debug") 80L else 150L
FA_Q_GRID <- 2:4
PGMM_Q_GRID <- 2:4
CFUST_Q_GRID <- 1:2
USTFA_Q_GRID <- 2:4
USTFA_NSTART <- if (RUN_MODE == "debug") 1L else 2L

# Ward.D2 forms an O(n^2) distance matrix. Skip it automatically on the full
# 13,910-row data to avoid a large memory spike.
WARD_MAX_N <- 6000L

# ----------------------------------------------------------------------------
# 2. PROJECT SETUP AND PACKAGES
# ----------------------------------------------------------------------------

if (!dir.exists(PROJECT_ROOT)) stop("PROJECT_ROOT does not exist: ", PROJECT_ROOT)
setwd(PROJECT_ROOT)
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)

RESULT_DIR <- file.path(PROJECT_ROOT, "results", "case2_gas_sensor_drift_concentrations")
DATA_DIR <- file.path(PROJECT_ROOT, "data", "case2_gas_sensor_drift_concentrations")
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
# 3. DOWNLOAD / LOAD UCI GAS SENSOR DRIFT DATA
# ----------------------------------------------------------------------------

GAS_ZIP_FILE <- file.path(DATA_DIR, "gas_sensor_drift_concentrations.zip")
GAS_RAW_DIR <- file.path(DATA_DIR, "uci_raw")
dir.create(GAS_RAW_DIR, recursive = TRUE, showWarnings = FALSE)

find_batch_files <- function(root = GAS_RAW_DIR) {
  ff <- list.files(root, pattern = "^batch([1-9]|10)\\.dat$",
                   recursive = TRUE, full.names = TRUE)
  if (!length(ff)) return(character(0))
  bid <- suppressWarnings(as.integer(sub("^batch([0-9]+)\\.dat$", "\\1", basename(ff))))
  ff <- ff[order(bid)]
  ff[!duplicated(basename(ff))]
}

obtain_gas_batches <- function() {
  cached <- find_batch_files()
  if (length(cached) == 10L && !FORCE_REFRESH_DATA) {
    message("Using cached UCI Gas Sensor batch files in: ", GAS_RAW_DIR)
    return(cached)
  }

  if (FORCE_REFRESH_DATA) {
    if (file.exists(GAS_ZIP_FILE)) unlink(GAS_ZIP_FILE)
    if (dir.exists(GAS_RAW_DIR)) unlink(GAS_RAW_DIR, recursive = TRUE, force = TRUE)
    dir.create(GAS_RAW_DIR, recursive = TRUE, showWarnings = FALSE)
  }

  ok <- FALSE
  last_error <- ""
  for (u in GAS_ZIP_URLS) {
    message("Downloading UCI Gas Sensor archive from: ", u)
    dl <- tryCatch({
      suppressWarnings(utils::download.file(u, GAS_ZIP_FILE, mode = "wb", quiet = FALSE))
      TRUE
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      FALSE
    })

    if (!isTRUE(dl) || !file.exists(GAS_ZIP_FILE) ||
        file.info(GAS_ZIP_FILE)$size < 1000000) next

    uz <- tryCatch({
      utils::unzip(GAS_ZIP_FILE, exdir = GAS_RAW_DIR)
      TRUE
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      FALSE
    })

    if (isTRUE(uz) && length(find_batch_files()) == 10L) {
      ok <- TRUE
      break
    }
  }

  if (!ok) {
    stop(
      "Could not obtain the ten UCI batch*.dat files automatically.\n",
      "Last download/unzip message: ", last_error, "\n",
      "Please download UCI dataset 270 (Gas Sensor Array Drift at Different ",
      "Concentrations), extract batch1.dat,...,batch10.dat under:\n",
      GAS_RAW_DIR
    )
  }

  find_batch_files()
}

parse_gas_batch <- function(path, batch_id) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  n <- length(lines)
  if (!n) stop("Empty batch file: ", path)

  gas_id <- integer(n)
  concentration <- numeric(n)
  Xb <- matrix(NA_real_, nrow = n, ncol = 128L)

  for (i in seq_len(n)) {
    tok <- strsplit(lines[i], "[[:space:]]+")[[1L]]
    if (length(tok) != 129L) {
      stop("Unexpected field count in ", basename(path), " line ", i,
           ": expected 129 tokens, found ", length(tok), ".")
    }

    head <- strsplit(tok[1L], ";", fixed = TRUE)[[1L]]
    if (length(head) != 2L) {
      stop("Could not parse gas_id;concentration in ", basename(path),
           " line ", i, ".")
    }
    gas_id[i] <- suppressWarnings(as.integer(head[1L]))
    concentration[i] <- suppressWarnings(as.numeric(head[2L]))

    fv <- strsplit(tok[-1L], ":", fixed = TRUE)
    idx <- suppressWarnings(as.integer(vapply(fv, `[`, character(1L), 1L)))
    val <- suppressWarnings(as.numeric(vapply(fv, `[`, character(1L), 2L)))

    if (length(idx) != 128L || anyNA(idx) || any(idx < 1L | idx > 128L) ||
        anyDuplicated(idx) > 0L || any(!is.finite(val))) {
      stop("Malformed feature:value fields in ", basename(path), " line ", i, ".")
    }
    Xb[i, idx] <- val
  }

  if (any(!gas_id %in% GAS_LEVELS)) {
    stop("Unexpected gas label in ", basename(path), ".")
  }
  if (any(!is.finite(concentration)) || any(concentration <= 0)) {
    stop("Invalid concentration value in ", basename(path), ".")
  }
  if (any(!is.finite(Xb))) {
    stop("Non-finite or missing sensor feature in ", basename(path), ".")
  }

  colnames(Xb) <- paste0("F", seq_len(ncol(Xb)))
  list(
    X = Xb,
    gas_id = gas_id,
    concentration = concentration,
    batch = rep.int(as.integer(batch_id), n)
  )
}

batch_files <- obtain_gas_batches()
batch_ids <- suppressWarnings(as.integer(
  sub("^batch([0-9]+)\\.dat$", "\\1", basename(batch_files))
))
if (!identical(sort(batch_ids), 1:10)) {
  stop("Expected exactly batch1.dat,...,batch10.dat after extraction.")
}

parsed_batches <- lapply(seq_along(batch_files), function(ii) {
  message("Parsing ", basename(batch_files[ii]), " ...")
  parse_gas_batch(batch_files[ii], batch_ids[ii])
})

X_all <- do.call(rbind, lapply(parsed_batches, `[[`, "X"))
gas_all <- unlist(lapply(parsed_batches, `[[`, "gas_id"), use.names = FALSE)
concentration_all <- unlist(lapply(parsed_batches, `[[`, "concentration"), use.names = FALSE)
batch_all <- unlist(lapply(parsed_batches, `[[`, "batch"), use.names = FALSE)

if (nrow(X_all) != 13910L || ncol(X_all) != 128L) {
  stop("Expected UCI dimension 13910 x 128, found ",
       nrow(X_all), " x ", ncol(X_all), ".")
}
if (length(gas_all) != nrow(X_all) || length(concentration_all) != nrow(X_all) ||
    length(batch_all) != nrow(X_all)) {
  stop("Metadata length mismatch after parsing UCI batch files.")
}

# Fixed pilot sample chosen using BATCH ONLY. Gas identity and concentration are
# not used to select rows. This preserves temporal coverage while keeping the
# first run substantially faster than the full 13,910-row analysis.
select_batch_stratified_indices <- function(batch, target_n, seed) {
  n <- length(batch)
  if (target_n >= n) return(seq_len(n))
  groups <- split(seq_len(n), batch)
  sizes <- vapply(groups, length, integer(1L))
  raw <- target_n * sizes / sum(sizes)
  alloc <- floor(raw)

  # Largest-remainder allocation gives exactly target_n rows while depending
  # only on temporal batch sizes, never on gas identity or concentration.
  remainder <- as.integer(target_n - sum(alloc))
  if (remainder > 0L) {
    ord <- order(raw - alloc, decreasing = TRUE)
    alloc[ord[seq_len(remainder)]] <- alloc[ord[seq_len(remainder)]] + 1L
  }
  if (any(alloc > sizes) || sum(alloc) != target_n) {
    stop("Internal pilot allocation error.")
  }

  set.seed(seed)
  idx <- unlist(Map(function(ids, m) sample(ids, size = m, replace = FALSE),
                    groups, as.integer(alloc)), use.names = FALSE)
  sort(idx)
}

if (identical(DATA_MODE, "pilot")) {
  analysis_index <- select_batch_stratified_indices(
    batch_all, min(PILOT_N, nrow(X_all)), MASTER_SEED + 270L
  )
} else {
  analysis_index <- seq_len(nrow(X_all))
}

X_raw_full <- X_all[analysis_index, , drop = FALSE]
gas_id <- gas_all[analysis_index]
concentration <- concentration_all[analysis_index]
batch_id <- batch_all[analysis_index]
source_row <- analysis_index

truth_factor <- factor(gas_id, levels = GAS_LEVELS, labels = GAS_NAMES)
truth <- as.integer(truth_factor)
truth_name <- as.character(truth_factor)

# No feature ranking or supervised screening. Remove only numerically constant
# columns in the fixed analysis population, then standardize without labels.
feature_sd <- apply(X_raw_full, 2L, stats::sd)
good_feature <- is.finite(feature_sd) & feature_sd > 1e-12
if (any(!good_feature)) {
  warning("Dropping ", sum(!good_feature),
          " numerically degenerate gas-sensor columns; no feature ranking was used.")
}
X_raw <- X_raw_full[, good_feature, drop = FALSE]
retained_feature_index <- which(good_feature)

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
                 file.path(DATA_DIR, "gas_sensor_feature_map.csv"), row.names = FALSE)

analysis_meta <- data.frame(
  SourceRow = source_row,
  Batch = batch_id,
  GasID = gas_id,
  Gas = truth_name,
  Concentration = concentration,
  stringsAsFactors = FALSE
)
utils::write.csv(analysis_meta,
                 file.path(DATA_DIR, "gas_sensor_analysis_metadata.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(analysis_meta, X_raw, check.names = FALSE),
  file.path(DATA_DIR, "gas_sensor_raw_analysis_matrix.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(analysis_meta, X, check.names = FALSE),
  file.path(DATA_DIR, "gas_sensor_standardized_analysis_matrix.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("CASE 2 DATA LOADED: UCI Gas Sensor Array Drift at Different Concentrations\n")
cat("Data mode:", DATA_MODE, "\n")
cat("n =", nrow(X), "; original p =", ncol(X_raw_full),
    "; analysis p =", ncol(X), "; K1 =", K1, "\n")
cat("Batch counts (batch used only for sampling/post-hoc interpretation):\n")
print(table(batch_id))
cat("True gas counts (labels used only for external evaluation):\n")
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
  "Dataset: UCI Gas Sensor Array Drift at Different Concentrations",
  "UCI dataset ID: 270",
  "DOI: 10.24432/C5MK6M",
  paste0("Primary source archive: ", GAS_ZIP_URLS[1L]),
  "UCI full dimension: 13910 observations x 128 sensor-response features",
  "The full dataset is organized into 10 temporal batches over 36 months.",
  "Gas classes: Ethanol, Ethylene, Ammonia, Acetaldehyde, Acetone, Toluene.",
  paste0("Data mode: ", DATA_MODE),
  paste0("Analysis n: ", nrow(X)),
  paste0("Pilot target n (ignored in full mode): ", PILOT_N),
  "Pilot rows, when used, are sampled within temporal batch only; gas labels,",
  "concentration values, ARI and MR are not used to choose the analysis rows.",
  paste0("Original p: ", ncol(X_raw_full)),
  paste0("Analysis p after only degenerate-column safeguard: ", ncol(X)),
  paste0("Degenerate columns removed: ", sum(!good_feature)),
  "No supervised feature selection was performed.",
  "All retained variables were standardized without gas labels.",
  "Gas identity is used only after fitting to compute first-layer ARI/MR.",
  "Batch and concentration are retained only for post-hoc pathway interpretation.",
  paste0("First-layer cluster count K1: ", K1),
  paste0("Second-layer candidate grid K2: ", paste(K2_GRID, collapse = ",")),
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

# Relabel arbitrary first-layer cluster IDs to gas labels with a Hungarian assignment.
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
      GasID = GAS_LEVELS[map],
      Gas = GAS_NAMES[map],
      stringsAsFactors = FALSE
    ),
    confusion = table(
      PredictedGas = factor(GAS_NAMES[mapped], levels = GAS_NAMES),
      TrueGas = factor(GAS_NAMES[truth_i], levels = GAS_NAMES)
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

# ----------------------------------------------------------------------------
# 6. ARCHITECTURE GRID + QUICK SANITY CHECKS
# ----------------------------------------------------------------------------

ARCH_GRID <- make_architecture_grid()
utils::write.csv(ARCH_GRID, file.path(RESULT_DIR, "architecture_grid.csv"), row.names = FALSE)
cat("\nArchitecture grid:\n")
print(ARCH_GRID)
cat("Planned deep fits =", nrow(ARCH_GRID) * N_DEEP_STARTS * 3L, "\n")

# Quick unsupervised sanity checks. Labels are used only to score after fitting.
cat("\nQuick unsupervised sanity checks:\n")
set.seed(MASTER_SEED)
km_check <- stats::kmeans(X, centers = K1, nstart = 50L, iter.max = 200L)$cluster
km_score <- score_partition(km_check)
cat(sprintf("  K-means: ARI=%.4f, MR=%.4f\n", km_score$ARI, km_score$MR))
if (nrow(X) <= WARD_MAX_N) {
  wd_check <- stats::cutree(stats::hclust(stats::dist(X), method = "ward.D2"), k = K1)
  wd_score <- score_partition(wd_check)
  cat(sprintf("  Ward.D2: ARI=%.4f, MR=%.4f\n", wd_score$ARI, wd_score$MR))
} else {
  cat("  Ward.D2: skipped in sanity check because n > WARD_MAX_N\n")
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
    active_set = 1L, M = MC_SIZE,
    warm_start = rdmm_warm,
    psi_floor = PSI_FLOOR, nu_bounds = NU_BOUNDS,
    min_iter = MIN_ITER, moving_window = MOVING_WINDOW,
    verbose = FALSE
  )
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

deep_fit_cache_file <- function(model, arch_id, start_id) {
  file.path(
    FIT_CACHE_DIR,
    paste0(model, "_", gsub("[^A-Za-z0-9_-]", "_", arch_id),
           "_start", start_id, ".rds")
  )
}

rdmm_cache_file <- function(arch_id, start_id) {
  deep_fit_cache_file("RDMM", arch_id, start_id)
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
          saveRDS(dg, deep_fit_cache_file("DGMM", arch$ArchID, ss))
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
          saveRDS(ds, deep_fit_cache_file("DStMM", arch$ArchID, ss))
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
# 9. POST-HOC INTERPRETATION OF THE BIC-SELECTED DStMM
# ----------------------------------------------------------------------------
# IMPORTANT: batch, concentration and gas labels enter only here, after all deep
# fitting and BIC architecture selection have been completed.

cramers_v <- function(x, y) {
  tab <- table(x, y)
  if (min(dim(tab)) <= 1L || sum(tab) <= 0L) return(NA_real_)
  cs <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  sqrt(as.numeric(cs$statistic) / (sum(tab) * (min(dim(tab)) - 1)))
}

write_table_matrix <- function(tab, file) {
  d <- as.data.frame.matrix(tab)
  d <- cbind(Row = rownames(d), d, check.names = FALSE)
  rownames(d) <- NULL
  utils::write.csv(d, file, row.names = FALSE)
}

posthoc_dstmm <- NULL
if (!is.null(best_deep_model) && any(best_deep_model$Model == "DStMM")) {
  ds_row <- best_deep_model[best_deep_model$Model == "DStMM", , drop = FALSE][1L, ]
  ds_file <- deep_fit_cache_file("DStMM", ds_row$ArchID, ds_row$Start)

  if (file.exists(ds_file)) {
    ds_best <- tryCatch(readRDS(ds_file), error = function(e) NULL)

    if (!is.null(ds_best) && !is.null(ds_best$s) && nrow(ds_best$s) == nrow(X)) {
      s_mat <- as.matrix(ds_best$s)
      s1 <- as.integer(s_mat[, 1L])
      s2 <- if (ncol(s_mat) >= 2L) as.integer(s_mat[, 2L]) else rep.int(1L, nrow(X))
      rel <- hungarian_relabel(s1, truth, groups = K1)
      mapped_gas <- GAS_NAMES[rel$mapped]

      utils::write.csv(
        rel$mapping,
        file.path(RESULT_DIR, "DStMM_Layer1_cluster_to_gas_mapping.csv"),
        row.names = FALSE
      )
      write_table_matrix(
        rel$confusion,
        file.path(RESULT_DIR, "DStMM_Layer1_confusion.csv")
      )

      posthoc_dstmm <- data.frame(
        SourceRow = source_row,
        TrueGas = truth_name,
        Batch = batch_id,
        Concentration = concentration,
        Layer1 = s1,
        Layer1MappedGas = mapped_gas,
        Layer2 = s2,
        Pathway = paste0(s1, "-", s2),
        stringsAsFactors = FALSE
      )
      utils::write.csv(
        posthoc_dstmm,
        file.path(RESULT_DIR, "DStMM_posthoc_observation_assignments.csv"),
        row.names = FALSE
      )

      write_table_matrix(
        table(posthoc_dstmm$Pathway, posthoc_dstmm$Batch),
        file.path(RESULT_DIR, "DStMM_Pathway_by_Batch.csv")
      )
      write_table_matrix(
        table(posthoc_dstmm$Layer2, posthoc_dstmm$Batch),
        file.path(RESULT_DIR, "DStMM_Layer2_by_Batch.csv")
      )
      write_table_matrix(
        table(posthoc_dstmm$Layer1MappedGas, posthoc_dstmm$Layer2),
        file.path(RESULT_DIR, "DStMM_MappedGas_by_Layer2.csv")
      )
      write_table_matrix(
        table(posthoc_dstmm$Pathway, posthoc_dstmm$TrueGas),
        file.path(RESULT_DIR, "DStMM_Pathway_by_TrueGas.csv")
      )

      # Effect-size diagnostics for whether the deeper allocation tracks temporal
      # drift. Report both overall and within fitted first-layer groups.
      assoc_rows <- list(
        data.frame(
          Scope = "Overall",
          Group = "All",
          N = nrow(posthoc_dstmm),
          CramersV_L2_Batch = cramers_v(posthoc_dstmm$Layer2, posthoc_dstmm$Batch),
          CramersV_Path_Batch = cramers_v(posthoc_dstmm$Pathway, posthoc_dstmm$Batch),
          stringsAsFactors = FALSE
        )
      )
      for (g in sort(unique(posthoc_dstmm$Layer1))) {
        dd <- posthoc_dstmm[posthoc_dstmm$Layer1 == g, , drop = FALSE]
        assoc_rows[[length(assoc_rows) + 1L]] <- data.frame(
          Scope = "Within fitted Layer1",
          Group = paste0("Layer1=", g),
          N = nrow(dd),
          CramersV_L2_Batch = cramers_v(dd$Layer2, dd$Batch),
          CramersV_Path_Batch = NA_real_,
          stringsAsFactors = FALSE
        )
      }
      assoc <- do.call(rbind, assoc_rows)
      utils::write.csv(
        assoc,
        file.path(RESULT_DIR, "DStMM_posthoc_batch_association.csv"),
        row.names = FALSE
      )

      # Concentration summaries are descriptive only; concentration was not used
      # during fitting or BIC selection.
      split_key <- interaction(posthoc_dstmm$Layer1MappedGas,
                               posthoc_dstmm$Layer2, drop = TRUE)
      conc_rows <- lapply(split(posthoc_dstmm, split_key), function(dd) {
        data.frame(
          MappedGas = dd$Layer1MappedGas[1L],
          Layer2 = dd$Layer2[1L],
          N = nrow(dd),
          MeanConcentration = mean(dd$Concentration),
          MedianConcentration = stats::median(dd$Concentration),
          Q1 = as.numeric(stats::quantile(dd$Concentration, 0.25, names = FALSE)),
          Q3 = as.numeric(stats::quantile(dd$Concentration, 0.75, names = FALSE)),
          stringsAsFactors = FALSE
        )
      })
      conc_summary <- do.call(rbind, conc_rows)
      rownames(conc_summary) <- NULL
      utils::write.csv(
        conc_summary,
        file.path(RESULT_DIR, "DStMM_concentration_by_MappedGas_Layer2.csv"),
        row.names = FALSE
      )

      writeLines(c(
        paste0("Selected DStMM architecture: ", ds_row$ArchID),
        paste0("Selected K2: ", ds_row$K2),
        paste0("Selected start: ", ds_row$Start),
        paste0("Data mode: ", DATA_MODE),
        paste0("n: ", nrow(X)),
        paste0("Overall Cramer's V(Layer2, Batch): ",
               signif(cramers_v(posthoc_dstmm$Layer2, posthoc_dstmm$Batch), 5)),
        paste0("Overall Cramer's V(Pathway, Batch): ",
               signif(cramers_v(posthoc_dstmm$Pathway, posthoc_dstmm$Batch), 5)),
        if (ds_row$K2 > 1L) {
          "K2 > 1: the BIC-selected DStMM contains genuine second-layer branching."
        } else {
          "K2 = 1: this data/run does not provide evidence for deeper mixture branching."
        }
      ), file.path(RESULT_DIR, "DStMM_posthoc_summary.txt"))
    } else {
      warning("Selected DStMM fit object is unavailable or has incompatible s matrix; post-hoc diagnostics skipped.")
    }
  } else {
    warning("Selected DStMM cache file not found; post-hoc diagnostics skipped: ", ds_file)
  }
}

# ----------------------------------------------------------------------------
# 10. BASELINES FOR TABLE 1
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
  if (inherits(km, "error")) add_row(baseline_failure("K-means", "Standardized gas-sensor features", "Distance", conditionMessage(km), el))
  else add_row(baseline_row("K-means", km, "Standardized gas-sensor features", "Distance", elapsed = el))
  cat("done\n")

  # Ward ---------------------------------------------------------------------
  cat("[2] Ward.D2 ... "); flush.console()
  if (nrow(X) > WARD_MAX_N) {
    add_row(baseline_failure(
      "Ward.D2", "Standardized gas-sensor features", "Hierarchical",
      paste0("Skipped because n=", nrow(X), " exceeds WARD_MAX_N=", WARD_MAX_N)
    ))
    cat("skipped (n too large for O(n^2) distance matrix)\n")
  } else {
    t0 <- proc.time()[[3L]]
    wd <- tryCatch(stats::cutree(stats::hclust(stats::dist(X), method = "ward.D2"), k = K1),
                   error = identity)
    el <- proc.time()[[3L]] - t0
    if (inherits(wd, "error")) add_row(baseline_failure("Ward.D2", "Standardized gas-sensor features", "Hierarchical", conditionMessage(wd), el))
    else add_row(baseline_row("Ward.D2", wd, "Standardized gas-sensor features", "Hierarchical", elapsed = el))
    cat("done\n")
  }

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
      if (inherits(pg, "error")) add_row(baseline_failure("PGMM", "Standardized gas-sensor features", "Gaussian factor mixture", conditionMessage(pg), el))
      else add_row(baseline_row("PGMM", pg$map, "Standardized gas-sensor features", "Gaussian factor mixture",
                                bic = if (!is.null(pg$bic)) suppressWarnings(max(as.numeric(pg$bic), na.rm = TRUE)) else NA_real_,
                                detail = paste0("best model=", pg$model, ", q=", pg$q), elapsed = el))
    } else add_row(baseline_failure("PGMM", "Standardized gas-sensor features", "Gaussian factor mixture", "Package pgmm not installed"))
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
        add_row(baseline_failure(method_name, "Standardized gas-sensor features",
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
          add_row(baseline_failure(method_name, "Standardized gas-sensor features",
                                   family_name, "No successful q fit", el))
        } else {
          # EMMIXmfa reports BIC as -2 logL + penalty; smaller is preferred.
          best_q <- ok$q[which.min(ok$BIC)]
          ff <- fits[[as.character(best_q)]]
          pred_ff <- if (!is.null(ff$clust)) ff$clust else tryCatch(
            stats::predict(ff, X), error = function(e) NULL
          )

          if (is.null(pred_ff)) {
            add_row(baseline_failure(method_name, "Standardized gas-sensor features",
                                     family_name, "Fitted model returned no clustering", el))
          } else {
            add_row(baseline_row(
              method_name, pred_ff, "Standardized gas-sensor features", family_name,
              bic = safe_num(ff$BIC), loglik = safe_num(ff$logL),
              detail = paste0("BIC-selected q=", best_q), elapsed = el
            ))
          }
        }
      }
      cat("done\n")
    }
  } else {
    add_row(baseline_failure("PGMM", "Standardized gas-sensor features", "Gaussian factor mixture", "Disabled: set RUN_FACTOR_BASELINES=TRUE for final run"))
    add_row(baseline_failure("MFA", "Standardized gas-sensor features", "Gaussian factor mixture", "Disabled: set RUN_FACTOR_BASELINES=TRUE for final run"))
    add_row(baseline_failure("MtFA", "Standardized gas-sensor features", "t factor mixture", "Disabled: set RUN_FACTOR_BASELINES=TRUE for final run"))
    add_row(baseline_failure("MCFA", "Standardized gas-sensor features", "Gaussian common-factor mixture", "Disabled: set RUN_FACTOR_BASELINES=TRUE for final run"))
    add_row(baseline_failure("MCtFA", "Standardized gas-sensor features", "t common-factor mixture", "Disabled: set RUN_FACTOR_BASELINES=TRUE for final run"))
  }

  # Skew-t baselines ---------------------------------------------------------
  if (RUN_SLOW_SKEWT_BASELINES) {
    # Direct skew-t FACTOR-ANALYZER competitor -------------------------------
    # Murray, Browne & McNicholas: mixture of unrestricted skew-t factor
    # analyzers (uMSTFA). This is fitted on the original standardized 128-D
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
          "uMSTFA (uskewFactors)", "Standardized gas-sensor features",
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
          "Standardized gas-sensor features",
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
        "uMSTFA (uskewFactors)", "Standardized gas-sensor features",
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
    add_row(baseline_failure("uMSTFA (uskewFactors)", "Standardized gas-sensor features", "unrestricted skew-t factor mixture", "Disabled: set RUN_SLOW_SKEWT_BASELINES=TRUE for final run"))
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
# 11. TABLE 1: METHOD COMPARISON
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
      Input = "Standardized gas-sensor features",
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
  if (RUN_MODE == "paper" && n_comp_ok < 8L) {
    warning(
      "Fewer than 8 comparison methods succeeded. Check comparison_method_status.csv ",
      "and install/repair the missing optional packages before using Table 1 in the paper."
    )
  }
}

# ----------------------------------------------------------------------------
# 12. FINAL CONSOLE SUMMARY
# ----------------------------------------------------------------------------

cat("\n\n============================================================\n")
cat("CASE 2 (GAS SENSOR DRIFT) FINISHED\n")
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
cat("  DStMM_posthoc_summary.txt (when selected DStMM cache is available)\n")
cat("  DStMM_Pathway_by_Batch.csv\n")
cat("  DStMM_concentration_by_MappedGas_Layer2.csv\n")
cat("\nCurrent DATA_MODE =", DATA_MODE, "; RUN_MODE =", RUN_MODE, "\n")
cat("Recommended workflow: debug+pilot first; then RUN_MODE='paper'.\n")
cat("If runtime is acceptable and final full-data evidence is desired, also set DATA_MODE='full'.\n")
cat("Slow factor/skew-t competitors run only in paper mode by default.\n")
cat("Direct skew-t factor competitor: uMSTFA via package 'uskewFactors'.\n")
