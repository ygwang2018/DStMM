rm(list = ls())
root <- normalizePath(if (file.exists("R/load_all.R")) "." else "..")
files <- c("experiment1_raw.csv", "experiment1_robustness_raw.csv", "experiment2_raw.csv")
for (ff in files) {
  path <- file.path(root, "results", ff)
  if (!file.exists(path)) next
  d <- read.csv(path, stringsAsFactors = FALSE)
  cat("\n", ff, "\n", sep = "")
  print(with(d, table(Model, Success, useNA = "ifany")))
  if (any(d$Model == "DGMM" & !d$Success)) {
    bad <- unique(d[d$Model == "DGMM" & !d$Success, "Error"])
    cat("DGMM errors:\n")
    print(bad)
  }
}
