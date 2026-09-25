required <- c("mvtnorm", "corpcor", "GIGrvg", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
message("Required packages are installed: ", paste(required, collapse = ", "))
