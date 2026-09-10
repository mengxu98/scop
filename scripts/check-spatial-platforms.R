args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply an output check directory")
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
result <- rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran", "--no-examples", "--no-tests"),
  error_on = "never", check_dir = args[1])
saveRDS(list(errors = result$errors, warnings = result$warnings, notes = result$notes),
  file.path(args[1], "check-summary.rds"))
if (length(result$errors)) quit(status = 1L)
