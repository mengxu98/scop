pkgload::load_all(".", quiet = TRUE, compile = FALSE)
res <- testthat::test_dir("tests/testthat", filter = "spatial|standard|seurat-fast",
  reporter = "summary", stop_on_failure = FALSE)
tab <- as.data.frame(res)
out <- commandArgs(trailingOnly = TRUE)
if (length(out)) utils::write.csv(tab[, !vapply(tab, is.list, logical(1)), drop = FALSE], out[1], row.names = FALSE)
if (any(tab$failed > 0 | tab$error)) quit(status = 1)
