scope <- new.env(parent = globalenv())
for (expression in parse(".github/pkgdown-build.R")) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
    identical(expression[[2L]], as.name("check_pages"))) {
    eval(expression, envir = scope)
  }
}
run <- function() {
  directory <- tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  valid <- file.path(directory, "valid.html")
  failed <- file.path(directory, "failed.html")
  writeLines('<html><span class="r-out">The error rate is 0</span></html>', valid)
  writeLines('<html><span class="r-err co">Error: backend execution failed</span></html>', failed)
  scope$check_pages(character())
  scope$check_pages(valid)
  error <- tryCatch(scope$check_pages(failed), error = identity)
  stopifnot(inherits(error, "error"), grepl("backend execution failed", conditionMessage(error), fixed = TRUE))
  error <- tryCatch(scope$check_pages(file.path(directory, "missing.html")), error = identity)
  stopifnot(inherits(error, "error"), grepl("Missing generated pages", conditionMessage(error), fixed = TRUE))
}
run()
cat("Website execution output checks passed.\n")
