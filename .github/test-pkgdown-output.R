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
driver <- normalizePath(".github/pkgdown-build.R")
check_build <- function() {
  directory <- tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  package <- file.path(directory, "fixture")
  library <- file.path(directory, "library")
  dir.create(package)
  dir.create(library)
  for (path in c("R", "man")) dir.create(file.path(package, path))
  writeLines(c(
    "Package: pkgdownfixture", "Version: 0.0.1", "Title: Website Build Fixture",
    "Description: Exercise executable examples and website build selection.",
    'Authors@R: person("Example", "Maintainer", email = "maintainer@example.org", role = c("aut", "cre"))',
    "License: CC0", "URL: https://example.org/"
  ), file.path(package, "DESCRIPTION"))
  writeLines("", file.path(package, "NAMESPACE"))
  writeLines(c("url: https://example.org/", "template:", "  bootstrap: 5"), file.path(package, "_pkgdown.yml"))
  for (topic in c("alpha", "beta")) {
    writeLines(c(
      paste0("\\name{", topic, "}"), paste0("\\alias{", topic, "}"),
      "\\title{Executable example}", "\\description{A real website build fixture.}",
      "\\examples{",
      sprintf('write("%s", file = Sys.getenv("PKGDOWN_TEST_COUNTER"), append = TRUE)', topic),
      "}"
    ), file.path(package, "man", paste0(topic, ".Rd")))
  }
  log <- file.path(directory, "build.log")
  status <- system2(file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", paste0("--library=", shQuote(library)), shQuote(package)),
    stdout = log, stderr = log
  )
  if (status != 0L) stop(paste(readLines(log), collapse = "\n"))
  counter <- file.path(directory, "executed.txt")
  plan <- list(
    full = TRUE, topics = c("alpha", "beta"), articles = character(),
    home = TRUE, news = TRUE, reference_index = TRUE, articles_index = TRUE,
    deleted = character(), source = "fixture"
  )
  render <- function() {
    jsonlite::write_json(plan, file.path(package, "pkgdown-plan.json"), auto_unbox = TRUE)
    callr::r(
      function(package, driver, counter) {
        setwd(package)
        Sys.setenv(PKGDOWN_TEST_COUNTER = counter)
        source(driver)
      },
      args = list(package, driver, counter), libpath = c(library, .libPaths()),
      stdout = log, stderr = log, user_profile = FALSE, system_profile = FALSE
    )
  }
  render()
  stopifnot(identical(readLines(counter), c("alpha", "beta")))
  stopifnot(file.exists(file.path(package, "docs/sitemap.xml")))
  stopifnot(file.exists(file.path(package, "docs/llms.txt")))
  stopifnot(file.exists(file.path(package, "docs/.pkgdown-build.json")))
  timings <- read.csv(file.path(package, "pkgdown-timings.csv"))
  stopifnot(all(c("reference/alpha", "reference/beta", "site-metadata", "full-site") %in% timings$page))
  plan$full <- plan$home <- plan$news <- plan$reference_index <- plan$articles_index <- FALSE
  plan$topics <- "beta"
  Sys.setFileTime(file.path(package, "docs/reference/alpha.html"), as.POSIXct("2000-01-01"))
  render()
  stopifnot(identical(readLines(counter), c("alpha", "beta", "beta")))
  plan$topics <- c("alpha", "beta")
  writeLines(c(
    "\\name{alpha}", "\\alias{alpha}", "\\title{Failure}", "\\description{Failure}",
    '\\examples{stop("backend execution failed")} '
  ), file.path(package, "man/alpha.Rd"))
  error <- tryCatch(render(), error = identity)
  stopifnot(inherits(error, "error"), grepl("backend execution failed", conditionMessage(error), fixed = TRUE))
  stopifnot(identical(readLines(counter), c("alpha", "beta", "beta")))
  timings <- read.csv(file.path(package, "pkgdown-timings.csv"))
  stopifnot("reference/alpha" %in% timings$page, !"reference/beta" %in% timings$page)
}
check_build()
cat("Website execution output checks passed.\n")
