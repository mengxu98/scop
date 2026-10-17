script <- normalizePath(".github/pkgdown-environment.R")
run <- function() {
  original <- getwd()
  directory <- tempfile()
  dir.create(directory)
  variables <- Sys.getenv(c("GITHUB_ENV", "GITHUB_OUTPUT"), unset = NA_character_)
  on.exit({
    setwd(original)
    unlink(directory, recursive = TRUE)
    for (name in names(variables)) {
      if (is.na(variables[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, stats::setNames(list(variables[[name]]), name))
    }
  })
  setwd(directory)
  dir.create("man")
  dir.create("vignettes")
  file.create("man/One.Rd", "man/Two.Rd", "vignettes/first.Rmd")
  Sys.setenv(GITHUB_ENV = file.path(directory, "environment"), GITHUB_OUTPUT = file.path(directory, "output"))
  for (scenario in c("unchanged", "R", "backend")) {
    plan <- list(
      full = FALSE, topics = I("One"), articles = character(), deleted = I("reference/Removed.html"),
      reasons = character(), published = "fixture", source = "current"
    )
    jsonlite::write_json(plan, "pkgdown-plan.json", auto_unbox = TRUE)
    scope <- new.env(parent = globalenv())
    scope$installed.packages <- function(...) {
      matrix(
        c("jsonlite", "2.0.0", "source-revision"),
        nrow = 1L,
        dimnames = list("jsonlite", c("Package", "Version", "RemoteSha"))
      )
    }
    previous <- list(
      R = if (scenario == "R") "0.0.0" else as.character(getRversion()),
      packages = list(jsonlite = if (scenario == "backend") "2.0.0:old-revision" else "2.0.0:source-revision")
    )
    scope$system2 <- function(...) jsonlite::toJSON(previous, auto_unbox = TRUE)
    sys.source(script, envir = scope)
    result <- jsonlite::read_json("pkgdown-plan.json")
    stopifnot(identical(result$full, scenario != "unchanged"))
    if (scenario != "unchanged") {
      stopifnot(identical(unlist(result$topics), c("One", "Two")))
      stopifnot(identical(result$articles, list("first")))
      stopifnot(identical(result$deleted, list("reference/Removed.html")))
    }
  }
}
run()
cat("Website environment invalidation checks passed.\n")
