script <- normalizePath(".github/pkgdown-dependencies.R")
run <- function() {
  original <- getwd()
  directory <- tempfile()
  dir.create(directory)
  on.exit({
    setwd(original)
    unlink(directory, recursive = TRUE)
  })
  setwd(directory)
  for (path in c("R", "man", "vignettes")) dir.create(path)
  writeLines(c(
    'RunToy <- function(x) { check_r("BackendPkg"); x[, 1] }',
    'Unrelated <- function() check_r("UnrelatedPkg")'
  ), "R/toy.R")
  writeLines(c(
    "\\name{toy}", "\\title{Toy}", "\\description{Toy}",
    "\\examples{", "RunToy(matrix(1))", "DirectPkg::fun()",
    "\\dontrun{IgnoredPkg::fun()}",
    "if (interactive()) InteractivePkg::fun()", "}"
  ), "man/toy.Rd")
  writeLines(c(
    "```{r setup}", "ArticlePkg::fun()", "```",
    "```{r external, eval=FALSE}", "ExternalPkg::fun()", "```"
  ), "vignettes/toy.Rmd")
  output <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
    env = c("PKGDOWN_TOPICS=toy", "PKGDOWN_ARTICLES=toy", "PKGDOWN_EXECUTE=true"), stdout = TRUE
  )
  stopifnot(is.null(attr(output, "status")))
  extra <- strsplit(sub("^extra=", "", output[startsWith(output, "extra=")]), ",", fixed = TRUE)[[1L]]
  stopifnot(all(c("BackendPkg", "DirectPkg", "ArticlePkg") %in% extra))
  stopifnot(!any(c("UnrelatedPkg", "IgnoredPkg", "InteractivePkg", "ExternalPkg") %in% extra))
}
run()
cat("Website dependency selection checks passed.\n")
