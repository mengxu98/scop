check_pages <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing generated pages: ", paste(missing, collapse = ", "))
  errors <- lapply(paths, function(path) {
    page <- xml2::read_html(path)
    messages <- xml2::xml_text(xml2::xml_find_all(page, "//*[contains(concat(' ', normalize-space(@class), ' '), ' r-err ')]"))
    if (length(messages)) paste(path, paste(unique(messages), collapse = "\n"), sep = ":\n") else character()
  })
  errors <- unlist(errors, use.names = FALSE)
  if (length(errors)) stop("Example execution failed:\n", paste(errors, collapse = "\n\n"))
}

plan <- jsonlite::read_json("pkgdown-plan.json", simplifyVector = TRUE)
timings <- data.frame(page = character(), seconds = numeric())
record <- function(name, expression) {
  start <- proc.time()[["elapsed"]]
  on.exit({
    timings <<- rbind(timings, data.frame(page = name, seconds = proc.time()[["elapsed"]] - start))
    utils::write.csv(timings, "pkgdown-timings.csv", row.names = FALSE)
  })
  force(expression)
}
capture.output(sessioninfo::session_info(), file = "pkgdown-session.txt")
installed <- installed.packages(fields = "RemoteSha")
installed <- installed[!duplicated(installed[, "Package"]), , drop = FALSE]
utils::write.table(installed[, c("Package", "Version", "RemoteSha")],
  file = "pkgdown-session.txt", append = TRUE, sep = "\t", row.names = FALSE
)
options(Ncpus = 2L, future.globals.maxSize = 4 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential", gc = TRUE)
pkg <- pkgdown::as_pkgdown(".")
build <- function() {
  if (plan$full) record("initialise", pkgdown::init_site(pkg))
  if (plan$home) record("home", pkgdown::build_home(pkg))
  if (plan$news) record("news", pkgdown::build_news(pkg, preview = FALSE))
  if (plan$reference_index) record("reference/index", pkgdown::build_reference_index(pkg))
  for (topic in plan$topics) {
    record(paste0("reference/", topic), {
      pkgdown::build_reference(pkg, topics = topic, examples = TRUE, run_dont_run = FALSE, devel = FALSE)
      check_pages(file.path(pkg$dst_path, "reference", paste0(topic, ".html")))
    })
  }
  if (plan$articles_index) record("articles/index", pkgdown::build_articles_index(pkg))
  for (article in plan$articles) {
    record(paste0("articles/", article), {
      pkgdown::build_article(article, pkg = pkg)
      check_pages(file.path(pkg$dst_path, "articles", paste0(article, ".html")))
    })
  }
  if (plan$full) {
    record("site-metadata", pkgdown::build_site(
      pkg,
      lazy = TRUE, examples = TRUE, run_dont_run = FALSE, install = FALSE,
      new_process = FALSE, preview = FALSE
    ))
  } else {
    record("search", pkgdown::build_search(pkg))
  }
}
record(if (plan$full) "full-site" else "incremental-site", build())
expected <- c(
  character(),
  if (length(plan$topics)) file.path("docs/reference", paste0(plan$topics, ".html")),
  if (length(plan$articles)) file.path("docs/articles", paste0(plan$articles, ".html"))
)
check_pages(expected)
if (length(plan$deleted) && any(file.exists(file.path("docs", plan$deleted)))) {
  stop("Deleted documentation pages remain in website output")
}
file.create("docs/.nojekyll")
installed <- installed.packages(fields = "RemoteSha")
installed <- installed[!duplicated(installed[, "Package"]), , drop = FALSE]
versions <- stats::setNames(
  paste(installed[, "Version"], ifelse(is.na(installed[, "RemoteSha"]), "", installed[, "RemoteSha"]), sep = ":"),
  installed[, "Package"]
)
previous <- if (file.exists("docs/.pkgdown-build.json")) {
  jsonlite::read_json("docs/.pkgdown-build.json")$packages
} else {
  list()
}
previous[names(versions)] <- as.list(versions)
jsonlite::write_json(
  list(source = plan$source, R = as.character(getRversion()), packages = previous),
  "docs/.pkgdown-build.json",
  auto_unbox = TRUE, pretty = TRUE
)
summary <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(summary)) {
  write(
    c(
      "", "| Page | Seconds |", "| --- | ---: |",
      sprintf("| %s | %.2f |", timings$page, timings$seconds)
    ),
    file = summary, append = TRUE
  )
}
