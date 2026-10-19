test_that("website verification rejects missing results and broken image paths", {
  root <- normalizePath(test_path("..", ".."), winslash = "/")
  script <- file.path(root, ".github", "pkgdown-build.R")
  skip_if_not(file.exists(script), "Website tooling is not included in installed packages")
  expressions <- parse(script)
  definition <- Filter(function(x) {
    is.call(x) && identical(x[[1L]], as.name("<-")) &&
      identical(x[[2L]], as.name("check_pages"))
  }, as.list(expressions))
  scope <- new.env(parent = globalenv())
  eval(definition[[1L]], scope)
  withr::local_dir(root)
  directory <- tempfile("pkgdown-")
  dir.create(directory)
  withr::defer(unlink(directory, recursive = TRUE))
  dir.create(file.path(directory, "reference"))
  page <- file.path(directory, "reference", "RunMonocle3.html")
  writeLines("<main><p>Code only</p></main>", page)
  expect_error(scope$check_pages(page), "Expected at least 3 result images; found 0")
  writeLines('<main><img src="missing.png"><img src="missing.png"><img src="missing.png"></main>', page)
  expect_error(scope$check_pages(page), "Missing image files")
  grDevices::png(file.path(dirname(page), "result.png"))
  graphics::plot(1:3)
  grDevices::dev.off()
  writeLines('<main><img src="result.png"><img src="result.png"><img src="result.png"></main>', page)
  expect_silent(scope$check_pages(page))
  writeLines('<main><img src="result.png"><img src="result.png"><img src="result.png"><pre class="r-err">backend failed</pre></main>', page)
  expect_error(scope$check_pages(page), "backend failed")
})

test_that("website dependency discovery includes qualified optional checks", {
  root <- normalizePath(test_path("..", ".."), winslash = "/")
  script <- file.path(root, ".github", "pkgdown-build.R")
  skip_if_not(file.exists(script), "Website tooling is not included in installed packages")
  withr::local_dir(root)
  withr::local_envvar(c(
    PKGDOWN_TOPICS = "RunSpatialBenchmark", PKGDOWN_ARTICLES = "",
    GITHUB_OUTPUT = "", PKGDOWN_EXECUTE = "false"
  ))
  output <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script), "dependencies"),
    stdout = TRUE, stderr = TRUE
  )
  expect_null(attr(output, "status"))
  packages <- strsplit(sub("^extra=", "", grep("^extra=", output, value = TRUE)), ",", fixed = TRUE)[[1L]]
  expect_true(all(c("Banksy", "BayesSpace", "github::lmweber/smoothclust") %in% packages))
})
