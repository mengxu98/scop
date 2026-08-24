test_that("kegg_release_version parses the classic /info/<org> Release line", {
  lines <- c(
    "hsa\tHomo sapiens (human)",
    "Release 119.0",
    "linked db\tpathway"
  )
  expect_identical(kegg_release_version(lines), "Release 119.0")
})

test_that("kegg_release_version falls back to the release-notes page", {
  testthat::skip_if_not_installed("httr")
  html <- paste0(
    "<html><head><title>KEGG Release Notes</title></head><body>",
    "<h3>KEGG Release Notes</h3>",
    "<p>Release 119.0, July 1, 2026</p>",
    "<p>Release 118.2, June 1, 2026</p>",
    "</body></html>"
  )
  testthat::local_mocked_bindings(
    GET = function(url, ...) html,
    content = identity,
    .package = "httr"
  )
  lines <- c(
    "hsa\tHomo sapiens (human)",
    "T01001\t20091 proteins, 16794 proteins with KOs"
  )
  expect_identical(kegg_release_version(lines), "Release 119.0")
})

test_that("kegg_release_version falls back to the retrieval date", {
  testthat::skip_if_not_installed("httr")
  testthat::local_mocked_bindings(
    GET = function(url, ...) stop("network down"),
    .package = "httr"
  )
  lines <- c("hsa\tHomo sapiens (human)")
  expect_identical(
    kegg_release_version(lines),
    paste0("Retrieved ", Sys.Date())
  )
})
