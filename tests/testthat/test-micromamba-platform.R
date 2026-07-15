test_that("micromamba_platform recognizes the Windows x86-64 architecture name", {
  testthat::local_mocked_bindings(
    is_windows = function() TRUE,
    .package = "scop"
  )

  platform <- getFromNamespace("micromamba_platform", "scop")("x86-64")
  expect_identical(platform, "win-64")
})
