test_that("startup metadata follows the loaded namespace in development", {
  withr::local_options(list(
    log_message.verbose = TRUE, scop_env_init = FALSE,
    scop_env_cache = getOption("scop_env_cache")
  ))
  # pkgload's source parent is not an installed library containing scop.
  hook <- getFromNamespace(".onAttach", "scop")
  nonexistent_library <- tempfile("not-an-installed-library-")
  text <- paste(capture_messages(hook(nonexistent_library, "scop")), collapse = "\n")
  expect_match(text, "Version:")
  expect_match(text, "update")
  expect_match(text, "Website:")
  expect_match(text, "https://mengxu98.github.io/scop/", fixed = TRUE)
  expect_length(capture_messages(suppressPackageStartupMessages(
    hook(nonexistent_library, "scop")
  )), 0L)
  withr::local_options(log_message.verbose = FALSE)
  expect_length(capture_messages(hook(nonexistent_library, "scop")), 0L)
})
