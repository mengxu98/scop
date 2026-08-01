test_that("runner_system2 supports an empty environment override", {
  output <- tempfile()
  error_output <- tempfile()
  status <- getFromNamespace("runner_system2", "scop")(
    command = file.path(R.home("bin"), "Rscript"),
    args = c("-e", shQuote("cat('runner-ok')")),
    env = character(),
    stdout = output,
    stderr = error_output
  )

  expect_identical(status, 0L)
  expect_identical(readLines(output, warn = FALSE), "runner-ok")
})

test_that("runner locks are exclusive and only their owner releases them", {
  lock_path <- tempfile("runner_lock_")
  acquire <- getFromNamespace("runner_acquire_lock", "scop")
  release <- getFromNamespace("runner_release_lock", "scop")
  lock <- acquire(lock_path, backend = "test backend")

  expect_error(
    acquire(lock_path, backend = "test backend"),
    "Another.*run"
  )
  expect_false(release(list(path = lock_path, token = "not-the-owner")))
  expect_true(file.exists(lock_path))
  expect_true(release(lock))
  expect_false(file.exists(lock_path))
})

test_that("runner JSON writes leave only a complete target", {
  path <- tempfile("runner_json_", fileext = ".json")
  write_json <- getFromNamespace("runner_write_json", "scop")
  read_json <- getFromNamespace("runner_read_json", "scop")

  write_json(list(value = "first"), path)
  write_json(list(value = "second"), path)

  expect_identical(read_json(path)$value, "second")
  leftovers <- list.files(
    dirname(path),
    pattern = paste0("^\\.", basename(path), "\\."),
    full.names = TRUE
  )
  expect_length(leftovers, 0L)
})

test_that("runner JSON writes report failed atomic replacements", {
  target <- tempfile("runner-json-directory-")
  dir.create(target)

  expect_error(
    getFromNamespace("runner_write_json", "scop")(
      list(value = "new"),
      target
    ),
    "Unable to atomically write"
  )
  expect_true(dir.exists(target))
})

test_that("runner log tails stay bounded and retain the final lines", {
  output <- tempfile()
  writeLines(
    c(paste("line", seq_len(5000L)), "", "final sentinel"),
    output
  )

  tail_lines <- getFromNamespace("runner_tail_lines", "scop")(
    output,
    max_lines = 20L,
    chunk_size = 37L
  )

  expect_length(tail_lines, 20L)
  expect_identical(tail_lines[[20L]], "final sentinel")
  expect_false(any(!nzchar(tail_lines)))
})
