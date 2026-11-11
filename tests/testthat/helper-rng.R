expect_rng_restored <- function(code, has_seed, kind) {
  withr::local_preserve_seed()
  original_kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(original_kind)))
  RNGkind(kind)
  set.seed(210)
  if (!has_seed) rm(".Random.seed", envir = .GlobalEnv)
  before <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  before_kind <- RNGkind()
  force(code)
  expect_identical(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), has_seed)
  expect_true(identical(get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE), before))
  expect_identical(RNGkind(), before_kind)
}
