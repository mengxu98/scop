test_that("tAge preprocessing names are normalized", {
  expect_identical(normalize_tage_model_preprocessing("scaled_diff"), "scaled_diff")
  expect_identical(
    normalize_tage_model_preprocessing(c("scaled_diff", "scaled_diff")),
    "scaled_diff"
  )
  expect_error(
    normalize_tage_model_preprocessing("scaled"),
    "model_preprocessing"
  )
})

test_that("tAge model filenames infer preprocessing names", {
  expect_identical(infer_tage_model_name("EN_Chronoage_scaleddiff.pkl"), "scaled_diff")
  expect_identical(infer_tage_model_name("EN_Chronoage_yugenediff.pkl"), "yugene_diff")
  expect_identical(infer_tage_model_name("EN_Chronoage_yugene.pkl"), "yugene")
})
