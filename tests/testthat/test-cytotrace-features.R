test_that("CytoTRACE2 model loading tolerates duplicated CSV row indexes", {
  data_dir <- tempfile()
  dir.create(data_dir)
  saveRDS(list(), file.path(data_dir, "model_parameters.rds"))
  writeLines(
    c(
      ",0",
      "0,A1bg",
      "0,A1cf",
      "1,A4galt"
    ),
    file.path(data_dir, "features_model_training_17.csv")
  )

  expect_identical(
    scop:::load_cytotrace2_data(data_dir, verbose = FALSE)$features,
    c("A1bg", "A1cf", "A4galt")
  )
})

test_that("CytoTRACE2 model loading rejects duplicated feature names", {
  data_dir <- tempfile()
  dir.create(data_dir)
  saveRDS(list(), file.path(data_dir, "model_parameters.rds"))
  writeLines(
    c(
      ",0",
      "0,A1bg",
      "1,A1bg"
    ),
    file.path(data_dir, "features_model_training_17.csv")
  )

  expect_error(
    scop:::load_cytotrace2_data(data_dir, verbose = FALSE),
    "duplicated feature names"
  )
})
