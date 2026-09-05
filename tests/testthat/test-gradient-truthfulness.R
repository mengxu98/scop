test_that("gradient flags independently control storage and variable features", {
  counts <- Matrix::Matrix(matrix(1:12, 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4))), sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts)
  other <- srt[["RNA"]]
  SeuratObject::Key(other) <- "other_"
  srt[["other"]] <- other
  srt$x <- 1:4
  srt$y <- c(1, 1, 2, 2)
  srt <- spatial_set_active_variable_features(srt, "other", "g3")
  srt@tools$SpatialGradientFeatures <- list(new = "unchanged")
  payload <- list(screening = data.frame(), significance = data.frame(),
    model_fits = data.frame(), top_variables = data.frame(variable = "g1"), parameters = data.frame())
  local_mocked_bindings(sgf_run_cpp_gradient = function(...) payload, .package = "scop")
  run <- function(...) RunSpatialGradientFeatures(srt, assay = "other", layer = "counts",
    variables = "g1", result_name = "new", verbose = FALSE, ...)
  for (empty in c(FALSE, TRUE)) {
    payload$top_variables <- data.frame(variable = if (empty) character() else "g1")
    for (save in c(FALSE, TRUE)) for (update in c(FALSE, TRUE)) {
      out <- run(store_results = save, set_variable_features = update)
      selected <- SeuratObject::VariableFeatures(out, assay = "other")
      selected <- as.character(selected[!is.na(selected)])
      expect_identical(selected,
        if (update) payload$top_variables$variable else "g3")
      expect_identical(out[["RNA"]], srt[["RNA"]])
      if (update && empty) expect_true(spatial_has_explicit_empty_variable_features(out, "other"))
      if (save) expect_identical(out@tools$SpatialGradientFeatures$new$top_variables,
        payload$top_variables)
      else expect_identical(out@tools, srt@tools)
    }
  }
  for (bad in list(NA, 1, logical(), c(TRUE, FALSE))) {
    expect_error(run(store_results = bad), "single non-missing logical")
    expect_error(run(set_variable_features = bad), "single non-missing logical")
  }
  payload$top_variables <- data.frame(variable = "absent")
  expect_error(run(store_results = FALSE, set_variable_features = TRUE), "absent from")
  payload$screening <- NULL
  expect_error(run(store_results = FALSE), "missing required")
  expect_identical(srt@tools$SpatialGradientFeatures, list(new = "unchanged"))
})

test_that("active variable feature helper supports legacy assays and ranked Assay5 selections", {
  counts <- Matrix::Matrix(matrix(1:12, 3,
    dimnames = list(paste0("g", 1:3), paste0("s", 1:4))), sparse = TRUE)
  for (assay in list(SeuratObject::CreateAssayObject(counts = counts),
                     SeuratObject::CreateAssay5Object(counts = counts))) {
    srt <- Seurat::CreateSeuratObject(assay)
    expect_error(spatial_set_active_variable_features(srt, "RNA", "absent"), "None of the features")
    srt <- spatial_set_active_variable_features(srt, "RNA", c("g2", "g1", "g2"))
    expect_identical(SeuratObject::VariableFeatures(srt), c("g2", "g1"))
    srt <- spatial_set_active_variable_features(srt, "RNA", character())
    expect_true(spatial_has_explicit_empty_variable_features(srt, "RNA"))
  }
})
