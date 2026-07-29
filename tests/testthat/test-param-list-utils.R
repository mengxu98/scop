test_that("validate_named_param_list accepts supported parameter containers", {
  expect_invisible(validate_named_param_list(list(), "params", require_list = TRUE))
  expect_invisible(
    validate_named_param_list(list(alpha = 1), "params", require_list = TRUE)
  )
  expect_invisible(validate_named_param_list(c(alpha = 1), "params"))
})

test_that("validate_named_param_list preserves type and name errors", {
  expect_error(
    validate_named_param_list(1, "params", require_list = TRUE),
    "params.*must be a list"
  )
  expect_error(
    validate_named_param_list(list(1), "params", require_list = TRUE),
    "params.*must contain named arguments only"
  )
  expect_error(
    validate_named_param_list(
      1,
      "params",
      require_list = TRUE,
      type_message = "must be a named list"
    ),
    "params.*must be a named list"
  )
})

test_that("shared named-list validation covers Giotto and standard workflows", {
  expect_invisible(validate_named_list(list(alpha = 1), "params"))
  expect_invisible(validate_named_list(list(), "params"))
  expect_error(validate_named_list(1, "params"), "params.*must be a named list")
  expect_error(validate_named_list(list(1), "params"), "params.*must be a named list")
})

test_that("scalar validation helpers preserve supported backend rules", {
  expect_invisible(validate_scalar_string("value", "name"))
  expect_invisible(
    validate_scalar_string(1, "name", require_character = FALSE)
  )
  expect_error(validate_scalar_string(1, "name"), "name.*single non-empty string")
  expect_error(
    validate_scalar_string(
      "",
      "name",
      message = "must be a non-empty character string"
    ),
    "name.*non-empty character string"
  )

  expect_invisible(validate_scalar_flag(TRUE, "flag"))
  expect_error(validate_scalar_flag(1, "flag"), "flag.*TRUE or FALSE")
  expect_error(validate_scalar_flag(NA, "flag"), "flag.*TRUE or FALSE")

  expect_identical(validate_scalar_integer(2L, "count"), 2L)
  expect_identical(
    validate_scalar_integer(
      0,
      "count",
      minimum = 0L,
      message = "must be a non-negative integer"
    ),
    0L
  )
  expect_error(validate_scalar_integer(2.8, "count"), "positive integer")
  expect_error(
    validate_scalar_integer(
      -1,
      "count",
      minimum = 0L,
      message = "must be a non-negative integer"
    ),
    "non-negative integer"
  )
})

test_that("shared Seurat validation preserves backend errors", {
  expect_invisible(validate_seurat_object(structure(list(), class = "Seurat")))
  expect_error(
    validate_seurat_object(list()),
    "srt.*must be a.*Seurat.*object"
  )
})

test_that("shared backend utilities preserve edge-case behavior", {
  expect_identical(validate_positive_integer(2.8, "cores"), 2L)
  expect_error(validate_positive_integer(0, "cores"), "positive integer")

  expect_identical(collapse_parameter_value(NULL), NA_character_)
  expect_identical(collapse_parameter_value(character()), "")
  expect_identical(collapse_parameter_value(c("a", "b")), "a,b")

  df <- data.frame(Gene_ID = c("1", "bad"), value = 1:2)
  expect_identical(
    pick_case_insensitive_column(df, c("gene_id", "gene")),
    "Gene_ID"
  )
  expect_null(pick_case_insensitive_column(df, "missing"))
  expect_equal(
    pick_numeric_column(df, c("missing", "Gene_ID")),
    c(1, NA_real_)
  )
  expect_true(all(is.na(pick_numeric_column(df, "missing"))))

  expect_identical(
    resource_ref("https://example.test/", "/a.rds"),
    "https://example.test/a.rds"
  )
})

test_that("resolve_reference_labels preserves missing and metadata errors", {
  reference <- SeuratObject::CreateSeuratObject(
    Matrix::Matrix(matrix(
      1:4,
      nrow = 2,
      ncol = 2,
      dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
    ), sparse = TRUE)
  )
  reference$label <- c("A", "B")

  expect_identical(
    as.character(resolve_reference_labels(reference, "label")),
    c("A", "B")
  )
  expect_error(
    resolve_reference_labels(reference),
    "reference_label.*single reference metadata column"
  )
  expect_error(
    resolve_reference_labels(reference, "missing"),
    "reference_label.*missing.*not present"
  )
})

test_that("resolve_common_features preserves feature order and filtering", {
  srt <- list(RNA = matrix(
    numeric(),
    nrow = 3,
    dimnames = list(c("gene3", "gene1", "gene2"), NULL)
  ))
  reference <- list(RNA = matrix(
    numeric(),
    nrow = 3,
    dimnames = list(c("gene2", "gene3", "gene4"), NULL)
  ))

  expect_identical(
    resolve_common_features(srt, reference, "RNA", "RNA"),
    c("gene3", "gene2")
  )
  expect_identical(
    resolve_common_features(
      srt,
      reference,
      "RNA",
      "RNA",
      features = c("gene2", "gene1")
    ),
    "gene2"
  )
})
