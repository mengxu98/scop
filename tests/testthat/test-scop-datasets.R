make_scop_dataset_fixture <- function() {
  base <- tempfile("scop-datasets-")
  collection_dir <- file.path(base, "Xenium")
  dir.create(collection_dir, recursive = TRUE)
  object <- list(name = "example", value = 1:3)
  rds_file <- file.path(collection_dir, "example.rds")
  saveRDS(object, rds_file)
  manifest <- data.frame(
    dataset = "example_dataset",
    file = "example.rds",
    object = "list",
    source = "unit-test",
    cells = NA_integer_,
    features = NA_integer_,
    assay = NA_character_,
    sha256 = unname(tools::sha256sum(rds_file)),
    size_bytes = file.info(rds_file)$size,
    description = "Example dataset fixture",
    stringsAsFactors = FALSE
  )
  utils::write.table(
    manifest,
    file = file.path(collection_dir, "manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  list(base = base, object = object)
}

test_that("ListScopDatasets reads a local datasets manifest", {
  fixture <- make_scop_dataset_fixture()
  manifest <- ListScopDatasets(
    collection = "Xenium",
    datasets_base_url = fixture$base
  )

  expect_equal(manifest$dataset, "example_dataset")
  expect_equal(manifest$file, "example.rds")
  expect_true("sha256" %in% colnames(manifest))
})

test_that("LoadScopDataset downloads, validates, caches, and loads an object", {
  fixture <- make_scop_dataset_fixture()
  cache_dir <- tempfile("scop-cache-")
  object <- LoadScopDataset(
    "example_dataset",
    collection = "Xenium",
    cache_dir = cache_dir,
    datasets_base_url = fixture$base,
    verbose = FALSE
  )
  path <- LoadScopDataset(
    "example_dataset",
    collection = "Xenium",
    cache_dir = cache_dir,
    datasets_base_url = fixture$base,
    return_path = TRUE,
    verbose = FALSE
  )

  expect_equal(object, fixture$object)
  expect_true(file.exists(path))
  expect_equal(basename(path), "example.rds")
})

test_that("LoadScopDataset rejects files that do not match the manifest", {
  fixture <- make_scop_dataset_fixture()
  rds_file <- file.path(fixture$base, "Xenium", "example.rds")
  saveRDS(list(name = "changed"), rds_file)

  expect_error(
    LoadScopDataset(
      "example_dataset",
      collection = "Xenium",
      cache_dir = tempfile("scop-cache-"),
      datasets_base_url = fixture$base,
      verbose = FALSE
    ),
    "failed size or sha256"
  )
})
