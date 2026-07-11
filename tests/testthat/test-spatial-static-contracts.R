spatial_dist_calls <- function(file) {
  parsed <- parse(file, keep.source = TRUE)
  tokens <- utils::getParseData(parsed)
  calls <- tokens[
    tokens$token == "SYMBOL_FUNCTION_CALL" & tokens$text == "dist",
    c("line1", "col1", "text"),
    drop = FALSE
  ]
  if (nrow(calls) == 0L) return(character())
  paste0(basename(file), ":", calls$line1)
}

spatial_test_package_root <- function() {
  normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
}

test_that("spatial registry covers the complete public surface", {
  registry <- scop:::spatial_method_registry()
  expect_equal(nrow(registry), length(unique(registry$method)))
  expect_gte(nrow(registry), 60L)
  expect_true(all(registry$method %in% getNamespaceExports("scop")))
  package_root <- spatial_test_package_root()
  expect_true(all(file.exists(file.path(package_root, "R", registry$implementation_files))))
  expect_false(any(
    registry$coordinate_requirement == "distance_sensitive" &
      registry$coordinate_space_current == "unknown"
  ))
  expect_true(all(registry$coordinate_space_current %in% c(
    "raw", "display", "legacy_display", "mixed", "none", "unknown"
  )))
  expect_true(all(registry$backend_requirement %in% c("all", "any")))
  backend_ids <- unique(unlist(strsplit(registry$backend_id, ";", fixed = TRUE)))
  expect_true(all(backend_ids %in% names(scop:::spatial_backend_registry())))
})

test_that("spatial discovery APIs share the registry contract", {
  methods <- ListSpatialMethods()
  expect_equal(methods$method, scop:::spatial_method_registry()$method)
  expect_true(all(c("analysis", "plot", "bridge", "workflow") %in% methods$kind))
  expect_true(all(ListSpatialMethods(kind = "bridge")$kind == "bridge"))
  expect_true(all(grepl("network", ListSpatialMethods(pattern = "network")$method, ignore.case = TRUE) |
    grepl("network", ListSpatialMethods(pattern = "network")$task, ignore.case = TRUE)))

  first <- suppressWarnings(SpatialBackendStatus(api_check = FALSE, refresh = TRUE))
  second <- suppressWarnings(SpatialBackendStatus(api_check = FALSE))
  expect_identical(first, second)
  expect_equal(SpatialBackendStatus(backend = "core")$availability, "available")
})

test_that("Giotto diagnostics and runtime routing share one symbol registry", {
  routing <- scop:::spatial_giotto_symbol_registry()
  backends <- scop:::spatial_backend_registry()
  expect_equal(
    backends$giotto$symbols,
    routing$symbol[routing$package == "Giotto" & routing$required]
  )
  expect_true(all(
    routing$symbol[routing$package == "GiottoClass" & routing$required] %in%
      backends$giotto_class$symbols
  ))
  expect_equal(
    routing$package[match("createGiottoObject", routing$symbol)],
    "GiottoClass"
  )
  expect_equal(
    routing$package[match("normalizeGiotto", routing$symbol)],
    "Giotto"
  )
  reverse_bridge <- SpatialBackendStatus(
    method = "giotto_to_srt",
    api_check = FALSE,
    refresh = TRUE
  )
  class_row <- reverse_bridge[reverse_bridge$backend_id == "giotto_class", , drop = FALSE]
  expect_match(class_row$required_symbols, scop:::spatial_giotto_converter_name())
})

test_that("SpatialResultInfo distinguishes ready empty partial and stale results", {
  counts <- matrix(
    c(3, 1, 0, 2, 0, 4, 1, 0, 2, 1, 3, 0),
    nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt@tools$SpatialNetwork <- list(
    method = "SpatialNetwork",
    active_graph = "empty_edges",
    graphs = list(
      empty_edges = list(
        nodes = data.frame(
          cell_id = colnames(srt), x = c(0, 1, 0, 1), y = c(0, 0, 1, 1),
          image = NA_character_
        ),
        edges = data.frame(from = character(), to = character(), distance = numeric(), weight = numeric()),
        parameters = list(method = "radius", radius = 0.01),
        source = list(image = NULL, coordinate_space = "raw")
      )
    )
  )
  info <- SpatialResultInfo(srt, tool_name = "SpatialNetwork")
  expect_equal(info$result_state, "ready")
  expect_equal(info$n_items, 1L)

  srt@tools$SpatialNetwork$graphs$empty_edges$nodes$cell_id[[1L]] <- "missing_cell"
  expect_equal(
    SpatialResultInfo(srt, tool_name = "SpatialNetwork")$result_state,
    "stale"
  )

  srt@tools$MistyR <- list(method = "MistyR", parameters = list(), results = list())
  expect_equal(nrow(SpatialResultInfo(srt, tool_name = "MistyR")), 0L)
  expect_equal(
    SpatialResultInfo(srt, tool_name = "MistyR", include_empty = TRUE)$result_state,
    "empty"
  )

  srt@tools$StatialKontextual <- list(table = data.frame(value = 1))
  expect_equal(
    SpatialResultInfo(srt, tool_name = "StatialKontextual")$result_state,
    "partial"
  )
})

test_that("no new dense distance calls enter sparse-required spatial paths", {
  registry <- scop:::spatial_method_registry()
  files <- unique(registry$implementation_files[registry$scalability == "sparse_required"])
  calls <- unlist(lapply(
    file.path(spatial_test_package_root(), "R", files),
    spatial_dist_calls
  ), use.names = FALSE)

  # These are the three migration targets tracked by the sparse graph Track.
  # Track C replaces this baseline with character() once the migrations land.
  expect_setequal(
    sub(":.*$", "", calls),
    c("RunSmoothClust.R", "RunSpatialNeighborhood.R", "RunSpatialVariableFeatures.R")
  )
})

test_that("new spatial contract files avoid forbidden optional backend calls", {
  files <- file.path(
    spatial_test_package_root(),
    c("R/SpatialRegistry.R", "R/RunSpatialNetwork.R", "R/SpatialCellPlot.R", "R/SpatialFrameworkConvert.R")
  )
  text <- unlist(lapply(files, readLines, warn = FALSE), use.names = FALSE)
  expect_false(any(grepl("\\b(library|require|requireNamespace)\\s*\\(", text, perl = TRUE)))
})
