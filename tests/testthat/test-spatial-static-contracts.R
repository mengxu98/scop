test_that("spatial registry covers the complete public surface", {
  registry <- scop:::spatial_method_registry()
  expect_equal(nrow(registry), length(unique(registry$method)))
  expect_gte(nrow(registry), 60L)
  expect_true(all(registry$method %in% getNamespaceExports("scop")))
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

test_that("spatial code does not bypass strict image resolution", {
  r_dir <- testthat::test_path("..", "..", "R")
  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  forbidden <- c(
    "image\\s*<-\\s*image\\s*%\\|\\|%\\s*images\\s*\\[",
    "GetTissueCoordinates\\([^)]*images\\s*\\[",
    "images\\s*\\[\\s*\\[?\\s*1L?\\s*\\]?\\s*\\]"
  )
  violations <- unlist(lapply(files, function(path) {
    lines <- readLines(path, warn = FALSE)
    hits <- unique(unlist(lapply(forbidden, grep, x = lines, perl = TRUE)))
    approved <- identical(basename(path), "SpatialCore.R") &
      trimws(lines[hits]) == "image <- images[[1L]]"
    hits <- hits[!approved]
    if (length(hits) == 0L) return(character())
    paste0(basename(path), ":", hits)
  }), use.names = FALSE)
  if (is.null(violations)) violations <- character()
  expect_identical(violations, character())
})

test_that("registered small analyses emit schema-v1 result families", {
  registry <- scop:::spatial_method_registry()
  target <- registry[
    registry$kind == "analysis" &
      registry$status == "stable" &
      !is.na(registry$tool_key) & nzchar(registry$tool_key) &
      !grepl("^RunSemla", registry$method),
    c("method", "task"),
    drop = FALSE
  ]
  collect_build_calls <- function(expr) {
    if (!is.call(expr)) return(list())
    found <- if (identical(expr[[1L]], as.name("spatial_result_build"))) list(expr) else list()
    for (i in seq_along(expr)[-1L]) {
      found <- c(found, collect_build_calls(expr[[i]]))
    }
    found
  }
  namespace <- asNamespace("scop")
  functions <- mget(ls(namespace, all.names = TRUE), namespace, inherits = FALSE)
  functions <- Filter(is.function, functions)
  calls <- unlist(lapply(functions, function(fun) collect_build_calls(body(fun))), recursive = FALSE)
  emitted <- do.call(rbind, lapply(calls, function(call) {
    args <- as.list(call)[-1L]
    result_type <- args[["result_type"]]
    provenance <- args[["provenance"]]
    if (!is.character(result_type) || !is.call(provenance) ||
      !identical(provenance[[1L]], as.name("list"))) return(NULL)
    provenance <- as.list(provenance)[-1L]
    producer <- provenance[["producer"]]
    if (!is.character(producer)) return(NULL)
    data.frame(producer = producer, result_type = result_type, stringsAsFactors = FALSE)
  }))
  expect_true(all(target$method %in% emitted$producer))
  actual <- emitted$result_type[match(target$method, emitted$producer)]
  expect_identical(unname(actual), unname(target$task))
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

  srt@tools$StatialKontextual <- list(
    table = data.frame(imageID = "slice1", test = "A:B", kontextual = 1, r = 10),
    parameters = list(coordinate_space = "raw", image = "slice1")
  )
  actual <- SpatialResultInfo(srt, tool_name = "StatialKontextual")
  expect_identical(actual$coordinate_space, "raw")
  expect_identical(actual$image, "slice1")

  srt@tools$StatialKontextual$cells <- c(colnames(srt), "missing_cell")
  stale <- SpatialResultInfo(srt, tool_name = "StatialKontextual")
  expect_identical(stale$result_state, "stale")
  expect_match(stale$empty_reason, "cells or nodes")
})
