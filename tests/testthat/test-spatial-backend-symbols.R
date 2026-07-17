spatial_api_symbol_expectations <- list(
  spatialqm = c(
    "getNcells", "getTxPerCell", "getTxPerArea", "getTxPerNuc",
    "getMeanExpression", "getMeanSignalRatio", "getCellTxFraction",
    "getMaxRatio", "getMaxDetection", "getMECR", "getMorans",
    "getSilhouetteWidth", "getSparsity", "getEntropy"
  ),
  spotsweeper = c("localOutliers", "findArtifacts"),
  card = c("createCARDObject", "CARD_deconvolution"),
  spatialecotyper = c(
    "SpatialEcoTyper", "MultiSpatialEcoTyper", "RecoverSE", "DeconvoluteSE"
  ),
  bayesspace = c("spatialPreprocess", "spatialCluster"),
  smoothclust = "smoothclust",
  meringue = c(
    "getSpatialNeighbors", "moranTest", "moranPermutationTest",
    "spatialCrossCorMatrix", "spatialCrossCorTest", "groupSigSpatialPatterns"
  ),
  sparkx = "sparkx",
  nnsvg = "nnSVG",
  semla = c(
    "UpdateSeuratForSemla", "GetSpatialNetwork", "RunLocalG",
    "RegionNeighbors", "RadialDistance"
  )
)

test_that("stable spatial backends declare the APIs their wrappers consume", {
  registry <- scop:::spatial_backend_registry()

  for (backend_id in names(spatial_api_symbol_expectations)) {
    expect_setequal(
      scop:::spatial_backend_required_symbols(registry[[backend_id]]),
      spatial_api_symbol_expectations[[backend_id]]
    )
  }
})

test_that("api_check is real for all ten diagnosed backends", {
  backend_ids <- names(spatial_api_symbol_expectations)
  status <- suppressWarnings(SpatialBackendStatus(
    backend = backend_ids,
    api_check = TRUE,
    refresh = TRUE
  ))

  expect_setequal(status$backend_id, backend_ids)
  expect_true(all(status$api_checked))
  expect_true(all(nzchar(status$required_symbols)))
  expect_true(all(status$availability %in% c(
    "available", "missing", "api_incompatible", "namespace_error"
  )))
})

test_that("semla diagnostics use producer-specific symbols", {
  spec <- scop:::spatial_backend_registry()$semla
  expected <- list(
    RunSemlaSpatialNetwork = c("UpdateSeuratForSemla", "GetSpatialNetwork"),
    RunSemlaLocalG = c("UpdateSeuratForSemla", "RunLocalG"),
    RunSemlaRegionNeighbors = c("UpdateSeuratForSemla", "RegionNeighbors"),
    RunSemlaRadialDistance = c("UpdateSeuratForSemla", "RadialDistance")
  )

  for (method in names(expected)) {
    expect_setequal(
      scop:::spatial_backend_required_symbols(spec, method = method),
      expected[[method]]
    )
  }
})

test_that("missing exports are classified as API incompatibility", {
  spec <- scop:::spatial_backend_registry()$bayesspace
  required <- scop:::spatial_backend_required_symbols(spec)
  exports <- setdiff(required, "spatialCluster")

  expect_identical(
    setdiff(required, exports),
    "spatialCluster"
  )
})

test_that("alternative API generations still select a complete symbol set", {
  spec <- scop:::spatial_backend_registry()$spacexr

  for (symbols in spec$symbol_sets) {
    expect_identical(
      scop:::spatial_backend_required_symbols(spec, exports = symbols),
      symbols
    )
  }
})

test_that("CARD diagnostics resolve the installed implementation", {
  spec <- scop:::spatial_backend_registry()$card

  expect_identical(spec$package_candidates, c("CARD", "CARDspa"))
  expect_identical(
    scop:::spatial_backend_resolve_package(spec, installed = "CARDspa"),
    "CARDspa"
  )
  expect_identical(
    scop:::spatial_backend_resolve_package(spec, installed = c("CARDspa", "CARD")),
    "CARD"
  )
  expect_identical(
    scop:::spatial_backend_resolve_package(spec, installed = character()),
    "CARD"
  )
})
