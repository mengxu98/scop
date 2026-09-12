#' @title Run SpatialEcoTyper spatial ecotype analysis
#'
#' @description
#' Discover, recover, or deconvolve spatial ecotypes using SpatialEcoTyper.
#'
#' @md
#' @inheritParams RunStandardWorkflow
#' @inheritParams thisutils::log_message
#' @inheritParams scop-params
#' @param srt Deprecated alias for `object`; supply exactly one of the two. It
#' will be removed in scop 1.0.0.
#' @param object A `Seurat` object. For `mode = "deconvolute"`, a numeric
#' expression matrix can also be supplied.
#' @param mode SpatialEcoTyper workflow. `"single"` runs single-sample de novo
#' discovery, `"multi"` runs conserved ecotype discovery across samples,
#' `"recover"` recovers pretrained SE labels, and `"deconvolute"` infers SE
#' abundances from bulk or spot-level expression.
#' @param assay Assay used for expression extraction. If `NULL`, the default
#' assay is used.
#' @param layer Assay layer used for expression extraction.
#' @param celltype.by Metadata column containing cell type annotations. Required
#' for `"single"`, `"multi"`, and `"recover"` unless `celltypes` is supplied for
#' `"recover"`.
#' @param sample.by Metadata column identifying samples for `mode = "multi"`.
#' @param x.by,y.by Metadata columns containing single-cell spatial coordinates.
#'   Used only when no image is present; image-backed discovery uses raw image
#'   coordinates through [SpatialCoordinates()].
#' @param image Optional image name for single discovery, or named sample-to-image
#'   map for multi discovery. Multi discovery resolves each sample independently.
#' @param dat Optional expression matrix used by `"recover"` or
#' `"deconvolute"`. If `NULL`, expression is extracted from `srt`.
#' @param celltypes Optional named vector of cell types passed to
#' `SpatialEcoTyper::RecoverSE()`.
#' @param features Optional feature vector used to subset the expression matrix.
#' @param outprefix Output prefix passed to `SpatialEcoTyper`. Use `NULL` to
#' avoid writing single-sample result files to the working directory.
#' @param outdir Output directory passed to multi-sample SpatialEcoTyper.
#' `NULL` creates a temporary directory.
#' @param radius Spatial neighborhood radius, in the same units as `x.by` and
#' `y.by`.
#' @param resolution Louvain clustering resolution used by `SpatialEcoTyper`.
#' @param nfeatures Number of variable features used by `SpatialEcoTyper`.
#' @param min.cts.per.region Minimum number of cell types required in a spatial
#' neighborhood.
#' @param npcs Number of principal components used for similarity networks.
#' @param min.cells Minimum number of cells or spatial meta-cells expressing a
#' feature.
#' @param min.features Minimum number of features detected in a cell or spatial
#' meta-cell.
#' @param iterations Number of similarity network fusion iterations.
#' @param minibatch Number of columns processed per mini-batch in SNF.
#' @param cores Number of CPU cores used by `SpatialEcoTyper`.
#' @param ncores Deprecated alias for `cores`; supply exactly one of the two. It
#' will be removed in scop 1.0.0.
#' @param grid.size Spatial grid size used to discretize coordinates.
#' @param filter.region.by.celltypes Optional cell types used to restrict spatial
#' neighborhoods.
#' @param k Number of spatial nearest neighbors used to construct spatial
#' meta-cells.
#' @param k.sn Number of nearest neighbors used to construct similarity networks.
#' @param dropcell Whether cells without spatial ecotype assignments are removed
#' from the returned `SpatialEcoTyper` metadata.
#' @param normalization.method,nmf_ranks,nrun.per.rank,min.coph,Region,downsample.by.region,subresolution,seed
#' Parameters passed to `SpatialEcoTyper::MultiSpatialEcoTyper()`.
#' @param scale Whether to scale expression for `"recover"` and
#' `"deconvolute"`.
#' @param Ws Pretrained basis matrices passed to `SpatialEcoTyper::RecoverSE()`.
#' @param ncell.per.run Number of cells processed per run by
#' `SpatialEcoTyper::RecoverSE()`.
#' @param min.score Minimum prediction score passed to
#' `SpatialEcoTyper::RecoverSE()`.
#' @param W Pretrained basis matrix passed to
#' `SpatialEcoTyper::DeconvoluteSE()`.
#' @param nsample.per.run Number of samples processed per run by
#' `SpatialEcoTyper::DeconvoluteSE()`.
#' @param sum2one Whether inferred SE abundances are normalized to sum to one.
#' @param prefix Prefix used for output metadata columns.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param store_results Whether to store raw results in `srt@tools`.
#' @param allow_partial Whether to allow missing SE labels for cells absent from
#' returned `SpatialEcoTyper` metadata. Default is `FALSE` to avoid silent
#' partial annotations.
#' @param ... Additional arguments passed to the selected SpatialEcoTyper
#' function.
#'
#' @details
#' Discovery requires cell-resolved spatial expression, cell-type labels,
#' and coordinates; multi-sample discovery also requires sample labels.
#' Recovery and deconvolution require the pretrained basis matrices `Ws`
#' and `W`, respectively.
#'
#' @return A `Seurat` object with SpatialEcoTyper results in metadata and raw
#' results stored in `srt@tools[[tool_name]]` when `store_results = TRUE`.
#' For matrix input with `mode = "deconvolute"`, the abundance matrix is
#' returned.
#' @seealso [SpatialSpotPlot()], [CellStatPlot()],
#' [SpatialEcoTyper tutorials](https://digitalcytometry.github.io/spatialecotyper/)
#' @export
#'
RunSpatialEcoTyper <- function(
  object,
  mode = c("single", "multi", "recover", "deconvolute"),
  assay = NULL,
  layer = "data",
  celltype.by = NULL,
  sample.by = NULL,
  x.by = "X",
  y.by = "Y",
  dat = NULL,
  celltypes = NULL,
  features = NULL,
  outprefix = NULL,
  outdir = NULL,
  radius = 50,
  resolution = 0.5,
  nfeatures = 300,
  min.cts.per.region = 2,
  npcs = 20,
  min.cells = 5,
  min.features = 10,
  iterations = 10,
  minibatch = 5000,
  cores = 4,
  ncores = NULL,
  grid.size = round(radius * 1.4),
  filter.region.by.celltypes = NULL,
  k = 20,
  k.sn = 50,
  dropcell = FALSE,
  normalization.method = "None",
  nmf_ranks = 10,
  nrun.per.rank = 30,
  min.coph = 0.95,
  Region = NULL,
  downsample.by.region = TRUE,
  subresolution = 30,
  seed = 1,
  scale = TRUE,
  Ws = NULL,
  ncell.per.run = 500,
  min.score = 0.6,
  W = NULL,
  nsample.per.run = 500,
  sum2one = TRUE,
  prefix = "SpatialEcoTyper",
  tool_name = "SpatialEcoTyper",
  store_results = TRUE,
  allow_partial = FALSE,
  verbose = TRUE,
  ...,
  image = NULL,
  srt = NULL
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  if (!is.null(ncores)) {
    .Deprecated(msg = paste0("`ncores` is deprecated; use `cores` instead. ",
      "It will be removed in scop 1.0.0."))
    cores <- ncores
  }
  mode <- match.arg(mode)

  has_seurat <- inherits(srt, "Seurat")
  if (!has_seurat && !identical(mode, "deconvolute")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object unless {.arg mode = 'deconvolute'}",
      message_type = "error"
    )
  }
  validate_scalar_string(prefix, "prefix", require_character = FALSE)
  validate_scalar_string(tool_name, "tool_name", require_character = FALSE)

  assay <- if (has_seurat) assay %||% SeuratObject::DefaultAssay(srt) else assay
  normdata <- spatialecotyper_get_data(
    srt = srt,
    dat = dat,
    assay = assay,
    layer = layer,
    features = features,
    require_seurat = !identical(mode, "deconvolute")
  )

  coordinate_input <- NULL
  if (identical(mode, "single")) {
    coordinate_input <- SpatialCoordinates(srt, image = image, coord.cols = c(x.by, y.by))
    cells <- intersect(colnames(normdata), coordinate_input$data$cell_id)
    if (length(cells) == 0L) {
      log_message("No expression cells match the selected image", message_type = "error")
    }
    if (!is.null(sample.by)) {
      spatialecotyper_check_meta_columns(srt[[]], sample.by)
      samples <- as.character(srt[[]][cells, sample.by])
      if (anyNA(samples) || any(!nzchar(samples)) || length(unique(samples)) != 1L) {
        log_message("Single discovery requires one sample; use mode = 'multi'", message_type = "error")
      }
    }
    normdata <- normdata[, cells, drop = FALSE]
    coordinate_input$data <- coordinate_input$data[cells, , drop = FALSE]
    coordinate_input$sources <- list(single = c(coordinate_input$source, list(transform = coordinate_input$transform)))
  } else if (identical(mode, "multi")) {
    validate_scalar_string(sample.by, "sample.by", require_character = FALSE)
    coordinate_input <- spatial_sample_coords(srt,
      sample.by = sample.by, image = image,
      coord.cols = c(x.by, y.by)
    )
  }
  check_r("digitalcytometry/SpatialEcoTyper", verbose = FALSE)

  if (identical(mode, "single")) {
    return(spatialecotyper_run_single(
      srt = srt,
      normdata = normdata,
      assay = assay,
      layer = layer,
      celltype.by = celltype.by,
      x.by = x.by,
      y.by = y.by,
      features = features,
      outprefix = outprefix,
      radius = radius,
      resolution = resolution,
      nfeatures = nfeatures,
      min.cts.per.region = min.cts.per.region,
      npcs = npcs,
      min.cells = min.cells,
      min.features = min.features,
      iterations = iterations,
      minibatch = minibatch,
      cores = cores,
      grid.size = grid.size,
      filter.region.by.celltypes = filter.region.by.celltypes,
      k = k,
      k.sn = k.sn,
      dropcell = dropcell,
      prefix = prefix,
      tool_name = tool_name,
      store_results = store_results,
      allow_partial = allow_partial,
      verbose = verbose,
      coordinate_input = coordinate_input,
      ...
    ))
  }

  if (identical(mode, "multi")) {
    return(spatialecotyper_run_multi(
      srt = srt,
      normdata = normdata,
      assay = assay,
      layer = layer,
      celltype.by = celltype.by,
      sample.by = sample.by,
      x.by = x.by,
      y.by = y.by,
      features = features,
      outdir = outdir,
      radius = radius,
      npcs = npcs,
      min.cells = min.cells,
      iterations = iterations,
      min.cts.per.region = min.cts.per.region,
      nfeatures = nfeatures,
      min.features = min.features,
      minibatch = minibatch,
      cores = cores,
      grid.size = grid.size,
      filter.region.by.celltypes = filter.region.by.celltypes,
      k = k,
      k.sn = k.sn,
      dropcell = dropcell,
      normalization.method = normalization.method,
      nmf_ranks = nmf_ranks,
      nrun.per.rank = nrun.per.rank,
      min.coph = min.coph,
      Region = Region,
      downsample.by.region = downsample.by.region,
      subresolution = subresolution,
      seed = seed,
      prefix = prefix,
      tool_name = tool_name,
      store_results = store_results,
      allow_partial = allow_partial,
      verbose = verbose,
      coordinate_input = coordinate_input,
      ...
    ))
  }

  if (identical(mode, "recover")) {
    return(spatialecotyper_run_recover(
      srt = srt,
      normdata = normdata,
      assay = assay,
      layer = layer,
      celltype.by = celltype.by,
      celltypes = celltypes,
      features = features,
      scale = scale,
      Ws = Ws,
      ncell.per.run = ncell.per.run,
      min.score = min.score,
      cores = cores,
      prefix = prefix,
      tool_name = tool_name,
      store_results = store_results,
      allow_partial = allow_partial,
      verbose = verbose,
      ...
    ))
  }

  spatialecotyper_run_deconvolute(
    srt = srt,
    normdata = normdata,
    assay = assay,
    layer = layer,
    features = features,
    scale = scale,
    W = W,
    nsample.per.run = nsample.per.run,
    sum2one = sum2one,
    cores = cores,
    prefix = prefix,
    tool_name = tool_name,
    store_results = store_results,
    verbose = verbose,
    ...
  )
}


spatialecotyper_run_single <- function(
  srt,
  normdata,
  assay,
  layer,
  celltype.by,
  x.by,
  y.by,
  features,
  outprefix,
  radius,
  resolution,
  nfeatures,
  min.cts.per.region,
  npcs,
  min.cells,
  min.features,
  iterations,
  minibatch,
  cores,
  grid.size,
  filter.region.by.celltypes,
  k,
  k.sn,
  dropcell,
  prefix,
  tool_name,
  store_results,
  allow_partial,
  verbose,
  coordinate_input,
  ...
) {
  metadata <- spatialecotyper_coordinate_metadata(
    srt = srt,
    cells = colnames(normdata),
    celltype.by = celltype.by,
    coords = coordinate_input$data
  )
  log_message(
    "Run {.pkg SpatialEcoTyper} single-sample discovery on {.val {ncol(normdata)}} cells and {.val {nrow(normdata)}} features",
    verbose = verbose
  )
  spatialecotyper_fun <- get_namespace_fun(
    "SpatialEcoTyper",
    "SpatialEcoTyper"
  )
  compatibility_env <- new.env(parent = environment(spatialecotyper_fun))
  for (symbol in c("summarise", "summarize", "across", "slice", "desc")) {
    compatibility_env[[symbol]] <- get_namespace_fun("dplyr", symbol)
  }
  compatibility_env[["where"]] <- get_namespace_fun("tidyselect", "where")
  compatibility_env[["%>%"]] <- get_namespace_fun("magrittr", "%>%")
  for (symbol in c("mostFrequent", "GetSpatialMetacells", "GetPCList")) {
    backend_fun <- get_namespace_fun("SpatialEcoTyper", symbol)
    if (!is.function(backend_fun)) {
      log_message(
        "The installed {.pkg SpatialEcoTyper} API is incompatible; missing function {.fn {symbol}}",
        message_type = "error"
      )
    }
    environment(backend_fun) <- compatibility_env
    compatibility_env[[symbol]] <- backend_fun
  }
  environment(spatialecotyper_fun) <- compatibility_env
  result <- spatialecotyper_fun(
    normdata = normdata,
    metadata = metadata,
    outprefix = outprefix,
    radius = radius,
    resolution = resolution,
    nfeatures = nfeatures,
    min.cts.per.region = min.cts.per.region,
    npcs = npcs,
    min.cells = min.cells,
    min.features = min.features,
    iterations = iterations,
    minibatch = minibatch,
    ncores = cores,
    grid.size = grid.size,
    filter.region.by.celltypes = filter.region.by.celltypes,
    k = k,
    k.sn = k.sn,
    dropcell = dropcell,
    ...
  )

  result_metadata <- spatialecotyper_extract_result_metadata(result)
  add_meta <- spatialecotyper_match_result_columns(
    result_metadata = result_metadata,
    cells = colnames(normdata),
    columns = "SE",
    names_out = paste0(prefix, "_SE"),
    allow_partial = allow_partial,
    verbose = verbose
  )
  srt <- Seurat::AddMetaData(srt, metadata = add_meta)

  srt <- spatialecotyper_store_tool(
    srt = srt,
    tool_name = tool_name,
    store_results = store_results,
    result = result,
    metadata = result_metadata,
    coordinate_input = coordinate_input,
    parameters = list(
      mode = "single",
      assay = assay,
      layer = layer,
      celltype.by = celltype.by,
      x.by = x.by,
      y.by = y.by,
      features = features,
      outprefix = outprefix,
      radius = radius,
      resolution = resolution,
      nfeatures = nfeatures,
      min.cts.per.region = min.cts.per.region,
      npcs = npcs,
      min.cells = min.cells,
      min.features = min.features,
      iterations = iterations,
      minibatch = minibatch,
      cores = cores,
      grid.size = grid.size,
      filter.region.by.celltypes = filter.region.by.celltypes,
      k = k,
      k.sn = k.sn,
      dropcell = dropcell,
      prefix = prefix,
      tool_name = tool_name,
      allow_partial = allow_partial
    )
  )
  log_message(
    "{.pkg SpatialEcoTyper} SE labels stored in metadata column {.val {paste0(prefix, '_SE')}}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

spatialecotyper_run_multi <- function(
  srt,
  normdata,
  assay,
  layer,
  celltype.by,
  sample.by,
  x.by,
  y.by,
  features,
  outdir,
  radius,
  npcs,
  min.cells,
  iterations,
  min.cts.per.region,
  nfeatures,
  min.features,
  minibatch,
  cores,
  grid.size,
  filter.region.by.celltypes,
  k,
  k.sn,
  dropcell,
  normalization.method,
  nmf_ranks,
  nrun.per.rank,
  min.coph,
  Region,
  downsample.by.region,
  subresolution,
  seed,
  prefix,
  tool_name,
  store_results,
  allow_partial,
  verbose,
  coordinate_input,
  ...
) {
  validate_scalar_string(sample.by, "sample.by", require_character = FALSE)
  if (!is.null(Region)) {
    validate_scalar_string(Region, "Region", require_character = FALSE)
  }
  meta <- srt[[]]
  spatialecotyper_check_meta_columns(
    meta = meta,
    cols = c(celltype.by, sample.by, Region)
  )
  samples <- as.character(meta[colnames(normdata), sample.by, drop = TRUE])
  if (any(is.na(samples) | !nzchar(samples))) {
    log_message(
      "{.arg sample.by} contains missing sample labels",
      message_type = "error"
    )
  }
  sample_levels <- unique(samples)
  data_list <- vector("list", length(sample_levels))
  metadata_list <- vector("list", length(sample_levels))
  names(data_list) <- sample_levels
  names(metadata_list) <- sample_levels
  for (sample_name in sample_levels) {
    sample_cells <- colnames(normdata)[samples == sample_name]
    data_list[[sample_name]] <- normdata[, sample_cells, drop = FALSE]
    metadata_list[[sample_name]] <- spatialecotyper_coordinate_metadata(
      srt = srt,
      cells = sample_cells,
      celltype.by = celltype.by,
      coords = coordinate_input$data
    )
    if (!is.null(Region)) {
      metadata_list[[sample_name]][[Region]] <- meta[sample_cells, Region, drop = TRUE]
    }
  }
  outdir <- outdir %||% {
    .inline0 <- prefix
    file.path(
      tempdir(),
      paste0(.inline0, "_", format(Sys.time(), "%Y%m%d%H%M%S"))
    )
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  log_message(
    "Run {.pkg SpatialEcoTyper} multi-sample discovery on {.val {length(data_list)}} samples",
    verbose = verbose
  )
  multi_fun <- get_namespace_fun("SpatialEcoTyper", "MultiSpatialEcoTyper")
  result <- multi_fun(
    data_list = data_list,
    metadata_list = metadata_list,
    outdir = outdir,
    normalization.method = normalization.method,
    nmf_ranks = nmf_ranks,
    nrun.per.rank = nrun.per.rank,
    min.coph = min.coph,
    radius = radius,
    min.cts.per.region = min.cts.per.region,
    nfeatures = nfeatures,
    min.features = min.features,
    Region = Region,
    downsample.by.region = downsample.by.region,
    subresolution = subresolution,
    minibatch = minibatch,
    ncores = cores,
    seed = seed,
    filter.region.by.celltypes = filter.region.by.celltypes,
    npcs = npcs,
    min.cells = min.cells,
    iterations = iterations,
    grid.size = grid.size,
    k = k,
    k.sn = k.sn,
    dropcell = dropcell,
    ...
  )
  result_metadata <- spatialecotyper_extract_multi_metadata(result, outdir)
  add_meta <- spatialecotyper_match_result_columns(
    result_metadata = result_metadata,
    cells = colnames(srt),
    columns = c("InitSE", "SE"),
    names_out = paste0(prefix, c("_InitSE", "_SE")),
    allow_partial = allow_partial,
    verbose = verbose
  )
  srt <- Seurat::AddMetaData(srt, metadata = add_meta)
  srt <- spatialecotyper_store_tool(
    srt = srt,
    tool_name = tool_name,
    store_results = store_results,
    result = result,
    metadata = result_metadata,
    coordinate_input = coordinate_input,
    parameters = list(
      mode = "multi",
      assay = assay,
      layer = layer,
      celltype.by = celltype.by,
      sample.by = sample.by,
      x.by = x.by,
      y.by = y.by,
      features = features,
      outdir = outdir,
      radius = radius,
      npcs = npcs,
      min.cells = min.cells,
      iterations = iterations,
      min.cts.per.region = min.cts.per.region,
      nfeatures = nfeatures,
      min.features = min.features,
      minibatch = minibatch,
      cores = cores,
      grid.size = grid.size,
      filter.region.by.celltypes = filter.region.by.celltypes,
      k = k,
      k.sn = k.sn,
      dropcell = dropcell,
      normalization.method = normalization.method,
      nmf_ranks = nmf_ranks,
      nrun.per.rank = nrun.per.rank,
      min.coph = min.coph,
      Region = Region,
      downsample.by.region = downsample.by.region,
      subresolution = subresolution,
      seed = seed,
      prefix = prefix,
      tool_name = tool_name,
      allow_partial = allow_partial
    )
  )
  log_message(
    "{.pkg SpatialEcoTyper} integrated SE labels stored in metadata columns {.val {colnames(add_meta)}}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

spatialecotyper_run_recover <- function(
  srt,
  normdata,
  assay,
  layer,
  celltype.by,
  celltypes,
  features,
  scale,
  Ws,
  ncell.per.run,
  min.score,
  cores,
  prefix,
  tool_name,
  store_results,
  allow_partial,
  verbose,
  ...
) {
  if (is.null(celltypes)) {
    validate_scalar_string(celltype.by, "celltype.by", require_character = FALSE)
    meta <- srt[[]]
    spatialecotyper_check_meta_columns(meta = meta, cols = celltype.by)
    celltypes <- meta[colnames(normdata), celltype.by, drop = TRUE]
  }
  celltypes <- spatialecotyper_prepare_celltypes(celltypes, colnames(normdata))
  log_message(
    "Recover {.pkg SpatialEcoTyper} SE labels for {.val {ncol(normdata)}} cells",
    verbose = verbose
  )
  recover_fun <- get_namespace_fun("SpatialEcoTyper", "RecoverSE")
  result <- recover_fun(
    dat = normdata,
    celltypes = celltypes,
    scale = scale,
    Ws = Ws,
    ncell.per.run = ncell.per.run,
    min.score = min.score,
    ncores = cores,
    ...
  )
  result_metadata <- spatialecotyper_extract_recover_metadata(result)
  add_meta <- spatialecotyper_match_result_columns(
    result_metadata = result_metadata,
    cells = colnames(srt),
    columns = c("InitSE", "SE", "PredScore"),
    names_out = paste0(prefix, c("_InitSE", "_SE", "_PredScore")),
    allow_partial = allow_partial,
    verbose = verbose
  )
  add_meta[[paste0(prefix, "_PredScore")]] <- suppressWarnings(
    as.numeric(add_meta[[paste0(prefix, "_PredScore")]])
  )
  srt <- Seurat::AddMetaData(srt, metadata = add_meta)
  srt <- spatialecotyper_store_tool(
    srt = srt,
    tool_name = tool_name,
    store_results = store_results,
    result = result,
    metadata = result_metadata,
    parameters = list(
      mode = "recover",
      assay = assay,
      layer = layer,
      celltype.by = celltype.by,
      features = features,
      scale = scale,
      Ws = Ws,
      ncell.per.run = ncell.per.run,
      min.score = min.score,
      cores = cores,
      prefix = prefix,
      tool_name = tool_name,
      allow_partial = allow_partial
    )
  )
  log_message(
    "{.pkg SpatialEcoTyper} recovered SE labels stored in metadata column {.val {paste0(prefix, '_SE')}}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

spatialecotyper_run_deconvolute <- function(
  srt,
  normdata,
  assay,
  layer,
  features,
  scale,
  W,
  nsample.per.run,
  sum2one,
  cores,
  prefix,
  tool_name,
  store_results,
  verbose,
  ...
) {
  log_message(
    "Infer {.pkg SpatialEcoTyper} SE abundance for {.val {ncol(normdata)}} samples",
    verbose = verbose
  )
  deconv_fun <- get_namespace_fun("SpatialEcoTyper", "DeconvoluteSE")
  result <- deconv_fun(
    dat = normdata,
    scale = scale,
    W = W,
    nsample.per.run = nsample.per.run,
    sum2one = sum2one,
    cores = cores,
    ...
  )
  abundance <- as.matrix(result)
  if (!inherits(srt, "Seurat")) {
    return(abundance)
  }
  add_meta <- spatialecotyper_abundance_metadata(
    abundance = abundance,
    cells = colnames(srt),
    prefix = prefix
  )
  srt <- Seurat::AddMetaData(srt, metadata = add_meta)
  srt <- spatialecotyper_store_tool(
    srt = srt,
    tool_name = tool_name,
    store_results = store_results,
    result = result,
    metadata = add_meta,
    parameters = list(
      mode = "deconvolute",
      assay = assay,
      layer = layer,
      features = features,
      scale = scale,
      W = W,
      nsample.per.run = nsample.per.run,
      sum2one = sum2one,
      cores = cores,
      prefix = prefix,
      tool_name = tool_name
    )
  )
  log_message(
    "{.pkg SpatialEcoTyper} abundance stored in metadata columns {.val {colnames(add_meta)}}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

spatialecotyper_get_data <- function(
  srt,
  dat = NULL,
  assay = NULL,
  layer = "data",
  features = NULL,
  require_seurat = TRUE
) {
  if (is.null(dat)) {
    if (!inherits(srt, "Seurat")) {
      if (require_seurat) {
        log_message(
          "{.arg srt} must be a {.cls Seurat} object",
          message_type = "error"
        )
      }
      dat <- srt
    } else {
      dat <- GetAssayData5(srt, assay = assay, layer = layer)
    }
  }
  if (is.data.frame(dat)) {
    dat <- as.matrix(dat)
  }
  if (is.null(rownames(dat)) || is.null(colnames(dat))) {
    log_message(
      "Expression data must contain feature and sample/cell names",
      message_type = "error"
    )
  }
  if (!is.null(features)) {
    features <- unique(as.character(features))
    missing_features <- setdiff(features, rownames(dat))
    features <- intersect(features, rownames(dat))
    if (length(features) == 0L) {
      log_message(
        "No requested {.arg features} are present in the selected assay/layer",
        message_type = "error"
      )
    }
    if (length(missing_features) > 0L) {
      log_message(
        "Ignoring {.val {length(missing_features)}} requested features absent from the selected assay/layer",
        message_type = "warning"
      )
    }
    dat <- dat[features, , drop = FALSE]
  }
  dat
}

spatialecotyper_check_meta_columns <- function(meta, cols) {
  missing_cols <- setdiff(cols, colnames(meta))
  if (length(missing_cols) > 0L) {
    log_message(
      "Missing metadata column{?s}: {.val {missing_cols}}",
      message_type = "error"
    )
  }
}

spatialecotyper_build_metadata <- function(meta, celltype.by, x.by, y.by) {
  x <- meta[[x.by]]
  y <- meta[[y.by]]
  celltype <- meta[[celltype.by]]
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.factor(y)) {
    y <- as.character(y)
  }

  metadata <- data.frame(
    X = suppressWarnings(as.numeric(x)),
    Y = suppressWarnings(as.numeric(y)),
    CellType = as.character(celltype),
    row.names = rownames(meta),
    stringsAsFactors = FALSE
  )
  invalid <- !is.finite(metadata$X) |
    !is.finite(metadata$Y) |
    is.na(metadata$CellType) |
    !nzchar(metadata$CellType)
  if (any(invalid)) {
    log_message(
      "SpatialEcoTyper metadata contains {.val {sum(invalid)}} invalid cell{?s}",
      message_type = "error"
    )
  }
  metadata
}

spatialecotyper_coordinate_metadata <- function(srt, cells, celltype.by, coords) {
  validate_scalar_string(celltype.by, "celltype.by", require_character = FALSE)
  spatialecotyper_check_meta_columns(srt[[]], celltype.by)
  if (!all(cells %in% rownames(coords))) {
    log_message("SpatialEcoTyper coordinates are missing expression cells", message_type = "error")
  }
  meta <- data.frame(
    X = coords[cells, "x"], Y = coords[cells, "y"],
    CellType = srt[[]][cells, celltype.by], row.names = cells
  )
  spatialecotyper_build_metadata(meta, "CellType", "X", "Y")
}

spatialecotyper_extract_result_metadata <- function(result) {
  if (
    !is.list(result) ||
      is.null(result$metadata) ||
      !is.data.frame(result$metadata)
  ) {
    log_message(
      "{.pkg SpatialEcoTyper} did not return a valid {.field metadata} data frame",
      message_type = "error"
    )
  }
  metadata <- result$metadata
  if (!"SE" %in% colnames(metadata)) {
    log_message(
      "{.pkg SpatialEcoTyper} did not return metadata with an {.field SE} column",
      message_type = "error"
    )
  }
  if (is.null(rownames(metadata)) || any(!nzchar(rownames(metadata)))) {
    log_message(
      "{.pkg SpatialEcoTyper} returned metadata without cell row names",
      message_type = "error"
    )
  }
  metadata
}

spatialecotyper_extract_multi_metadata <- function(result, outdir) {
  if (is.data.frame(result)) {
    metadata <- result
  } else if (is.list(result) && is.data.frame(result$metadata)) {
    metadata <- result$metadata
  } else if (is.list(result) && is.data.frame(result$metadatas)) {
    metadata <- result$metadatas
  } else {
    rds_file <- file.path(outdir, "MultiSE_metadata_final.rds")
    tsv_file <- file.path(outdir, "MultiSE_metadata_final.tsv")
    if (file.exists(rds_file)) {
      metadata <- readRDS(rds_file)
    } else if (file.exists(tsv_file)) {
      metadata <- utils::read.delim(tsv_file, check.names = FALSE)
    } else {
      log_message(
        "{.pkg SpatialEcoTyper} did not return integrated metadata and no MultiSE metadata file was found",
        message_type = "error"
      )
    }
  }
  required <- c("InitSE", "SE")
  missing_cols <- setdiff(required, colnames(metadata))
  if (length(missing_cols) > 0L) {
    log_message(
      "{.pkg SpatialEcoTyper} integrated metadata is missing column{?s}: {.val {missing_cols}}",
      message_type = "error"
    )
  }
  metadata
}

spatialecotyper_extract_recover_metadata <- function(result) {
  if (!is.data.frame(result)) {
    log_message(
      "{.pkg SpatialEcoTyper} did not return a valid recovery metadata data frame",
      message_type = "error"
    )
  }
  required <- c("InitSE", "SE", "PredScore")
  missing_cols <- setdiff(required, colnames(result))
  if (length(missing_cols) > 0L) {
    log_message(
      "{.pkg SpatialEcoTyper} recovery metadata is missing column{?s}: {.val {missing_cols}}",
      message_type = "error"
    )
  }
  result
}

spatialecotyper_result_ids <- function(result_metadata) {
  if ("CID" %in% colnames(result_metadata)) {
    ids <- as.character(result_metadata$CID)
  } else {
    ids <- rownames(result_metadata)
  }
  if (is.null(ids) || any(is.na(ids) | !nzchar(ids))) {
    log_message(
      "{.pkg SpatialEcoTyper} returned metadata without usable cell IDs",
      message_type = "error"
    )
  }
  ids
}

spatialecotyper_match_result_columns <- function(
  result_metadata,
  cells,
  columns,
  names_out,
  allow_partial = FALSE,
  verbose = TRUE
) {
  missing_cols <- setdiff(columns, colnames(result_metadata))
  if (length(missing_cols) > 0L) {
    log_message(
      "{.pkg SpatialEcoTyper} metadata is missing column{?s}: {.val {missing_cols}}",
      message_type = "error"
    )
  }
  ids <- spatialecotyper_result_ids(result_metadata)
  result_metadata <- result_metadata[!duplicated(ids), , drop = FALSE]
  rownames(result_metadata) <- ids[!duplicated(ids)]
  out <- data.frame(row.names = cells, check.names = FALSE)
  common <- intersect(cells, rownames(result_metadata))
  missing_cells <- setdiff(cells, rownames(result_metadata))
  if (length(missing_cells) > 0L && !isTRUE(allow_partial)) {
    log_message(
      "{.pkg SpatialEcoTyper} returned no SE result for {.val {length(missing_cells)}} cell{?s}. Set {.arg allow_partial = TRUE} to keep missing values as {.val NA}.",
      message_type = "error"
    )
  }
  if (length(missing_cells) > 0L) {
    log_message(
      "{.pkg SpatialEcoTyper} returned no SE result for {.val {length(missing_cells)}} cell{?s}; storing {.val NA}",
      message_type = "warning",
      verbose = verbose
    )
  }
  for (i in seq_along(columns)) {
    values <- stats::setNames(rep(NA, length(cells)), cells)
    values[common] <- result_metadata[common, columns[[i]], drop = TRUE]
    missing_values <- names(values)[is.na(values) | !nzchar(as.character(values))]
    if (length(missing_values) > 0L && !isTRUE(allow_partial)) {
      log_message(
        "{.pkg SpatialEcoTyper} returned missing values for {.val {length(missing_values)}} cell{?s} in {.field {columns[[i]]}}. Set {.arg allow_partial = TRUE} to keep missing values as {.val NA}.",
        message_type = "error"
      )
    }
    out[[names_out[[i]]]] <- values
  }
  out
}

spatialecotyper_prepare_celltypes <- function(celltypes, cells) {
  if (is.null(names(celltypes))) {
    if (length(celltypes) != length(cells)) {
      log_message(
        "{.arg celltypes} must have one value per expression column when unnamed",
        message_type = "error"
      )
    }
    names(celltypes) <- cells
  } else {
    if (anyDuplicated(names(celltypes)) > 0L) {
      log_message(
        "{.arg celltypes} names must be unique",
        message_type = "error"
      )
    }
    missing_cells <- setdiff(cells, names(celltypes))
    if (length(missing_cells) > 0L) {
      log_message(
        "{.arg celltypes} is missing {.val {length(missing_cells)}} expression column name{?s}",
        message_type = "error"
      )
    }
  }
  celltypes <- celltypes[cells]
  if (any(is.na(celltypes) | !nzchar(as.character(celltypes)))) {
    log_message(
      "{.arg celltypes} contains missing cell type labels",
      message_type = "error"
    )
  }
  as.character(celltypes)
}

spatialecotyper_abundance_metadata <- function(abundance, cells, prefix) {
  if (is.null(rownames(abundance)) || is.null(colnames(abundance))) {
    log_message(
      "{.pkg SpatialEcoTyper} abundance matrix must contain row and column names",
      message_type = "error"
    )
  }
  if (all(cells %in% colnames(abundance))) {
    abundance <- t(abundance[, cells, drop = FALSE])
  } else if (all(cells %in% rownames(abundance))) {
    abundance <- abundance[cells, , drop = FALSE]
  } else {
    log_message(
      "{.pkg SpatialEcoTyper} abundance matrix sample names do not match Seurat cells",
      message_type = "error"
    )
  }
  se_names <- make.names(colnames(abundance), unique = TRUE)
  colnames(abundance) <- paste0(prefix, "_Abundance_", se_names)
  add_meta <- as.data.frame(abundance, check.names = FALSE)
  if (ncol(add_meta) > 0L) {
    dominant <- colnames(abundance)[max.col(as.matrix(abundance), ties.method = "first")]
    dominant <- sub(paste0("^", prefix, "_Abundance_"), "", dominant)
    add_meta[[paste0(prefix, "_DominantSE")]] <- dominant
  }
  add_meta
}

spatialecotyper_store_tool <- function(
  srt,
  tool_name,
  store_results,
  result,
  metadata,
  parameters,
  coordinate_input = NULL
) {
  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      result = result,
      metadata = metadata,
      parameters = parameters
    )
    if (!is.null(coordinate_input)) {
      srt@tools[[tool_name]]$coordinates <- coordinate_input$data
      srt@tools[[tool_name]]$source <- list(
        coordinate_space = "raw", samples = coordinate_input$sources,
        coordinate_contract_version = .spatial_coordinate_contract_version
      )
      srt@tools[[tool_name]] <- spatial_tag_coordinate_contract(srt@tools[[tool_name]])
    }
  }
  srt
}
