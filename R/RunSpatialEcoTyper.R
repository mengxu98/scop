#' @title Run SpatialEcoTyper spatial ecotype analysis
#'
#' @description
#' Run SpatialEcoTyper workflows through the optional `SpatialEcoTyper`
#' package and write spatial ecotype labels or abundances back to a `Seurat`
#' object when possible.
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object. For `mode = "deconvolute"`, a numeric
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
#' @param ncores Number of CPU cores used by `SpatialEcoTyper`.
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
#' @return A `Seurat` object with SpatialEcoTyper results in metadata and raw
#' results stored in `srt@tools[[tool_name]]` when `store_results = TRUE`.
#' For matrix input with `mode = "deconvolute"`, the abundance matrix is
#' returned.
#' @concept spatial-producer
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' srt <- visium_human_pancreas_sub
#' srt$CellType <- srt$coda_label
#' srt$SpatialEcoTyper_SE <- ifelse(srt$x > stats::median(srt$x), "SE1", "SE2")
#' srt$sample <- ifelse(srt$y > stats::median(srt$y), "slice_a", "slice_b")
#'
#' SpatialEcoTyperSpatialPlot(
#'   srt,
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#' SpatialEcoTyperCompositionPlot(
#'   srt,
#'   group.by = "CellType",
#'   sample.by = "sample",
#'   position = "fill"
#' )
#'
#' srt <- RunSpatialEcoTyper(
#'   srt,
#'   celltype.by = "CellType",
#'   x.by = "x",
#'   y.by = "y",
#'   nfeatures = 100,
#'   ncores = 1,
#'   verbose = FALSE
#' )
#'
#' srt <- RunSpatialEcoTyper(
#'   srt,
#'   mode = "multi",
#'   celltype.by = "CellType",
#'   sample.by = "sample",
#'   x.by = "x",
#'   y.by = "y",
#'   nfeatures = 100,
#'   ncores = 1,
#'   verbose = FALSE
#' )
RunSpatialEcoTyper <- function(
  srt,
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
  ncores = 4,
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
  ...
) {
  mode <- match.arg(mode)
  check_r("digitalcytometry/SpatialEcoTyper", verbose = FALSE)

  has_seurat <- inherits(srt, "Seurat")
  if (!has_seurat && !identical(mode, "deconvolute")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object unless {.arg mode = 'deconvolute'}",
      message_type = "error"
    )
  }
  spatialecotyper_assert_scalar_string(prefix, "prefix")
  spatialecotyper_assert_scalar_string(tool_name, "tool_name")

  assay <- if (has_seurat) assay %||% SeuratObject::DefaultAssay(srt) else assay
  normdata <- spatialecotyper_get_data(
    srt = srt,
    dat = dat,
    assay = assay,
    layer = layer,
    features = features,
    require_seurat = !identical(mode, "deconvolute")
  )

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
      ncores = ncores,
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
      ncores = ncores,
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
      ncores = ncores,
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
    ncores = ncores,
    prefix = prefix,
    tool_name = tool_name,
    store_results = store_results,
    verbose = verbose,
    ...
  )
}

#' @title SpatialEcoTyper spatial plot
#'
#' @description
#' Plot SpatialEcoTyper labels on spatial coordinates using scop's spatial
#' plotting style.
#'
#' @md
#' @inheritParams SpatialSpotPlot
#' @param group.by Metadata column containing SpatialEcoTyper labels.
#' @param x.by,y.by Metadata coordinate columns used when no image coordinates
#' are available.
#' @param palette,palcolor Palette passed to `palette_colors()`.
#' @param ... Additional arguments passed to [SpatialSpotPlot()].
#'
#' @return A `ggplot`, `patchwork`, or list of `ggplot` objects.
#'
#' @examples
#' counts <- matrix(
#'   c(3, 0, 1, 2, 0, 4, 1, 0, 2, 1, 3, 0),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#' rownames(counts) <- c("EPCAM", "COL1A1", "PTPRC")
#' colnames(counts) <- paste0("spot", 1:4)
#' srt <- Seurat::CreateSeuratObject(counts)
#' srt$X <- c(0, 1, 0, 1)
#' srt$Y <- c(0, 0, 1, 1)
#' srt$SpatialEcoTyper_SE <- c("SE1", "SE1", "SE2", "SE2")
#' srt$CellType <- c("Epithelial", "Fibroblast", "Immune", "Epithelial")
#'
#' SpatialEcoTyperSpatialPlot(
#'   srt,
#'   overlay_image = FALSE,
#'   coord.cols = c("X", "Y"),
#'   pt.size = 4
#' )
#' @export
SpatialEcoTyperSpatialPlot <- function(
  srt,
  group.by = "SpatialEcoTyper_SE",
  x.by = "X",
  y.by = "Y",
  image = NULL,
  overlay_image = TRUE,
  coord.cols = c(x.by, y.by),
  palette = "Spectral",
  palcolor = NULL,
  theme_use = "theme_scop",
  ...
) {
  SpatialSpotPlot(
    srt = srt,
    group.by = group.by,
    image = image,
    overlay_image = overlay_image,
    coord.cols = coord.cols,
    palette = palette,
    palcolor = palcolor,
    theme_use = theme_use,
    ...
  )
}

#' @title SpatialEcoTyper composition plot
#'
#' @description
#' Plot SpatialEcoTyper composition across cell types, samples, or other
#' metadata groups using scop's default plotting theme and palettes.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object.
#' @param se.by Metadata column containing SpatialEcoTyper labels.
#' @param group.by Metadata column used as the composition group.
#' @param sample.by Optional metadata column used for faceting.
#' @param position Bar position. `"fill"` shows fractions and `"stack"` shows
#' counts.
#' @param palette,palcolor Palette passed to `palette_colors()`.
#' @param legend.position Legend position.
#' @param theme_use Theme function name. Default is `"theme_scop"`.
#' @param theme_args Additional arguments passed to the theme function.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' counts <- matrix(
#'   c(3, 0, 1, 2, 0, 4, 1, 0, 2, 1, 3, 0),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#' rownames(counts) <- c("EPCAM", "COL1A1", "PTPRC")
#' colnames(counts) <- paste0("spot", 1:4)
#' srt <- Seurat::CreateSeuratObject(counts)
#' srt$SpatialEcoTyper_SE <- c("SE1", "SE1", "SE2", "SE2")
#' srt$CellType <- c("Epithelial", "Fibroblast", "Immune", "Epithelial")
#' srt$sample <- c("slice1", "slice1", "slice2", "slice2")
#'
#' SpatialEcoTyperCompositionPlot(
#'   srt,
#'   group.by = "CellType",
#'   sample.by = "sample",
#'   position = "fill"
#' )
#' @export
SpatialEcoTyperCompositionPlot <- function(
  srt,
  se.by = "SpatialEcoTyper_SE",
  group.by = "CellType",
  sample.by = NULL,
  position = c("fill", "stack"),
  palette = "Spectral",
  palcolor = NULL,
  legend.position = "right",
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  position <- match.arg(position)
  meta <- srt[[]]
  needed_cols <- c(se.by, group.by, sample.by)
  missing_cols <- setdiff(needed_cols, colnames(meta))
  if (length(missing_cols) > 0L) {
    log_message(
      "Missing metadata column{?s}: {.val {missing_cols}}",
      message_type = "error"
    )
  }
  dat <- data.frame(
    SE = meta[[se.by]],
    Group = meta[[group.by]],
    Sample = if (is.null(sample.by)) "All" else meta[[sample.by]],
    stringsAsFactors = FALSE
  )
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  if (nrow(dat) == 0L) {
    log_message("No complete SpatialEcoTyper composition records to plot", message_type = "error")
  }
  dat$SE <- factor(dat$SE, levels = unique(as.character(dat$SE)))
  dat$Group <- factor(dat$Group, levels = unique(as.character(dat$Group)))
  dat$Sample <- factor(dat$Sample, levels = unique(as.character(dat$Sample)))

  plot_dat <- as.data.frame(
    stats::xtabs(~ SE + Group + Sample, data = dat),
    stringsAsFactors = FALSE
  )
  colnames(plot_dat) <- c("SE", "Group", "Sample", "n")
  plot_dat <- plot_dat[plot_dat$n > 0, , drop = FALSE]
  fill_cols <- palette_colors(
    x = levels(dat$Group),
    palette = palette,
    palcolor = palcolor
  )
  p <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(x = .data$SE, y = .data$n, fill = .data$Group)
  ) +
    ggplot2::geom_col(width = 0.78, position = position) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = fill_cols, drop = FALSE) +
    ggplot2::labs(
      x = se.by,
      y = if (identical(position, "fill")) "Fraction" else "Cell count",
      fill = group.by
    ) +
    spatialecotyper_plot_theme(
      theme_use = theme_use,
      theme_args = theme_args
    ) +
    ggplot2::theme(legend.position = legend.position)
  if (identical(position, "fill")) {
    p <- p + ggplot2::scale_y_continuous(labels = scales::percent)
  }
  if (!is.null(sample.by) && length(levels(dat$Sample)) > 1L) {
    p <- p + ggplot2::facet_wrap(~ Sample, scales = "free_y")
  }
  p
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
  ncores,
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
  ...
) {
  metadata <- spatialecotyper_seurat_metadata(
    srt = srt,
    cells = colnames(normdata),
    celltype.by = celltype.by,
    x.by = x.by,
    y.by = y.by
  )
  log_message(
    "Run {.pkg SpatialEcoTyper} single-sample discovery on {.val {ncol(normdata)}} cells and {.val {nrow(normdata)}} features",
    verbose = verbose
  )
  spatialecotyper_fun <- get_namespace_fun(
    "SpatialEcoTyper",
    "SpatialEcoTyper"
  )
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
    ncores = ncores,
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
    cells = colnames(srt),
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
      ncores = ncores,
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
  ncores,
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
  ...
) {
  spatialecotyper_assert_scalar_string(sample.by, "sample.by")
  if (!is.null(Region)) {
    spatialecotyper_assert_scalar_string(Region, "Region")
  }
  meta <- srt[[]]
  spatialecotyper_check_meta_columns(
    meta = meta,
    cols = c(celltype.by, sample.by, x.by, y.by, Region)
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
    metadata_list[[sample_name]] <- spatialecotyper_build_metadata(
      meta = meta[sample_cells, , drop = FALSE],
      celltype.by = celltype.by,
      x.by = x.by,
      y.by = y.by
    )
    if (!is.null(Region)) {
      metadata_list[[sample_name]][[Region]] <- meta[sample_cells, Region, drop = TRUE]
    }
  }
  outdir <- outdir %||% spatialecotyper_temp_outdir(prefix)
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
    ncores = ncores,
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
      ncores = ncores,
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
  ncores,
  prefix,
  tool_name,
  store_results,
  allow_partial,
  verbose,
  ...
) {
  if (is.null(celltypes)) {
    spatialecotyper_assert_scalar_string(celltype.by, "celltype.by")
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
    ncores = ncores,
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
      ncores = ncores,
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
  ncores,
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
    ncores = ncores,
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
      ncores = ncores,
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

spatialecotyper_seurat_metadata <- function(
  srt,
  cells,
  celltype.by,
  x.by,
  y.by
) {
  spatialecotyper_assert_scalar_string(celltype.by, "celltype.by")
  spatialecotyper_assert_scalar_string(x.by, "x.by")
  spatialecotyper_assert_scalar_string(y.by, "y.by")
  meta <- srt[[]]
  spatialecotyper_check_meta_columns(meta = meta, cols = c(celltype.by, x.by, y.by))
  missing_cells <- setdiff(cells, rownames(meta))
  if (length(missing_cells) > 0L) {
    log_message(
      "Metadata is missing {.val {length(missing_cells)}} cell{?s} from the selected assay/layer",
      message_type = "error"
    )
  }
  spatialecotyper_build_metadata(
    meta = meta[cells, , drop = FALSE],
    celltype.by = celltype.by,
    x.by = x.by,
    y.by = y.by
  )
}

spatialecotyper_assert_scalar_string <- function(x, arg) {
  if (
    is.null(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x)
  ) {
    log_message(
      "{.arg {arg}} must be a single non-empty string",
      message_type = "error"
    )
  }
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
  invalid <- is.na(metadata$X) |
    is.na(metadata$Y) |
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
  parameters
) {
  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      result = result,
      metadata = metadata,
      parameters = parameters
    )
    srt@tools[[tool_name]] <- spatial_result_build(
      bundle = srt@tools[[tool_name]],
      method = "SpatialEcoTyper",
      result_type = "ecotype",
      provenance = list(producer = "RunSpatialEcoTyper", backend_id = "spatialecotyper")
    )
  }
  srt
}

spatialecotyper_temp_outdir <- function(prefix) {
  file.path(
    tempdir(),
    paste0(prefix, "_", format(Sys.time(), "%Y%m%d%H%M%S"))
  )
}

spatialecotyper_plot_theme <- function(theme_use = "theme_scop", theme_args = list()) {
  if (identical(theme_use, "theme_scop")) {
    return(do.call(theme_scop, theme_args))
  }
  do.call(theme_use, theme_args)
}
