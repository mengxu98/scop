#' @title Run multi-sample spatial integration
#'
#' @description
#' Integrate multi-slice or multi-sample spatial transcriptomics data with an
#' optional spatial backend and store standardized embeddings, domains, and
#' aligned coordinates in a `Seurat` object.
#'
#' @md
#' @inheritParams RunSpatialVariableFeatures
#' @inheritParams SpatialSpotPlot
#' @param object A merged spatial `Seurat` object or a list of spatial `Seurat`
#' objects.
#' @param method Spatial integration backend.
#' @param sample.by Metadata column identifying samples for a merged `Seurat`
#' object. For list input, list names are copied into this column.
#' @param reduction.name Name of the integrated embedding reduction. If `NULL`,
#' a method-specific name is used.
#' @param cluster_colname Metadata column used for spatial domain labels. If
#' `NULL`, a method-specific name is used.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param ... Additional backend-specific arguments.
#'
#' @return A `Seurat` object with spatial integration results stored in
#' metadata, reductions, and `srt@tools[[tool_name]]`.
#' @export
#'
#' @examples
#' \dontrun{
#' srt <- RunSpatialIntegration(
#'   object = spatial,
#'   method = "PRECAST",
#'   sample.by = "sample",
#'   assay = "Spatial",
#'   coord.cols = c("col", "row")
#' )
#'
#' SpatialIntegrationPlot(srt, plot_type = "spatial")
#' SpatialIntegrationPlot(srt, plot_type = "embedding")
#' }
RunSpatialIntegration <- function(
  object,
  method = c("PRECAST", "BASS", "SpatialMNN"),
  sample.by = NULL,
  assay = NULL,
  layer = "counts",
  coord.cols = c("col", "row"),
  features = NULL,
  image = NULL,
  reduction.name = NULL,
  cluster_colname = NULL,
  tool_name = "SpatialIntegration",
  store_results = TRUE,
  verbose = TRUE,
  ...
) {
  method <- match.arg(method)
  spatial_integration_assert_scalar_string(tool_name, "tool_name")
  sample.by <- spatial_integration_resolve_sample_by(object, sample.by)
  reduction.name <- reduction.name %||% paste0("SpatialIntegration_", method)
  cluster_colname <- cluster_colname %||% paste0("SpatialIntegration_", method, "_domain")
  spatial_integration_assert_scalar_string(reduction.name, "reduction.name")
  spatial_integration_assert_scalar_string(cluster_colname, "cluster_colname")

  input <- spatial_integration_prepare_input(
    object = object,
    sample.by = sample.by,
    assay = assay,
    layer = layer,
    features = features,
    image = image,
    coord.cols = coord.cols
  )

  log_message(
    "Run {.pkg {method}} spatial integration with {.val {length(input$samples)}} sample{?s} and {.val {length(input$cells)}} spot{?s}",
    verbose = verbose
  )
  backend <- spatial_integration_run_backend(
    method = method,
    input = input,
    verbose = verbose,
    ...
  )
  result <- spatial_integration_standardize_result(
    backend = backend,
    input = input,
    method = method
  )
  srt <- spatial_integration_apply_result(
    srt = input$srt,
    result = result,
    method = method,
    sample.by = sample.by,
    assay = input$assay,
    layer = layer,
    image = image,
    coord.cols = coord.cols,
    reduction.name = reduction.name,
    cluster_colname = cluster_colname,
    tool_name = tool_name,
    store_results = store_results
  )
  log_message(
    "{.pkg {method}} spatial integration results stored in {.code srt@tools[[{tool_name}]]}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

#' @title Plot spatial integration results
#'
#' @description
#' Visualize standardized results produced by [RunSpatialIntegration()].
#'
#' @md
#' @inheritParams SpatialSpotPlot
#' @param srt A `Seurat` object containing spatial integration results.
#' @param method Stored integration method. If `NULL`, the active method stored
#' in `srt@tools[[tool_name]]` is used.
#' @param plot_type Plot type: `"spatial"`, `"embedding"`, `"alignment"`, or
#' `"composition"`.
#' @param group.by Metadata column used for coloring. Defaults to the stored
#' spatial domain column.
#' @param sample.by Metadata column used for facets or composition grouping.
#' Defaults to the stored sample column.
#' @param reduction Reduction used for embedding plots. Defaults to the stored
#' integration reduction.
#' @param use_aligned Whether spatial plots should use aligned coordinates when
#' available.
#' @param tool_name Name of the `srt@tools` entry created by
#' [RunSpatialIntegration()].
#' @param combine Whether to combine plots when delegated plotting returns a
#' list.
#' @param ... Additional arguments passed to [SpatialSpotPlot()] or
#' [CellDimPlot()].
#'
#' @return A `ggplot`, patchwork object, or list of plots.
#' @export
SpatialIntegrationPlot <- function(
  srt,
  method = NULL,
  plot_type = c("spatial", "embedding", "alignment", "composition"),
  group.by = NULL,
  sample.by = NULL,
  reduction = NULL,
  cluster_colname = NULL,
  coord.cols = c("col", "row"),
  use_aligned = FALSE,
  tool_name = "SpatialIntegration",
  combine = TRUE,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  plot_type <- match.arg(plot_type)
  bundle <- spatial_integration_get_bundle(srt, tool_name = tool_name)
  method <- method %||% bundle$active_method
  method_bundle <- spatial_integration_get_method_bundle(bundle, method)
  parameters <- method_bundle$parameters %||% bundle$parameters %||% list()
  group.by <- group.by %||% cluster_colname %||% parameters$cluster_colname
  sample.by <- sample.by %||% parameters$sample.by
  reduction <- reduction %||% parameters$reduction.name
  if (is.null(group.by) || !group.by %in% colnames(srt@meta.data)) {
    log_message(
      "A valid {.arg group.by} metadata column is required for spatial integration plotting",
      message_type = "error"
    )
  }

  if (identical(plot_type, "spatial")) {
    coord_use <- coord.cols
    if (isTRUE(use_aligned)) {
      coord_use <- parameters$aligned_coord_cols
      if (is.null(coord_use) || !all(coord_use %in% colnames(srt@meta.data))) {
        log_message(
          "Aligned coordinates are not available for {.arg use_aligned = TRUE}",
          message_type = "error"
        )
      }
    }
    return(SpatialSpotPlot(
      srt = srt,
      group.by = group.by,
      split.by = sample.by,
      coord.cols = coord_use,
      combine = combine,
      palette = palette,
      palcolor = palcolor,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    ))
  }

  if (identical(plot_type, "embedding")) {
    if (is.null(reduction) || !reduction %in% SeuratObject::Reductions(srt)) {
      log_message(
        "A valid integrated {.arg reduction} is required for embedding plots",
        message_type = "error"
      )
    }
    return(CellDimPlot(
      srt = srt,
      group.by = group.by,
      reduction = reduction,
      split.by = sample.by,
      combine = combine,
      palette = palette,
      palcolor = palcolor,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    ))
  }

  if (identical(plot_type, "composition")) {
    return(spatial_integration_composition_plot(
      srt = srt,
      group.by = group.by,
      sample.by = sample.by,
      palette = palette,
      palcolor = palcolor,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }

  spatial_integration_alignment_plot(
    srt = srt,
    group.by = group.by,
    sample.by = sample.by,
    parameters = parameters,
    palette = palette,
    palcolor = palcolor,
    theme_use = theme_use,
    theme_args = theme_args
  )
}

spatial_integration_prepare_input <- function(
  object,
  sample.by,
  assay,
  layer,
  features,
  image,
  coord.cols
) {
  if (inherits(object, "Seurat")) {
    if (!sample.by %in% colnames(object@meta.data)) {
      log_message(
        "{.arg sample.by} {.val {sample.by}} is not present in {.arg object}",
        message_type = "error"
      )
    }
    samples <- as.character(object@meta.data[[sample.by]])
    if (any(is.na(samples) | !nzchar(samples))) {
      log_message(
        "{.arg sample.by} contains missing sample labels",
        message_type = "error"
      )
    }
    srt <- object
    srt_list <- Seurat::SplitObject(object, split.by = sample.by)
    names(srt_list) <- names(srt_list) %||% unique(samples)
  } else {
    srt_list <- spatial_integration_as_list(object, sample.by = sample.by)
    srt <- spatial_integration_merge_list(srt_list, sample.by = sample.by)
  }
  samples <- names(srt_list)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  features_use <- spatial_integration_features(
    features = features,
    expr = expr,
    srt_list = srt_list,
    assay = assay
  )
  coords <- spatial_dim_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    overlay_image = FALSE
  )$data
  cells <- intersect(colnames(srt), colnames(expr))
  cells <- intersect(cells, rownames(coords))
  coords <- coords[cells, , drop = FALSE]
  keep_coords <- is.finite(coords$x) & is.finite(coords$y)
  cells <- cells[keep_coords]
  coords <- coords[cells, , drop = FALSE]
  if (length(cells) < 3L) {
    log_message(
      "At least three spatial spots with finite coordinates are required",
      message_type = "error"
    )
  }
  expr <- expr[features_use, cells, drop = FALSE]
  keep_features <- Matrix::rowSums(expr != 0) > 0
  keep_cells <- Matrix::colSums(expr != 0) > 0
  expr <- expr[keep_features, keep_cells, drop = FALSE]
  coords <- coords[colnames(expr), , drop = FALSE]
  cells <- colnames(expr)
  if (nrow(expr) == 0L || ncol(expr) == 0L) {
    log_message(
      "No non-zero features or spots remain for spatial integration",
      message_type = "error"
    )
  }
  expr <- spatial_integration_sparse_matrix(expr)
  split_cells <- split(cells, srt@meta.data[cells, sample.by, drop = TRUE])
  list(
    srt = srt,
    srt_list = spatial_integration_split_merged(
      srt = srt,
      cells_by_sample = split_cells
    ),
    expr = expr,
    expr_list = lapply(split_cells, function(x) expr[, x, drop = FALSE]),
    coords = coords,
    coords_list = lapply(split_cells, function(x) coords[x, , drop = FALSE]),
    samples = names(split_cells),
    cells = cells,
    features = rownames(expr),
    assay = assay,
    sample.by = sample.by
  )
}

spatial_integration_as_list <- function(object, sample.by) {
  if (inherits(object, "Seurat")) {
    if (!sample.by %in% colnames(object@meta.data)) {
      log_message(
        "{.arg sample.by} {.val {sample.by}} is not present in {.arg object}",
        message_type = "error"
      )
    }
    samples <- as.character(object@meta.data[[sample.by]])
    if (any(is.na(samples) | !nzchar(samples))) {
      log_message(
        "{.arg sample.by} contains missing sample labels",
        message_type = "error"
      )
    }
    out <- Seurat::SplitObject(object, split.by = sample.by)
    names(out) <- names(out) %||% unique(samples)
    return(out)
  }
  if (!is.list(object) || length(object) == 0L) {
    log_message(
      "{.arg object} must be a {.cls Seurat} object or a non-empty list of {.cls Seurat} objects",
      message_type = "error"
    )
  }
  is_seurat <- vapply(object, inherits, logical(1), what = "Seurat")
  if (!all(is_seurat)) {
    log_message(
      "Every element of {.arg object} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  nms <- names(object)
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    nms <- paste0("sample", seq_along(object))
  }
  nms <- make.unique(as.character(nms), sep = "_")
  names(object) <- nms
  for (nm in nms) {
    object[[nm]]@meta.data[[sample.by]] <- nm
  }
  object
}

spatial_integration_merge_list <- function(srt_list, sample.by) {
  if (length(srt_list) == 1L) {
    srt <- srt_list[[1L]]
    srt@meta.data[[sample.by]] <- names(srt_list)[1L]
    return(srt)
  }
  merge(
    x = srt_list[[1L]],
    y = srt_list[-1L],
    add.cell.ids = names(srt_list),
    merge.data = TRUE
  )
}

spatial_integration_split_merged <- function(srt, cells_by_sample) {
  lapply(cells_by_sample, function(cells) {
    srt[, cells]
  })
}

spatial_integration_features <- function(features, expr, srt_list, assay) {
  common <- Reduce(
    intersect,
    lapply(srt_list, function(x) rownames(x[[assay]]))
  )
  if (!is.null(features)) {
    common <- intersect(unique(features), common)
  }
  common <- intersect(common, rownames(expr))
  if (length(common) == 0L) {
    log_message(
      "No requested {.arg features} are shared across spatial samples",
      message_type = "error"
    )
  }
  common
}

spatial_integration_sparse_matrix <- function(mat) {
  if (!inherits(mat, "Matrix")) {
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
  }
  if (!inherits(mat, "dgCMatrix")) {
    mat <- methods::as(mat, "dgCMatrix")
  }
  mat@x[!is.finite(mat@x)] <- 0
  Matrix::drop0(mat)
}

spatial_integration_run_backend <- function(method, input, verbose = TRUE, ...) {
  params <- list(...)
  switch(method,
    PRECAST = spatial_integration_run_precast(input, params, verbose = verbose),
    BASS = spatial_integration_run_bass(input, params, verbose = verbose),
    SpatialMNN = spatial_integration_run_spatialmnn(input, params, verbose = verbose)
  )
}

spatial_integration_run_precast <- function(input, params, verbose = TRUE) {
  check_r("PRECAST", verbose = FALSE)
  create_fun <- get_namespace_fun("PRECAST", "CreatePRECASTObject")
  adj_fun <- get_namespace_fun("PRECAST", "AddAdjList")
  par_fun <- get_namespace_fun("PRECAST", "AddParSetting")
  run_fun <- get_namespace_fun("PRECAST", "PRECAST")
  select_fun <- get_namespace_fun("PRECAST", "SelectModel")
  obj <- spatial_integration_call(
    create_fun,
    c(
      list(seuList = input$srt_list, project = "scop_spatial_integration"),
      params$create_params %||% list()
    )
  )
  obj <- spatial_integration_call(
    adj_fun,
    c(list(PRECASTObj = obj), params$adj_params %||% list())
  )
  obj <- spatial_integration_call(
    par_fun,
    c(list(PRECASTObj = obj), params$par_params %||% list())
  )
  obj <- spatial_integration_call(
    run_fun,
    c(list(PRECASTObj = obj), params$run_params %||% list())
  )
  obj <- spatial_integration_call(
    select_fun,
    c(list(PRECASTObj = obj), params$select_params %||% list())
  )
  spatial_integration_extract_backend(
    raw_result = obj,
    method = "PRECAST",
    input = input
  )
}

spatial_integration_run_bass <- function(input, params, verbose = TRUE) {
  check_r("BASS", verbose = FALSE)
  create_fun <- get_namespace_fun("BASS", "createBASSObject")
  preprocess_fun <- get_namespace_fun("BASS", "BASS.preprocess")
  run_fun <- get_namespace_fun("BASS", "BASS.run")
  post_fun <- tryCatch(get_namespace_fun("BASS", "BASS.postprocess"), error = function(e) NULL)
  bass <- spatial_integration_call(
    create_fun,
    c(
      list(
        X = input$expr_list,
        xy = input$coords_list,
        C = params$C %||% params$n_domains %||% 7,
        R = params$R %||% length(input$samples)
      ),
      params$create_params %||% list()
    )
  )
  bass <- spatial_integration_call(
    preprocess_fun,
    c(list(BASS = bass), params$preprocess_params %||% list())
  )
  bass <- spatial_integration_call(run_fun, c(list(BASS = bass), params$run_params %||% list()))
  if (is.function(post_fun)) {
    bass <- spatial_integration_call(post_fun, c(list(BASS = bass), params$postprocess_params %||% list()))
  }
  spatial_integration_extract_backend(
    raw_result = bass,
    method = "BASS",
    input = input
  )
}

spatial_integration_run_spatialmnn <- function(input, params, verbose = TRUE) {
  pkg <- if (requireNamespace("SpatialMNN", quietly = TRUE)) "SpatialMNN" else "spatialMNN"
  check_r(pkg, verbose = FALSE)
  run_fun <- get_namespace_fun(pkg, "spatialMNN")
  res <- spatial_integration_call(
    run_fun,
    c(list(seu_ls = input$srt_list), params$run_params %||% params)
  )
  spatial_integration_extract_backend(
    raw_result = res,
    method = "SpatialMNN",
    input = input
  )
}

spatial_integration_extract_backend <- function(raw_result, method, input) {
  if (is.list(raw_result) && any(c("embedding", "domains", "aligned_coords") %in% names(raw_result))) {
    raw_result$raw_result <- raw_result$raw_result %||% raw_result
    return(raw_result)
  }
  if (identical(method, "BASS")) {
    domains <- tryCatch(raw_result@results$z, error = function(e) NULL)
    if (is.matrix(domains)) {
      domains <- domains[, 1L]
    }
    if (!is.null(domains)) {
      domains <- as.character(domains)
      names(domains) <- input$cells[seq_along(domains)]
    }
    return(list(domains = domains, raw_result = raw_result))
  }
  if (identical(method, "SpatialMNN") && is.list(raw_result)) {
    merged <- spatial_integration_extract_spatialmnn_list(raw_result, input)
    merged$raw_result <- raw_result
    return(merged)
  }
  list(raw_result = raw_result)
}

spatial_integration_extract_spatialmnn_list <- function(raw_result, input) {
  cells <- character()
  domains <- character()
  for (obj in raw_result) {
    if (!inherits(obj, "Seurat")) {
      next
    }
    domain_col <- grep("cluster|domain|niche|mnn", colnames(obj@meta.data), ignore.case = TRUE, value = TRUE)[1L]
    if (is.na(domain_col)) {
      next
    }
    cells <- c(cells, colnames(obj))
    domains <- c(domains, as.character(obj@meta.data[[domain_col]]))
  }
  if (length(cells) > 0L) {
    names(domains) <- cells
  } else {
    domains <- NULL
  }
  list(domains = domains)
}

spatial_integration_standardize_result <- function(backend, input, method) {
  embedding <- spatial_integration_standardize_embedding(backend$embedding, input$cells)
  domains <- spatial_integration_standardize_named_vector(backend$domains, input$cells)
  aligned_coords <- spatial_integration_standardize_coords(
    coords = backend$aligned_coords,
    cells = input$cells
  )
  if (is.null(embedding) && is.null(domains) && is.null(aligned_coords)) {
    log_message(
      "{.pkg {method}} did not return embeddings, spatial domains, or aligned coordinates that can be matched to Seurat cells",
      message_type = "error"
    )
  }
  list(
    embedding = embedding,
    domains = domains,
    aligned_coords = aligned_coords,
    features = input$features,
    raw_result = backend$raw_result %||% backend
  )
}

spatial_integration_standardize_embedding <- function(embedding, cells) {
  if (is.null(embedding)) {
    return(NULL)
  }
  embedding <- as.matrix(embedding)
  if (is.null(rownames(embedding))) {
    if (nrow(embedding) != length(cells)) {
      log_message(
        "Backend embedding must have row names or one row per Seurat cell",
        message_type = "error"
      )
    }
    rownames(embedding) <- cells
  }
  cells_use <- cells[cells %in% rownames(embedding)]
  if (length(cells_use) == 0L) {
    log_message(
      "Backend embedding rows could not be matched to Seurat cells",
      message_type = "error"
    )
  }
  embedding <- embedding[cells_use, , drop = FALSE]
  storage.mode(embedding) <- "double"
  if (is.null(colnames(embedding))) {
    colnames(embedding) <- paste0("dim", seq_len(ncol(embedding)))
  }
  embedding
}

spatial_integration_standardize_named_vector <- function(x, cells) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x) || is.matrix(x)) {
    x <- x[, 1L, drop = TRUE]
  }
  x <- as.character(x)
  if (is.null(names(x))) {
    if (length(x) != length(cells)) {
      log_message(
        "Backend domain labels must be named or have one value per Seurat cell",
        message_type = "error"
      )
    }
    names(x) <- cells
  }
  cells_use <- cells[cells %in% names(x)]
  if (length(cells_use) == 0L) {
    log_message(
      "Backend domain labels could not be matched to Seurat cells",
      message_type = "error"
    )
  }
  x[cells_use]
}

spatial_integration_standardize_coords <- function(coords, cells) {
  if (is.null(coords)) {
    return(NULL)
  }
  coords <- as.data.frame(coords, check.names = FALSE)
  if (ncol(coords) < 2L) {
    log_message(
      "Backend aligned coordinates must contain at least two columns",
      message_type = "error"
    )
  }
  if (is.null(rownames(coords))) {
    if (nrow(coords) != length(cells)) {
      log_message(
        "Backend aligned coordinates must have row names or one row per Seurat cell",
        message_type = "error"
      )
    }
    rownames(coords) <- cells
  }
  cells_use <- cells[cells %in% rownames(coords)]
  if (length(cells_use) == 0L) {
    log_message(
      "Backend aligned coordinates could not be matched to Seurat cells",
      message_type = "error"
    )
  }
  out <- coords[cells_use, seq_len(2L), drop = FALSE]
  colnames(out) <- c("x", "y")
  out$x <- as.numeric(out$x)
  out$y <- as.numeric(out$y)
  out
}

spatial_integration_apply_result <- function(
  srt,
  result,
  method,
  sample.by,
  assay,
  layer,
  image,
  coord.cols,
  reduction.name,
  cluster_colname,
  tool_name,
  store_results
) {
  all_cells <- colnames(srt)
  if (!is.null(result$domains)) {
    domain_full <- rep(NA_character_, length(all_cells))
    names(domain_full) <- all_cells
    domain_full[names(result$domains)] <- as.character(result$domains)
    srt[[cluster_colname]] <- domain_full
  }

  aligned_coord_cols <- NULL
  if (!is.null(result$aligned_coords)) {
    aligned_coord_cols <- paste0(reduction.name, c("_aligned_x", "_aligned_y"))
    x_full <- rep(NA_real_, length(all_cells))
    y_full <- rep(NA_real_, length(all_cells))
    names(x_full) <- names(y_full) <- all_cells
    x_full[rownames(result$aligned_coords)] <- result$aligned_coords$x
    y_full[rownames(result$aligned_coords)] <- result$aligned_coords$y
    srt[[aligned_coord_cols[1L]]] <- x_full
    srt[[aligned_coord_cols[2L]]] <- y_full
  }

  if (!is.null(result$embedding)) {
    key <- spatial_integration_reduction_key(reduction.name)
    embedding_full <- matrix(
      NA_real_,
      nrow = length(all_cells),
      ncol = ncol(result$embedding),
      dimnames = list(all_cells, paste0(key, seq_len(ncol(result$embedding))))
    )
    embedding_full[rownames(result$embedding), ] <- result$embedding
    srt[[reduction.name]] <- SeuratObject::CreateDimReducObject(
      embeddings = embedding_full,
      key = key,
      assay = assay
    )
  }

  parameters <- list(
    method = method,
    sample.by = sample.by,
    assay = assay,
    layer = layer,
    image = image,
    coord.cols = coord.cols,
    reduction.name = reduction.name,
    cluster_colname = cluster_colname,
    aligned_coord_cols = aligned_coord_cols,
    tool_name = tool_name
  )
  method_bundle <- list(
    embedding = result$embedding,
    domains = result$domains,
    aligned_coords = result$aligned_coords,
    raw_result = result$raw_result,
    parameters = parameters
  )
  if (isTRUE(store_results)) {
    old <- srt@tools[[tool_name]] %||% list()
    methods_store <- old$methods %||% list()
    methods_store[[method]] <- method_bundle
    srt@tools[[tool_name]] <- list(
      active_method = method,
      methods = methods_store,
      parameters = parameters,
    features = result$features,
      samples = unique(as.character(srt@meta.data[[sample.by]])),
      cells = all_cells
    )
  }
  srt
}

spatial_integration_get_bundle <- function(srt, tool_name) {
  bundle <- srt@tools[[tool_name]]
  if (is.null(bundle)) {
    log_message(
      "Cannot find spatial integration results in {.code srt@tools[[{tool_name}]]}",
      message_type = "error"
    )
  }
  bundle
}

spatial_integration_get_method_bundle <- function(bundle, method) {
  method_bundle <- bundle$methods[[method]]
  if (is.null(method_bundle)) {
    log_message(
      "Cannot find spatial integration method {.val {method}} in stored results",
      message_type = "error"
    )
  }
  method_bundle
}

spatial_integration_composition_plot <- function(
  srt,
  group.by,
  sample.by,
  palette,
  palcolor,
  theme_use,
  theme_args
) {
  if (is.null(sample.by) || !sample.by %in% colnames(srt@meta.data)) {
    log_message(
      "A valid {.arg sample.by} metadata column is required for composition plots",
      message_type = "error"
    )
  }
  dat <- srt@meta.data[, c(sample.by, group.by), drop = FALSE]
  dat <- dat[!is.na(dat[[sample.by]]) & !is.na(dat[[group.by]]), , drop = FALSE]
  colnames(dat) <- c("sample", "domain")
  if (nrow(dat) == 0L) {
    log_message(
      "No non-missing sample/domain labels are available for composition plotting",
      message_type = "error"
    )
  }
  tab <- as.data.frame(table(dat$sample, dat$domain), stringsAsFactors = FALSE)
  colnames(tab) <- c("sample", "domain", "count")
  tab <- tab[tab$count > 0, , drop = FALSE]
  tab$fraction <- ave(tab$count, tab$sample, FUN = function(x) x / sum(x))
  cols <- palette_colors(unique(as.character(tab$domain)), palette = palette, palcolor = palcolor)
  ggplot2::ggplot(tab, ggplot2::aes(x = .data$sample, y = .data$fraction, fill = .data$domain)) +
    ggplot2::geom_col(width = 0.75, color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = sample.by, y = "Fraction", fill = group.by) +
    spatial_integration_plot_theme(theme_use, theme_args)
}

spatial_integration_alignment_plot <- function(
  srt,
  group.by,
  sample.by,
  parameters,
  palette,
  palcolor,
  theme_use,
  theme_args
) {
  aligned_cols <- parameters$aligned_coord_cols
  raw_cols <- parameters$coord.cols
  if (
    is.null(aligned_cols) ||
      !all(aligned_cols %in% colnames(srt@meta.data)) ||
      is.null(raw_cols) ||
      !all(raw_cols %in% colnames(srt@meta.data))
  ) {
    log_message(
      "Alignment plots require stored raw and aligned coordinate columns",
      message_type = "error"
    )
  }
  meta <- srt@meta.data
  dat_raw <- data.frame(
    x = meta[[raw_cols[1L]]],
    y = meta[[raw_cols[2L]]],
    group = meta[[group.by]],
    sample = if (!is.null(sample.by) && sample.by %in% colnames(meta)) meta[[sample.by]] else "All",
    coordinate = "Raw",
    stringsAsFactors = FALSE
  )
  dat_aligned <- data.frame(
    x = meta[[aligned_cols[1L]]],
    y = meta[[aligned_cols[2L]]],
    group = meta[[group.by]],
    sample = if (!is.null(sample.by) && sample.by %in% colnames(meta)) meta[[sample.by]] else "All",
    coordinate = "Aligned",
    stringsAsFactors = FALSE
  )
  dat <- rbind(dat_raw, dat_aligned)
  dat <- dat[is.finite(dat$x) & is.finite(dat$y), , drop = FALSE]
  if (nrow(dat) == 0L) {
    log_message(
      "No finite raw/aligned coordinates are available for plotting",
      message_type = "error"
    )
  }
  dat$group <- factor(as.character(dat$group), levels = unique(as.character(dat$group)))
  cols <- palette_colors(levels(dat$group), palette = palette, palcolor = palcolor)
  ggplot2::ggplot(dat, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$group)) +
    ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.1, size = 1.2, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = cols, na.value = "grey80") +
    ggplot2::facet_grid(.data$sample ~ .data$coordinate) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL, fill = group.by) +
    spatial_integration_plot_theme(theme_use, theme_args)
}

spatial_integration_plot_theme <- function(theme_use = "theme_scop", theme_args = list()) {
  if (identical(theme_use, "theme_scop")) {
    return(do.call(theme_scop, theme_args))
  }
  if (is.character(theme_use)) {
    return(do.call(get(theme_use, mode = "function"), theme_args))
  }
  do.call(theme_use, theme_args)
}

spatial_integration_call <- function(fun, args) {
  args <- args[!vapply(args, is.null, logical(1))]
  fmls <- tryCatch(names(formals(fun)), error = function(e) NULL)
  if (!is.null(fmls) && !"..." %in% fmls) {
    args <- args[names(args) %in% fmls]
  }
  do.call(fun, args)
}

spatial_integration_resolve_sample_by <- function(object, sample.by) {
  if (!is.null(sample.by)) {
    spatial_integration_assert_scalar_string(sample.by, "sample.by")
    return(sample.by)
  }
  if (inherits(object, "Seurat")) {
    log_message(
      "{.arg sample.by} is required when {.arg object} is a merged {.cls Seurat} object",
      message_type = "error"
    )
  }
  ".scop_spatial_sample"
}

spatial_integration_reduction_key <- function(reduction.name) {
  key <- gsub("[^A-Za-z0-9]", "", reduction.name)
  if (!nzchar(key)) {
    key <- "SpatialIntegration"
  }
  paste0(key, "_")
}

spatial_integration_assert_scalar_string <- function(x, arg) {
  if (is.null(x) || !is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message(
      "{.arg {arg}} must be a single non-empty string",
      message_type = "error"
    )
  }
  invisible(TRUE)
}
