new_giotto2 <- function(
  giotto,
  source = list(),
  results = list(),
  active = NULL,
  history = list(),
  parameters = list()
) {
  structure(
    list(
      giotto = giotto,
      source = source,
      results = results,
      active = active,
      history = history,
      parameters = parameters
    ),
    class = c("giotto2", "list")
  )
}

#' @export
print.giotto2 <- function(x, ...) {
  cat("<giotto2>\n")
  cat("  cells: ", length(x$source$cells %||% character()), "\n", sep = "")
  cat("  features: ", length(x$source$features %||% character()), "\n", sep = "")
  cat("  results: ", paste(names(x$results), collapse = ", "), "\n", sep = "")
  cat("  active: ", x$active %||% "<none>", "\n", sep = "")
  invisible(x)
}

giotto_validate_scop_object <- function(x, allow_seurat = FALSE) {
  if (inherits(x, "giotto2")) {
    return(invisible(TRUE))
  }
  if (isTRUE(allow_seurat) && inherits(x, "Seurat")) {
    return(invisible(TRUE))
  }
  log_message(
    "{.arg x} must be a {.cls giotto2} object",
    message_type = "error"
  )
}

giotto_validate_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    log_message(
      "{.arg {arg}} must be TRUE or FALSE",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

giotto_append_history <- function(x, step, parameters = list()) {
  entry <- list(
    step = step,
    parameters = parameters,
    time = Sys.time()
  )
  x$history[[length(x$history) + 1L]] <- entry
  x
}

giotto_update_result <- function(x, name, result, active = TRUE) {
  x$results[[name]] <- result
  if (isTRUE(active)) {
    x$active <- name
  }
  x
}

giotto_do_call <- function(name, args) {
  fun <- giotto_get_fun(name)
  assign(name, fun, envir = environment())
  assign(name, fun, envir = parent.frame())
  for (env in sys.frames()) {
    try(assign(name, fun, envir = env), silent = TRUE)
  }
  global_env <- globalenv()
  global_exists <- exists(name, envir = global_env, inherits = FALSE)
  global_old <- if (global_exists) get(name, envir = global_env, inherits = FALSE) else NULL
  assign(name, fun, envir = global_env)
  on.exit({
    if (global_exists) {
      assign(name, global_old, envir = global_env)
    } else if (exists(name, envir = global_env, inherits = FALSE)) {
      rm(list = name, envir = global_env)
    }
  }, add = TRUE)
  fmls <- giotto_formal_names(fun)
  if (!is.null(fmls) && !"..." %in% fmls) {
    dropped <- setdiff(names(args), fmls)
    if (length(dropped) > 0L) {
      log_message(
        "Unsupported Giotto argument(s): {.val {dropped}}",
        message_type = "error"
      )
    }
  }
  do.call(name, args, envir = asNamespace("Giotto"))
}

#' @title Convert Seurat to an internal Giotto workflow object
#'
#' @description
#' Create a `giotto2` object from a Seurat object. The converter is
#' SCT-aware: raw counts remain the default Giotto input, while SCT normalized
#' values are optionally added as an extra Giotto expression layer. The input
#' Seurat object is not modified.
#'
#' @inheritParams RunGiottoCluster
#' @param sct.assay Name of the SCT assay.
#' @param use_sct How to handle SCT data. `"auto"` keeps counts as the main
#' Giotto expression and records SCT availability. `"none"` ignores SCT.
#' `"normalized"` adds SCT normalized values as an additional expression layer.
#' @param use_official Whether to try `Giotto::seuratToGiottoV5()` before
#' falling back to the scop-controlled converter.
#'
#' @return A `giotto2` object.
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' g <- structure(
#'   list(
#'     giotto = list(
#'       umap = cbind(
#'         UMAP_1 = as.numeric(scale(spatial$x)),
#'         UMAP_2 = as.numeric(scale(spatial$y))
#'       )
#'     ),
#'     source = list(
#'       cells = colnames(spatial),
#'       features = rownames(spatial),
#'       coordinates = data.frame(
#'         cell_ID = colnames(spatial),
#'         sdimx = spatial$x,
#'         sdimy = spatial$y
#'       )
#'     ),
#'     results = list(
#'       cluster = list(
#'         table = data.frame(
#'           cluster = paste0("cluster_", (seq_len(ncol(spatial)) - 1) %% 3 + 1),
#'           row.names = colnames(spatial)
#'         )
#'       ),
#'       spatial_network = list(
#'         table = data.frame(
#'           from = colnames(spatial)[1:8],
#'           to = colnames(spatial)[2:9]
#'         )
#'       )
#'     ),
#'     active = "cluster"
#'   ),
#'   class = c("giotto2", "list")
#' )
#'
#' GiottoPlot(g, plot_type = "cluster")
#' GiottoPlot(g, plot_type = "network")
#'
#' g <- SeuratToScopGiotto(
#'   spatial,
#'   assay = "Spatial",
#'   layer = "counts",
#'   coord.cols = c("x", "y"),
#'   verbose = FALSE
#' )
#' @export
SeuratToScopGiotto <- function(
  srt,
  assay = NULL,
  layer = "counts",
  sct.assay = "SCT",
  use_sct = c("auto", "none", "normalized"),
  image = NULL,
  coord.cols = c("x", "y"),
  features = NULL,
  conversion_params = list(),
  use_official = TRUE,
  verbose = TRUE,
  seed = 11
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  use_sct <- match.arg(use_sct)
  giotto_validate_named_list(conversion_params, "conversion_params")
  giotto_require(verbose = verbose)
  giotto_old_options <- giotto_disable_r_only_options()
  on.exit(options(giotto_old_options), add = TRUE)

  assays <- SeuratObject::Assays(srt)
  assay <- giotto_pick_seurat_assay(
    srt = srt,
    assay = assay,
    sct.assay = sct.assay
  )
  if (!assay %in% assays) {
    log_message("{.arg assay} {.val {assay}} is not present in {.cls Seurat}", message_type = "error")
  }

  input <- giotto_prepare_input(
    srt = srt,
    assay = assay,
    layer = layer,
    features = features,
    image = image,
    coord.cols = coord.cols
  )

  official <- NULL
  if (isTRUE(use_official)) {
    official <- giotto_try_seurat_to_giotto_v5(
      srt = srt,
      image = image,
      verbose = verbose
    )
  }
  gobject <- official %||% giotto_create_object(
    input = input,
    layer = layer,
    conversion_params = conversion_params,
    verbose = verbose
  )

  sct_info <- giotto_add_sct_expression(
    gobject = gobject,
    srt = srt,
    cells = input$cells,
    features = input$features,
    sct.assay = sct.assay,
    use_sct = use_sct,
    verbose = verbose
  )
  gobject <- sct_info$gobject

  out <- new_giotto2(
    giotto = gobject,
    source = list(
      type = "Seurat",
      assay = assay,
      layer = layer,
      image = image,
      coord.cols = coord.cols,
      cells = input$cells,
      features = input$features,
      coordinates = input$spatial_locs,
      metadata = input$metadata,
      sct_assay = sct.assay,
      sct_mode = use_sct,
      sct_added_layer = sct_info$layer,
      used_official_converter = !is.null(official)
    ),
    results = list(),
    active = NULL,
    history = list(),
    parameters = list(
      spat_unit = "cell",
      feat_type = "rna",
      expression_values = giotto_expression_values(layer),
      seed = seed
    )
  )
  out <- giotto_append_history(
    out,
    step = "SeuratToScopGiotto",
    parameters = list(
      assay = assay,
      layer = layer,
      image = image,
      coord.cols = coord.cols,
      use_sct = use_sct,
      used_official_converter = !is.null(official)
    )
  )
  giotto_validate_runtime_object(out, verbose = verbose)
  out
}

giotto_pick_seurat_assay <- function(srt, assay = NULL, sct.assay = "SCT") {
  assays <- SeuratObject::Assays(srt)
  if (!is.null(assay)) {
    return(assay)
  }
  default_assay <- SeuratObject::DefaultAssay(srt)
  if (
    identical(default_assay, sct.assay) &&
      sct.assay %in% assays &&
      any(assays != sct.assay)
  ) {
    if ("RNA" %in% assays) {
      return("RNA")
    }
    return(assays[assays != sct.assay][1L])
  }
  default_assay
}

giotto_try_seurat_to_giotto_v5 <- function(srt, image = NULL, verbose = TRUE) {
  fun <- tryCatch(giotto_get_fun("seuratToGiottoV5"), error = function(e) NULL)
  if (!is.function(fun)) {
    return(NULL)
  }
  args <- list(sobject = srt, verbose = FALSE)
  if (!is.null(image)) {
    args$spatial_assay <- image
  }
  tryCatch(
    giotto_suppress_known_warnings(giotto_call(fun, args)),
    error = function(e) {
      log_message(
        "Giotto official Seurat converter failed; using scop fallback converter",
        message_type = "warning",
        verbose = verbose
      )
      NULL
    }
  )
}

giotto_add_sct_expression <- function(
  gobject,
  srt,
  cells,
  features,
  sct.assay = "SCT",
  use_sct = "auto",
  verbose = TRUE
) {
  out <- list(gobject = gobject, layer = NULL)
  if (identical(use_sct, "none") || !sct.assay %in% SeuratObject::Assays(srt)) {
    return(out)
  }
  if (!identical(use_sct, "normalized")) {
    return(out)
  }
  sct_data <- tryCatch(
    GetAssayData5(srt, assay = sct.assay, layer = "data"),
    error = function(e) NULL
  )
  if (is.null(sct_data)) {
    log_message(
      "SCT normalized data were not found; continuing without an SCT Giotto expression layer",
      message_type = "warning",
      verbose = verbose
    )
    return(out)
  }
  common_features <- intersect(features, rownames(sct_data))
  common_cells <- intersect(cells, colnames(sct_data))
  if (length(common_features) == 0L || length(common_cells) == 0L) {
    return(out)
  }
  createExprObj <- giotto_get_fun("createExprObj")
  setExpression <- giotto_get_fun("setExpression")
  expr_obj <- giotto_suppress_known_warnings(
    giotto_call(
      createExprObj,
      list(
        expression_data = sct_data[common_features, common_cells, drop = FALSE],
        name = "sct_normalized",
        spat_unit = "cell",
        feat_type = "rna",
        provenance = "cell"
      )
    )
  )
  out$gobject <- giotto_call(
    setExpression,
    list(gobject = gobject, x = expr_obj, verbose = FALSE)
  )
  out$layer <- "sct_normalized"
  out
}

giotto_validate_runtime_object <- function(x, verbose = TRUE) {
  meta <- tryCatch(
    giotto_get_cell_metadata(x$giotto),
    error = function(e) NULL
  )
  if (is.null(meta) || nrow(meta) == 0L) {
    log_message(
      "Converted Giotto object does not expose cell metadata",
      message_type = "error",
      verbose = verbose
    )
  }
  invisible(TRUE)
}

#' @title Run a Giotto workflow
#'
#' @description
#' Run basic or full Giotto analysis on a `giotto2` object. Seurat input is
#' converted first with [SeuratToScopGiotto()].
#'
#' @param x A `giotto2` or Seurat object.
#' @param steps `"basic"` runs preprocessing, PCA/UMAP, nearest-network
#' clustering, and spatial network construction. `"full"` additionally runs
#' spatial genes, spatial modules, optional cell proximity, and HMRF.
#' @param group.by Metadata column used for cell proximity enrichment.
#' @param return_seurat Whether to return a Seurat object when `x` is Seurat.
#' If `FALSE`, returns the internal `giotto2` workflow object.
#' @param store_results Whether to store the internal Giotto workflow object in
#' `srt@tools[[tool_name]]` when returning Seurat.
#' @param tool_name Name used to store the Giotto workflow object in `srt@tools`.
#' @param verbose Whether to print progress messages.
#' @param seed Random seed for reproducible Giotto calls.
#' @param ... Passed to [SeuratToScopGiotto()] when `x` is Seurat.
#'
#' @return A Seurat object by default for Seurat input, otherwise a `giotto2`
#' workflow object.
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' g <- structure(
#'   list(
#'     source = list(
#'       cells = colnames(spatial),
#'       features = rownames(spatial),
#'       coordinates = data.frame(
#'         cell_ID = colnames(spatial),
#'         sdimx = spatial$x,
#'         sdimy = spatial$y
#'       )
#'     ),
#'     results = list(
#'       cluster = list(
#'         table = data.frame(
#'           cluster = paste0("cluster_", (seq_len(ncol(spatial)) - 1) %% 3 + 1),
#'           row.names = colnames(spatial)
#'         )
#'       )
#'     ),
#'     active = "cluster"
#'   ),
#'   class = c("giotto2", "list")
#' )
#' GiottoPlot(g, plot_type = "cluster")
#'
#' g <- RunGiottoWorkflow(
#'   spatial,
#'   steps = "basic",
#'   assay = "Spatial",
#'   layer = "counts",
#'   coord.cols = c("x", "y"),
#'   return_seurat = FALSE,
#'   verbose = FALSE
#' )
#' @export
RunGiottoWorkflow <- function(
  x,
  steps = c("basic", "full"),
  group.by = NULL,
  return_seurat = inherits(x, "Seurat"),
  store_results = TRUE,
  tool_name = "Giotto",
  verbose = TRUE,
  seed = 11,
  ...
) {
  steps <- match.arg(steps)
  giotto_validate_flag(return_seurat, "return_seurat")
  giotto_validate_flag(store_results, "store_results")
  giotto_validate_scalar_string(tool_name, "tool_name")
  if (inherits(x, "Seurat")) {
    srt <- x
    dots <- list(...)
    x <- do.call(
      SeuratToScopGiotto,
      c(list(srt = srt, verbose = verbose, seed = seed), dots)
    )
    cluster_formals <- names(formals(RunGiottoCluster))
    cluster_args <- dots[names(dots) %in% cluster_formals]
    cluster_result <- do.call(
      RunGiottoCluster,
      utils::modifyList(
        list(
          srt = srt,
          features = x$source$features,
          image = x$source$image,
          coord.cols = x$source$coord.cols,
          store_giotto = FALSE,
          verbose = verbose,
          seed = seed
        ),
        cluster_args
      )
    )
    x$giotto <- cluster_result$giotto
    x <- giotto_update_result(
      x,
      name = "cluster",
      result = list(
        table = data.frame(
          cell = rownames(cluster_result$clusters),
          cluster = cluster_result$clusters$cluster,
          row.names = rownames(cluster_result$clusters),
          stringsAsFactors = FALSE
        ),
        vector = cluster_result$cluster_vector,
        metadata = cluster_result$giotto_metadata,
        parameters = cluster_result$parameters
      )
    )
    x <- giotto_append_history(
      x,
      step = "RunGiottoWorkflow:cluster",
      parameters = cluster_result$parameters
    )
    x$parameters$pca_name <- cluster_result$parameters$preprocess_params$name %||% "pca"
    x <- GiottoReduce(
      x,
      reduction = "umap",
      dims = cluster_result$parameters$dims %||% 1:20,
      verbose = verbose,
      seed = seed
    )
    x <- GiottoSpatialNetwork(x, network_name = "Delaunay_network", verbose = verbose)
    if (identical(steps, "full")) {
      spatial_genes_result <- do.call(
        RunGiottoSpatialGenes,
        utils::modifyList(
          list(
            srt = srt,
            features = x$source$features,
            image = x$source$image,
            coord.cols = x$source$coord.cols,
            store_giotto = FALSE,
            verbose = verbose,
            seed = seed
          ),
          dots[names(dots) %in% names(formals(RunGiottoSpatialGenes))]
        )
      )
      x <- giotto_update_result(
        x,
        name = "spatial_genes",
        result = list(
          table = spatial_genes_result$results,
          top_features = spatial_genes_result$top_features,
          raw = spatial_genes_result$raw_result,
          parameters = spatial_genes_result$parameters
        )
      )
      spatial_modules_result <- tryCatch(
        do.call(
          RunGiottoSpatialModules,
          utils::modifyList(
            list(
              srt = srt,
              features = x$results$spatial_genes$top_features %||% x$source$features,
              k = min(10L, max(2L, length(unique(x$results$spatial_genes$top_features %||% x$source$features)) - 1L)),
              image = x$source$image,
              coord.cols = x$source$coord.cols,
              store_giotto = FALSE,
              verbose = verbose,
              seed = seed
            ),
            dots[names(dots) %in% names(formals(RunGiottoSpatialModules))]
          )
        ),
        error = function(e) {
          log_message(
            "Giotto spatial modules failed during full workflow; continuing with completed results",
            message_type = "warning",
            verbose = verbose
          )
          NULL
        }
      )
      if (!is.null(spatial_modules_result)) {
        x <- giotto_update_result(
          x,
          name = "spatial_modules",
          result = list(
            spatial_cor = spatial_modules_result$spatial_cor,
            modules = spatial_modules_result$modules,
            module_tables = spatial_modules_result$module_tables,
            parameters = spatial_modules_result$parameters
          )
        )
      }
      if (!is.null(group.by)) {
        proximity_result <- tryCatch(
          RunGiottoCellProximity(
            srt = srt,
            group.by = group.by,
            image = x$source$image,
            coord.cols = x$source$coord.cols,
            store_giotto = FALSE,
            verbose = verbose,
            seed = seed
          ),
          error = function(e) {
            log_message(
              "Giotto cell proximity failed during full workflow; continuing with completed results",
              message_type = "warning",
              verbose = verbose
            )
            NULL
          }
        )
        if (!is.null(proximity_result)) {
          x <- giotto_update_result(
            x,
            name = "cell_proximity",
            result = list(
              table = proximity_result$enrichment,
              raw = proximity_result$raw_result,
              parameters = proximity_result$parameters
            )
          )
        }
      }
      x <- tryCatch(
        GiottoHMRF(x, verbose = verbose, seed = seed),
        error = function(e) {
          log_message(
            "Giotto HMRF failed during full workflow; returning completed non-HMRF results",
            message_type = "warning",
            verbose = verbose
          )
          x
        }
      )
    }
    if (isTRUE(return_seurat)) {
      result <- if (!is.null(x$results$hmrf)) "hmrf" else "cluster"
      srt <- AddGiottoToSeurat(
        srt = srt,
        x = x,
        result = result,
        tool_name = tool_name,
        store_result = store_results
      )
      return(srt)
    }
    return(x)
  }
  giotto_validate_scop_object(x)
  x <- GiottoPreprocess(x, verbose = verbose, seed = seed)
  x <- GiottoReduce(x, reduction = "pca", verbose = verbose, seed = seed)
  x <- GiottoReduce(x, reduction = "umap", verbose = verbose, seed = seed)
  x <- GiottoCluster(x, verbose = verbose, seed = seed)
  x <- GiottoSpatialNetwork(x, network_name = "Delaunay_network", verbose = verbose)
  if (identical(steps, "full")) {
    x <- GiottoSpatialGenes(x, verbose = verbose, seed = seed)
    x <- GiottoSpatialModules(x, verbose = verbose, seed = seed)
    if (!is.null(group.by)) {
      x <- GiottoCellProximity(x, group.by = group.by, verbose = verbose, seed = seed)
    }
    x <- GiottoHMRF(x, verbose = verbose, seed = seed)
  }
  x
}

#' @title Preprocess an internal Giotto workflow object
#'
#' @param x A `giotto2` workflow object.
#' @param filter_params Additional parameters reserved for future filtering.
#' @param norm_params Additional parameters passed to `Giotto::normalizeGiotto()`.
#' @param stat_params Additional parameters passed to `Giotto::addStatistics()`.
#' @param hvf_params Additional parameters passed to `Giotto::calculateHVF()`.
#' @param verbose Whether to print progress messages.
#' @param seed Random seed for reproducible Giotto calls.
#'
#' @return A `giotto2` workflow object.
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' g <- structure(
#'   list(
#'     source = list(
#'       cells = colnames(spatial),
#'       features = rownames(spatial)[1:100],
#'       coordinates = data.frame(cell_ID = colnames(spatial), sdimx = spatial$x, sdimy = spatial$y)
#'     ),
#'     results = list(),
#'     active = NULL
#'   ),
#'   class = c("giotto2", "list")
#' )
#' GiottoPlot(g, plot_type = "spatial")
#'
#' g <- SeuratToScopGiotto(spatial, assay = "Spatial", coord.cols = c("x", "y"), verbose = FALSE)
#' g <- GiottoPreprocess(g, verbose = FALSE)
#' @export
GiottoPreprocess <- function(
  x,
  filter_params = list(),
  norm_params = list(),
  stat_params = list(),
  hvf_params = list(),
  verbose = TRUE,
  seed = 11
) {
  giotto_validate_scop_object(x)
  giotto_validate_named_list(filter_params, "filter_params")
  giotto_validate_named_list(norm_params, "norm_params")
  giotto_validate_named_list(stat_params, "stat_params")
  giotto_validate_named_list(hvf_params, "hvf_params")
  normalizeGiotto <- giotto_get_fun("normalizeGiotto")
  normalized <- FALSE
  x$giotto <- tryCatch(
    {
      out <- giotto_do_call(
        "normalizeGiotto",
        giotto_merge_args(
          list(
            gobject = x$giotto,
            verbose = FALSE
          ),
          norm_params,
          arg_name = "norm_params"
        )
      )
      normalized <<- TRUE
      out
    },
    error = function(e) {
      log_message(
        "Giotto normalization failed; continuing with existing expression values",
        message_type = "warning",
        verbose = verbose
      )
      x$giotto
    }
  )
  addStatistics <- tryCatch(giotto_get_fun("addStatistics"), error = function(e) NULL)
  if (is.function(addStatistics)) {
    x$giotto <- tryCatch(
      giotto_do_call(
        "addStatistics",
        giotto_merge_args(
          list(
            gobject = x$giotto,
            spat_unit = x$parameters$spat_unit %||% "cell",
            feat_type = x$parameters$feat_type %||% "rna",
            return_gobject = TRUE,
            verbose = verbose
          ),
          stat_params,
          arg_name = "stat_params"
        )
      ),
      error = function(e) {
        log_message(
          "Giotto statistics calculation failed; continuing without added statistics",
          message_type = "warning",
          verbose = verbose
        )
        x$giotto
      }
    )
  }
  calculateHVF <- tryCatch(giotto_get_fun("calculateHVF"), error = function(e) NULL)
  if (is.function(calculateHVF)) {
    x$giotto <- tryCatch(
      giotto_do_call(
        "calculateHVF",
        giotto_merge_args(
          list(
            gobject = x$giotto,
            spat_unit = x$parameters$spat_unit %||% "cell",
            feat_type = x$parameters$feat_type %||% "rna",
            expression_values = if (isTRUE(normalized)) "normalized" else (x$parameters$expression_values %||% "raw"),
            return_gobject = TRUE,
            show_plot = FALSE,
            set_seed = TRUE,
            seed_number = seed,
            verbose = verbose
          ),
          hvf_params,
          arg_name = "hvf_params",
          reserved = c("gobject", "return_gobject")
        )
      ),
      error = function(e) {
        log_message(
          "Giotto HVF calculation failed; continuing without HVF annotations",
          message_type = "warning",
          verbose = verbose
        )
        x$giotto
      }
    )
  }
  x$parameters$expression_values <- if (isTRUE(normalized)) {
    "normalized"
  } else {
    x$parameters$expression_values %||% "raw"
  }
  giotto_append_history(
    x,
    step = "GiottoPreprocess",
    parameters = list(
      filter_params = filter_params,
      norm_params = norm_params,
      stat_params = stat_params,
      hvf_params = hvf_params,
      seed = seed
    )
  )
}

#' @title Run Giotto dimensional reduction
#'
#' @param x A `giotto2` workflow object.
#' @param reduction Dimensional reduction to run.
#' @param dims Dimensions to use.
#' @param name Name for the Giotto reduction.
#' @param features Features used for the reduction.
#' @param params Additional parameters passed to the Giotto reduction function.
#' @param verbose Whether to print progress messages.
#' @param seed Random seed for reproducible Giotto calls.
#'
#' @return A `giotto2` workflow object.
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' embedding <- cbind(
#'   UMAP_1 = as.numeric(scale(spatial$x)),
#'   UMAP_2 = as.numeric(scale(spatial$y))
#' )
#' rownames(embedding) <- colnames(spatial)
#' g <- structure(
#'   list(
#'     giotto = list(umap = embedding),
#'     source = list(cells = colnames(spatial), features = rownames(spatial)),
#'     results = list(
#'       cluster = list(
#'         table = data.frame(
#'           cell = colnames(spatial),
#'           cluster = spatial$coda_label,
#'           row.names = colnames(spatial)
#'         )
#'       )
#'     ),
#'     parameters = list(umap_name = "umap")
#'   ),
#'   class = c("giotto2", "list")
#' )
#' GiottoPlot(g, plot_type = "dim")
#'
#'   g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
#'   g <- GiottoPreprocess(g)
#'   g <- GiottoReduce(g, reduction = "pca", dims = 1:10)
#'   g <- GiottoReduce(g, reduction = "umap", dims = 1:10)
#' @export
GiottoReduce <- function(
  x,
  reduction = c("pca", "umap"),
  dims = 1:20,
  name = NULL,
  features = NULL,
  params = list(),
  verbose = TRUE,
  seed = 11
) {
  giotto_validate_scop_object(x)
  reduction <- match.arg(reduction)
  giotto_validate_named_list(params, "params")
  cells_n <- length(x$source$cells %||% character())
  features_n <- length(x$source$features %||% character())
  dims <- unique(as.integer(dims))
  dims <- dims[is.finite(dims) & dims > 0L]
  dims <- dims[dims <= max(1L, min(cells_n, features_n) - 1L)]
  if (length(dims) == 0L) {
    dims <- 1L
  }
  features <- features %||% x$source$features
  if (identical(reduction, "pca")) {
    name <- name %||% "pca"
    x$giotto <- giotto_do_call(
      "runPCA",
      giotto_merge_args(
        list(
          gobject = x$giotto,
          spat_unit = x$parameters$spat_unit %||% "cell",
          feat_type = x$parameters$feat_type %||% "rna",
          expression_values = x$parameters$expression_values %||% "normalized",
          feats_to_use = features,
          name = name,
          ncp = max(dims),
          center = TRUE,
          scale_unit = TRUE,
          set_seed = TRUE,
          seed_number = seed,
          verbose = verbose
        ),
        params,
        arg_name = "params",
        reserved = c("gobject")
      )
    )
    x$parameters$pca_name <- name
  } else {
    name <- name %||% "umap"
    x$giotto <- giotto_do_call(
      "runUMAP",
      giotto_merge_args(
        list(
          gobject = x$giotto,
          spat_unit = x$parameters$spat_unit %||% "cell",
          feat_type = x$parameters$feat_type %||% "rna",
          dim_reduction_to_use = "pca",
          dim_reduction_name = x$parameters$pca_name %||% "pca",
          dimensions_to_use = dims,
          name = name,
          n_neighbors = max(2L, min(15L, cells_n - 1L)),
          set_seed = TRUE,
          seed_number = seed,
          verbose = verbose
        ),
        params,
        arg_name = "params",
        reserved = c("gobject")
      )
    )
    x$parameters$umap_name <- name
  }
  giotto_append_history(
    x,
    step = paste0("GiottoReduce:", reduction),
    parameters = list(dims = dims, name = name, seed = seed)
  )
}

#' @title Run Giotto nearest-network clustering
#'
#' @param x A `giotto2` workflow object.
#' @param method Giotto clustering method.
#' @param dims Dimensions used to build the nearest-neighbor network.
#' @param k Number of nearest neighbors.
#' @param resolution Clustering resolution.
#' @param network_name Name for the Giotto nearest-neighbor network.
#' @param cluster_name Name for the Giotto cluster result.
#' @param params Additional parameters passed to the Giotto clustering function.
#' @param verbose Whether to print progress messages.
#' @param seed Random seed for reproducible Giotto calls.
#'
#' @return A `giotto2` workflow object.
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' coords <- data.frame(
#'   cell_ID = colnames(spatial),
#'   sdimx = spatial$x,
#'   sdimy = spatial$y,
#'   row.names = colnames(spatial)
#' )
#' g <- structure(
#'   list(
#'     source = list(cells = colnames(spatial), coordinates = coords),
#'     results = list(
#'       cluster = list(
#'         table = data.frame(
#'           cell = colnames(spatial),
#'           cluster = spatial$coda_label,
#'           row.names = colnames(spatial)
#'         )
#'       )
#'     ),
#'     parameters = list(k = 8, resolution = 0.4)
#'   ),
#'   class = c("giotto2", "list")
#' )
#' GiottoPlot(g, plot_type = "cluster")
#'
#'   g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
#'   g <- GiottoPreprocess(g)
#'   g <- GiottoReduce(g, reduction = "pca", dims = 1:10)
#'   g <- GiottoCluster(g, dims = 1:10, k = 8, resolution = 0.4)
#' @export
GiottoCluster <- function(
  x,
  method = c("leiden", "louvain"),
  dims = 1:20,
  k = 20,
  resolution = 1,
  network_name = "scop_NN",
  cluster_name = NULL,
  params = list(),
  verbose = TRUE,
  seed = 11
) {
  giotto_validate_scop_object(x)
  method <- match.arg(method)
  giotto_validate_named_list(params, "params")
  cells_n <- length(x$source$cells %||% character())
  k <- min(as.integer(k), max(1L, cells_n - 1L))
  dims <- unique(as.integer(dims))
  dims <- dims[is.finite(dims) & dims > 0L]
  dims <- dims[dims <= max(1L, cells_n - 1L)]
  x$giotto <- giotto_do_call(
    "createNearestNetwork",
    list(
      gobject = x$giotto,
      spat_unit = x$parameters$spat_unit %||% "cell",
      feat_type = x$parameters$feat_type %||% "rna",
      dim_reduction_to_use = "pca",
      dim_reduction_name = x$parameters$pca_name %||% "pca",
      dimensions_to_use = dims,
      k = k,
      name = network_name
    )
  )
  cluster_name <- cluster_name %||% paste0(method, "_clus")
  cluster_fun_name <- if (identical(method, "leiden")) "doLeidenCluster" else "doLouvainCluster"
  x$giotto <- giotto_do_call(
    cluster_fun_name,
    giotto_merge_args(
      list(
        gobject = x$giotto,
        spat_unit = x$parameters$spat_unit %||% "cell",
        feat_type = x$parameters$feat_type %||% "rna",
        network_name = network_name,
        nn_network_to_use = "sNN",
        resolution = resolution,
        name = cluster_name,
        return_gobject = TRUE,
        set_seed = TRUE,
        seed_number = seed
      ),
      params,
      arg_name = "params",
      reserved = c("gobject", "return_gobject")
    )
  )
  meta <- giotto_get_cell_metadata(x$giotto)
  clusters <- giotto_extract_clusters(
    giotto_metadata = meta,
    cells = x$source$cells,
    cluster_name = cluster_name,
    method = method
  )
  x <- giotto_update_result(
    x,
    name = "cluster",
    result = list(
      table = data.frame(
        cell = names(clusters),
        cluster = as.character(clusters),
        row.names = names(clusters),
        stringsAsFactors = FALSE
      ),
      vector = clusters,
      metadata = meta,
      parameters = list(
        method = method,
        dims = dims,
        k = k,
        resolution = resolution,
        network_name = network_name,
        cluster_name = cluster_name
      )
    )
  )
  giotto_append_history(
    x,
    step = "GiottoCluster",
    parameters = list(method = method, dims = dims, k = k, resolution = resolution, seed = seed)
  )
}

#' @title Create a Giotto spatial network
#'
#' @param x A `giotto2` workflow object.
#' @param network_method Spatial network method.
#' @param network_name Name for the Giotto spatial network.
#' @param params Additional parameters passed to `Giotto::createSpatialNetwork()`.
#' @param verbose Whether to print progress messages.
#'
#' @return A `giotto2` workflow object.
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' coords <- data.frame(
#'   cell_ID = colnames(spatial),
#'   sdimx = spatial$x,
#'   sdimy = spatial$y,
#'   row.names = colnames(spatial)
#' )
#' edges <- data.frame(
#'   from = colnames(spatial)[1:79],
#'   to = colnames(spatial)[2:80]
#' )
#' g <- structure(
#'   list(
#'     source = list(cells = colnames(spatial), coordinates = coords),
#'     results = list(spatial_network = list(name = "Delaunay_network", table = edges))
#'   ),
#'   class = c("giotto2", "list")
#' )
#' GiottoPlot(g, plot_type = "network")
#'
#'   g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
#'   g <- GiottoSpatialNetwork(g, network_method = "Delaunay")
#' @export
GiottoSpatialNetwork <- function(
  x,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  params = list(),
  verbose = TRUE
) {
  giotto_validate_scop_object(x)
  network_method <- match.arg(network_method)
  giotto_validate_named_list(params, "params")
  network_name <- network_name %||% params[["name"]] %||% paste0(network_method, "_network")
  spatial_network <- giotto_create_spatial_network(
    gobject = x$giotto,
    network_name = network_name,
    network_method = network_method,
    network_params = params,
    verbose = verbose
  )
  x$giotto <- spatial_network$gobject
  network_table <- giotto_get_spatial_network_table(
    gobject = x$giotto,
    network_name = spatial_network$network_name,
    spat_unit = x$parameters$spat_unit %||% "cell",
    verbose = verbose
  )
  x <- giotto_update_result(
    x,
    name = "spatial_network",
    result = list(
      name = spatial_network$network_name,
      method = network_method,
      table = network_table,
      parameters = params
    ),
    active = FALSE
  )
  giotto_append_history(
    x,
    step = "GiottoSpatialNetwork",
    parameters = list(network_method = network_method, network_name = spatial_network$network_name)
  )
}

giotto_get_spatial_network_table <- function(gobject, network_name, spat_unit = "cell", verbose = TRUE) {
  table <- tryCatch(
    giotto_do_call(
      "getSpatialNetwork",
      list(
        gobject = gobject,
        spat_unit = spat_unit,
        name = network_name,
        output = "networkDT",
        verbose = verbose
      )
    ),
    error = function(e) NULL
  )
  if (is.null(table)) {
    return(NULL)
  }
  as.data.frame(table, stringsAsFactors = FALSE)
}

#' @title Run Giotto spatial gene detection
#'
#' @param x A `giotto2` workflow object.
#' @param features Features to test.
#' @param network_method Spatial network method.
#' @param network_name Name for the Giotto spatial network.
#' @param bin_method Binarization method passed to `Giotto::binSpect()`.
#' @param top_n Number of top spatial genes to store.
#' @param params Additional parameters passed to `Giotto::binSpect()`.
#' @param verbose Whether to print progress messages.
#' @param seed Random seed for reproducible Giotto calls.
#'
#' @return A `giotto2` workflow object.
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' coords <- data.frame(
#'   cell_ID = colnames(spatial),
#'   sdimx = spatial$x,
#'   sdimy = spatial$y,
#'   row.names = colnames(spatial)
#' )
#' spatial_gene_table <- data.frame(
#'   feat_ID = rownames(spatial)[1:8],
#'   spatGeneRank = seq_len(8),
#'   adj.p.value = seq(0.001, 0.04, length.out = 8)
#' )
#' g <- structure(
#'   list(
#'     source = list(cells = colnames(spatial), coordinates = coords),
#'     results = list(
#'       spatial_genes = list(
#'         table = spatial_gene_table,
#'         top_features = spatial_gene_table$feat_ID
#'       )
#'     )
#'   ),
#'   class = c("giotto2", "list")
#' )
#' GiottoPlot(g, plot_type = "spatial_genes", top_n = 6)
#'
#'   g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
#'   g <- GiottoPreprocess(g)
#'   g <- GiottoSpatialNetwork(g)
#'   g <- GiottoSpatialGenes(g, features = rownames(spatial)[1:50], top_n = 10)
#' @export
GiottoSpatialGenes <- function(
  x,
  features = NULL,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  bin_method = c("kmeans", "rank"),
  top_n = 100,
  params = list(),
  verbose = TRUE,
  seed = 11
) {
  giotto_validate_scop_object(x)
  network_method <- match.arg(network_method)
  bin_method <- match.arg(bin_method)
  features <- features %||% x$source$features
  network_info <- giotto_ensure_spatial_network(
    x = x,
    network_method = network_method,
    network_name = network_name,
    verbose = verbose
  )
  x <- network_info$x
  network_name <- network_info$network_name
  raw <- giotto_do_call(
    "binSpect",
    giotto_merge_args(
      list(
        gobject = x$giotto,
        spat_unit = x$parameters$spat_unit %||% "cell",
        feat_type = x$parameters$feat_type %||% "rna",
        bin_method = bin_method,
        expression_values = x$parameters$expression_values %||% "normalized",
        subset_feats = features,
        spatial_network_name = network_name,
        do_parallel = FALSE,
        cores = 1,
        seed = seed,
        verbose = verbose
      ),
      params,
      arg_name = "params"
    )
  )
  table <- giotto_result_to_data_frame(raw)
  top_features <- giotto_top_features(table, n = as.integer(top_n))
  x <- giotto_update_result(
    x,
    name = "spatial_genes",
    result = list(
      table = table,
      top_features = top_features,
      raw = raw,
      parameters = list(network_name = network_name, bin_method = bin_method, top_n = top_n)
    )
  )
  giotto_append_history(
    x,
    step = "GiottoSpatialGenes",
    parameters = list(network_name = network_name, bin_method = bin_method, top_n = top_n, seed = seed)
  )
}

#' @title Run Giotto spatial co-expression modules
#'
#' @param x A `giotto2` workflow object.
#' @param features Features to test.
#' @param network_method Spatial network method.
#' @param network_name Name for the Giotto spatial network.
#' @param cor_method Correlation method used by Giotto.
#' @param k Number of spatial co-expression modules.
#' @param detect_params Additional parameters passed to
#' `Giotto::detectSpatialCorFeats()`.
#' @param cluster_params Additional parameters passed to
#' `Giotto::clusterSpatialCorFeats()`.
#' @param verbose Whether to print progress messages.
#' @param seed Random seed for reproducible Giotto calls.
#'
#' @return A `giotto2` workflow object.
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' features <- rownames(spatial)[1:6]
#' module_table <- expand.grid(
#'   feat_ID = features,
#'   variable = paste0("module_", 1:3),
#'   stringsAsFactors = FALSE
#' )
#' module_table$spat_cor <- seq(-0.7, 0.8, length.out = nrow(module_table))
#' g <- structure(
#'   list(
#'     source = list(cells = colnames(spatial), features = features),
#'     results = list(
#'       spatial_modules = list(
#'         module_tables = list(result.cor_DT = module_table),
#'         features = features
#'       )
#'     )
#'   ),
#'   class = c("giotto2", "list")
#' )
#' GiottoPlot(g, plot_type = "spatial_modules", top_n = 6)
#'
#'   g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
#'   g <- GiottoPreprocess(g)
#'   g <- GiottoSpatialNetwork(g)
#'   g <- GiottoSpatialModules(g, features = rownames(spatial)[1:50], k = 3)
#' @export
GiottoSpatialModules <- function(
  x,
  features = NULL,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  cor_method = c("pearson", "spearman", "kendall"),
  k = 10,
  detect_params = list(),
  cluster_params = list(),
  verbose = TRUE,
  seed = 11
) {
  giotto_validate_scop_object(x)
  network_method <- match.arg(network_method)
  cor_method <- match.arg(cor_method)
  giotto_validate_named_list(detect_params, "detect_params")
  giotto_validate_named_list(cluster_params, "cluster_params")
  features <- features %||% x$results$spatial_genes$top_features %||% x$source$features
  k <- min(as.integer(k), max(2L, length(unique(features)) - 1L))
  network_info <- giotto_ensure_spatial_network(
    x = x,
    network_method = network_method,
    network_name = network_name,
    verbose = verbose
  )
  x <- network_info$x
  network_name <- network_info$network_name
  spat_cor <- giotto_do_call(
    "detectSpatialCorFeats",
    giotto_merge_args(
      list(
        gobject = x$giotto,
        spat_unit = x$parameters$spat_unit %||% "cell",
        feat_type = x$parameters$feat_type %||% "rna",
        method = "network",
        expression_values = x$parameters$expression_values %||% "normalized",
        subset_feats = features,
        spatial_network_name = network_name,
        cor_method = cor_method
      ),
      detect_params,
      arg_name = "detect_params"
    )
  )
  modules <- giotto_do_call(
    "clusterSpatialCorFeats",
    giotto_merge_args(
      list(
        spatCorObject = spat_cor,
        name = cluster_params[["name"]] %||% "spat_clus",
        k = as.integer(k),
        return_obj = TRUE
      ),
      cluster_params,
      arg_name = "cluster_params",
      reserved = c("spatCorObject", "return_obj")
    )
  )
  x <- giotto_update_result(
    x,
    name = "spatial_modules",
    result = list(
      spatial_cor = spat_cor,
      modules = modules,
      module_tables = giotto_collect_data_frames(modules),
      parameters = list(network_name = network_name, cor_method = cor_method, k = k)
    )
  )
  giotto_append_history(
    x,
    step = "GiottoSpatialModules",
    parameters = list(network_name = network_name, cor_method = cor_method, k = k, seed = seed)
  )
}

#' @title Run Giotto cell proximity enrichment
#'
#' @param x A `giotto2` workflow object.
#' @param group.by Metadata column containing cell or spot groups.
#' @param network_method Spatial network method.
#' @param network_name Name for the Giotto spatial network.
#' @param number_of_simulations Number of label simulations used by Giotto.
#' @param adjust_method Multiple-testing correction method.
#' @param params Additional parameters passed to `Giotto::cellProximityEnrichment()`.
#' @param verbose Whether to print progress messages.
#' @param seed Random seed for reproducible Giotto calls.
#'
#' @return A `giotto2` workflow object.
#' @examples
#' proximity <- data.frame(
#'   group_1 = c("Ductal", "Ductal", "Endocrine", "Stromal"),
#'   group_2 = c("Endocrine", "Stromal", "Stromal", "Ductal"),
#'   enrichment = c(1.6, 0.8, 1.3, 0.7),
#'   p.adj = c(0.01, 0.08, 0.03, 0.12)
#' )
#' g <- structure(
#'   list(
#'     results = list(cell_proximity = list(table = proximity)),
#'     parameters = list(network_method = "Delaunay", number_of_simulations = 100)
#'   ),
#'   class = c("giotto2", "list")
#' )
#' GiottoPlot(g, plot_type = "cell_proximity")
#'
#'   data(visium_human_pancreas_sub)
#'   spatial <- visium_human_pancreas_sub
#'   g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
#'   g <- GiottoSpatialNetwork(g)
#'   g <- GiottoCellProximity(g, group.by = "coda_label", number_of_simulations = 100)
#' @export
GiottoCellProximity <- function(
  x,
  group.by,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  number_of_simulations = 1000,
  adjust_method = "fdr",
  params = list(),
  verbose = TRUE,
  seed = 11
) {
  giotto_validate_scop_object(x)
  giotto_validate_scalar_string(group.by, "group.by")
  network_method <- match.arg(network_method)
  adjust_method <- match.arg(adjust_method, c("none", "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel", "BY"))
  network_info <- giotto_ensure_spatial_network(
    x = x,
    network_method = network_method,
    network_name = network_name,
    verbose = verbose
  )
  x <- network_info$x
  network_name <- network_info$network_name
  raw <- giotto_do_call(
    "cellProximityEnrichment",
    giotto_merge_args(
      list(
        gobject = x$giotto,
        spat_unit = x$parameters$spat_unit %||% "cell",
        feat_type = x$parameters$feat_type %||% "rna",
        spatial_network_name = network_name,
        cluster_column = group.by,
        number_of_simulations = as.integer(number_of_simulations),
        adjust_method = adjust_method,
        set_seed = TRUE,
        seed_number = seed
      ),
      params,
      arg_name = "params"
    )
  )
  table <- giotto_extract_proximity_enrichment(raw)
  x <- giotto_update_result(
    x,
    name = "cell_proximity",
    result = list(
      table = table,
      raw = raw,
      parameters = list(group.by = group.by, network_name = network_name, number_of_simulations = number_of_simulations)
    )
  )
  giotto_append_history(
    x,
    step = "GiottoCellProximity",
    parameters = list(group.by = group.by, network_name = network_name, number_of_simulations = number_of_simulations, seed = seed)
  )
}

#' @title Run Giotto HMRF spatial domains
#'
#' @param x A `giotto2` workflow object.
#' @param spatial_genes Spatial genes used by Giotto HMRF.
#' @param network_name Name for the Giotto spatial network.
#' @param k Number of HMRF domains.
#' @param betas HMRF beta values.
#' @param hmrf_name Name for the HMRF result.
#' @param params Additional parameters passed to `Giotto::doHMRF()`.
#' @param verbose Whether to print progress messages.
#' @param seed Random seed for reproducible Giotto calls.
#'
#' @return A `giotto2` workflow object.
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' coords <- data.frame(
#'   cell_ID = colnames(spatial),
#'   sdimx = spatial$x,
#'   sdimy = spatial$y,
#'   row.names = colnames(spatial)
#' )
#' hmrf_meta <- data.frame(
#'   cell_ID = colnames(spatial),
#'   scop_HMRF_k4_b10 = paste0("domain_", as.integer(factor(spatial$coda_label)) %% 4 + 1),
#'   row.names = colnames(spatial)
#' )
#' g <- structure(
#'   list(
#'     source = list(cells = colnames(spatial), coordinates = coords),
#'     results = list(
#'       hmrf = list(
#'         table = hmrf_meta["scop_HMRF_k4_b10"],
#'         metadata = hmrf_meta
#'       )
#'     )
#'   ),
#'   class = c("giotto2", "list")
#' )
#' GiottoPlot(g, plot_type = "hmrf")
#'
#'   g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
#'   g <- GiottoPreprocess(g)
#'   g <- GiottoSpatialNetwork(g)
#'   g <- GiottoHMRF(g, spatial_genes = rownames(spatial)[1:30], k = 4)
#' @export
GiottoHMRF <- function(
  x,
  spatial_genes = NULL,
  network_name = "Delaunay_full",
  k = 20,
  betas = c(0, 10, 20),
  hmrf_name = "scop_HMRF",
  params = list(),
  verbose = TRUE,
  seed = 11
) {
  giotto_validate_scop_object(x)
  giotto_validate_named_list(params, "params")
  network_info <- giotto_ensure_spatial_network(
    x = x,
    network_method = "Delaunay",
    network_name = network_name,
    verbose = verbose
  )
  x <- network_info$x
  spatial_genes <- spatial_genes %||% x$results$spatial_genes$top_features %||% utils::head(x$source$features, 50)
  raw <- giotto_do_call(
    "doHMRF",
    giotto_merge_args(
      list(
        gobject = x$giotto,
        spat_unit = x$parameters$spat_unit %||% "cell",
        feat_type = x$parameters$feat_type %||% "rna",
        expression_values = x$parameters$expression_values %||% "normalized",
        spatial_network_name = network_info$network_name,
        spatial_genes = spatial_genes,
        seed = seed,
        name = hmrf_name,
        k = as.integer(k),
        betas = betas,
        overwrite_output = TRUE
      ),
      params,
      arg_name = "params"
    )
  )
  x$giotto <- giotto_do_call(
    "addHMRF",
    list(
      gobject = x$giotto,
      spat_unit = x$parameters$spat_unit %||% "cell",
      feat_type = x$parameters$feat_type %||% "rna",
      HMRFoutput = raw,
      k = as.integer(k),
      betas_to_add = betas,
      hmrf_name = hmrf_name
    )
  )
  meta <- giotto_get_cell_metadata(x$giotto)
  hmrf_cols <- grep(hmrf_name, colnames(meta), value = TRUE, fixed = TRUE)
  x <- giotto_update_result(
    x,
    name = "hmrf",
    result = list(
      table = meta[, hmrf_cols, drop = FALSE],
      raw = raw,
      metadata = meta,
      parameters = list(network_name = network_info$network_name, spatial_genes = spatial_genes, k = k, betas = betas, hmrf_name = hmrf_name)
    )
  )
  giotto_append_history(
    x,
    step = "GiottoHMRF",
    parameters = list(network_name = network_info$network_name, k = k, betas = betas, hmrf_name = hmrf_name, seed = seed)
  )
}

giotto_ensure_spatial_network <- function(x, network_method = "Delaunay", network_name = NULL, verbose = TRUE) {
  if (
    !is.null(network_name) &&
      !is.null(x$results$spatial_network$name) &&
      identical(network_name, x$results$spatial_network$name)
  ) {
    return(list(x = x, network_name = network_name))
  }
  if (is.null(network_name) && !is.null(x$results$spatial_network$name)) {
    return(list(x = x, network_name = x$results$spatial_network$name))
  }
  x <- GiottoSpatialNetwork(
    x,
    network_method = network_method,
    network_name = network_name,
    verbose = verbose
  )
  list(x = x, network_name = x$results$spatial_network$name)
}

#' @title Add Giotto results back to Seurat
#'
#' @param srt A Seurat object.
#' @param x A `giotto2` workflow object.
#' @param result Giotto result to copy back.
#' @param name Metadata column name to write. If `NULL`, a default name is used.
#' @param tool_name Name used to store the Giotto workflow object in
#' `srt@tools`.
#' @param store_result Whether to store the Giotto workflow object in
#' `srt@tools[[tool_name]]`.
#'
#' @return A Seurat object.
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' g <- structure(
#'   list(
#'     source = list(cells = colnames(spatial)),
#'     results = list(
#'       cluster = list(
#'         table = data.frame(
#'           cell = colnames(spatial),
#'           cluster = spatial$coda_label,
#'           row.names = colnames(spatial)
#'         )
#'       )
#'     )
#'   ),
#'   class = c("giotto2", "list")
#' )
#' spatial <- AddGiottoToSeurat(
#'   spatial,
#'   g,
#'   result = "cluster",
#'   name = "Giotto_cluster",
#'   store_result = FALSE
#' )
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "Giotto_cluster",
#'   plot_type = "point",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#' @export
AddGiottoToSeurat <- function(
  srt,
  x,
  result = c("cluster", "hmrf"),
  name = NULL,
  tool_name = "Giotto",
  store_result = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  giotto_validate_scop_object(x)
  result <- match.arg(result)
  if (identical(result, "cluster")) {
    table <- x$results$cluster$table
    if (is.null(table) || !"cluster" %in% colnames(table)) {
      log_message("No Giotto cluster result is available", message_type = "error")
    }
    col <- name %||% "Giotto_cluster"
    cells <- intersect(rownames(table), colnames(srt))
    srt@meta.data[[col]] <- NA_character_
    srt@meta.data[cells, col] <- as.character(table[cells, "cluster"])
  } else {
    meta <- x$results$hmrf$metadata
    hmrf_cols <- colnames(x$results$hmrf$table %||% data.frame())
    if (is.null(meta) || length(hmrf_cols) == 0L) {
      log_message("No Giotto HMRF result is available", message_type = "error")
    }
    col <- name %||% hmrf_cols[1L]
    ids <- giotto_metadata_cell_ids(meta, cells = x$source$cells)
    values <- meta[[hmrf_cols[1L]]]
    names(values) <- ids
    cells <- intersect(names(values), colnames(srt))
    srt@meta.data[[col]] <- NA
    srt@meta.data[cells, col] <- values[cells]
  }
  if (isTRUE(store_result)) {
    srt@tools[[tool_name]] <- x
  }
  srt
}

#' @export
GiottoPlot.giotto2 <- function(
  x,
  plot_type = c("spatial", "cluster", "dim", "network", "spatial_genes", "spatial_modules", "cell_proximity", "hmrf"),
  result = NULL,
  feature = NULL,
  top_n = 20,
  palette = "Chinese",
  palcolor = NULL,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  pt.size = NULL,
  pt.alpha = 0.95,
  stroke = 0.08,
  bg_color = "grey25",
  edge.color = "grey70",
  edge.alpha = 0.45,
  edge.linewidth = 0.25,
  legend.position = "right",
  theme_use = "theme_scop",
  theme_args = list(),
  title = NULL,
  subtitle = NULL,
  ...
) {
  giotto_validate_scop_object(x)
  plot_type <- match.arg(plot_type)
  if (plot_type %in% c("spatial", "cluster", "hmrf")) {
    dat <- giotto_plot_spatial_data(x, plot_type = plot_type, result = result)
    color_col <- if (plot_type == "spatial") NULL else ".value"
    return(giotto_scop_spatial_points(
      dat = dat,
      color_col = color_col,
      palette = palette,
      palcolor = palcolor,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      stroke = stroke,
      bg_color = bg_color,
      legend.position = legend.position,
      theme_use = theme_use,
      theme_args = theme_args,
      title = title %||% giotto_plot_default_title(plot_type),
      subtitle = subtitle
    ))
  }
  if (identical(plot_type, "dim")) {
    dat <- giotto_dim_plot_data(x, result = result)
    return(giotto_scop_spatial_points(
      dat = dat,
      color_col = ".value",
      x_col = "dim1",
      y_col = "dim2",
      palette = palette,
      palcolor = palcolor,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      stroke = stroke,
      bg_color = bg_color,
      legend.position = legend.position,
      theme_use = theme_use,
      theme_args = theme_args,
      title = title %||% "Giotto dimension reduction",
      subtitle = subtitle
    ))
  }
  if (identical(plot_type, "spatial_genes")) {
    return(giotto_plot_spatial_genes_scop(x, top_n = top_n, palette = palette, palcolor = palcolor, theme_use = theme_use, theme_args = theme_args, title = title, subtitle = subtitle))
  }
  if (identical(plot_type, "cell_proximity")) {
    return(giotto_plot_proximity_scop(x, heatmap_palette = heatmap_palette, heatmap_palcolor = heatmap_palcolor, theme_use = theme_use, theme_args = theme_args, title = title, subtitle = subtitle))
  }
  if (identical(plot_type, "spatial_modules")) {
    return(giotto_plot_modules_scop(x, top_n = top_n, heatmap_palette = heatmap_palette, heatmap_palcolor = heatmap_palcolor, theme_use = theme_use, theme_args = theme_args, title = title, subtitle = subtitle))
  }
  if (identical(plot_type, "network")) {
    network_dat <- giotto_network_plot_data(x, result = result)
    return(giotto_scop_network(
      nodes = network_dat$nodes,
      edges = network_dat$edges,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      bg_color = bg_color,
      edge.color = edge.color,
      edge.alpha = edge.alpha,
      edge.linewidth = edge.linewidth,
      legend.position = legend.position,
      theme_use = theme_use,
      theme_args = theme_args,
      title = title %||% "Giotto spatial network",
      subtitle = subtitle %||% "Edges from Giotto spatial network"
    ))
  }
}

#' @export
plot.giotto2 <- function(x, y = NULL, ...) {
  GiottoPlot(x, ...)
}

giotto_plot_spatial_data <- function(x, plot_type = "spatial", result = NULL) {
  coords <- x$source$coordinates
  if (is.null(coords) || nrow(coords) == 0L) {
    log_message("No spatial coordinates are available in {.cls giotto2}", message_type = "error")
  }
  dat <- data.frame(
    cell = coords$cell_ID %||% rownames(coords),
    x = coords$sdimx %||% coords$x,
    y = coords$sdimy %||% coords$y,
    stringsAsFactors = FALSE
  )
  if (identical(plot_type, "cluster")) {
    table <- x$results$cluster$table
    if (is.null(table)) {
      log_message("No Giotto cluster result is available", message_type = "error")
    }
    values <- table$cluster
    names(values) <- rownames(table)
    dat$.value <- as.character(values[dat$cell])
  } else if (identical(plot_type, "hmrf")) {
    meta <- x$results$hmrf$metadata
    hmrf_table <- x$results$hmrf$table
    if (is.null(meta) || is.null(hmrf_table) || ncol(hmrf_table) == 0L) {
      log_message("No Giotto HMRF result is available", message_type = "error")
    }
    col <- result %||% colnames(hmrf_table)[1L]
    ids <- giotto_metadata_cell_ids(meta, cells = x$source$cells)
    values <- meta[[col]]
    names(values) <- ids
    dat$.value <- as.character(values[dat$cell])
  }
  dat
}

giotto_network_plot_data <- function(x, result = NULL) {
  nodes <- giotto_plot_spatial_data(x, plot_type = "spatial", result = result)
  net <- x$results$spatial_network$table
  if (is.null(net) || nrow(net) == 0L) {
    return(list(nodes = nodes, edges = NULL))
  }
  net <- as.data.frame(net, stringsAsFactors = FALSE)
  if (all(c("sdimx_begin", "sdimy_begin", "sdimx_end", "sdimy_end") %in% colnames(net))) {
    edges <- data.frame(
      x = net$sdimx_begin,
      y = net$sdimy_begin,
      xend = net$sdimx_end,
      yend = net$sdimy_end,
      stringsAsFactors = FALSE
    )
  } else if (all(c("from", "to") %in% colnames(net))) {
    coords <- nodes[, c("cell", "x", "y")]
    from <- coords[match(net$from, coords$cell), , drop = FALSE]
    to <- coords[match(net$to, coords$cell), , drop = FALSE]
    edges <- data.frame(
      x = from$x,
      y = from$y,
      xend = to$x,
      yend = to$y,
      stringsAsFactors = FALSE
    )
  } else {
    edges <- NULL
  }
  if (!is.null(edges)) {
    edges <- edges[stats::complete.cases(edges), , drop = FALSE]
  }
  list(nodes = nodes, edges = edges)
}

giotto_scop_spatial_points <- function(
  dat,
  color_col = NULL,
  x_col = "x",
  y_col = "y",
  palette = "Chinese",
  palcolor = NULL,
  pt.size = NULL,
  pt.alpha = 0.95,
  stroke = 0.08,
  bg_color = "grey25",
  legend.position = "right",
  theme_use = "theme_scop",
  theme_args = list(),
  title = NULL,
  subtitle = NULL
) {
  pt.size <- pt.size %||% 1.8
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]))
  if (is.null(color_col)) {
    p <- p + ggplot2::geom_point(color = bg_color, size = pt.size, alpha = pt.alpha)
  } else {
    dat[[color_col]] <- factor(dat[[color_col]])
    cols <- palette_colors(levels(dat[[color_col]]), palette = palette, palcolor = palcolor)
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[color_col]])) +
      ggplot2::geom_point(shape = 21, color = bg_color, stroke = stroke, size = pt.size, alpha = pt.alpha) +
      ggplot2::scale_fill_manual(values = cols, na.value = "grey85")
  }
  p +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL, fill = NULL) +
    giotto_plot_theme(theme_use, theme_args) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = legend.position,
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey35")
    )
}

giotto_scop_network <- function(
  nodes,
  edges = NULL,
  pt.size = NULL,
  pt.alpha = 0.95,
  bg_color = "grey25",
  edge.color = "grey70",
  edge.alpha = 0.45,
  edge.linewidth = 0.25,
  legend.position = "right",
  theme_use = "theme_scop",
  theme_args = list(),
  title = NULL,
  subtitle = NULL
) {
  pt.size <- pt.size %||% 1.8
  p <- ggplot2::ggplot(nodes, ggplot2::aes(x = .data$x, y = .data$y))
  if (!is.null(edges) && nrow(edges) > 0L) {
    p <- p + ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      inherit.aes = FALSE,
      color = edge.color,
      alpha = edge.alpha,
      linewidth = edge.linewidth,
      lineend = "round"
    )
  }
  p +
    ggplot2::geom_point(color = bg_color, size = pt.size, alpha = pt.alpha) +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    giotto_plot_theme(theme_use, theme_args) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = legend.position,
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey35")
    )
}

giotto_plot_default_title <- function(plot_type) {
  switch(plot_type,
    spatial = "Giotto spatial coordinates",
    cluster = "Giotto clusters",
    hmrf = "Giotto HMRF domains",
    "Giotto plot"
  )
}

giotto_dim_plot_data <- function(x, result = NULL) {
  reduction <- result %||% x$parameters$umap_name %||% "umap"
  dat <- giotto_find_reduction_data(x$giotto, cells = x$source$cells, reduction = reduction)
  if (is.null(dat)) {
    log_message("No Giotto dimension reduction coordinates are available for plotting", message_type = "error")
  }
  table <- x$results$cluster$table
  dat$.value <- if (!is.null(table)) {
    values <- table$cluster
    names(values) <- rownames(table)
    as.character(values[dat$cell])
  } else {
    "Giotto"
  }
  dat
}

giotto_find_reduction_data <- function(x, cells, reduction = "umap") {
  accessor_dat <- giotto_get_reduction_matrix(
    gobject = x,
    cells = cells,
    reduction = reduction
  )
  if (!is.null(accessor_dat)) {
    return(accessor_dat)
  }
  seen <- new.env(parent = emptyenv())
  recurse <- function(obj) {
    key <- tryCatch(paste0(class(obj)[1], ":", utils::object.size(obj)), error = function(e) NULL)
    if (!is.null(key) && exists(key, envir = seen, inherits = FALSE)) {
      return(NULL)
    }
    if (!is.null(key)) {
      assign(key, TRUE, envir = seen)
    }
    if (is.matrix(obj) || inherits(obj, "data.frame")) {
      df <- as.data.frame(obj)
      rn <- rownames(df)
      if (!is.null(rn) && length(intersect(rn, cells)) > 0L && ncol(df) >= 2L) {
        num_cols <- colnames(df)[vapply(df, is.numeric, logical(1))]
        if (length(num_cols) >= 2L) {
          out <- data.frame(
            cell = rn,
            dim1 = df[[num_cols[1L]]],
            dim2 = df[[num_cols[2L]]],
            stringsAsFactors = FALSE
          )
          return(out[out$cell %in% cells, , drop = FALSE])
        }
      }
    }
    if (is.list(obj) || methods::is(obj, "S4")) {
      vals <- if (is.list(obj)) obj else lapply(methods::slotNames(obj), function(nm) methods::slot(obj, nm))
      nms <- names(vals) %||% character(length(vals))
      ord <- order(!grepl(reduction, nms, ignore.case = TRUE))
      for (i in ord) {
        hit <- tryCatch(recurse(vals[[i]]), error = function(e) NULL)
        if (!is.null(hit)) {
          return(hit)
        }
      }
    }
    NULL
  }
  recurse(x)
}

giotto_get_reduction_matrix <- function(gobject, cells, reduction = "umap") {
  method <- if (grepl("umap", reduction, ignore.case = TRUE)) {
    "umap"
  } else if (grepl("pca", reduction, ignore.case = TRUE)) {
    "pca"
  } else {
    reduction
  }
  candidates <- unique(c(reduction, method, NA_character_))
  for (candidate in candidates) {
    mat <- tryCatch(
      giotto_do_call(
        "getDimReduction",
        list(
          gobject = gobject,
          spat_unit = "cell",
          feat_type = "rna",
          reduction = "cells",
          reduction_method = method,
          name = if (is.na(candidate)) NULL else candidate,
          output = "matrix"
        )
      ),
      error = function(e) NULL
    )
    if (!is.null(mat) && ncol(mat) >= 2L && nrow(mat) > 0L) {
      df <- as.data.frame(mat)
      rn <- rownames(df)
      if (is.null(rn)) {
        rn <- cells[seq_len(min(length(cells), nrow(df)))]
      }
      num_cols <- colnames(df)[vapply(df, is.numeric, logical(1))]
      if (length(num_cols) >= 2L) {
        out <- data.frame(
          cell = rn,
          dim1 = df[[num_cols[1L]]],
          dim2 = df[[num_cols[2L]]],
          stringsAsFactors = FALSE
        )
        return(out[out$cell %in% cells, , drop = FALSE])
      }
    }
  }
  NULL
}

giotto_plot_spatial_genes_scop <- function(x, top_n = 20, palette = "Chinese", palcolor = NULL, theme_use = "theme_scop", theme_args = list(), title = NULL, subtitle = NULL) {
  table <- x$results$spatial_genes$table
  if (is.null(table) || nrow(table) == 0L) {
    log_message("No Giotto spatial gene result is available", message_type = "error")
  }
  table <- utils::head(table, as.integer(top_n))
  feature_col <- giotto_plot_pick_col(table, c("feats", "feat_ID", "gene", "feature", "features"))
  score_col <- giotto_plot_pick_numeric_col(table, c("score", "spatGeneRank", "estimate", "statistic", "p.value", "adj.p.value"))
  table$.feature <- factor(as.character(table[[feature_col]]), levels = rev(unique(as.character(table[[feature_col]]))))
  table$.score <- table[[score_col]]
  if (grepl("p[._]?value|p\\.adj|adj", score_col, ignore.case = TRUE)) {
    table$.score <- -log10(pmax(table$.score, .Machine$double.xmin))
  }
  fill_col <- palette_colors(palette = palette, palcolor = palcolor, n = 1L)[[1L]]
  ggplot2::ggplot(table, ggplot2::aes(x = .data$.score, y = .data$.feature)) +
    ggplot2::geom_col(fill = fill_col, width = 0.72) +
    ggplot2::labs(
      title = title %||% "Giotto spatial genes",
      subtitle = subtitle %||% paste0("Top ", nrow(table), " features"),
      x = giotto_plot_pretty_label(score_col),
      y = NULL
    ) +
    giotto_plot_theme(theme_use, theme_args) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

giotto_plot_proximity_scop <- function(x, heatmap_palette = "RdBu", heatmap_palcolor = NULL, theme_use = "theme_scop", theme_args = list(), title = NULL, subtitle = NULL) {
  table <- x$results$cell_proximity$table
  if (is.null(table) || nrow(table) == 0L) {
    log_message("No Giotto cell proximity result is available", message_type = "error")
  }
  pseudo <- giotto_result(
    result_type = "cell_proximity",
    giotto = x$giotto,
    enrichment = table,
    parameters = x$results$cell_proximity$parameters %||% list()
  )
  GiottoPlot.giotto2_cell_proximity(
    pseudo,
    heatmap_palette = heatmap_palette,
    heatmap_palcolor = heatmap_palcolor,
    theme_use = theme_use,
    theme_args = theme_args,
    title = title %||% "Giotto cell proximity enrichment",
    subtitle = subtitle
  )
}

giotto_plot_modules_scop <- function(x, top_n = 20, heatmap_palette = "RdBu", heatmap_palcolor = NULL, theme_use = "theme_scop", theme_args = list(), title = NULL, subtitle = NULL) {
  pseudo <- giotto_result(
    result_type = "spatial_modules",
    giotto = x$giotto,
    module_tables = x$results$spatial_modules$module_tables,
    features = x$source$features,
    parameters = x$results$spatial_modules$parameters %||% list()
  )
  GiottoPlot.giotto2_spatial_modules(
    pseudo,
    top_n = top_n,
    heatmap_palette = heatmap_palette,
    heatmap_palcolor = heatmap_palcolor,
    theme_use = theme_use,
    theme_args = theme_args,
    title = title %||% "Giotto spatial co-expression",
    subtitle = subtitle %||% "Spatial correlation among top features"
  )
}
