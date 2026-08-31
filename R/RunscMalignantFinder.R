#' @title Run scMalignantFinder malignant cell identification
#'
#' @description
#' Run the Python package `scMalignantFinder` on a Seurat or AnnData object and
#' append malignant-cell predictions to Seurat metadata. The pretrained model
#' files are not bundled with `scop`; provide a directory containing
#' `model.joblib` and `ordered_feature.tsv` through `pretrain_dir`.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A Seurat object.
#' @param adata Optional Python AnnData object.
#' @param h5ad Optional path to an `.h5ad` file.
#' @param assay Assay used when `srt` is supplied. Default is `"RNA"`.
#' @param layer Layer used when `srt` is supplied. Default is `"counts"`.
#' @param cells Optional cells to run. If supplied with `srt`, results are
#' appended to these cells and other cells receive `NA`.
#' @param pretrain_dir Directory containing pretrained `scMalignantFinder`
#' model files: `model.joblib` and `ordered_feature.tsv`.
#' @param train_h5ad_path Optional training `.h5ad` file used when training a
#' model from scratch.
#' @param feature_path Optional feature file used when training from
#' scratch.
#' @param model_method Model used when training from scratch. One of
#' `"LogisticRegression"`, `"RandomForest"`, or `"XGBoost"`.
#' @param norm_type Passed to `scMalignantFinder`. Use `TRUE` for raw counts
#' that should be library-size normalized; use `FALSE` for already normalized
#' input. If `NULL`, defaults to `TRUE` only for Seurat counts input.
#' @param use_raw Whether to use `adata.raw.X` when available.
#' @param n_thread Number of threads used by `scMalignantFinder`.
#' @param prefix Optional prefix for output metadata columns. Default preserves
#' the original `scMalignantFinder` column names.
#' @param return_seurat Whether to return a Seurat object when `srt` is
#' supplied. If `FALSE`, returns a data frame of predictions.
#'
#' @return A Seurat object with `scMalignantFinder_prediction` and
#' `malignancy_probability` added, or a data frame when `return_seurat = FALSE`.
#'
#' @references
#' Yu Q, Li YY, Chen Y. scMalignantFinder distinguishes malignant cells in
#' single-cell and spatial transcriptomics by leveraging cancer signatures.
#' Communications Biology. 2025.
#'
#' @seealso [srt_to_adata], [RunscMalignantRegion], [RunscMalignantStates]
#'
#' @export
#'
#' @examplesIf FALSE
#' data(pancreas_sub)
#' pancreas_sub <- RunscMalignantFinder(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "counts",
#'   pretrain_dir = "path/to/pretrained_model"
#' )
#' CellDimPlot(pancreas_sub, group.by = "malignancy_probability")
RunscMalignantFinder <- function(
  srt = NULL,
  adata = NULL,
  h5ad = NULL,
  assay = "RNA",
  layer = "counts",
  cells = NULL,
  pretrain_dir = NULL,
  train_h5ad_path = NULL,
  feature_path = NULL,
  model_method = c("LogisticRegression", "RandomForest", "XGBoost"),
  norm_type = NULL,
  use_raw = FALSE,
  n_thread = 1,
  prefix = "",
  return_seurat = !is.null(srt),
  verbose = TRUE
) {
  model_method <- match.arg(model_method)
  scmf_check_one_input(srt = srt, adata = adata, h5ad = h5ad)
  norm_type <- scmf_resolve_norm_type(norm_type, srt = srt, layer = layer)
  h5ad <- scmf_expand_path(h5ad)
  pretrain_dir <- scmf_expand_path(pretrain_dir)
  train_h5ad_path <- scmf_expand_path(train_h5ad_path)
  feature_path <- scmf_expand_path(feature_path)
  if (!is.null(cells) && is.null(srt)) {
    log_message(
      "{.arg cells} is only supported when {.arg srt} is supplied",
      message_type = "error"
    )
  }
  if (is.null(pretrain_dir) && (is.null(train_h5ad_path) || is.null(feature_path))) {
    log_message(
      "Provide either {.arg pretrain_dir} or both {.arg train_h5ad_path} and {.arg feature_path}",
      message_type = "error"
    )
  }

  log_message(
    "Running {.pkg scMalignantFinder} malignant cell identification...",
    message_type = "running",
    verbose = verbose
  )
  scmf_prepare_python(verbose = verbose)
  if (is.null(pretrain_dir) && identical(model_method, "XGBoost")) {
    check_python("xgboost", verbose = verbose)
  }

  if (!is.null(srt)) {
    if (!inherits(srt, "Seurat")) {
      log_message(
        "{.arg srt} must be a {.cls Seurat} object",
        message_type = "error"
      )
    }
    cells <- scmf_cells(srt, cells)
    srt_input <- if (length(cells) == ncol(srt)) {
      srt
    } else {
      subset(srt, cells = cells)
    }
    adata <- if (isTRUE(use_raw)) {
      srt_to_adata(
        srt = srt_input,
        assay_x = assay,
        layer_x = layer,
        assay_y = character(0),
        verbose = verbose
      )
    } else {
      python_anndata_from_matrix(
        GetAssayData5(srt_input, assay = assay, layer = layer)
      )
    }
  }

  test_input <- adata %||% h5ad
  scmf <- import_scmalignantfinder(convert = TRUE)
  obs <- scmf$run_scmalignantfinder(
    test_input = test_input,
    pretrain_dir = pretrain_dir,
    train_h5ad_path = train_h5ad_path,
    feature_path = feature_path,
    model_method = model_method,
    norm_type = norm_type,
    use_raw = use_raw,
    n_thread = as.integer(n_thread),
    return_obs = TRUE,
    verbose = verbose
  )
  obs <- as.data.frame(obs)

  if (!isTRUE(return_seurat)) {
    return(obs)
  }
  if (is.null(srt)) {
    log_message(
      "{.arg return_seurat = TRUE} requires {.arg srt}",
      message_type = "error"
    )
  }

  source_cols <- c("scMalignantFinder_prediction", "malignancy_probability")
  output_cols <- {
    .inline0 <- source_cols
    .inline1 <- prefix
    stats::setNames(
      if (nzchar(.inline1)) paste0(.inline1, .inline0) else .inline0,
      .inline0
    )
  }
  srt <- scmf_append_obs_to_srt(
    srt = srt,
    obs = obs,
    source_cols = source_cols,
    output_cols = output_cols,
    verbose = verbose
  )
  srt@tools$scMalignantFinder <- list(
    method = "scMalignantFinder",
    task = "malignant_cell_identification",
    columns = unname(output_cols),
    parameters = list(
      assay = assay,
      layer = layer,
      cells = cells,
      pretrain_dir = pretrain_dir,
      train_h5ad_path = train_h5ad_path,
      feature_path = feature_path,
      model_method = model_method,
      norm_type = norm_type,
      use_raw = use_raw,
      n_thread = as.integer(n_thread)
    )
  )
  srt
}

#' @title Run scMalignantFinder malignant spatial region identification
#'
#' @description
#' Use `scMalignantFinder` spatial utilities to score malignant signatures and
#' infer malignant spatial regions. For Seurat input, provide `spatial.cols`
#' when spatial-neighborhood refinement is requested and the converted AnnData
#' object does not already contain spatial coordinates.
#'
#' @md
#' @inheritParams RunscMalignantFinder
#' @param signature_gmt Path to the malignant signature `.gmt` file, such as
#' `sc_malignant_deg.gmt` from the scMalignantFinder resources.
#' @param features Features in `adata.obs` used for region clustering.
#' If `NULL`, uses `malignancy_probability` and `Malignant_up`, plus
#' `image_score` when `image = TRUE`.
#' @param nclus Number of clusters used by the spatial region model.
#' @param define_feature Feature used to select the malignant region cluster.
#' @param spatial_nn Whether to refine labels by spatial neighbors.
#' @param spatial.cols Optional two metadata columns used as spatial
#' coordinates for Seurat input.
#' @param spatial_key Key written to `adata.obsm` for spatial coordinates.
#' @param image Whether to call `scMalignantFinder.spatial.image_cal`.
#' Default is `FALSE` because Seurat-to-AnnData conversion does not generally
#' carry histology images.
#' @param backend Region-inference backend. `"cpp"` uses a compiled sparse
#' AUCell scorer, Ward clustering, and spatial-neighbor refinement for Seurat
#' input when `image = FALSE`. `"python"` retains the official
#' scMalignantFinder/Squidpy path and is required for image features, AnnData,
#' and h5ad input.
#'
#' @return A Seurat object with malignant region metadata, or a data frame when
#' `return_seurat = FALSE`.
#'
#' @export
#'
#' @examplesIf FALSE
#' srt <- RunscMalignantRegion(
#'   srt,
#'   signature_gmt = "path/to/sc_malignant_deg.gmt",
#'   spatial.cols = c("x", "y")
#' )
RunscMalignantRegion <- function(
  srt = NULL,
  adata = NULL,
  h5ad = NULL,
  assay = "RNA",
  layer = "counts",
  cells = NULL,
  signature_gmt,
  features = NULL,
  nclus = 3,
  define_feature = "Malignant_up",
  spatial_nn = TRUE,
  spatial.cols = NULL,
  spatial_key = "spatial",
  image = FALSE,
  norm_type = NULL,
  prefix = "scMalignantFinder_",
  return_seurat = !is.null(srt),
  verbose = TRUE,
  backend = c("cpp", "python")
) {
  backend_missing <- missing(backend)
  scmf_check_one_input(srt = srt, adata = adata, h5ad = h5ad)
  backend <- match.arg(backend)
  if (is.null(srt) || isTRUE(image)) {
    if (backend_missing) {
      backend <- "python"
    } else if (identical(backend, "cpp")) {
      reason <- if (isTRUE(image)) "image feature extraction" else "AnnData or h5ad input"
      log_message(
        "{.arg backend = 'cpp'} does not support {reason}; use {.arg backend = 'python'}",
        message_type = "error"
      )
    }
  }
  norm_type <- scmf_resolve_norm_type(norm_type, srt = srt, layer = layer)
  h5ad <- scmf_expand_path(h5ad)
  signature_gmt <- scmf_expand_path(signature_gmt)
  if (missing(signature_gmt) || is.null(signature_gmt) || !nzchar(signature_gmt)) {
    log_message("{.arg signature_gmt} is required", message_type = "error")
  }
  if (!file.exists(path.expand(signature_gmt))) {
    log_message(
      "{.file {signature_gmt}} does not exist",
      message_type = "error"
    )
  }
  if (!is.null(cells) && is.null(srt)) {
    log_message(
      "{.arg cells} is only supported when {.arg srt} is supplied",
      message_type = "error"
    )
  }
  if (is.null(features)) {
    features <- c("malignancy_probability", "Malignant_up")
    if (isTRUE(image)) {
      features <- c(features, "image_score")
    }
  }
  if (!is.null(srt) && !is.null(spatial.cols)) {
    if (!inherits(srt, "Seurat")) {
      log_message(
        "{.arg srt} must be a {.cls Seurat} object",
        message_type = "error"
      )
    }
    precheck_cells <- scmf_cells(srt, cells)
    precheck_srt <- if (length(precheck_cells) == ncol(srt)) {
      srt
    } else {
      subset(srt, cells = precheck_cells)
    }
    scmf_get_spatial_coordinates(precheck_srt, spatial.cols)
  }

  spatial_coordinates <- NULL
  if (!is.null(srt)) {
    if (!inherits(srt, "Seurat")) {
      log_message(
        "{.arg srt} must be a {.cls Seurat} object",
        message_type = "error"
      )
    }
    cells <- scmf_cells(srt, cells)
    srt_input <- if (length(cells) == ncol(srt)) {
      srt
    } else {
      subset(srt, cells = cells)
    }
    spatial_coordinates <- scmf_get_spatial_coordinates(srt_input, spatial.cols)
  }

  log_message(
    "Running {.pkg scMalignantFinder} malignant region identification with the {.val {backend}} backend...",
    message_type = "running",
    verbose = verbose
  )
  if (identical(backend, "cpp")) {
    expr <- GetAssayData5(
      srt_input,
      assay = assay,
      layer = layer
    )
    obs <- scmf_native_region(
      expr = expr,
      metadata = srt_input[[]],
      signature_gmt = signature_gmt,
      features = features,
      nclus = nclus,
      define_feature = define_feature,
      spatial_nn = spatial_nn,
      spatial_coordinates = spatial_coordinates
    )
  } else {
    scmf_prepare_python(verbose = verbose)
    if (!is.null(srt)) {
      adata <- srt_to_adata(
        srt = srt_input,
        assay_x = assay,
        layer_x = layer,
        assay_y = character(0),
        verbose = verbose
      )
    }
    test_input <- adata %||% h5ad
    scmf <- import_scmalignantfinder(convert = TRUE)
    obs <- scmf$run_scmalignant_region(
      test_input = test_input,
      signature_gmt = signature_gmt,
      features = features,
      nclus = as.integer(nclus),
      define_feature = define_feature,
      spatial_nn = spatial_nn,
      spatial_coordinates = spatial_coordinates,
      spatial_key = spatial_key,
      image = image,
      norm_type = norm_type,
      return_obs = TRUE,
      verbose = verbose
    )
    obs <- as.data.frame(obs)
  }

  if (!isTRUE(return_seurat)) {
    return(obs)
  }
  if (is.null(srt)) {
    log_message(
      "{.arg return_seurat = TRUE} requires {.arg srt}",
      message_type = "error"
    )
  }

  source_cols <- intersect(
    c("cluster", "region_prediction", "Malignant_up", "image_score", features),
    colnames(obs)
  )
  output_cols <- scmf_region_output_names(source_cols, prefix = prefix)
  srt <- scmf_append_obs_to_srt(
    srt = srt,
    obs = obs,
    source_cols = source_cols,
    output_cols = output_cols,
    verbose = verbose
  )
  srt@tools$scMalignantFinder_region <- list(
    method = "scMalignantFinder",
    task = "malignant_region_identification",
    columns = unname(output_cols),
    parameters = list(
      assay = assay,
      layer = layer,
      cells = cells,
      signature_gmt = signature_gmt,
      features = features,
      nclus = as.integer(nclus),
      define_feature = define_feature,
      spatial_nn = spatial_nn,
      spatial.cols = spatial.cols,
      spatial_key = spatial_key,
      image = image,
      norm_type = norm_type,
      backend = backend
    )
  )
  srt
}

#' @title Run scMalignantFinder cancer cell state scoring
#'
#' @description
#' Score cancer cell state gene sets with `scMalignantFinder` AUCell utilities
#' and append the resulting activity scores to Seurat metadata.
#'
#' @md
#' @inheritParams RunscMalignantFinder
#' @param gene_sets Path to a `.gmt` file containing cancer cell state gene
#' sets, such as `Malignant_MPs.Gavish_2023.gmt` from the scMalignantFinder
#' resources.
#' @param backend State-scoring backend. `"cpp"` uses a compiled sparse
#' AUCell implementation for Seurat input. `"python"` retains the official
#' scMalignantFinder path and is used for AnnData or h5ad input.
#'
#' @return A Seurat object with cancer-state AUCell scores, or a data frame when
#' `return_seurat = FALSE`.
#'
#' @export
#'
#' @examplesIf FALSE
#' srt <- RunscMalignantStates(
#'   srt,
#'   gene_sets = "path/to/Malignant_MPs.Gavish_2023.gmt"
#' )
RunscMalignantStates <- function(
  srt = NULL,
  adata = NULL,
  h5ad = NULL,
  assay = "RNA",
  layer = "counts",
  cells = NULL,
  gene_sets,
  norm_type = NULL,
  prefix = "scMalignantState_",
  return_seurat = !is.null(srt),
  verbose = TRUE,
  backend = c("cpp", "python")
) {
  backend_missing <- missing(backend)
  scmf_check_one_input(srt = srt, adata = adata, h5ad = h5ad)
  backend <- match.arg(backend)
  if (is.null(srt)) {
    if (backend_missing) {
      backend <- "python"
    } else if (identical(backend, "cpp")) {
      log_message(
        "{.arg backend = 'cpp'} requires {.arg srt}; use {.arg backend = 'python'} for AnnData or h5ad input",
        message_type = "error"
      )
    }
  }
  norm_type <- scmf_resolve_norm_type(norm_type, srt = srt, layer = layer)
  h5ad <- scmf_expand_path(h5ad)
  gene_sets <- scmf_expand_path(gene_sets)
  if (missing(gene_sets) || is.null(gene_sets) || !nzchar(gene_sets)) {
    log_message("{.arg gene_sets} is required", message_type = "error")
  }
  if (!file.exists(path.expand(gene_sets))) {
    log_message("{.file {gene_sets}} does not exist", message_type = "error")
  }
  if (!is.null(cells) && is.null(srt)) {
    log_message(
      "{.arg cells} is only supported when {.arg srt} is supplied",
      message_type = "error"
    )
  }

  log_message(
    "Running {.pkg scMalignantFinder} cancer cell state scoring...",
    message_type = "running",
    verbose = verbose
  )
  if (!is.null(srt)) {
    if (!inherits(srt, "Seurat")) {
      log_message(
        "{.arg srt} must be a {.cls Seurat} object",
        message_type = "error"
      )
    }
    cells <- scmf_cells(srt, cells)
  }

  if (identical(backend, "cpp")) {
    expr <- GetAssayData5(
      srt,
      assay = assay,
      layer = layer
    )[, cells, drop = FALSE]
    state_sets <- scmf_adapt_gene_sets(
      scmf_read_gene_sets(gene_sets),
      rownames(expr)
    )
    obs <- as.data.frame(run_aucell_scores(
      expr_counts = expr,
      gene_sets = state_sets,
      strategy = "topk",
      algorithm = "ctxcore",
      seed = 2L,
      tie_method = "numpy",
      auc_threshold = max(1 / nrow(expr), scmf_auc_threshold(expr)),
      normalize_by_signature_max = TRUE
    ))
  } else {
    scmf_prepare_python(verbose = verbose)
    if (!is.null(srt)) {
      srt_input <- if (length(cells) == ncol(srt)) {
        srt
      } else {
        subset(srt, cells = cells)
      }
      adata <- srt_to_adata(
        srt = srt_input,
        assay_x = assay,
        layer_x = layer,
        assay_y = character(0),
        verbose = verbose
      )
    }

    test_input <- adata %||% h5ad
    scmf <- import_scmalignantfinder(convert = TRUE)
    obs <- scmf$run_scmalignant_states(
      test_input = test_input,
      gene_sets = gene_sets,
      norm_type = norm_type,
      return_obs = TRUE,
      verbose = verbose
    )
    obs <- as.data.frame(obs)
  }

  if (!isTRUE(return_seurat)) {
    return(obs)
  }
  if (is.null(srt)) {
    log_message(
      "{.arg return_seurat = TRUE} requires {.arg srt}",
      message_type = "error"
    )
  }

  source_cols <- colnames(obs)
  output_cols <- stats::setNames(
    paste0(prefix, make.unique(make.names(source_cols))),
    source_cols
  )
  srt <- scmf_append_obs_to_srt(
    srt = srt,
    obs = obs,
    source_cols = source_cols,
    output_cols = output_cols,
    verbose = verbose
  )
  srt@tools$scMalignantFinder_states <- list(
    method = "scMalignantFinder",
    task = "cancer_cell_state_scoring",
    columns = unname(output_cols),
    parameters = list(
      assay = assay,
      layer = layer,
      cells = cells,
      gene_sets = gene_sets,
      norm_type = norm_type,
      backend = backend
    )
  )
  srt
}

scmf_read_gene_sets <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parsed <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3L) {
      return(NULL)
    }
    list(
      name = fields[[1L]],
      genes = unique(fields[-c(1L, 2L)])
    )
  })
  parsed <- Filter(Negate(is.null), parsed)
  gene_sets <- lapply(parsed, `[[`, "genes")
  names(gene_sets) <- make.unique(vapply(parsed, `[[`, character(1), "name"))
  gene_sets <- lapply(gene_sets, function(genes) {
    genes[!is.na(genes) & nzchar(genes)]
  })
  gene_sets <- gene_sets[lengths(gene_sets) > 0L]
  if (length(gene_sets) == 0L) {
    log_message(
      "{.arg gene_sets} contains no non-empty GMT records",
      message_type = "error"
    )
  }
  gene_sets
}

scmf_adapt_gene_sets <- function(gene_sets, feature_names, min_fraction = 0.5) {
  adapted <- lapply(gene_sets, function(genes) {
    genes <- unique(as.character(genes))
    overlap <- intersect(genes, feature_names)
    if (length(genes) > 0L && length(overlap) / length(genes) >= min_fraction) {
      overlap
    } else {
      NULL
    }
  })
  adapted <- Filter(Negate(is.null), adapted)
  if (length(adapted) == 0L) {
    log_message(
      "No signatures retain at least 50% gene overlap with the expression matrix",
      message_type = "error"
    )
  }
  adapted
}

scmf_auc_threshold <- function(expr) {
  detected <- if (inherits(expr, "sparseMatrix")) {
    Matrix::colSums(expr != 0)
  } else {
    colSums(expr != 0)
  }
  as.numeric(stats::quantile(
    detected,
    probs = 0.01,
    names = FALSE,
    type = 7
  )) / nrow(expr)
}

scmf_native_region <- function(
  expr,
  metadata,
  signature_gmt,
  features,
  nclus,
  define_feature,
  spatial_nn,
  spatial_coordinates
) {
  expr <- gene_set_scoring_to_dgC(expr)
  gene_sets <- scmf_adapt_gene_sets(
    scmf_read_gene_sets(signature_gmt),
    rownames(expr)
  )
  auc_threshold <- max(1 / nrow(expr), scmf_auc_threshold(expr))
  auc_scores <- run_aucell_scores(
    expr_counts = expr,
    gene_sets = gene_sets,
    strategy = "topk",
    algorithm = "ctxcore",
    seed = 2L,
    tie_method = "numpy",
    auc_threshold = auc_threshold,
    normalize_by_signature_max = TRUE
  )
  obs <- as.data.frame(metadata, stringsAsFactors = FALSE)
  obs <- obs[colnames(expr), , drop = FALSE]
  for (score_name in colnames(auc_scores)) {
    obs[[score_name]] <- auc_scores[, score_name]
  }
  missing_features <- setdiff(unique(c(features, define_feature)), colnames(obs))
  if (length(missing_features) > 0L) {
    log_message(
      "Region feature(s) not found after AUCell scoring: {.val {missing_features}}",
      message_type = "error"
    )
  }
  feature_matrix <- as.matrix(obs[, features, drop = FALSE])
  storage.mode(feature_matrix) <- "double"
  if (ncol(feature_matrix) < 2L || any(!is.finite(feature_matrix))) {
    log_message(
      "{.arg features} must contain at least two finite numeric columns",
      message_type = "error"
    )
  }
  nclus <- suppressWarnings(as.integer(nclus)[1L])
  if (is.na(nclus) || nclus < 1L || nclus > nrow(feature_matrix)) {
    log_message(
      "{.arg nclus} must be between 1 and the number of selected cells",
      message_type = "error"
    )
  }
  cluster <- if (nclus == 1L) {
    rep.int(0L, nrow(feature_matrix))
  } else {
    as.integer(stats::cutree(
      stats::hclust(stats::dist(feature_matrix), method = "ward.D2"),
      k = nclus
    )) - 1L
  }
  cluster_means <- tapply(obs[[define_feature]], cluster, mean)
  malignant_cluster <- as.integer(names(cluster_means)[which.max(cluster_means)])
  prediction <- ifelse(cluster == malignant_cluster, "Malignant", "Normal")

  if (isTRUE(spatial_nn)) {
    if (is.null(spatial_coordinates)) {
      log_message(
        "{.arg spatial_nn = TRUE} with the C++ backend requires {.arg spatial.cols}",
        message_type = "error"
      )
    }
    if (nrow(spatial_coordinates) >= 2L) {
      coords <- data.frame(
        cell_id = rownames(obs),
        x = spatial_coordinates[, 1L],
        y = spatial_coordinates[, 2L]
      )
      graph <- spatial_graph_compute(
        coords = coords,
        method = "knn",
        k = min(6L, nrow(coords) - 1L),
        directed = TRUE,
        weight = "binary"
      )
      raw_prediction <- prediction
      split_neighbors <- split(graph$edges$to, graph$edges$from)
      for (cell_i in seq_along(prediction)) {
        neighbors <- split_neighbors[[as.character(cell_i)]]
        if (length(neighbors) < 2L) {
          next
        }
        neighbor_labels <- raw_prediction[neighbors]
        counts <- table(factor(
          neighbor_labels,
          levels = c("Malignant", "Normal")
        ))
        winner <- which.max(counts)
        if (counts[[winner]] / length(neighbor_labels) > 0.5) {
          prediction[[cell_i]] <- names(counts)[[winner]]
        }
      }
    }
  }
  obs$cluster <- factor(cluster)
  obs$region_prediction <- factor(
    prediction,
    levels = c("Normal", "Malignant")
  )
  keep <- unique(c(
    "cluster",
    "region_prediction",
    colnames(auc_scores),
    features
  ))
  obs[, keep, drop = FALSE]
}

scmf_prepare_python <- function(verbose = TRUE) {
  explicit_python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  python_initialized <- isTRUE(reticulate::py_available(initialize = FALSE))
  if ((nzchar(explicit_python) || python_initialized) &&
    isTRUE(scmf_python_classifier_available())) {
    return(invisible(TRUE))
  }
  envname <- "scmalignantfinder_env"
  PrepareEnv(
    envname = envname,
    modules = "scmalignantfinder",
    verbose = verbose
  )
  ok <- check_python(
    "scMalignantFinder",
    envname = envname,
    pip = TRUE,
    verbose = verbose
  )
  if (isFALSE(ok) || !isTRUE(scmf_python_classifier_available())) {
    log_message(
      "Failed to locate a usable {.pkg scMalignantFinder} classifier module in the active Python environment.",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

scmf_python_classifier_available <- function() {
  isTRUE(tryCatch(
    reticulate::py_module_available("scMalignantFinder.classifier"),
    error = function(e) FALSE
  ))
}

import_scmalignantfinder <- function(convert = TRUE) {
  scop_python_import("scmalignantfinder", convert = convert)
}

scmf_check_one_input <- function(srt = NULL, adata = NULL, h5ad = NULL) {
  n_input <- sum(!vapply(list(srt, adata, h5ad), is.null, logical(1)))
  if (n_input != 1L) {
    log_message(
      "Provide exactly one of {.arg srt}, {.arg adata}, or {.arg h5ad}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

scmf_resolve_norm_type <- function(norm_type, srt = NULL, layer = NULL) {
  if (!is.null(norm_type)) {
    return(isTRUE(norm_type))
  }
  !is.null(srt) && identical(layer, "counts")
}

scmf_expand_path <- function(path) {
  if (is.null(path)) {
    return(NULL)
  }
  path.expand(path)
}

scmf_get_spatial_coordinates <- function(srt, spatial.cols = NULL) {
  if (is.null(spatial.cols)) {
    return(NULL)
  }
  if (length(spatial.cols) != 2L || anyNA(spatial.cols) || any(!nzchar(spatial.cols))) {
    log_message("{.arg spatial.cols} must contain exactly two metadata column names", message_type = "error")
  }
  meta <- srt[[]]
  missing_spatial <- setdiff(spatial.cols, colnames(meta))
  if (length(missing_spatial) > 0) {
    log_message(
      "{.arg spatial.cols} not found in metadata: {.val {missing_spatial}}",
      message_type = "error"
    )
  }
  spatial_df <- meta[, spatial.cols, drop = FALSE]
  is_numeric <- vapply(spatial_df, is.numeric, logical(1))
  if (!all(is_numeric)) {
    log_message("{.arg spatial.cols} must be finite numeric metadata columns", message_type = "error")
  }
  spatial_coordinates <- as.matrix(spatial_df)
  if (any(!is.finite(spatial_coordinates))) {
    log_message("{.arg spatial.cols} must be finite numeric metadata columns", message_type = "error")
  }
  storage.mode(spatial_coordinates) <- "double"
  spatial_coordinates
}

scmf_cells <- function(srt, cells = NULL) {
  if (is.null(cells)) {
    return(colnames(srt))
  }
  cells <- as.character(cells)
  missing_cells <- setdiff(cells, colnames(srt))
  if (length(missing_cells) > 0) {
    log_message(
      "{.arg cells} contains cells not found in {.arg srt}: {.val {missing_cells}}",
      message_type = "error"
    )
  }
  cells
}


scmf_region_output_names <- function(source_cols, prefix = "scMalignantFinder_") {
  base_names <- source_cols
  base_names[base_names == "cluster"] <- "region_cluster"
  base_names[base_names == "region_prediction"] <- "region_prediction"
  stats::setNames(paste0(prefix, make.names(base_names)), source_cols)
}

scmf_append_obs_to_srt <- function(
  srt,
  obs,
  source_cols,
  output_cols,
  verbose = TRUE
) {
  obs <- as.data.frame(obs)
  missing_cols <- setdiff(source_cols, colnames(obs))
  if (length(missing_cols) > 0) {
    log_message(
      "{.pkg scMalignantFinder} did not return expected column(s): {.val {missing_cols}}",
      message_type = "error"
    )
  }
  if (is.null(rownames(obs)) || anyNA(rownames(obs)) || any(!nzchar(rownames(obs)))) {
    log_message(
      "{.pkg scMalignantFinder} results must have cell barcodes as row names",
      message_type = "error"
    )
  }

  matched_cells <- intersect(colnames(srt), rownames(obs))
  if (length(matched_cells) == 0) {
    log_message(
      "Unable to align {.pkg scMalignantFinder} results to Seurat cells",
      message_type = "error"
    )
  }
  ignored_cells <- setdiff(rownames(obs), colnames(srt))
  if (length(ignored_cells) > 0) {
    log_message(
      "Ignoring {.val {length(ignored_cells)}} returned cell{?s} not present in {.arg srt}",
      message_type = "warning",
      verbose = verbose
    )
  }

  for (col in source_cols) {
    col_values <- obs[[col]]
    if (is.factor(col_values)) {
      col_values <- as.character(col_values)
    }
    names(col_values) <- rownames(obs)
    values <- scmf_empty_vector(col_values, ncol(srt))
    names(values) <- colnames(srt)
    values[matched_cells] <- col_values[matched_cells]
    srt[[unname(output_cols[[col]])]] <- values
  }
  srt
}

scmf_empty_vector <- function(example, n) {
  if (is.integer(example)) {
    return(rep(NA_integer_, n))
  }
  if (is.numeric(example)) {
    return(rep(NA_real_, n))
  }
  if (is.logical(example)) {
    return(rep(NA, n))
  }
  rep(NA_character_, n)
}
