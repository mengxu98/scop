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
  if (!((nzchar(Sys.getenv("RETICULATE_PYTHON")) || reticulate::py_available(initialize = FALSE)) &&
    isTRUE(scmf_python_classifier_available()))) {
    PrepareEnv(modules = "scanpy", verbose = verbose)
    if (!isTRUE(scmf_python_classifier_available())) {
      ok <- check_python("scMalignantFinder", pip = TRUE, verbose = verbose)
      if (isFALSE(ok) || !isTRUE(scmf_python_classifier_available())) {
        scmf_install_python_github(verbose = verbose)
      }
    }
    if (!isTRUE(scmf_python_classifier_available())) {
      log_message(
        "Failed to install or locate a usable {.pkg scMalignantFinder} classifier module. Install it manually in the active {.pkg scop} Python environment.",
        message_type = "error"
      )
    }
  }
  if (is.null(pretrain_dir) && identical(model_method, "XGBoost")) {
    scmf_check_xgboost_python(verbose = verbose)
  }

  if (!is.null(srt)) {
    scmf_assert_seurat(srt)
    cells <- scmf_cells(srt, cells)
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
  output_cols <- scmf_prediction_output_names(source_cols, prefix = prefix)
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
  verbose = TRUE
) {
  scmf_check_one_input(srt = srt, adata = adata, h5ad = h5ad)
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
    scmf_assert_seurat(srt)
    precheck_cells <- scmf_cells(srt, cells)
    precheck_srt <- if (length(precheck_cells) == ncol(srt)) {
      srt
    } else {
      subset(srt, cells = precheck_cells)
    }
    scmf_get_spatial_coordinates(precheck_srt, spatial.cols)
  }

  log_message(
    "Running {.pkg scMalignantFinder} malignant region identification...",
    message_type = "running",
    verbose = verbose
  )
  if (!((nzchar(Sys.getenv("RETICULATE_PYTHON")) || reticulate::py_available(initialize = FALSE)) &&
    isTRUE(scmf_python_classifier_available()))) {
    PrepareEnv(modules = "scanpy", verbose = verbose)
    if (!isTRUE(scmf_python_classifier_available())) {
      ok <- check_python("scMalignantFinder", pip = TRUE, verbose = verbose)
      if (isFALSE(ok) || !isTRUE(scmf_python_classifier_available())) {
        scmf_install_python_github(verbose = verbose)
      }
    }
    if (!isTRUE(scmf_python_classifier_available())) {
      log_message(
        "Failed to install or locate a usable {.pkg scMalignantFinder} classifier module. Install it manually in the active {.pkg scop} Python environment.",
        message_type = "error"
      )
    }
  }

  spatial_coordinates <- NULL
  if (!is.null(srt)) {
    scmf_assert_seurat(srt)
    cells <- scmf_cells(srt, cells)
    srt_input <- if (length(cells) == ncol(srt)) {
      srt
    } else {
      subset(srt, cells = cells)
    }
    spatial_coordinates <- scmf_get_spatial_coordinates(srt_input, spatial.cols)
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
      norm_type = norm_type
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
  verbose = TRUE
) {
  scmf_check_one_input(srt = srt, adata = adata, h5ad = h5ad)
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
  if (!((nzchar(Sys.getenv("RETICULATE_PYTHON")) || reticulate::py_available(initialize = FALSE)) &&
    isTRUE(scmf_python_classifier_available()))) {
    PrepareEnv(modules = "scanpy", verbose = verbose)
    if (!isTRUE(scmf_python_classifier_available())) {
      ok <- check_python("scMalignantFinder", pip = TRUE, verbose = verbose)
      if (isFALSE(ok) || !isTRUE(scmf_python_classifier_available())) {
        scmf_install_python_github(verbose = verbose)
      }
    }
    if (!isTRUE(scmf_python_classifier_available())) {
      log_message(
        "Failed to install or locate a usable {.pkg scMalignantFinder} classifier module. Install it manually in the active {.pkg scop} Python environment.",
        message_type = "error"
      )
    }
  }

  if (!is.null(srt)) {
    scmf_assert_seurat(srt)
    cells <- scmf_cells(srt, cells)
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
      norm_type = norm_type
    )
  )
  srt
}

scmf_install_python_github <- function(verbose = TRUE) {
  python <- tryCatch(reticulate::py_config()$python, error = function(e) "")
  if (!nzchar(python)) {
    python <- Sys.which("python3")
  }
  if (!nzchar(python)) {
    return(invisible(FALSE))
  }
  log_message(
    "Installing {.pkg scMalignantFinder} from GitHub with active Python...",
    message_type = "running",
    verbose = verbose
  )
  out <- tryCatch(
    system2(
      python,
      c(
        "-m", "pip", "install",
        "--no-cache-dir",
        "git+https://github.com/Jonyyqn/scMalignantFinder.git"
      ),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) structure(conditionMessage(e), status = 1L)
  )
  status <- attr(out, "status") %||% 0L
  if (!identical(status, 0L)) {
    log_message(
      "GitHub installation of {.pkg scMalignantFinder} failed: {.val {paste(out, collapse = '; ')}}",
      message_type = "warning",
      verbose = verbose
    )
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

scmf_python_classifier_available <- function() {
  code <- paste(
    "import importlib.util, pathlib, sys",
    "spec = importlib.util.find_spec('scMalignantFinder')",
    "assert spec is not None and spec.submodule_search_locations",
    "pkg_dir = pathlib.Path(list(spec.submodule_search_locations)[0])",
    "module_path = pkg_dir / 'classifier.py'",
    "assert module_path.exists()",
    "module_spec = importlib.util.spec_from_file_location('_scop_scmf_classifier_check', module_path)",
    "module = importlib.util.module_from_spec(module_spec)",
    "sys.modules['_scop_scmf_classifier_check'] = module",
    "module_spec.loader.exec_module(module)",
    "assert hasattr(module, 'scMalignantFinder')",
    sep = "\n"
  )
  isTRUE(tryCatch({
    reticulate::py_run_string(code, local = TRUE)
    TRUE
  }, error = function(e) FALSE))
}

scmf_check_xgboost_python <- function(verbose = TRUE) {
  ok <- check_python("xgboost", pip = TRUE, verbose = verbose)
  if (isFALSE(ok)) {
    log_message(
      "{.pkg xgboost} is required when training {.pkg scMalignantFinder} with {.val XGBoost}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

import_scmalignantfinder <- function(convert = TRUE) {
  reticulate::import_from_path(
    "scmalignantfinder",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = convert
  )
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

scmf_assert_seurat <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  invisible(TRUE)
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

scmf_prediction_output_names <- function(source_cols, prefix = "") {
  stats::setNames(
    if (nzchar(prefix)) paste0(prefix, source_cols) else source_cols,
    source_cols
  )
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
